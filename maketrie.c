#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>
#include <omp.h>
#include <ctype.h>
#include <time.h>

#include <assert.h>

#ifndef FALSE
#define TRUE (0==0)
#define FALSE (!TRUE)
#endif

static FILE *duplicates = NULL;
static FILE *rejects = NULL;
static FILE *read_file = NULL;
static FILE *read_index = NULL;
static FILE *trie_file = NULL;
static FILE *sorted_and_unique_reads = NULL;

#define GBs 4L

typedef unsigned long long EDGE;
typedef unsigned long long INDEX;

#define ENDS_WORD (1ULL<<63ULL)
#define EDGE_MASK (ENDS_WORD-1UL)

typedef struct cell
{
   EDGE edge[5];
} CELL;
CELL *trie_cell;

static INDEX MAX_SIZE = ((INDEX) 0L);

#define ROOT_CELL ((INDEX)1L)

static INDEX last_used_edge = ROOT_CELL;

static CELL empty;

#define _A_ 0
#define _C_ 1
#define _G_ 2
#define _T_ 3
#define _N_ 4

char *trt = "ACGTN";

static long freq[256];
static long letters = 0L;
static int seq = 0, dups = 0;

#define MAX_LINE 1024

static long length[MAX_LINE];
static int read_length = 0;

#define CORES_PER_NODE 16ULL

static int mpirank = 0, mpisize = 1;

static long long TASKS_PER_NODE, PROCESSORS_PER_NODE;

#define TAG_DATA 1
#define TAG_ACK 2

#define TAG_SEND_RAW_MEM 3
#define TAG_READ_READ 4
#define TAG_WRITE_READ 5
#define TAG_EXIT_PROGRAM 6

#define TAG_ADD_READ 7
#define TAG_GET_NEXT_FREE_EDGE 8
#define TAG_OUTPUT_DUPINFO 9
#define TAG_OUTPUT_READ 10
#define TAG_WALK_AND_PRINT_TRIE_INTERNAL 13
#define TAG_DUMP_TRIE 14

static long long CHUNKBITS, CHUNKSIZE, CHUNKMASK;

static void shut_down_other_nodes (void)
{
   int target_rank;
   long int value = 0L;
   MPI_Status status;

   for (target_rank = 1; target_rank < mpisize; target_rank++) {
      if (target_rank != mpirank) {
         MPI_Send (&value, 1, MPI_LONG, target_rank, TAG_EXIT_PROGRAM,
                   MPI_COMM_WORLD);
         MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &status);
      }
   }
}

static void send_bytes (int dest, void *memp, long bytes)
{
   MPI_Send (memp, bytes, MPI_BYTE, dest, TAG_SEND_RAW_MEM, MPI_COMM_WORLD);
}

static void remote_walk_and_print_trie_internal (int target_rank, char *s,
                                                 EDGE edge, int len)
{
   MPI_Status status;
   long value = 0L;
   long stringlength;

   stringlength = strlen (s) + 1;
   MPI_Send (&value, 1, MPI_LONG, target_rank,
             TAG_WALK_AND_PRINT_TRIE_INTERNAL, MPI_COMM_WORLD);
   MPI_Send (&stringlength, 1, MPI_LONG, target_rank, TAG_DATA,
             MPI_COMM_WORLD);
   send_bytes (target_rank, s, stringlength);
   MPI_Send (&edge, 1, MPI_LONG, target_rank, TAG_DATA, MPI_COMM_WORLD);
   MPI_Send (&len, 1, MPI_INT, target_rank, TAG_DATA, MPI_COMM_WORLD);
   MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
}

static void remote_dump_trie (int target_rank, char *filename)
{
   MPI_Status status;
   long value = 0L;
   long stringlength;

   stringlength = strlen (filename) + 1;
   MPI_Send (&value, 1, MPI_LONG, target_rank, TAG_DUMP_TRIE, MPI_COMM_WORLD);
   MPI_Send (&stringlength, 1, MPI_LONG, target_rank, TAG_DATA,
             MPI_COMM_WORLD);
   send_bytes (target_rank, filename, stringlength);
   MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
}

static void remote_setread (int target_rank, INDEX index, CELL value)
{
   long int dummy;
   MPI_Status status;

   assert (target_rank > mpirank);
   MPI_Send (&dummy, 1, MPI_LONG, target_rank, TAG_WRITE_READ,
             MPI_COMM_WORLD);
   MPI_Send (&index, 1, MPI_LONG_LONG, target_rank, TAG_DATA, MPI_COMM_WORLD);
   assert (sizeof (value.edge[0]) == sizeof (long));
   MPI_Send (&value.edge, sizeof (value.edge) / sizeof (value.edge[0]),
             MPI_LONG, target_rank, TAG_DATA, MPI_COMM_WORLD);
   MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
}

static void remote_getread (int target_rank, INDEX index, CELL * valuep)
{
   MPI_Status status;
   long dummy;

   MPI_Send (&dummy, 1, MPI_LONG, target_rank, TAG_READ_READ, MPI_COMM_WORLD);
   MPI_Send (&index, 1, MPI_LONG_LONG, target_rank, TAG_DATA, MPI_COMM_WORLD);
   assert (sizeof (valuep->edge[0]) == sizeof (long));
   MPI_Recv (valuep, sizeof (valuep->edge) / sizeof (valuep->edge[0]),
             MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
}

static void setread (INDEX index, CELL value)
{
int target_rank;

   target_rank = index >> CHUNKBITS;
   if (target_rank == mpirank) {
      int i;

      for (i = 0; i < 5; i++) trie_cell[index & CHUNKMASK].edge[i] = value.edge[i];
   } else {
      if (target_rank < mpirank) {
         fprintf (stderr,
                  "PROGRAM BUG: Node %d requested access to trie_cell[%lld] on node %d"
                  " - assert that we never feed backwards...\n",
                  mpirank, index, target_rank);
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
      if (target_rank >= mpisize) {
         fprintf (stderr,
                  "ERROR: array bounds exceeded!  Requested access to trie_cell[%lld]\n",
                  index);
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
      remote_setread (target_rank, index, value);
   }
}

static void getread (INDEX index, CELL * valuep)
{
int target_rank;

   target_rank = index >> CHUNKBITS;
   if (target_rank == mpirank) {
      int i;

      for (i = 0; i < 5; i++) valuep->edge[i] = trie_cell[index & CHUNKMASK].edge[i];
   } else {
      if (target_rank < mpirank) {
         fprintf (stderr,
                  "PROGRAM BUG: Node %d requested access to trie_cell[%lld] on node %d"
                  " - assert that we never feed backwards...\n",
                  mpirank, index, target_rank);
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
      if (target_rank >= mpisize) {
         fprintf (stderr,
                  "ERROR: array bounds exceeded!  Requested access to trie_cell[%lld]\n",
                  index);
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
      remote_getread (target_rank, index, valuep);
   }
}

static int add_read (char *s, EDGE edge, long read_number, int len);
static INDEX get_next_free_edge (void);

static int local_add_read (char *s, EDGE edge, long read_number, int len)
{
   int c;

   assert ((edge >> CHUNKBITS) == mpirank);
   assert ((*s != '\0') && (*s != '\n') && (*s != '\r'));

   c = *s++;

   letters++;
   freq[c]++;

   if (c == 'A') c = _A_;
   else if (c == 'C') c = _C_;
   else if (c == 'G') c = _G_;
   else if (c == 'T') c = _T_;
   else c = _N_;

   if ((*s == '\0') || (*s == '\n') || (*s == '\r')) {
      if (trie_cell[edge & CHUNKMASK].edge[c] & ENDS_WORD) {
         long original_read = trie_cell[edge & CHUNKMASK].edge[c] & EDGE_MASK;

         fprintf (duplicates, "%ld:0 %ld\n", original_read, read_number);
         if (ferror (duplicates)) {
            fprintf (stderr,
                     "\n\n************* add_read() (duplicates) failed, %s\n",
                     strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
         dups++;
      } else {
         assert ((trie_cell[edge & CHUNKMASK].edge[c] & EDGE_MASK) == 0);
         trie_cell[edge & CHUNKMASK].edge[c] = (ENDS_WORD | read_number);
      }

      return len + 1;
   }

   if (trie_cell[edge & CHUNKMASK].edge[c] == 0LL) {

      INDEX new_edge = get_next_free_edge ();

      setread (new_edge, empty);

      trie_cell[edge & CHUNKMASK].edge[c] = new_edge;
      if (new_edge >= MAX_SIZE) {
         fprintf (stderr,
                  "Ran out of free edges after %d reads (last_used_edge = %lld, MAX_SIZE = %lld)\n",
                  seq, new_edge, MAX_SIZE);
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
   } else {

   }
   return add_read (s, trie_cell[edge & CHUNKMASK].edge[c], read_number,
                    len + 1);

}
static INDEX remote_get_next_free_edge (target_rank)
{
   MPI_Status status;
   long value = 0L;
   INDEX free_edge;

   assert (target_rank != mpirank);
   MPI_Send (&value, 1, MPI_LONG, target_rank, TAG_GET_NEXT_FREE_EDGE,
             MPI_COMM_WORLD);
   assert (sizeof (INDEX) == sizeof (long long));
   MPI_Recv (&free_edge, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);

   return free_edge;
}

static int remote_add_read (long target_rank, char *s, long edge,
                            long read_number, int len)
{
   MPI_Status status;
   long value = 0L;
   long stringlength;

   assert (target_rank != mpirank);
   stringlength = strlen (s) + 1;
   MPI_Send (&value, 1, MPI_LONG, target_rank, TAG_ADD_READ, MPI_COMM_WORLD);
   MPI_Send (&stringlength, 1, MPI_LONG, target_rank, TAG_DATA,
             MPI_COMM_WORLD);
   send_bytes (target_rank, s, stringlength);
   assert (sizeof (EDGE) == sizeof (long));
   MPI_Send (&edge, 1, MPI_LONG, target_rank, TAG_DATA, MPI_COMM_WORLD);
   assert (sizeof (read_number) == sizeof (long));
   MPI_Send (&read_number, 1, MPI_LONG, target_rank, TAG_DATA,
             MPI_COMM_WORLD);
   assert (sizeof (len) == sizeof (int));
   MPI_Send (&len, 1, MPI_INT, target_rank, TAG_DATA, MPI_COMM_WORLD);
   MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);

   return (int) value;
}

static void remote_output_read (long target_rank, char *s, EDGE readindex)
{
   MPI_Status status;
   long value = 0L;
   long stringlength;

   assert (target_rank != mpirank);
   stringlength = strlen (s) + 1;
   MPI_Send (&value, 1, MPI_LONG, target_rank, TAG_OUTPUT_READ,
             MPI_COMM_WORLD);
   MPI_Send (&stringlength, 1, MPI_LONG, target_rank, TAG_DATA,
             MPI_COMM_WORLD);
   send_bytes (target_rank, s, stringlength);
   MPI_Send (&readindex, 1, MPI_LONG_LONG, target_rank, TAG_DATA,
             MPI_COMM_WORLD);
   MPI_Recv (&value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
}

static void accept_get_next_free_edge (int caller)
{
   INDEX new_edge;

   new_edge = get_next_free_edge ();
   MPI_Send (&new_edge, 1, MPI_LONG_LONG, caller, 0, MPI_COMM_WORLD);
}

static void accept_add_read (int myrank, long value, MPI_Status status)
{
   char *s;
   EDGE edge;
   long read_number, stringlength;
   int len;
   int caller;

   MPI_Recv (&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   caller = status.MPI_SOURCE;

   s = malloc (stringlength);
   if (s == NULL) fprintf (stderr, "CRAP!  Failed to malloc... why\?\?\?\n");
   MPI_Recv (s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   assert (sizeof (EDGE) == sizeof (long));
   MPI_Recv (&edge, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
   assert (sizeof (read_number) == sizeof (long));
   MPI_Recv (&read_number, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   assert (sizeof (len) == sizeof (int));
   MPI_Recv (&len, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);

   {
      long long int target_rank = (long long int) edge >> CHUNKBITS;
      assert (target_rank == mpirank);
   }

   value = add_read (s, edge, read_number, len);
   free (s);
   s = NULL;

   MPI_Send (&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);
}

static void walk_and_print_trie_internal (char *s, EDGE edge, int len);
static void dump_trie (char *filename);
static void accept_walk_and_print_trie_internal (int myrank, long value,
                                                 MPI_Status status)
{
   char s[MAX_LINE];
   EDGE edge;
   long stringlength;
   int len;
   int caller;

   MPI_Recv (&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   caller = status.MPI_SOURCE;
   MPI_Recv (s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   MPI_Recv (&edge, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
   MPI_Recv (&len, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
   walk_and_print_trie_internal (s, edge, len);
   MPI_Send (&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);
}

static void accept_dump_trie (int myrank, long value, MPI_Status status)
{
   char filename[MAX_LINE];
   long stringlength;
   int caller;

   MPI_Recv (&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   caller = status.MPI_SOURCE;
   MPI_Recv (filename, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   dump_trie (filename);
   MPI_Send (&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);
}

static void output_read (char *s, EDGE readindex);
static void accept_output_read (int myrank, long value, MPI_Status status)
{
   char *s;
   EDGE readindex;
   long stringlength;
   int caller;

   MPI_Recv (&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   caller = status.MPI_SOURCE;

   s = malloc (stringlength);
   if (s == NULL) fprintf (stderr, "CRAP!  Failed to malloc... why\?\?\?\n");
   MPI_Recv (s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
   MPI_Recv (&readindex, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

   output_read (s, readindex);

   free (s);
   s = NULL;

   MPI_Send (&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);
}

static void output_read (char *s, EDGE readindex)
{
   time_t curtime;

   if (mpirank == mpisize - 1) {
      static int printed = 0;

      fprintf (sorted_and_unique_reads, "%s %12lld\n", s, readindex);
      if (ferror (sorted_and_unique_reads)) {
         fprintf (stderr, "\n\n************* output_read() failed, %s\n",
                  strerror (errno));
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
      printed++;
      if ((printed % 1000000) == 0) {
         time (&curtime);
         fprintf (stderr, "%d unique and sorted reads written back at %s",
                  printed, ctime (&curtime));
      }
   } else {
      remote_output_read (mpisize - 1, s, readindex);
   }
}

static INDEX get_next_free_edge (void)
{
   if (((last_used_edge + (INDEX) 1) >> CHUNKBITS) != mpirank) {
      INDEX edge;
      static int next_guy = 1;
      static int init = FALSE;

      if (!init) {
         next_guy = mpirank + 1;
         init = TRUE;
         if (mpirank == mpisize) {
            fprintf (stderr,
                     "ERROR: not enough RAM for this input file (%d * %lld cells used)."
                     "  Try resubmtting with some more processors.\n",
                     mpisize, CHUNKSIZE);
         }
      }
      edge = remote_get_next_free_edge (next_guy);
      next_guy = (edge >> CHUNKBITS);
      return edge;
   } else
      return ++last_used_edge;
}

static int add_read (char *s, EDGE edge, long read_number, int len)
{
   int len2;
   long long int target_rank = (long long) edge >> CHUNKBITS;

   if (len == 0) {
      assert (edge == ROOT_CELL);
      if (read_number > EDGE_MASK) {
         fprintf (stderr, "maketrie: too many READs! (%ld)  Limit is %lld\n",
                  read_number, EDGE_MASK);
         assert (read_number <= EDGE_MASK);
      }
      seq++;
   }

   if (target_rank == mpirank) {
      len2 = local_add_read (s, edge, read_number, len);
      if (len == 0) length[len2]++;
      return len2;
   } else {
      return remote_add_read ((long) target_rank, s, edge, read_number, len);
   }

}

static void walk_and_print_trie_internal (char *s, EDGE edge, int len)
{
   int i;
   int target_rank = edge >> CHUNKBITS;

   if (target_rank != mpirank) {
      remote_walk_and_print_trie_internal (target_rank, s, edge, len);
      return;
   }

   s[len + 1] = '\0';
   for (i = 0; i < 5; i++) {
      s[len] = trt[i];
      if (trie_cell[edge & CHUNKMASK].edge[i] & ENDS_WORD) {
         output_read (s, trie_cell[edge & CHUNKMASK].edge[i] & EDGE_MASK);
      } else if (trie_cell[edge & CHUNKMASK].edge[i]) {
         walk_and_print_trie_internal (s, trie_cell[edge & CHUNKMASK].edge[i],
                                       len + 1);
      }
   }

}

static void dump_trie (char *filename)
{
   time_t curtime;
   int rc;

   fprintf (stderr, "maketrie[%d]: ", mpirank);
   if (mpirank == 0) {
      trie_file = fopen (filename, "w");
      fprintf (stderr, "Writing");
   } else {
      trie_file = fopen (filename, "a");
      fprintf (stderr, "Appending");
   }
   time (&curtime);
   fprintf (stderr, " to dumped trie %s at %s\n", filename, ctime (&curtime));
   if (trie_file == NULL) {
      fprintf (stderr, "maketrie[%d]: Cannot save trie to %s - %s\n", mpirank,
               filename, strerror (errno));
      shut_down_other_nodes ();
      MPI_Finalize ();
      exit (EXIT_FAILURE);
   }
   if (last_used_edge == (mpirank * CHUNKSIZE - 1)) return;

   fwrite (trie_cell, last_used_edge + 1 - (mpirank * CHUNKSIZE),
           sizeof (CELL), trie_file);
   if (ferror (trie_file)) {
      fprintf (stderr, "maketrie[%d]: Error saving trie to %s - %s\n",
               mpirank, filename, strerror (errno));
      shut_down_other_nodes ();
      MPI_Finalize ();
      exit (EXIT_FAILURE);
   }
   time (&curtime);
   fprintf (stderr, "Written at %s\n", ctime (&curtime));
   if (trie_file) {
      rc = fclose (trie_file);
      trie_file = NULL;
      if (rc == EOF) {
         fprintf (stderr, "maketrie[%d]: Error saving trie to %s - %s\n",
                  mpirank, filename, strerror (errno));
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
   }

   if (last_used_edge == ((mpirank + 1) * CHUNKSIZE - 1)) {
      fprintf (stderr, "Not done.  Asking next rank to continue...\n");

      if (mpirank != mpisize - 1) remote_dump_trie (mpirank + 1, filename);
   }
}

static void walk_and_print_trie (void)
{
   char s[MAX_LINE];
   time_t curtime;

   time (&curtime);
   fprintf (stderr, "Printing sorted reads at %s", ctime (&curtime));
   walk_and_print_trie_internal (s, ROOT_CELL, 0);
   if (mpirank == mpisize - 1) {
      if (sorted_and_unique_reads) {
         int rc = fclose (sorted_and_unique_reads);
         sorted_and_unique_reads = NULL;
         if (rc == EOF) {
            fprintf (stderr,
                     "maketrie[%d]: Error closing sorted output - %s\n",
                     mpirank, strerror (errno));
         }
      }
   }
   time (&curtime);
   fprintf (stderr, "Printing sorted reads complete at %s", ctime (&curtime));
}

int main (int argc, char **argv)
{
   time_t curtime;
   char fname[1024];
   char line[MAX_LINE];
   int i, c, rc, lineno = 1, number_of_lengths = 0;
   long read_number = 0;
   int namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];

   time (&curtime);
   fprintf (stderr, "Program started at %s", ctime (&curtime));

   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &mpisize);
   MPI_Comm_rank (MPI_COMM_WORLD, &mpirank);
   fprintf (stderr, "I am rank %d of world size %d\n", mpirank, mpisize);
   MPI_Get_processor_name (processor_name, &namelen);
   if (processor_name && strchr (processor_name, '.')) *strchr (processor_name, '.') = '\0';

   if ((mpirank == 0) && (argc > 2)) {
      fprintf (stderr, "warning: extra parameter %s ignored...\n", argv[2]);
   }

   if (argc >= 2) {

      if (mpirank == 0) {
         read_file = fopen (argv[1], "r");
         if (read_file == NULL) {
            fprintf (stderr, "maketrie: cannot open input \"%s\" - %s\n",
                     argv[1], strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
         fprintf (stderr, "Input: %s\n", argv[1]);
      }

      assert (strlen (argv[1]) + 10 < 1024);

      sprintf (fname, "%s-dups-%05d", argv[1], mpirank);
      duplicates = fopen (fname, "w");
      if (duplicates == NULL) {
         fprintf (stderr,
                  "maketrie on rank %d: cannot open output \"%s\" - %s\n",
                  mpirank, fname, strerror (errno));
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      }
      fprintf (stderr, "Output: %s\n", fname);

      if (mpirank == mpisize - 1) {
         sprintf (fname, "%s-sorted", argv[1]);
         sorted_and_unique_reads = fopen (fname, "w");
         if (sorted_and_unique_reads == NULL) {
            fprintf (stderr,
                     "maketrie[%d]: cannot create sorted output \"%s\" - %s\n",
                     mpirank, fname, strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
         fprintf (stderr, "Sorted READ Output: %s\n", fname);

         sprintf (fname, "%s-rejects", argv[1]);
         rejects = fopen (fname, "w");
         if (rejects == NULL) {
            fprintf (stderr,
                     "maketrie[%d]: cannot create reject file \"%s\" - %s\n",
                     mpirank, fname, strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
         fprintf (stderr, "Reject Output: %s\n", fname);
      }

      if (mpirank == 0) {
         sprintf (fname, "%s-index", argv[1]);
         read_index = fopen (fname, "w");
         if (read_index == NULL) {
            fprintf (stderr, "maketrie: cannot create index \"%s\" - %s\n",
                     fname, strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
         fprintf (stderr, "READ Index Output: %s\n", fname);
      }

   } else {
      if (mpirank == 0) fprintf (stderr, "syntax: maketrie input.fastq\n");
      MPI_Finalize ();
      exit (EXIT_FAILURE);
   }

   CHUNKBITS = (((INDEX) sizeof (INDEX)) * 8ULL) - 1ULL;
   {
      FILE *meminfo, *cpuinfo;
      long long memsize;
      static char line[1024];

      PROCESSORS_PER_NODE = 0ULL;
      cpuinfo = fopen ("/proc/cpuinfo", "r");
      if (cpuinfo) {
         for (;;) {
            char *s = fgets (line, 1023, cpuinfo);
            if (s == NULL) break;
            if (strncmp (line, "processor", strlen ("processor")) == 0) {
               PROCESSORS_PER_NODE++;
            }
         }
         fprintf (stderr, "Discovered %lld processors per node\n",
                  PROCESSORS_PER_NODE);
      }
      if (PROCESSORS_PER_NODE == 0ULL) PROCESSORS_PER_NODE = CORES_PER_NODE;

      TASKS_PER_NODE = PROCESSORS_PER_NODE / (long long) omp_get_max_threads ();
      fprintf (stderr,
               "Node %s, Rank %d, and running %lld ranks on this node.   <-------------------------------\n",
               processor_name, mpirank, TASKS_PER_NODE);

      meminfo = fopen ("/proc/meminfo", "r");
      if (meminfo) for (;;) {
            int count;
            char *s = fgets (line, 1023, meminfo);

            if (s == NULL) break;

            count = sscanf (line, "MemTotal:     %lld kB", &memsize);
            if (count == 1) {

               memsize *= 1024ULL;
               memsize /= TASKS_PER_NODE;
               memsize /= (long long) sizeof (CELL);

               CHUNKBITS = 1ULL;
               for (;;) {
                  if ((1ULL << CHUNKBITS) >= memsize) break;
                  CHUNKBITS++;
               }
               CHUNKBITS -= 1ULL;
               fprintf (stderr,
                        "rounding down memsize to %lldM cells per core"
                        " (%lld bits), ie %lldM cells per node\n",
                        (1ULL << CHUNKBITS) >> 24ULL, CHUNKBITS,
                        ((1ULL << CHUNKBITS) * TASKS_PER_NODE) >> 24ULL);
               break;
            } else {
               CHUNKBITS = CHUNKBITS >> 1ULL;
            }
         }
   }
   CHUNKSIZE = (1ULL << CHUNKBITS);
   CHUNKMASK = (CHUNKSIZE - 1ULL);

   while (CHUNKBITS >= 16ULL) {
      fprintf (stderr,
               "Node %d: trying calloc of %lld cells of %d bytes each.\n",
               mpirank, CHUNKSIZE, (int) sizeof (CELL));
      trie_cell = calloc (CHUNKSIZE, sizeof (CELL));
      if (trie_cell == NULL) {
         CHUNKBITS -= 1ULL;
         CHUNKSIZE = (1ULL << CHUNKBITS);
         CHUNKMASK = (CHUNKSIZE - 1ULL);
      } else
         break;
   }

#ifdef MULTINODE_DEBUG100

   free (trie_cell);
   CHUNKBITS = 8;
   CHUNKSIZE = (1ULL << CHUNKBITS);
   CHUNKMASK = (CHUNKSIZE - 1ULL);
   trie_cell = calloc (CHUNKSIZE, sizeof (CELL));
#endif

#ifdef MULTINODE_DEBUG1K

   free (trie_cell);
   CHUNKBITS = 9;
   CHUNKSIZE = (1ULL << CHUNKBITS);
   CHUNKMASK = (CHUNKSIZE - 1ULL);
   trie_cell = calloc (CHUNKSIZE, sizeof (CELL));
#endif

   if (trie_cell == NULL) {
      fprintf (stderr,
               "rpctest: rank %d unable to allocate array of %lld longs\n",
               mpirank, CHUNKSIZE);
      shut_down_other_nodes ();
      MPI_Finalize ();
      exit (EXIT_FAILURE);
   } else {
      fprintf (stderr, "Node %d: allocated %lld-item long array\n", mpirank,
               CHUNKSIZE);
   }

   fprintf (stderr,
            "node %d: using %dM-items.  Launching listener now.\n",
            mpirank, (int) (CHUNKSIZE >> 24ULL)
    );

   MAX_SIZE = CHUNKSIZE * mpisize;
   fprintf (stderr, "setting MAX_SIZE to %lld (%lld * %d)\n", MAX_SIZE,
            CHUNKSIZE, mpisize);

   for (i = 0; i < 256; i++) freq[i] = 0;
   for (i = 0; i < MAX_LINE; i++) length[i] = 0;
   for (i = 0; i < 5; i++) empty.edge[i] = 0LL;

   for (i = 0; i < 5; i++) trie_cell[ROOT_CELL].edge[i] = 0LL;
#ifdef TWONODE_DEBUG
   last_used_edge = CHUNKSIZE - 100ULL;
#endif
   fprintf (stderr,
            "last_used_edge: %lld,  CHUNKSIZE: %lld,  MAX_SIZE: %lld\n",
            last_used_edge, CHUNKSIZE, MAX_SIZE);

   if (mpirank == 0) {

      fprintf (stderr,
               "\nCombined system is using %lldM trie edges distributed across %d ranks\n\n",
               (long long) mpisize * (CHUNKSIZE >> 24ULL), mpisize);

      for (;;) {
         int len;
         off_t read_start;
         char *s;

         read_start = ftello (read_file);
         if (read_index) fwrite (&read_start, sizeof (off_t), 1, read_index);

         s = fgets (line, MAX_LINE, read_file);
         if (s == NULL) break;

         lineno++;
         fgets (line, MAX_LINE, read_file);
         s = strchr (line, '\n');
         if (s) *s = '\0';
         len = add_read (line, ROOT_CELL, read_number++, 0);

         lineno++;
         fgets (line, MAX_LINE, read_file);
         if (line[0] != '+') {
            fprintf (stderr, "Input data format error in READ file line %d\n",
                     lineno);
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }

         lineno++;
         fgets (line, MAX_LINE, read_file);
         lineno++;
         if ((read_number % 1000000) == 0) {
            time_t curtime;
            time (&curtime);
            fprintf (stderr, "%ld READs loaded at %s", read_number,
                     ctime (&curtime));
         }
         if (read_number == 0x7FFFFFFF) {
            fprintf (stderr,
                     "maketrie: an assumption was wrong.  We have an input file "
                     "with more than %d READs.  Code fix needed.\n",
                     0x7FFFFFFF);
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
      }
      if (duplicates) {
         rc = fclose (duplicates);
         duplicates = NULL;
         if (rc == EOF) {
            fprintf (stderr, "maketrie: error closing %s-duplicates - %s\n",
                     argv[1], strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
      }

      if (read_file) {
         rc = fclose (read_file);
         read_file = NULL;
         if (rc == EOF) {
            fprintf (stderr, "maketrie: error closing %s - %s\n", argv[1],
                     strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
      }

      if (read_index) {
         rc = fclose (read_index);
         read_index = NULL;
         if (rc == EOF) {
            fprintf (stderr, "maketrie: error closing %s-index - %s\n",
                     argv[1], strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
      }

      fprintf (stderr,
               "\nread trie built using %lld nodes (%0.0f%% of capacity)\n",
               last_used_edge, 100.0 * last_used_edge / MAX_SIZE);
      fprintf (stderr,
               "\nTotal of %d reads indexed and sorted, including %d (%0.0f%%) duplicates"
               " (dup count is temporarily inaccurate when using multiple nodes)\n",
               seq, dups, dups * 100.0 / seq);
      fprintf (stderr, "\nFrequencies:\n");
      for (c = 0; c < 256; c++) if (freq[c]) fprintf (stderr, "   %c  %ld\n", c, freq[c]);
      for (i = 0; i < MAX_LINE; i++) {
         if (length[i]) {
            number_of_lengths++;
            read_length = i;
         }
      }
      fprintf (stderr, "\n");
      if (number_of_lengths == 0) {
         fprintf (stderr, "Error: No READs found!  Bad input file?\n");
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      } else if (number_of_lengths != 1) {
         fprintf (stderr,
                  "Error: this code does not handle READs of differing lengths\n");
         fprintf (stderr, "\nWe found READs of lengths:\n");
         for (i = 0; i < MAX_LINE; i++) {
            if (length[i]) {
               fprintf (stderr, "     %d  (%ld)", i, length[i]);
            }
         }
         fprintf (stderr, "\n");
         fprintf (stderr,
                  "\nPlease clean the data first with a program like 'fastqc'.\n\n");
         shut_down_other_nodes ();
         MPI_Finalize ();
         exit (EXIT_FAILURE);
      } else {
         fprintf (stderr, "READ length: %d\n", read_length);
      }
      fprintf (stderr, "\n");
      fflush (stderr);

      walk_and_print_trie ();
      sprintf (fname, "%s-edges", argv[1]);
      dump_trie (fname);

      if (rejects) {
         rc = fclose (rejects);
         rejects = NULL;
         if (rc == EOF) {
            fprintf (stderr,
                     "maketrie: cannot close reject file \"%s\" - %s\n",
                     argv[1], strerror (errno));
            shut_down_other_nodes ();
            MPI_Finalize ();
            exit (EXIT_FAILURE);
         }
      }

      time (&curtime);
      fprintf (stderr, "Program complete at %s", ctime (&curtime));
      shut_down_other_nodes ();
      free (trie_cell);
      trie_cell = NULL;

   } else {
      INDEX index;
      CELL read_value;
      long longvalue;
      MPI_Status status;

      for (i = 0; i < 5; i++) empty.edge[i] = 0LL;
      last_used_edge = mpirank * CHUNKSIZE - 1;

      for (;;) {
         int caller;

         MPI_Recv (&longvalue, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &status);
         caller = status.MPI_SOURCE;

         if (status.MPI_TAG == TAG_READ_READ) {
            MPI_Recv (&index, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                      MPI_COMM_WORLD, &status);
            getread (index, &read_value);
            MPI_Send (&read_value.edge,
                      sizeof (read_value.edge) / sizeof (read_value.edge[0]),
                      MPI_LONG, caller, 0, MPI_COMM_WORLD);

         } else if (status.MPI_TAG == TAG_WRITE_READ) {
            MPI_Recv (&index, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                      MPI_COMM_WORLD, &status);
            assert (sizeof (read_value.edge[0]) == sizeof (long));
            MPI_Recv (&read_value.edge,
                      sizeof (read_value.edge) / sizeof (read_value.edge[0]),
                      MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                      &status);
            setread (index, read_value);
            MPI_Send (&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);

         } else if (status.MPI_TAG == TAG_ADD_READ) {
            accept_add_read (mpirank, longvalue, status);

         } else if (status.MPI_TAG == TAG_GET_NEXT_FREE_EDGE) {
            accept_get_next_free_edge (caller);
            MPI_Send (&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);

         } else if (status.MPI_TAG == TAG_OUTPUT_READ) {
            accept_output_read (mpirank, longvalue, status);

         } else if (status.MPI_TAG == TAG_WALK_AND_PRINT_TRIE_INTERNAL) {
            accept_walk_and_print_trie_internal (mpirank, longvalue, status);

         } else if (status.MPI_TAG == TAG_DUMP_TRIE) {
            accept_dump_trie (mpirank, longvalue, status);

         } else if (status.MPI_TAG == TAG_EXIT_PROGRAM) {
            fprintf (stderr, "Node %d asked to exit\n", mpirank);
            MPI_Send (&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);
            fprintf (stderr, "Node %d exit acknowleged\n", mpirank);
            break;

         } else {

         }
      }
      fprintf (stderr,
               "Node %d exiting cleanly.  local base = %lld,"
               "  last_used_edge = %lld,  local maximum = %lld\n",
               mpirank, mpirank * CHUNKSIZE, last_used_edge,
               (mpirank + 1) * CHUNKSIZE);
      if (duplicates) {
         rc = fclose (duplicates);
         duplicates = NULL;
         if (rc == EOF) {
            fprintf (stderr, "maketrie: error closing %s-duplicates - %s\n",
                     argv[1], strerror (errno));
         }
      }

      if (trie_cell != NULL) free (trie_cell);
      trie_cell = NULL;

   }

   MPI_Finalize ();

   exit (EXIT_SUCCESS);
   return EXIT_FAILURE;
}
