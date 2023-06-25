#ifdef AMOS_OVERLAPS
// User-tweakable parameters:
#define MIN_OVERLAP 14
                           // Overlaps of < MIN_OVERLAP letters are not interesting.
                           // 14 was by observation from a sample of approximately 15M reads
                           // Smaller samples may need a smaller minimum overlap & vice-versa.
                           // Later this cutoff will be determined dynamically.
#define MAX_OVERLAPS 8
                           // Don't print more than MAX_OVERLAPS items for any one
                           // index position on a single read.  Arbitrary.
#else
#define MIN_OVERLAP 1      // Not relevant when outputting trie nodes instead of list of reads
#define MAX_OVERLAPS 999999
#endif

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>
#include <omp.h>
#include <ctype.h>
#include <time.h>  // for info only

#define _XOPEN_SOURCE 500
#include <unistd.h>
ssize_t pread(int fd, void *buf, size_t count, off_t offset);

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#ifndef FALSE
#define TRUE (0==0)
#define FALSE (!TRUE)
#endif

static FILE *overlaps = NULL;
static FILE *read_file_sorted = NULL;
static int memory_mapped = FALSE;
static int trie_file_fd = -1;

//#define GBs 4L

typedef unsigned long long EDGE;
typedef unsigned long long INDEX;
#define ENDS_WORD (1ULL<<63ULL)
#define EDGE_MASK (ENDS_WORD-1UL)

typedef struct cell {
  EDGE edge[5];
} CELL;
CELL *trie_cell;

#define ROOT_CELL ((INDEX)1L)

static INDEX last_used_edge = ROOT_CELL;

#define _A_ 0
#define _C_ 1
#define _G_ 2
#define _T_ 3
#define _N_ 4

char *trt = "ACGTN"; // translate table back from the above

#define MAX_LINE 1024

static int read_length = 0; // This is the length we found in the sorted-read file.


#define CORES_PER_NODE 16ULL

static int mpirank = 0, mpisize = 1, cluster_base = 0, cluster_size = 1;

static long long TASKS_PER_NODE, PROCESSORS_PER_NODE;

#define TAG_DATA 1
#define TAG_ACK 2

#define TAG_SEND_RAW_MEM 3
#define TAG_EXIT_PROGRAM 6

#define TAG_LOCATE_OVERLAPS 11
#define TAG_PRINT_OVERLAPS 12

ssize_t retrying_pread(int filedes, void *buffer, size_t size, off_t offset)
{
  ssize_t requested = size;
  ssize_t rc;
  for (;;) {
    //fprintf(stderr, "Node %d: Load 0x%lx bytes to offset 0x%lx at %p\n", mpirank, size, offset, buffer);
    rc = pread(filedes, buffer, size, offset);
    if (rc <= 0) return rc;
    buffer = (char *)buffer + rc; size -= rc; offset += rc;
    if (size == 0L) return requested;
    //fprintf(stderr, "Node %d: Got 0x%lx bytes - fetching 0x%lx more...\n", mpirank,rc, size);
  }
}

static long long CHUNKBITS, CHUNKSIZE, CHUNKMASK;

static void shut_down_other_nodes(void) // within this cluster-group
{
  int target_rank;
  long int value = 0L;
  //MPI_Status status;

  for (target_rank = cluster_base+1; target_rank < cluster_base+cluster_size; target_rank++) {
    if (target_rank != mpirank) { // normally called only by cluster_base except if error
      MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_EXIT_PROGRAM, // The command.
	   MPI_COMM_WORLD);/* always use this */
#ifdef NEVER
      // Get acknowlegement that it has been written (to synchronise)
      MPI_Recv(&value,/* message buffer - not used */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   /*rank,*/     MPI_ANY_SOURCE,/* receive from any sender */
	   /* TAG_ACK,*/ MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
#endif
    }
  }
}

static void send_bytes(int dest, void *memp, long bytes) { // just a type-checked macro
  MPI_Send(memp, bytes, MPI_BYTE, dest, TAG_SEND_RAW_MEM, MPI_COMM_WORLD);
}

static void locate_overlaps(char *s, EDGE edge, long read_number, int matching_offset);
static void print_overlaps(EDGE edge, long read_number, int matching_offset, int *number_printed);
static void local_locate_overlaps(char *s, EDGE edge, long read_number, int matching_offset)
{
  int c;

  for (;;) {
    //fprintf(stderr, "locate_overlaps(\"%s\", %llx, %ld, %d)\n", s, edge, read_number, matching_offset);
    c = *s++;

    if (c == 'A') c = _A_; // architecture-dependent whether this or a lookup table is faster..
    else if (c == 'C') c = _C_;
    else if (c == 'G') c = _G_;
    else if (c == 'T') c = _T_;
    else c = _N_; // some other char

    edge = trie_cell[edge&CHUNKMASK].edge[c] & EDGE_MASK;
    if (edge == 0LL) return; // no matches down this path

    // edge may now be stored on a different processor so do *NOT* access trie_cell[edge&CHUNKMASK]
    //  except with a 'safe' procedure
  
    if ((*s == '\0') || (*s == '\n') || (*s == '\r')) {
      int print_count = 0;
      // this character matched and was the last letter in the string.
      print_overlaps(edge, read_number, matching_offset, &print_count); // edge points to read_number
      return;
    }

    //fprintf(stderr, "recurse: locate_overlaps(\"%s\", %llx, %ld, %d)\n", s, edge, read_number,
    //        matching_offset);
    if ((edge >> CHUNKBITS) != (mpirank%cluster_size)) {
      locate_overlaps(/* modified */ s, edge, read_number, matching_offset);
      return;
    } // else optimise tail recursion by going round the loop again.
  }
}

#ifdef AMOS_OVERLAPS
static void local_print_overlaps(EDGE edge, long read_number, int matching_offset, int *number_printed)
{
  int i;

  //fprintf(stderr, "Node %d: print_overlaps(%lld, %ld, %d, %d)\n", mpirank, edge, read_number,
  //        matching_offset, *number_printed);

  for (i = 0; i < 5; i++) {
    if ((*number_printed) >= MAX_OVERLAPS) return; // Enough!
    if (trie_cell[edge & CHUNKMASK].edge[i]&ENDS_WORD) {
      // Do we want to include self-overlaps?
      // i.e. where trie_cell[edge & CHUNKMASK].edge[i]&EDGE_MASK == read_number ???
      fprintf(overlaps, "{OVL\nadj:N\nrds:%ld,%lld\nscr:0\nahg:%d\nbhg:%d\n}\n",
              1+read_number, // Hopefully I have these two in the right order now...
              1+(trie_cell[edge & CHUNKMASK].edge[i]&EDGE_MASK),
              matching_offset, matching_offset);
              // NOTE: I've been numbering reads from 0 up like a computer scientist should;
              // apparently bioinformaticists count from 1 up like engineers do :-(  Hence "+1" above.
      *number_printed = (1 + (*number_printed));
      // should we warn if we're truncating excessive overlaps?
      if (ferror(overlaps)) {
        fprintf(stderr, "\n\n************* print_overlaps() failed, %s\n", strerror(errno));
        shut_down_other_nodes();
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    } else if (trie_cell[edge & CHUNKMASK].edge[i]) {
      // not final letter, and this letter is present with more to follow
      // recurse, *safely*, to locate all leaf nodes
      if ((*number_printed) < MAX_OVERLAPS) {
        print_overlaps(trie_cell[edge & CHUNKMASK].edge[i], read_number, matching_offset,
                       number_printed);
      }
    }
  }
}
#endif

static int remote_locate_overlaps(long target_rank, char *s, long edge,
                                  long read_number, int matching_offset)
{  // pass to another node
  MPI_Status status;
  long value = 0L; // generic up-front parameter
  long stringlength;

//#pragma omp critical
  {
  stringlength = strlen(s)+1;

  // COMMAND CODE
  //fprintf(stderr, "sending TAG_LOCATE_OVERLAPS\n");
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_LOCATE_OVERLAPS,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // SEND PARAMETERS - later optimise by sending one param via 'value' above
  //fprintf(stderr, "sending stringlength=%ld\n", stringlength);
  MPI_Send(&stringlength,/* message buffer - address part */
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  //fprintf(stderr, "Sending string %s\n", s);
  send_bytes(target_rank, s, stringlength);

  //fprintf(stderr, "sending edge %llx\n", edge);
  MPI_Send(&edge,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&read_number,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&matching_offset,
	   1,/* one data item */
	   MPI_INT,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // Get acknowlegement that it has been written (to synchronise)
  MPI_Recv(&value,/* message buffer - used for len function result */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
}
  // If len were not needed we could fire & forget, by removing the Recv above...
  return (int)value;

}

#ifdef AMOS_OVERLAPS
static int remote_print_overlaps(long target_rank, long edge, long read_number,
                                 int matching_offset, int *number_printed)
{ // pass to another node
  MPI_Status status;
  long value = 0L; // generic up-front parameter

//#pragma omp critical
  {
  // COMMAND CODE
  value = (long) *number_printed;
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_PRINT_OVERLAPS,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // SEND PARAMETERS - later optimise by sending one param via 'value' above

  MPI_Send(&edge,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&read_number,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&matching_offset,
	   1,/* one data item */
	   MPI_INT,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // Get acknowlegement that it has been written (to synchronise)
  MPI_Recv(&value,/* message buffer - used for len function result */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
  *number_printed = (int)value;
  }
  return (int)value;

}
#endif

static void accept_locate_overlaps(int myrank, long value, MPI_Status status)
{ // receive request from RPC mechanism
  char s[MAX_LINE];
  EDGE edge;
  long read_number, stringlength;
  int matching_offset;
  int caller;

  // RECEIVE PARAMETERS (possible optional first parameter passed in as 'value')
  //fprintf(stderr, "Receiving stringlen\n");
  MPI_Recv(&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  //fprintf(stderr, "Received stringlen=%ld\n", stringlength);
  caller = status.MPI_SOURCE;

  //fprintf(stderr, "Receiving string\n");
  MPI_Recv(s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  //fprintf(stderr, "Received string=%s\n", s);

  MPI_Recv(&edge, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  //fprintf(stderr, "received  edge %llx\n", edge);

  MPI_Recv(&read_number, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  MPI_Recv(&matching_offset, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  locate_overlaps(s, edge, read_number, matching_offset);

  MPI_Send(&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege and return result

}

static void accept_print_overlaps(int myrank, long value, MPI_Status status)
{ // receive request from RPC mechanism
  EDGE edge;
  long read_number;
  int matching_offset;
  int number_printed = 0;
  int caller;

  // RECEIVE PARAMETERS (possible optional first parameter passed in as 'value')
  caller = status.MPI_SOURCE;
  number_printed = (int)value;

  MPI_Recv(&edge, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  MPI_Recv(&read_number, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  MPI_Recv(&matching_offset, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  print_overlaps(edge, read_number, matching_offset, &number_printed);

  // Acknowlege and return result
  value=(long)number_printed; MPI_Send(&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD);

}

static void locate_overlaps(char *s, EDGE edge, long read_number, int matching_offset)
{
  long long int target_rank = (long long)edge >> CHUNKBITS;

  //fprintf(stderr, "locate_overlaps(\"%s\", %llx)  target_rank=%lld  mpirank=%d\n", s,
  //        edge, target_rank, mpirank);
  if (target_rank == (mpirank%cluster_size)) {
    local_locate_overlaps(s, edge, read_number, matching_offset);
  } else {
    remote_locate_overlaps((long)(target_rank+cluster_base), s, edge, read_number, matching_offset);
  }

  return;
}

static void print_overlaps(EDGE edge, long read_number, int matching_offset,
                           /* COPY-IN/COPY-OUT: */ int *number_printed)
{

  // This is called by 'locate_overlaps'.  It can either print out our internal format which
  // is a read_number, offset, and node_id, or it can walk the trie starting with node_id
  // and output actual overlaps, eg AMOS "OVL" records.

#ifdef AMOS_OVERLAPS
  long long int target_rank = (long long)edge >> CHUNKBITS;

  // A tree-walk is necessary to find the leaves
  if (target_rank == (mpirank%cluster_size)) {
    local_print_overlaps(edge, read_number, matching_offset, number_printed);
  } else {
    remote_print_overlaps(target_rank+cluster_base, edge, read_number, matching_offset, number_printed);
  }
#else
  // no walk needed - just print the common node - this is actually more useful,
  // but is not a paradigm used by current assemblers...
  // There is no need to heed "MIN_OVERLAP" or "MAX_OVERLAPS" when all we're printing
  // is one node for all overlaps of a certain length.  Those tweaks are only useful
  // when we walk the trie at this node and generate a large list of actual overlaps.
  fprintf(overlaps, "%ld:%d @%lld\n", read_number, matching_offset, edge);
#endif
  return;
}

int main(int argc, char **argv)
{
  /*  (local declarations) */
  off_t segment_size;
  time_t curtime;
  char fname[1024];
  //char line[MAX_LINE];
  int rc; // i, c
  long read_number = 0;
  int namelen;
  int required=MPI_THREAD_SERIALIZED; // Required level of MPI threading support
  /* Each thread will call MPI routines, but these calls will be coordinated
     to occur only one at a time within a process.
  */
  int provided;                       // Provided level of MPI threading support
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  

  // Initialize MPI with threading
  MPI_Init_thread(&argc, &argv, required, &provided);

  // Determine the MPI rank, number of processes, and processor name
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  //fprintf(stderr, "I am rank %d of world size %d\n", mpirank, mpisize);
  MPI_Get_processor_name(processor_name, &namelen);
  if (processor_name && strchr(processor_name, '.')) *strchr(processor_name, '.') = '\0';

  time(&curtime); if (mpirank == 0) fprintf(stderr, "Program started at %s", ctime(&curtime));

  if (provided < required) {

    // Insufficient support, degrade to 1 thread and warn the user

    if (mpirank == 0) {
      printf("Warning:  This MPI implementation provides insufficient"
	     " threading support.\n");
    }
    omp_set_num_threads(1);
  }

  /*  Command-line parameter handling and file opening */
  if ((mpirank == 0) && (argc > 2)) {
    fprintf(stderr, "warning: extra parameter %s ignored...\n", argv[2]);
  }

  if (argc >= 2) { // command-line handling could be cleaner.

    /*  Open all files. */

#ifdef AMOS_OVERLAPS
    sprintf(fname, "%s-ovl-%05d.afg", argv[1], mpirank);
#else
    sprintf(fname, "%s-%05d.ovl", argv[1], mpirank);
#endif
    overlaps = fopen(fname, "w");
    if (overlaps == NULL) {
      fprintf(stderr, "findoverlaps[%d]: cannot create overlap output \"%s\" - %s\n",
              mpirank, fname, strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    //fprintf(stderr, "Overlap Output[%d]: %s\n", mpirank, fname);

  } else {
    if (mpirank == 0) fprintf(stderr, "syntax: findoverlaps input.fastq\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  /* OS-dependent (Linux) code below lets us grab as much memory as possible in a single chunk. */

  // DO NOT throw this away and attempt to build the trie using individual cells claimed from malloc.
  // Sequential memory allocation is *critical* to the speed of the program.

  CHUNKBITS = (((INDEX)sizeof(INDEX)) * 8ULL) - 1ULL;
  {
    FILE *meminfo, *cpuinfo;
    long long memsize;
    static char line[1024 /*MAX_LINE*/];
    PROCESSORS_PER_NODE = 0ULL;
    cpuinfo = fopen("/proc/cpuinfo", "r");
    if (cpuinfo) {
      for (;;) {
        char *s = fgets(line, 1023 /*MAX_LINE*/, cpuinfo);
        if (s == NULL) break;
        if (strncmp(line, "processor", strlen("processor")) == 0) {
          PROCESSORS_PER_NODE++;
        }
      }
      if (mpirank == 0) fprintf(stderr, "Discovered %lld processors per node\n", PROCESSORS_PER_NODE);
    }
    if (PROCESSORS_PER_NODE == 0ULL) PROCESSORS_PER_NODE = CORES_PER_NODE;

    //    TASKS_PER_NODE = PROCESSORS_PER_NODE/(long long)omp_get_max_threads();
    //    fprintf(stderr, "Node %s, Rank %d, and running %lld ranks (%lld/%lld) on this node.\n",
    //	    processor_name, mpirank, TASKS_PER_NODE, (long long)PROCESSORS_PER_NODE, (long long)omp_get_max_threads());

    TASKS_PER_NODE = 1;

    meminfo = fopen("/proc/meminfo", "r");
    if (meminfo) for (;;) {
      int count;
      char *s = fgets(line, 1023 /*MAX_LINE*/, meminfo);
      if (s == NULL) break;
      //fprintf(stderr, "checking %s", line);
      count = sscanf(line, "MemTotal:     %lld kB", &memsize);
      if (count == 1) {
        fprintf(stderr, "got memsize %lld kB\n", memsize);
        memsize *= 1024ULL; // re-scale to bytes
        memsize /= TASKS_PER_NODE; // availability per processor
        memsize /= (long long)sizeof(CELL); // convert to cell count

        CHUNKBITS = 1ULL;
        for (;;) {
          if ((1ULL << CHUNKBITS) >= memsize) break;
          CHUNKBITS++;
	}
        CHUNKBITS -= 1ULL;
        if (mpirank == 0) fprintf(stderr,
        "Rounding down memsize to %lldM cells per core (%lld bits), ie %lldM cells per node\n",
                (1ULL<<CHUNKBITS)>>24ULL,
                CHUNKBITS,
                ((1ULL<<CHUNKBITS)*TASKS_PER_NODE)>>24ULL);
	break;
      } else {
        CHUNKBITS = CHUNKBITS>>1ULL;
      }
    }
  }
  CHUNKSIZE = (1ULL<<CHUNKBITS); /* index to local array is 0..CHUNKSIZE */
  CHUNKMASK = (CHUNKSIZE-1ULL);

  while (CHUNKBITS >= 16ULL) { // absolute minimum, no point in going below this!
    fprintf(stderr, "Node %d: trying calloc of %lld cells of %d bytes each.\n",
            mpirank, CHUNKSIZE, (int)sizeof(CELL));
    trie_cell = calloc(CHUNKSIZE, sizeof(CELL)); // trie_cell[0..CHUNKSIZE]
    if (trie_cell == NULL) {
      CHUNKBITS -= 1ULL;
      CHUNKSIZE = (1ULL<<CHUNKBITS);
      CHUNKMASK = (CHUNKSIZE-1ULL);
    } else break; // malloc successful
  }

  if (trie_cell == NULL) {
    fprintf(stderr,
            "findoverlaps: rank %d unable to allocate array of %lld longs\n", mpirank, CHUNKSIZE);
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    fprintf(stderr, "Node %d: allocated %lld-item long array\n", mpirank, CHUNKSIZE);
  }

  fprintf(stderr,
          "node %d: using %dM-items.  Launching listener now.\n",
          mpirank,
          (int)(CHUNKSIZE >> 24ULL)
         );

  {
    off_t file_length;

    sprintf(fname, "%s-edges", argv[1]);
    trie_file_fd = open(fname, O_RDONLY);
    if (trie_file_fd < 0) {
      fprintf(stderr, "findoverlaps: cannot access trie file %s - %s\n", fname, strerror(errno));
      exit(EXIT_FAILURE);
    }
    file_length = lseek(trie_file_fd, (off_t)0LL, SEEK_END);
    last_used_edge = file_length/sizeof(CELL)-1LL;

    cluster_size = (last_used_edge + (CHUNKSIZE-1)) / CHUNKSIZE;

    if (mpisize < cluster_size) {
      if (mpirank == 0) fprintf(stderr, "ERROR: We need %d or more compute nodes to be allocated - we only have %d\n",
                                cluster_size, mpisize);
      MPI_Finalize(); // release this processor back to the pool.  I have no use for it.
      exit(0);
    }

    if ((mpisize / cluster_size) > 1) {
      if (mpirank == 0) fprintf(stderr, "We expect a parallel speedup by a factor of %d\n",
              mpisize / cluster_size);
    }

    if (mpirank >= ((mpisize / cluster_size) * cluster_size)) {
      fprintf(stderr,
              "*** WARNING: Node %d is not needed (last required node is %d) - releasing it...\n",
              mpirank, (mpisize / cluster_size) * cluster_size - 1);
      if (overlaps) fclose(overlaps); overlaps = NULL; // Would prefer never to have opened it...
      MPI_Finalize(); // release this processor back to the pool.  I have no use for it.
      exit(0);
    }

    cluster_base = (mpirank / cluster_size) * cluster_size;

    fprintf(stderr, "Node %d: cluster is %d..%d\n", mpirank, cluster_base, cluster_base+cluster_size-1);

    if ((mpirank%cluster_size) < cluster_size-1) {
      // LOAD 'CHUNKSIZE'
      segment_size = CHUNKSIZE * sizeof(CELL);
      fprintf(stderr, "Node %d: Loading full sized chunk. (%lld)\n", mpirank, (long long int)segment_size);
    } else if ((mpirank%cluster_size) == cluster_size-1) {
      // LOAD the remaining data <= CHUNKSIZE at the end of the array
      segment_size = ((last_used_edge & CHUNKMASK)+1) * sizeof(CELL);
      fprintf(stderr, "Node %d: Loading remainder of last chunk. (%lld)\n", mpirank, (long long int)segment_size);
    } else {
      // (mpirank%cluster_size) >= cluster_size - can't happen
      //assert((mpirank%cluster_size) < cluster_size);
      exit(1);
    }

    //fprintf(stderr, "Node %d: mapping %ld bytes of %s at offset %ld\n",
    //        mpirank, segment_size, argv[1], (mpirank%cluster_size)*CHUNKSIZE*sizeof(CELL));


    // MEMORY MAPPING HAS BEEN REMOVED DUE TO SLOW SPEED:
    if (trie_cell != NULL) free(trie_cell);
    trie_cell = NULL;//mmap(NULL, (size_t)segment_size, PROT_READ, MAP_PRIVATE/*SHARED*/,
                     //     trie_file_fd, (off_t)(mpirank%cluster_size)*(off_t)CHUNKSIZE*sizeof(CELL));


    if ((trie_cell == NULL) || (trie_cell == (void *)-1)) {
      ssize_t rc;
      //fprintf(stderr, "findoverlaps[%d]: failed to map %s - %s - loading into RAM instead\n",
      //         mpirank, fname, strerror(errno));
      trie_cell = malloc(segment_size);
      rc = retrying_pread(trie_file_fd, trie_cell, (size_t)segment_size,
                          (off_t)(mpirank%cluster_size)*(off_t)CHUNKSIZE*sizeof(CELL));
      if (rc != segment_size) {
	fprintf(stderr,
                "findoverlaps[%d]: failed to fetch %ld bytes from offset 0x%llx on file %d, rc = %ld\n",
		mpirank,segment_size, (mpirank%cluster_size)*CHUNKSIZE, trie_file_fd, rc);
	exit(1);
      }
    } else {
      memory_mapped = TRUE;
      // We have enough physical memory to back this VM, but by mapping rather than loading, it should
      // start up more quickly.  (Turns out it slows down overall... removed...)
      //fprintf(stderr, "findoverlaps[%d]: loaded %ld bytes at %p\n", mpirank, segment_size, trie_cell);
    }
    //fprintf(stderr, "Node %d: last_used_edge: %lld,  CHUNKSIZE: %lld\n",
    //        mpirank, last_used_edge, CHUNKSIZE);
  }
  

  if ((mpirank%cluster_size) == 0) {
    // --------------------------- MAIN BODY OF CODE ON PRIMARY PROCESSORS ----------------------------
    if (mpirank == 0) fprintf(stderr,
            "\nCombined system is using %lldM trie edges distributed across %d nodes.\n\n",
            (long long)mpisize * (CHUNKSIZE >> 24ULL), mpisize
           );
    if (mpirank == 0) fprintf(stderr,
            "This is enough for %d copies of the database...\n\n",
            mpisize / cluster_size
           );


    // having built our key data structure, now we find the overlaps

    /*  PHASE TWO: DO AN ALL-TO-ALL COMPARISON BETWEEN THE READS */

    if (mpirank == 0) if (mpirank == 0) fprintf(stderr,
            "\nLocating overlaps by comparing each unique READ"
            " against the trie of all unique reads...\n\n");

    // IF WE ARE RUNNING MORE THAN ONE GROUP OF NODES, 'BASE' PROCESSOR OF EACH GROUP NEEDS
    // TO READ DIFFERENT SECTIONS OF THIS FILE.  EITHER INTERLEAVING LINES, OR BY LARGER BLOCKS.

    sprintf(fname, "%s-sorted", argv[1]);
    read_number = 0;
    read_file_sorted = fopen(fname, "r");
    if (read_file_sorted == NULL) {
      fprintf(stderr, "findoverlaps[%d]: cannot reopen input \"%s\" - %s\n",
              mpirank, argv[1], strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    // This loop is executed on the base node of each group.  This doesn't give perfectly balanced
    // processing, but it is fairly close.  Rather than split the data up into chunks, we simply
    // have each of the <N> processing groups skip all but 1/N of the lines from the input file.

    time(&curtime); if (mpirank == 0) fprintf(stderr, "\nStarting comparisons at %s", ctime(&curtime));

    for (;;) {  // THIS IS THE "EMBARASSINGLY PARALLEL" MAIN LOOP, RUNNING ON MULTIPLE NODES.
                // If we needed to, we could make this a null-process in a polling loop.
                // But Master/Slave is working well enough for now...
      char line[MAX_LINE];
      int len, c;
      char *s;
      INDEX original_read_number;

      s = line;
      for (;;) {
        c = fgetc(read_file_sorted);
        if ((c == EOF) || ferror(read_file_sorted)) break;
        *s++ = c;
        if (c == '\n') break;
      }
      if ((c == EOF) || ferror(read_file_sorted)) break;

      *--s = '\0';

      if (read_length == 0) {
        s = line;
        while (*s++ != ' ') read_length += 1;
      }

      read_number++;

      if (strlen(line) != read_length+13) continue;

      s = line+read_length+1;
      while (*s == ' ') s++;
      if (!isdigit(*s)) continue;

      original_read_number = atoll(s);

      line[read_length] = '\0';
      s = line; // The loop below skips the first letter of the target -
                //  we already handled dups elsewhere

      // Adding the pragma below brings the time taken to process fake100m
      // down from about 1hr to about 40m.  Not a huge win but not bad for 
      // a single line of meta-code.  The results are identical except for
      // being in a slightly different order.  The parallelism is only happening
      // on the master rank - we're not accepting multiple requests yet
      // on the slave ranks.

      // HOWEVER!!!! when I added parallel MPI support (multiple compute groups),
      // this broke.  Had to comment it out...  have not yet retested this
      // subsequent to getting the MPI parallelism to work. (which may be enough)

      // just process every <n>th line when we have <n> compute groups.
      if ((read_number%(mpisize/cluster_size)) == (cluster_base/cluster_size)) {

//#pragma omp parallel for
        for (len = read_length-1; len >= MIN_OVERLAP; len--) {
          // (MIN_OVERLAP may not be needed if AMOS_OVERLAPS is not defined.)

    	  locate_overlaps(s-len+read_length, ROOT_CELL, original_read_number,
            /* OFFSET 1..37 */ read_length-len);   // The meat in the sandwich...

        }

      }

      if ((read_number % 1000000) == 0) {
        time(&curtime);
        if (mpirank == 0) fprintf(stderr, "%ld READs processed for overlaps at %s", read_number,
                ctime(&curtime));
      }
    }

    if (read_file_sorted) {
      rc = fclose(/* input */ read_file_sorted); read_file_sorted = NULL;
      if (rc == EOF) {
        fprintf(stderr, "findoverlaps[%d]: cannot close input \"%s\" - %s\n", 
                mpirank, argv[1], strerror(errno));
        shut_down_other_nodes();
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }

    time(&curtime); fprintf(stderr, "Program group %d of %d complete at %s",
                             mpirank/cluster_size, mpisize/cluster_size, ctime(&curtime));
    shut_down_other_nodes();
    fflush(overlaps);

  } else {
    /*  ALL OTHER PROCESSORS RUN AN ACCEPT/DISPATCH LOOP FOR REMOTE PROCEDURE CALLS... */
    // (currently single-threaded.  An MPI guru could write this better.)

    // Also it is VERY IMPORTANT to note that Node 0 does *not* run a copy of
    // this dispatcher.  (because it is single-threaded too.)  Would be nice if
    // that option were available, eg to allow callbacks.

    // One optimisation trick planned for the far future...  if A calls B and B calls C,
    // then C could pass the result back directly to A, allowing B to return to its
    // displatch loop and accept more commands.  Which will be fine if it's always A
    // that initiates certain types of call... - just need to be careful with
    // specifying 'MPI_ANY_SOURCE' rather than 'caller'.

    //INDEX index;
    //CELL read_value;
    long longvalue;
    MPI_Status status;

    //time(&curtime); fprintf(stderr, "\nSlave %d starting at %s", mpirank, ctime(&curtime));
    // ultimately would like to run 16 listeners in parallel...
    for (;;) {
      int caller, done = FALSE;

      // Accept an RPC request.  Get the command code.
      MPI_Recv(&longvalue, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      caller = status.MPI_SOURCE;

      /*
       * Check the tag of the received message and dispatch accordingly
       */
      if (status.MPI_TAG == TAG_LOCATE_OVERLAPS) {  // CALL locate_overlaps()
	accept_locate_overlaps(mpirank, longvalue, status);

      } else if (status.MPI_TAG == TAG_PRINT_OVERLAPS) {  // CALL print_overlaps()
	accept_print_overlaps(mpirank, longvalue, status);

      } else if (status.MPI_TAG == TAG_EXIT_PROGRAM) {  // Kill clients
        //fprintf(stderr, "Node %d asked to exit\n", mpirank);
        //MPI_Send(&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege
        //fprintf(stderr, "Node %d exit acknowleged\n", mpirank);
        done = TRUE;

      } else {
        // UNKNOWN - CODING ERROR?
      }

      if (done) break;
    }
    //fprintf(stderr, "Node %d exiting cleanly.  local base = %lld\n",
    //	    mpirank, (mpirank%cluster_size) * CHUNKSIZE);

    if (overlaps) {
      rc = fclose(/* output */overlaps); overlaps = NULL;
      if (rc == EOF) {
        fprintf(stderr, "findoverlaps[%d]: Error closing %s-ovl-%05d.afg - %s\n",
                mpirank, argv[1], mpirank, strerror(errno));
      }
    }

  }

  if (memory_mapped) {
    fprintf(stderr, "findoverlaps[%d]: unmapping %s-edges (%ld bytes)\n",
             mpirank, argv[1], segment_size);
    rc = munmap(trie_cell, segment_size);
    fprintf(stderr, "findoverlaps[%d]: Error unmapping %s-edges - %s\n",
            mpirank, argv[1], strerror(errno));
  } else if (trie_cell != NULL) free(trie_cell);
  trie_cell = NULL;

  MPI_Finalize();

  exit(EXIT_SUCCESS);
  return EXIT_FAILURE;
}
