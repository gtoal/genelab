/*
    nearmatch ~gtoal/genelab/data/40kreads-schliesky.fastq AAACCAGCAGATCCAGCACCAACGACGACGACATCAGTCTCAGCATAAGTGATCATATCCGTCATGTACCTTCTCGTCATCTCACGGGACACGATCGATTC
 */

// Currently implemented with MMAP, unfortunately it only works for small files... when the index is 24Gb,
// mmap reports that it can't allocate memory.  What it really means is that there is no slot large enough
// in the virtual address space where the file could be connected.

// Only workaround is load as much as possible in RAM, and use a function to get an entry for elements
// outside the array...  To do... (be careful to perform minimal disk I/O)

#define _FILE_OFFSET_BITS 64

typedef unsigned long long EDGE;
typedef unsigned long long INDEX;
#define ENDS_WORD (1ULL<<63ULL)
#define EDGE_MASK (ENDS_WORD-1ULL)

#define MAX_LINE 1024
#define MIN_OVERLAP 13
#define ALLOWED_ERRORS 3 /* Temp */

//static INDEX MAX_SIZE = ((INDEX)0L);

#define ROOT_CELL ((INDEX)1L)
// Node 0 is unused, 0 is needed as a terminator.

// originally I used 'next_free_edge' and it was exclusive.  Have now changed that
// to 'last_used_edge' which is inclusive.  getting a free edge is now a function
// call that we can make atomic with an OMP PRAGMA.


typedef struct node {
  // ACGT maps to 0 1 2 3, anything else such as 'N' is 4
  EDGE edge[5];
} CELL;
CELL *trie_cell;

static INDEX last_used_edge; // inclusive
static INDEX contig_number = 1;

#define _A_ 0
#define _C_ 1
#define _G_ 2
#define _T_ 3
#define _N_ 4

char *trt = "ACGTN"; // translate table back from the above

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>

#define _XOPEN_SOURCE 500
#include <unistd.h>
ssize_t pread(int fd, void *buf, size_t count, off_t offset);

#ifndef FALSE
#define TRUE (0==0)
#define FALSE (0!=0)
#endif

static int freq[256];
static int printed = FALSE;
static long long *read_sequence_no_to_file_offset, contig_size;
static FILE *read_file, *gene_file, *afg_contig_file, *afg_tle_file, *afg_reads_file;
static int trie_fd = -1, index_fd = -1;
static off_t file_length;
static char gene_file_name[MAX_LINE];
static char trie_file_name[MAX_LINE];
static char index_file_name[MAX_LINE];

static char *stringat(long long textp)
{
  char line[MAX_LINE], qual[MAX_LINE];
  char *s;
  (void)fseek(read_file, textp, SEEK_SET);
  s = fgets(line, MAX_LINE, read_file);
  fgets(line, MAX_LINE, read_file);
  s = strchr(line, '\n'); if (s) *s = '\0'; // trim trailing newline
  fgets(qual, MAX_LINE, read_file);
  s = fgets(qual, MAX_LINE, read_file);
  s = strchr(qual, '\n'); if (s) *s = '\0'; // trim trailing newline
  strcat(line, ";");
  strcat(line, qual);
  return strdup(line);
}

void print_match(INDEX edge)
{
  long long location;
  char *s, *q;

  if (read_sequence_no_to_file_offset == NULL) {
    size_t rc;
    rc = pread(index_fd, &location, sizeof(long long), edge*sizeof(long long));
    if (rc != sizeof(long long)) {
      fprintf(stderr, "overlap: failed to fetch %ld bytes from offset 0x%llx on file %d, rc = %ld\n", sizeof(long long), (long long)edge*sizeof(long long), index_fd, rc); exit(1);
      exit(1);
    }
  } else {
    location = read_sequence_no_to_file_offset[edge];
  }

  s = stringat(location);q = strchr(s, ';'); if (q) *q++ = '\0';
  fprintf(stdout, "%s (read #%lld)\n", s, edge);
  free(s); s = NULL;
  printed++;
}

void print_remaining_trie(INDEX trie_index, int trie_fd, int index_fd, int offset)
{
  char *s;
  int e, sp;
  INDEX edge;
  CELL tmp;
  CELL *this;

  //fprintf(stderr, "> print_remaining_trie(%lld)\n", trie_index);
  if (trie_cell && (trie_index <= last_used_edge)) {
    this = &trie_cell[trie_index];
  } else {
    ssize_t rc = pread(trie_fd, &tmp, sizeof(CELL), trie_index*sizeof(CELL));
    if (rc != sizeof(CELL)) {
      fprintf(stderr, "overlap: failed to fetch %ld bytes from offset 0x%llx on file %d, rc = %ld\n", sizeof(CELL), (long long)trie_index*sizeof(CELL), trie_fd, rc); exit(1);
      exit(1);
    }
    this = &tmp;
  }

  for (e = _A_; e <= _N_; e++) {
    edge = this->edge[e]&EDGE_MASK;
    if (this->edge[e]&ENDS_WORD) {
      if (edge != 0LL) print_match(edge);
    } else {
      if (edge != 0LL) print_remaining_trie(edge, trie_fd, index_fd, offset);
    }    
  }
  //fprintf(stderr, "< print_remaining_trie(%lld)\n", trie_index);
}

static CELL *fetch_trie_cell(INDEX idx) {
  static CELL tmp; // NOT THREAD SAFE.  Only one call at a time.
  ssize_t rc;

  //fprintf(stderr, "Fetch: %lld\n", idx);
  rc = pread(trie_fd, &tmp, sizeof(CELL), idx*sizeof(CELL));
  if (rc != sizeof(CELL)) {
    fprintf(stderr, "nearmatch: failed to fetch %d bytes from offset 0x%llx on file %d, rc = %d\n",
            (int)sizeof(CELL), idx*sizeof(CELL), trie_fd, (int)rc); exit(1);
    exit(1);
  }
  return &tmp;
}

void lookup_read(INDEX trie_index, char *s, int actual_errors, int allowed_errors)
{
  int c;
  INDEX edge;
  CELL *this;
  //fprintf(stderr, "> lookup_read(%lld, \"%s\")\n", trie_index, s);

  this = ((trie_cell && (trie_index <= last_used_edge)) ? &trie_cell[trie_index] : fetch_trie_cell(trie_index));

  c = *s;

  if (c == '\0') {
    // This is a prefix of a string or strings at trie_cell[trie_index]
    //fprintf(stderr, "warning: target string is shorter than the reads in this database\n");
    print_remaining_trie(trie_index, trie_fd, index_fd, strlen(s));
    //fprintf(stderr, "<2lookup_read(%lld, \"%s\")\n", trie_index, s);
    return;
  }

  //fprintf(stderr, "Match: %c @%lld\n", c, trie_index);
  if (c == 'A') c = _A_;
  else if (c == 'C') c = _C_;
  else if (c == 'G') c = _G_;
  else if (c == 'T') c = _T_;
  else if (c == 'N') c = _N_;
  else {
    fprintf(stderr, "nearmatch: bad character '%c' at %s\n", c, s);
    c = _N_; // some other char
  }


  if ((this->edge[c]&ENDS_WORD) && (*(s+1) != '\0')) {
    fprintf(stderr, "warning: target string is longer than the reads in this database - ignoringing the excess at the end: %s\n", s+1);
  }

  // Apart from printing multiple copies of one read when the target string contains Ns, this works...

  edge = this->edge[c]&EDGE_MASK;
  if (edge) { // Simple literal match...
    if (this->edge[c]&ENDS_WORD) {
      //fprintf(stderr, "3: print_match(%lld)\n", edge);
      print_match(edge);
    } else {
      //fprintf(stderr, "3: lookup_read(%lld, \"%s\") this=%p\n", edge, s+1, this);
      lookup_read(edge, s+1, actual_errors, allowed_errors); // Matching an 'N' doesn't count against allowed errors
      this = ((trie_cell && (trie_index <= last_used_edge)) ? &trie_cell[trie_index] : fetch_trie_cell(trie_index)); // Restore
    }
  }

  if (c != _N_) { // avoid double matching
    int e;
    edge = this->edge[_N_]&EDGE_MASK;
    if (edge) { // Given letter in target string matches an N in the database
      if (this->edge[_N_]&ENDS_WORD) {
        //fprintf(stderr, "2: print_match(%lld) c=%d this=%p\n", edge, c, this);
	print_match(edge);
      } else {
        //fprintf(stderr, "2: lookup_read(%lld, \"%s\")\n", edge, s+1);
        lookup_read(edge, s+1, actual_errors, allowed_errors); // Matching an 'N' doesn't count against allowed errors
        this = ((trie_cell && (trie_index <= last_used_edge)) ? &trie_cell[trie_index] : fetch_trie_cell(trie_index)); // Restore
      }
    }

    // NOW DO ERROR CORRECTION
    if (actual_errors < allowed_errors) {
      for (e = _A_; e < _N_; e++) {
        edge = this->edge[e]&EDGE_MASK;
        if ((c != e) && edge) {
          if (this->edge[e]&ENDS_WORD) {
	    //fprintf(stderr, "1: print_match(%lld)\n", edge);
            print_match(edge);
          } else {
	    //fprintf(stderr, "1: lookup_read(%lld, \"%s\")\n", edge, s+1);
            lookup_read(edge, s+1, actual_errors+1, allowed_errors); // Matching an 'N' doesn't count against allowed errors
            this = ((trie_cell && (trie_index <= last_used_edge)) ? &trie_cell[trie_index] : fetch_trie_cell(trie_index)); // Restore
          }
        }
      }
    }
  } else { // An N in the target string should match any letter in the database
    int e;
    // This handles an 'N' in the given target string.
    for (e = _A_; e < _N_; e++) { // if c was _N_, it already matched the N in first test above
      edge = this->edge[e]&EDGE_MASK;
      if ((c != e) && edge) {
        if (this->edge[e]&ENDS_WORD) {
	  //fprintf(stderr, "1: print_match(%lld)\n", edge);
          print_match(edge);
        } else {
	  //fprintf(stderr, "1: lookup_read(%lld, \"%s\")\n", edge, s+1);
          lookup_read(edge, s+1, actual_errors, allowed_errors); // Matching an 'N' doesn't count against allowed errors
          this = ((trie_cell && (trie_index <= last_used_edge)) ? &trie_cell[trie_index] : fetch_trie_cell(trie_index)); // Restore
        }
      }
    }
  }
  //fprintf(stderr, "< lookup_read(%lld, \"%s\")\n", trie_index, s);
}

int main(int argc, char **argv)
{
  char *target, *target_tail;
  long long trie_index;
  int loops, indent; //, len
  char most_frequent;

  if (argc != 3) {
    fprintf(stderr, "syntax: nearmatch file.fastq ACTUAL_READ\n");
    exit(EXIT_FAILURE);
  }

  read_file = fopen(argv[1], "r");
  if (read_file == NULL) {
    fprintf(stderr, "nearmatch: cannot access fastq file %s - %s\n", argv[1], strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(trie_file_name, "%s-edges", argv[1]);
  trie_fd = open(trie_file_name, O_RDONLY);
  if (trie_fd < 0) {
    fprintf(stderr, "nearmatch: cannot access trie file %s - %s\n", trie_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(index_file_name, "%s-index", argv[1]);
  index_fd = open(index_file_name, O_RDONLY);
  if (index_fd < 0) {
    fprintf(stderr, "nearmatch: cannot access trie file %s - %s\n", index_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  file_length = lseek(trie_fd, (off_t)0LL, SEEK_END);
  //fprintf(stderr, "trie: %lld entries\n", (long long)file_length/sizeof(CELL)-1LL);
  trie_cell = mmap(NULL, (size_t)file_length, PROT_READ, MAP_PRIVATE/*SHARED*/, trie_fd, (off_t)0LL);
  if ((trie_cell == NULL) || (trie_cell == (void *)-1)) {
    //fprintf(stderr, "nearmatch: failed to map %s - %s - so using direct access instead\n", trie_file_name, strerror(errno));
    trie_cell = NULL; 
    last_used_edge = -1; // force all accesses via disk.
    // later can load as much as possible into ram, thuis caching the majority of accesses.
  }

  file_length = lseek(index_fd, (off_t)0LL, SEEK_END);
  last_used_edge = (long long)file_length/sizeof(long long)-1LL;
  //fprintf(stderr, "index: %lld entries (last_used_edge = %lld)\n", last_used_edge-1, last_used_edge);
  read_sequence_no_to_file_offset = mmap(NULL, (size_t)file_length, PROT_READ, MAP_SHARED, index_fd, (off_t)0LL);
  if ((read_sequence_no_to_file_offset == NULL) || (read_sequence_no_to_file_offset == (void *)-1)) {
    fprintf(stderr, "nearmatch: failed to map %s - %s - so using direct access instead\n", trie_file_name, strerror(errno));
    read_sequence_no_to_file_offset = NULL;
  } else {
    //fprintf(stderr, "Successfully mmap'd sequence_number_to_file_offset[%d]\n", (int)((size_t)file_length/sizeof(INDEX))-1);
  }

  target = strdup(argv[2]);

  lookup_read(ROOT_CELL, target, 0, ALLOWED_ERRORS); // could set allowed_errors to a percentage of the read length
  exit(EXIT_SUCCESS);
  return EXIT_FAILURE;
}
