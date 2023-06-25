/*
    makeafg 256seq.fastq 256seq.fastq-NEWCONTIG
 */

// Take a contig as a simple text string, and the file of reads that it was extracted from, and generate a .afg
// matching up all the reads which match 100% (or with a small allowed number of errors) so that it can be viewed
// in 'tablet' etc.

// NOT YET WRITTEN.  CURRENLT THIS IS THE CODE FROM glocate.c

#define _FILE_OFFSET_BITS 64

typedef unsigned long long EDGE;
typedef unsigned long long INDEX;
#define ENDS_WORD (1ULL<<63ULL)
#define EDGE_MASK (ENDS_WORD-1ULL)

#define MAX_LINE 1024

#define MIN_OVERLAP 13

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

static int freq[256];

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

void RED(FILE *f, long long read_id, char *seq, char *qlt);
void TLE(long long read_id, char *seq, int overlap_len, long long offset);
void walk_trie(INDEX trie_index, int trie_fd, int index_fd, int offset)
{
  char *s;
  int e, sp;
  INDEX edge;
  CELL tmp;
  CELL *this;

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

  for (e = 0; e < 5; e++) {
    edge = this->edge[e]&EDGE_MASK;
    if (this->edge[e]&ENDS_WORD) {
      long long location;
      char *q;

      if (read_sequence_no_to_file_offset == NULL) {
	ssize_t rc;
        rc = pread(index_fd, &location, sizeof(long long), edge*sizeof(long long));
        if (rc != sizeof(long long)) {
          fprintf(stderr, "overlap: failed to fetch %ld bytes from offset 0x%llx on file %d, rc = %ld\n", sizeof(long long), (long long)edge*sizeof(long long), index_fd, rc); exit(1);
          exit(1);
        }
      } else {
	location = read_sequence_no_to_file_offset[edge];
      }

      s = stringat(location);q = strchr(s, ';'); *q++ = '\0';
      for (sp = 0; sp < offset; sp++) fputc(' ', stdout);
      fprintf(stdout, "%s (read #%lld)", s, edge);
      fprintf(stdout, " %c", s[strlen(s)-offset]);
      freq[(int)(s[strlen(s)-offset])] += 1;
      RED(afg_reads_file, edge, s, q);
      TLE(edge, s, offset, contig_size);
      free(s); s = NULL;
      fputc('\n', stdout);
    } else {
      if (edge != 0LL) walk_trie(edge, trie_fd, index_fd, offset);
    }    
  }
}

static CELL *fetch_trie_cell(INDEX idx) {
  static CELL tmp; // NOT THREAD SAFE.  Only one call at a time.
  ssize_t rc;

  rc = pread(trie_fd, &tmp, sizeof(CELL), idx*sizeof(CELL));
  if (rc != sizeof(CELL)) {
    fprintf(stderr, "glocate: failed to fetch %d bytes from offset 0x%llx on file %d, rc = %d\n",
            (int)sizeof(CELL), idx*sizeof(CELL), trie_fd, (int)rc); exit(1);
    exit(1);
  }
  return &tmp;
}

long long lookup_read(INDEX trie_index, char *s)
{
  int c;
  INDEX edge;
  CELL *this;

  if (trie_cell && (trie_index <= last_used_edge)) {
    this = &trie_cell[trie_index];
  } else {
    this = fetch_trie_cell(trie_index);
  }

  c = *s;

  if (c == '\0') {
    //fprintf(stderr, "This is a prefix of a string or strings at trie_cell[%lld]\n", trie_index);
    return trie_index;
  }

  //fprintf(stderr, "Match: %c @%lld\n", c, trie_index);
  if (c == 'A') c = _A_;
  else if (c == 'C') c = _C_;
  else if (c == 'G') c = _G_;
  else if (c == 'T') c = _T_;
  else if (c == 'N') c = _N_;
  else {
    fprintf(stderr, "glocate: bad character '%c' at %s\n", c, s);
    c = _N_; // some other char
  }

  edge = this->edge[c]&EDGE_MASK;

  if (edge == 0LL) return 0LL;

  if (this->edge[c]&ENDS_WORD) {

    if (*(s+1) != '\0') {
      fprintf(stderr, "warning: target string is longer than the reads in this database - excess is: %s\n", s+1);
    }
    //    assert(*s == '\0');
    return edge;

  } else if (!(this->edge[c]&ENDS_WORD)) {
    return lookup_read(edge, s+1);
  }
  return 0LL; // NOT FOUND.

}

// Eliminate duplicates and sort into numerical order.  Doesn't scale well due to external sort program,
// but as long as we're just exploring individual contigs and not assembling the whole genome, this should
// work fine.  And if we were assembling the whole genome, we'd generate the reads.afg by translating our
// fastq file with no need to eliminate duplicates

// sort -u < 256seq.fastq-GTGTCCATTATGTTGAACAAGAACAGTTCCATCTCCTTTT-reads.afg |sort -n -k2 -t:|tr ' ' '\n'> ZZZ;mv ZZZ 256seq.fastq-GTGTCCATTATGTTGAACAAGAACAGTTCCATCTCCTTTT-reads.afg

void RED(FILE *f, long long read_id, char *seq, char *qlt)
{ // just the raw reads
  fprintf(f, "{RED ");
  fprintf(f, "iid:%lld ", read_id);
  fprintf(f, "eid:%lld ", read_id);
  fprintf(f, "seq: %s . ", seq);
  fprintf(f, "qlt: %s . ", qlt);
  fprintf(f, "}\n");

}

void CTG_begin(void)
{ // Contig_t: Sequence_t
  static int next = 0;
  next++;

  fprintf(afg_contig_file, "{CTG\n");
  fprintf(afg_contig_file, "iid:%d\n", next);
  fprintf(afg_contig_file, "eid:%d-0\n", next);
  fprintf(afg_contig_file, "seq:\n");
  // output the contig here
}

void CTG_end(char *firstq,long long contig_length) {
  long long i;
  fprintf(afg_contig_file, "\n.\n");
  fprintf(afg_contig_file, "qlt:\n");
  // output the quality info here.  We'll fake it.
  for (i = 0LL-strlen(firstq); i < contig_length; i++) {
    fputc(64, afg_contig_file);
    if (((i % 80LL) == 79LL) && (i+1LL != contig_length)) fputc('\n', afg_contig_file);
  }
  fprintf(afg_contig_file, "\n.\n");
  // now insert all the tiling info that was writen to the .tle file ...
  fprintf(afg_contig_file, "}\n");
}

// sort -n -t: -k2 < 256seq.fastq-GTGTCCATTATGTTGAACAAGAACAGTTCCATCTCCTTTT-tle.afg > ZZZ;mv ZZZ  256seq.fastq-GTGTCCATTATGTTGAACAAGAACAGTTCCATCTCCTTTT-tle.afg
// and pass through only the highest overlap for each read
void TLE(long long read_id, char *seq, int overlap_len, long long offset) {
  fprintf(afg_tle_file, "{TLE ");
  fprintf(afg_tle_file, "src:%lld ", read_id);
  fprintf(afg_tle_file, "off:%lld ", offset+overlap_len);
  fprintf(afg_tle_file, "clr:0,%03d ", (int)(strlen(seq)-overlap_len));
  fprintf(afg_tle_file, "}\n");
}

int main(int argc, char **argv)
{
  char *target, *target_tail;
  long long trie_index;
  int loops, indent; //, len
  char most_frequent;

  if (argc != 3) {
    fprintf(stderr, "syntax: glocate file.fastq ACTUAL_READ\n");
    exit(EXIT_FAILURE);
  }

  read_file = fopen(argv[1], "r");
  if (read_file == NULL) {
    fprintf(stderr, "glocate: cannot access fastq file %s - %s\n", argv[1], strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(gene_file_name, "%s-%s", argv[1], argv[2]);
  gene_file = fopen(gene_file_name, "w");
  if (gene_file == NULL) {
    fprintf(stderr, "glocate: cannot write contig file %s - %s\n", gene_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(gene_file_name, "%s-%s-reads.afg", argv[1], argv[2]);
  afg_reads_file = fopen(gene_file_name, "w");
  if (afg_reads_file == NULL) {
    fprintf(stderr, "glocate: cannot write afg reads file %s - %s\n", gene_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(gene_file_name, "%s-%s-contig.afg", argv[1], argv[2]);
  afg_contig_file = fopen(gene_file_name, "w");
  if (afg_contig_file == NULL) {
    fprintf(stderr, "glocate: cannot write afg contig file %s - %s\n", gene_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(gene_file_name, "%s-%s-tle.afg", argv[1], argv[2]);
  afg_tle_file = fopen(gene_file_name, "w");
  if (afg_tle_file == NULL) {
    fprintf(stderr, "glocate: cannot write afg contig file %s - %s\n", gene_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }
  // the three afg files are later concatenated into a single afg file suitable for use in 'tablet'

  sprintf(trie_file_name, "%s-edges", argv[1]);
  trie_fd = open(trie_file_name, O_RDONLY);
  if (trie_fd < 0) {
    fprintf(stderr, "glocate: cannot access trie file %s - %s\n", trie_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  sprintf(index_file_name, "%s-index", argv[1]);
  index_fd = open(index_file_name, O_RDONLY);
  if (index_fd < 0) {
    fprintf(stderr, "glocate: cannot access trie file %s - %s\n", index_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  file_length = lseek(trie_fd, (off_t)0LL, SEEK_END);
  //fprintf(stderr, "trie: %lld entries\n", (long long)file_length/sizeof(CELL)-1LL);
  trie_cell = mmap(NULL, (size_t)file_length, PROT_READ, MAP_PRIVATE/*SHARED*/, trie_fd, (off_t)0LL);
  if ((trie_cell == NULL) || (trie_cell == (void *)-1)) {
    //fprintf(stderr, "glocate: failed to map %s - %s - so using direct access instead\n", trie_file_name, strerror(errno));
    trie_cell = NULL; 
    last_used_edge = -1; // force all accesses via disk.
    // later can load as much as possible into ram, thuis caching the majority of accesses.
  }

  file_length = lseek(index_fd, (off_t)0LL, SEEK_END);
  last_used_edge = (long long)file_length/sizeof(long long)-1LL;
  //fprintf(stderr, "index: %lld entries (last_used_edge = %lld)\n", last_used_edge-1, last_used_edge);
  read_sequence_no_to_file_offset = mmap(NULL, (size_t)file_length, PROT_READ, MAP_SHARED, index_fd, (off_t)0LL);
  if ((read_sequence_no_to_file_offset == NULL) || (read_sequence_no_to_file_offset == (void *)-1)) {
    fprintf(stderr, "glocate: failed to map %s - %s - so using direct access instead\n", trie_file_name, strerror(errno));
    read_sequence_no_to_file_offset = NULL;
  } else {
    //fprintf(stderr, "Successfully mmap'd sequence_number_to_file_offset[%d]\n", (int)((size_t)file_length/sizeof(INDEX))-1);
  }

  target = strdup(argv[2]);
  fprintf(stdout, "%s\n", target);
  fprintf(gene_file, "%s", target);
  contig_size = 0LL;
  CTG_begin();
  fprintf(afg_contig_file, "%s", target);
  loops = 0;
  for (;;) {
    target_tail = target;
    freq['C'] = freq['G'] = freq['A'] = freq['T'] = freq['N'] = 0;
    indent = 0;
    trie_index = lookup_read(ROOT_CELL, target_tail);
    if (trie_index != 0ULL) {
      // output root read in various forms
      char *s, *q;
      //int sp;
      INDEX edge = trie_index;
      long long location;

      if (read_sequence_no_to_file_offset == NULL) {
	ssize_t rc;
        rc = pread(index_fd, &location, sizeof(long long), edge*sizeof(long long));
        if (rc != sizeof(long long)) {
          fprintf(stderr, "overlap: failed to fetch %ld bytes from offset 0x%llx on file %d, rc = %ld\n", sizeof(long long), (long long)edge*sizeof(long long), index_fd, rc); exit(1);
          exit(1);
        }
      } else {
	location = read_sequence_no_to_file_offset[edge];
      }

      s = stringat(location);q = strchr(s, ';'); *q++ = '\0';
      fprintf(stdout, "%s (read #%lld)", s, edge);
      RED(afg_reads_file, edge, s, q);
      TLE(edge, s, 0, contig_size);
      free(s); s = NULL;
      fputc('\n', stdout);
    }
    for (;;) { // we only look one character ahead
      target_tail += 1;
      loops += 1; indent += 1;
      //    next_char = 0;
      trie_index = lookup_read(ROOT_CELL, target_tail);
      if (trie_index == 0ULL) continue;
      walk_trie(trie_index, trie_fd, index_fd, indent);
      if (strlen(target_tail) < 16) break;
    }

    // this is crude.  What we need here is to identify buckets of frequencies, taking into
    // account the random sampling rate.  We need to be able to identify when there is a fork
    // in the road, so that the user can make a manual decision which branch to follow.
    most_frequent = 'N';
    if ((freq['C'] >= freq['G']) && (freq['C'] >= freq['A']) && (freq['C'] >= freq['T'])) most_frequent = 'C';
    if ((freq['G'] >= freq['C']) && (freq['G'] >= freq['A']) && (freq['G'] >= freq['T'])) most_frequent = 'G';
    if ((freq['A'] >= freq['C']) && (freq['A'] >= freq['G']) && (freq['A'] >= freq['T'])) most_frequent = 'A';
    if ((freq['T'] >= freq['C']) && (freq['T'] >= freq['G']) && (freq['T'] >= freq['A'])) most_frequent = 'T';

    {
    int slop_factor = freq[(int)most_frequent]/8;
    if (most_frequent == 'C') {
      if ((freq['G']+slop_factor >= freq['C']) || (freq['A']+slop_factor >= freq['C']) || (freq['T']+slop_factor >= freq['C'])) { // (typo fixed! - rerun last results!)
        most_frequent = 'N';
      }
    } else if (most_frequent == 'G') {
      if ((freq['C']+slop_factor >= freq['G']) || (freq['A']+slop_factor >= freq['G']) || (freq['T']+slop_factor >= freq['G'])) {
        most_frequent = 'N';
      }
    } else if (most_frequent == 'A') {
      if ((freq['C']+slop_factor >= freq['A']) || (freq['G']+slop_factor >= freq['A']) || (freq['T']+slop_factor >= freq['A'])) {
        most_frequent = 'N';
      }
    } else if (most_frequent == 'T') {
      if ((freq['C']+slop_factor >= freq['T']) || (freq['G']+slop_factor >= freq['T']) || (freq['A']+slop_factor >= freq['T'])) {
        most_frequent = 'N';
      }
    }
    if (most_frequent == 'C') {
      if (freq['G']+freq['A']+freq['T'] > freq[(int)most_frequent]-slop_factor) {
        most_frequent = 'N';
      }
    } else if (most_frequent == 'G') {
      if (freq['C']+freq['A']+freq['T'] > freq[(int)most_frequent]-slop_factor) {
        most_frequent = 'N';
      }
    } else if (most_frequent == 'A') {
      if (freq['C']+freq['G']+freq['T'] > freq[(int)most_frequent]-slop_factor) {
        most_frequent = 'N';
      }
    } else if (most_frequent == 'T') {
      if (freq['C']+freq['G']+freq['A'] > freq[(int)most_frequent]-slop_factor) {
        most_frequent = 'N';
      }
    }
    }

    memmove(target, target+1, strlen(target+1));
    target[strlen(target)-1] = most_frequent;
    if (most_frequent == 'N') break;
    fputc(most_frequent, gene_file); fflush(gene_file);
    fputc(most_frequent, afg_contig_file);  if ((contig_size % 80LL) == 79LL) fputc('\n', afg_contig_file); fflush(afg_contig_file); contig_size += 1LL;



    { // Warn if there is any other significant branching choice, even if it is not as common as the primary one...
    fprintf(stdout, "#%c  c: %d  g: %d  a: %d  t: %d @%d\n%s\n",
	    most_frequent, freq['C'], freq['G'], freq['A'], freq['T'], contig_size, target);
    }
  }

  CTG_end(target, contig_size);
  fprintf(stdout, "#N  c: %d  g: %d  a: %d  t: %d\n", freq['C'], freq['G'], freq['A'], freq['T']);
  fprintf(stderr, "Exited with target = \"%s\", c: %d  g: %d  a: %d  t: %d\n", target, freq['C'], freq['G'], freq['A'], freq['T']);
  exit(EXIT_SUCCESS);
  return EXIT_FAILURE;
}