/*< This is the source code for a program that builds a fast data structure
    containing all the reads from a gene sequencing system.  Click on comments
    like this one to drill down to the code that implements the main algorithm. */

/*< This source file is a complete program and no other files are required to
    duplicate/test this work.  Download maketrie.c from the same web directory
    for a compilable copy.

    It compiles with something like:

        <tt>mpicc -openmp -O -o maketrie -Wall maketrie.c</tt>

    and can be used with your own data in fastq format by:

        <tt>./maketrie file.fastq</tt>

    The code will run on a single node but does need to be compiled on an MPI
    system.  It currently uses MPI to scale the amount of RAM used to the problem
    size, however the more memory that each processor has, the more efficiently
    the code will run.  A run of an earlier version of this code (which determined
    overlaps as well as building the trie) took 24hrs on 8 regular compute nodes,
    but ran in 4 hrs on a single compute node which had 8 times the normal amount
    of RAM.

    This version of the code is deliberately simple and does not use any
    optimisations that would obfuscate the algorithm, in order to serve as a
    blueprint for others to include this algorithm their gene assembly pipelines.

    Graham Toal &lt;<tt>gtoal@gtoal.com</tt>&gt;
    May 2013.
 */
/*
      This code is intended to solve one part of de novo gene assembly.

   Currently it solves a simplified version of the problem described as
   follows:

     * Generate a gene (a string of length 'G') consisting of the letters CGAT
       chosen randomly.  G is a large number between tens of thousands and
       billions.  For this test, we do *not* read the complementary string.
       (Add that later. Shouldn't add significantly more complexity to
        overlap detection, but will double the amount of data to be handled)
       Neither do we yet handle repeats (eg take a 4k substring and insert
       elsewhere, optionally with about 1 in 10 random letter changes to
       simulate evolutionary convergence - not relevant to overlap detection,
       more significant for contig reassembly)

     * Sample that string at random locations and extract 'K' substrings
       of length 'L'.  (These are called 'Reads' or 'Short Reads'). 
       L is fixed for any one extract.
       We're using L=37 although this code is independent of the size
       of L.  The output READs are *not* tagged with the location they
       were found at in the gene.  (see accompanying program genelab.c for
       the software that creates this dummy data).  We generate K READs
       such that there should be sufficient samples to cover G and that
       all READs have a high probability of overlapping with a nearby
       READ.

     * read in the READs and build a trie in memory containing all the READs.
       This trie is likely to be larger than the ram capacity of a single
       compute node.  (A secondary goal is to see how large a gene can
       be handled on a single compute node)

     * trie-building is not a parallel activity but we can make use of
       the memory of multiple nodes to keep the entire trie in RAM with
       relatively little inter-processor communications overhead. Early
       test results of trie-construction for G=100,000,000, K=40,000,000,
       and L=37 are 5 minutes on a single node; for G=200,000,000,
       K=80,000,000, and L=37 we require 3 nodes and the process takes 15
       minutes.  There is an overhead for using more than one node but it is
       not significant.

    *  A second program (findoverlaps) then determines all the overlaps
       between READs, and outputs a graph of these overlaps to a file.
       A subsequent phase (preferably written by someone else) will extract
       the lowest-cost path through that graph which covers the majority of
       READs.  Locating the overlaps is now the most expensive part of this
       code, but fortunately it is embarassingly parallel!  By grouping nodes
       in subgroups of sufficient nodes to hold a full copy of the trie, we
       can have multiple copies of these subgroups, each of which can handle
       another subset of N.

   (duplicate READs are found on the fly as they are added to the trie.
   We only record the first occurance in the trie, writing the duplicate
   information out to file at the point that the duplication is detected.
   Unlike the later overlap detection, this process can happen on any
   of the compute nodes, so there are multiple output files containing
   subsets of the results, which need to be externally concatenated after
   the run completes.  This isn't a high overhead (it's just "cat") but
   is a hidden cost when calculating the effectiveness of this algorithm)

   At this point in development, this code builds the trie as above.
   Once the trie is built correctly and confirmed, I will add the code
   to calculate overlaps.  This is actually a relatively low cost function
   whose complexity is K * L, ie effectively linear.  It can be parallelised
   trivially in L (using OMP) and with a little more effort in K (using MPI),
   ie overlap determination is embarassingly parallel.

   To determine an overlap, take a single READ and walk the trie from
   the first letter in the READ.  If the trie-walk has not prematurely
   terminated before hitting the end of the READ, there is an overlap.
   The overlap can be defined by storing the trie index at the furthest
   point that was reached, and later by walking what is left in that 
   subtree from that stored index, you will reach the end of all the
   matching READs.  An identifier for the READ is stored at the end
   of the trie entry for that READ.

   Then repeat the last processes starting at the second, third, etc
   letter in the READ to locate overlaps starting at those offsets.
   This can either be done in a loop of size L or by L processes in
   parallel.

   The cost of a single comparison is O(L) (reducing from L-1 to 2 (or
   'MIN_OVERLAP' which will be larger - explanation elsewhere) as the
   overlap length shortens.  An overlap of 1 is not useful so we will
   probably have a cut-off point later in terms of a minimal overlap -
   not defined by length but by how many other READs overlap at this
   offset...  A graph of offset vs #overlap_count is also potentially
   useful to determine critical sections of overlap - if overlaps
   were random, it would be a logarithmic curve. If there are flat
   sections, that tells us about a point of inflexion which may be
   useful in graph re-assembly later)

   Because trie lookups do not modify the trie, they are safe to execute
   in parallel. (OMP, not MPI)

   The final output of this code will be a file containing a complete
   overlap graph which can be used to reconstitute G.  Although I have
   some ideas on how this might be done, I am going to restrict this
   project to just the overlap calculation part to keep it tractable.

   I do understand that the subset of the problem that this code solves
   is not exactly the same as the real-life problem, but I'm throwing
   it out there to people working in this area to see if the two domains
   are sufficiently close that this can be adapted into a real system.

   Please supply feedback including requests for any key components that
   would make it applicable to real life gene assembly to:

   Graham Toal <gtoal@gtoal.com>
   (started in mid April 2013, working reliably by May 19th, 2013)

   This code was tested using four nodes on the Texas Advanced Computer
   Center's "Stampede" cluster.

   PS This program has been algorithmically optimised but not low-level
   optimised.  A lot of operations can be pipelined and the low-level
   code could be improved.  However it is sufficiently fast that I
   haven't felt the need so far, and it is a lot easier to understand
   the way it is written compared to what I would have to do to it to
   get perhaps a factor of 2 speedup at best.
 */


/*< <em>(Current 'To do' list.)</em> */
/*> TO DO:
      7) Work out why there's a lost 15 minutes at the start of the
         matching.  Is it file-system related?  Could we go faster
         by walking the dawg than by re-reading the index file???
         [Only at TACC - doesn't happen on UTPA cluster]
     11) Store data in the terminal nodes, more than just a pointer to
         fastq entry for that read.  Could store dup count, quality data,
         ultimately pointers to contig data, ie make this a full graph
         suitable for input to a gene assembly graph algorithm.
         Note that to preserve the purity of the space allocation (ie
         the almost-free 'last_used_node++'), we must make sure that
         any data we want to attach to the leaves *must* be made to fit
         into 40 bytes, whether it fits or not.

         NOTE: 5 long longs * 8 bytes * 8 bits = 320 bits.  Even using
         3 bits per letter (to allow us to encode 'N's), that is 106 letters
         we could store in a leaf!

         So if storing the actual read at the leaf is helpful, we
         could esily do it.

         (And of course with only minimal effort we could use >1 leaf node,
         if we really can't cope with the Procustean solution)
     14) Because I used 0L as a special NULL pointer for trie edges (ROOT_NODE being 1L),
         I may have screwed up and made read #0 impossible to access.  I think
         I should renumber the reads from a base of 1 up.  Serendipitously,
         it appears from various documentation that I've since read, that
         reads are conventially numbered from 1 upwards anyway!  So it won't
         hurt to do a mass renaming.
     15) BIG ONE HERE!: if we re-build the trie from the sorted data, adding
         from L=1 to L=K using iterative deepening, we can guarantee to place
         the first <N> chars of each string in the first memory block.  Since
         the majority of searches of matches are failures, this will speed up
         the failures tremendously!  And the beauty of it is that this is an
         optional post-processing step that will speed up the data structure
         but not change the contents in any significant way that affects
         other code.  Well worth doing if a saved trie is going to be accessed
         repeatedly, especially with the code that runs on a single machine
         and caches the first chunk of the trie_cell[] array...
     16) In the short term:
         give all available memory if only using 1 node.
         give all avaiable memory to last node if using > 1 nodes - remember to
         tweak chunksize to match so that it doesn't try to pass calls on to
         a non-existant downstream node
         In the long term:
         Use modulo arithmetic rather than binary-and to separate by chunk size,
         and remove constraint that chunks must be powers of two.
 */
/*>*//*>*/

/*< Declarations */

/*< '<tt>#include</tt>' header files. */
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>
#include <omp.h>
#include <ctype.h>
#include <time.h>  // for info only

// #define NDEBUG 1 // Set this if you want to remove the overhead of the assert tests.
#include <assert.h>

#ifndef FALSE
#define TRUE (0==0)
#define FALSE (!TRUE)
#endif
/*>*/

/*< Input/Output Files. */
static FILE *duplicates = NULL; // only rank 0 can write to this file.
static FILE *rejects = NULL; // only rank 0 can write to this file.
static FILE *read_file = NULL;
static FILE *read_index = NULL;
static FILE *trie_file = NULL;
static FILE *sorted_and_unique_reads = NULL; // only rank (mpisize-1) can write to this file.
/*>*/

/* This code goes to some effort to work out at runtime how much RAM is
   available, and it allocates as much as it can up front for use by the
   main data structure which is an array of structs.

   In the event that this code is ported to a system where the dynamic
   memory determination doesn't work, we do have a way to fall back on
   a hard-coded default, which in the distribution is set to the relatively
   low value of 4Gb below, in the hope that the code 'just works' the
   first time someone compiles it.  If this affects you and you're unable
   to improve the dynamic memory-sizing code to understand your system,
   increase the number below to suit your hardware... */
#define GBs 4L

/*< <em>(A note on 64-bit architecture)</em> */
// When I first wrote this code I tried to make it optional whether you
// used 32 bit or 64 bit data, to allow for a version for people on older
// processors.  However at this point I doubt that anyone is trying to do
// gene assembly on 32 bit architectures and if they were, they wouldn't
// be able to handle very large data, so I intend to simplify the code now
// and assume 64-bit "long long" data throughout.

typedef unsigned long long EDGE;
typedef unsigned long long INDEX;
#define ENDS_WORD (1ULL<<63ULL)
#define EDGE_MASK (ENDS_WORD-1UL)
/*>*/

/*< THE MAIN DATA STRUCTURE */
// the structure below is a 5-way branching trie.  Each edge contains a
// pointer (actually an index) to other edges in the array, and packed
// into each edge is a bit to denote 'end of word' (ie the last letter
// in a READ).  To simplify the structure I have removed earlier support
// for multiple-length strings.  Now all strings must be the same
// length.  The actual length is not important, just that every string
// added to te data structure is the same length as every other.

// We could have saved some space by using a 4-way trie rather than a
// 5-way trie, but I included the 'N' code in reads in the hope of
// future-proofing the code.  It's not much of an overhead.

// If I occasionally refer to this structure in the code below as the 'DAWG'
// that is because it is a variation of a structure used heavily in spelling
// checkers/correctors and in word games such as boggle, scrabble, etc.
// The mnemonic stands for "Directed Acyclic Word Graph".
// In this case we stay with the functionally compatible trie and do not
// perform the tail-compression that would convert the trie into a dag.

// (This is forced on us by the fact that we tag each READ at the end of
// its trie with the sequence number of the READ so that we can build
// the connectivity/overlap graph, thus the tails are no longer identical)

// Note: although not done yet, the use of a DAWG allows for cheap detection
// of 1-letter errors.  (Actually n-letters for any n, but I don't know
// enough about the problem domain to know if n>1 is helpful?)

// See also...
// http://www.slideshare.net/agbiotec/overview-of-genome-assembly-algorithms
// for some background on the whole process


// A dynamically-allocated array [0..MAX_SIZE-1] of these trie edge cells forms the main data structure.
// The array is partitioned across all available compute nodes.  We run one task ("PE") per
// compute-node and give it all available memory.  The code will work if allocated to one
// task per core instead, but engenders an unnecessary inter-node communication overhead.

typedef struct cell {
  // ACGT maps to 0 1 2 3, anything else such as 'N' is 4
  EDGE edge[5];
} CELL;
CELL *trie_cell;

static INDEX MAX_SIZE = ((INDEX)0L);

#define ROOT_CELL ((INDEX)1L)
// Node 0 is unused, 0 is needed as a terminator.

// originally I used 'next_free_edge' and it was exclusive.  Have now changed that
// to 'last_used_edge' which is inclusive.  getting a free edge is now a function
// call that we can make atomic with an OMP PRAGMA.
static INDEX last_used_edge = ROOT_CELL; // should this be 0L instead since it is inclusive?  ROOT_CELL is special, but where do I initialise it?

static CELL empty; // A constant cell, containing 0LL all nodes.

// trie_cell[edge].edge[_G_] will have ENDS_WORD clear and will point to 0 if there is no G at that position
//                                have ENDS_WORD clear and will point to another trie_cell if it is a G but not the last letter in the string
//                                have ENDS_WORD set if G is the last letter in the string, and it will point to the read number


/*< <em>('To do' notes on the data structure)</em> */
// It would be more consistent if this were shared among all nodes.
// the fact that it is not is only OK because only the last node is
// likely to be the target of adding a new edge, but I should code
//an explicit test for this just to be sure...

// this output confirms that we need to feed back updates of last_used_edge to all nodes :-(  (200m test, 3 nodes, 64 bit):

//  maketrie.721908.out:Node 2 exiting cleanly.  local base = 1073741824,  last_used_edge = 1602285502,  local maximum = 1610612736
//  maketrie.721908.out:Node 1 exiting cleanly.  local base = 536870912,  last_used_edge = 1096686229,  local maximum = 1073741824

// *OR* need to restructure add_read so that a call to the next node does the increment implicitly, leaving our
// own last_used_edge permanently stuck at local maximum.  eg   the_new_edge = add_read_at_last_used_edge(...)
/*>*//*>*/


/*< (minor housekeeping: translate gene letters ACGTN to internal codes 0..4 on input and output */
#define _A_ 0
#define _C_ 1
#define _G_ 2
#define _T_ 3
#define _N_ 4

char *trt = "ACGTN"; // translate table back from the above
/*>*/

/*< (Statistics gathering) */
// Actually, the stats need to be gathered across multiple nodes, which is not
// currently implemented.  This data isn't used anywhere to we could do away with
// this altogether, or when all the important code is written we could revisit
// this and implement it properly using an RPC call to communicate the data
// across each of the nodes.
static long freq[256];
static long letters = 0L;
static int seq = 0, dups = 0; // should seq be a long int?  Do we get that many reads?
/*>*/

/*< Input of reads.  "fastq" format etc. */
// There is a lbrary for reading fastq format which I may use later, but the
// format is so trivial I currently handle it ad hoc on the fly.
// We do want to optimise the I/O heavily as we are using such large
// files that I/O overhead becomes noticable (at least, if we are
// successful in getting the actual processing to be so fast that the
// I/O is no longer swamped by the processing times).

#define MAX_LINE 1024 // Longest input line.  We're looking at lines of
                      // 37, 60, 70, 101 etc - maybe even in the 400's
                      // in some cases, so 1K should be overkill. 
                      // I want to avoid all use of the heap if possible
                      // in order to make memory usage predictable.

// We check the length of all READs and raise an error if they are not all the same.
static long length[MAX_LINE];
static int read_length = 0; // This is the length we found in the READ file.
/*>*/

/*< MPI support. */
//------------------------------------------------------------------------
// This hard-wired default shouldn't be needed now unless we
// fail to grope the peripherals correctly...

// thumper: #define CORES_PER_NODE 12ULL
// stampede:
#define CORES_PER_NODE 16ULL

// MPI support:
static int mpirank = 0, mpisize = 1;

// Used in calculation to detemine how much memory each task on a node
// should allocate for itself.  If we have more than one task per node,
// the first one is in danger of claiming all available memory and the
// subsequent ones will fail to claim memory, resulting in an error (or
// worse, possibly a crash)

static long long TASKS_PER_NODE, PROCESSORS_PER_NODE;
// for example with 16 processors per node we can either run 16 tasks
// each with 1/16 of the available memory, or 1 task using all the memory
/*>*/

/*< The code is structured around a simulated Remote Procedure Call mechanism,
    implemented on top of MPI.  Rank 0 remains an imperative main program
    but all the other ranks run MPI listeners to handle procedure calls
    originating in rank 0. */

// All MPI calls made by the program must follow the same template.

// 1) send the command code.  This may include one optional 'long' data item.  No Ack.
// 2) Send any more parameters.  Send multiple items in a single call when possible.
// 3) Receive results
// 4) if no results expected, receive empty 'done' acknowlegement

// NOTE.  Callers should be single-threaded *FOR NOW* - don't ever send two requests in parallel
// however receivers can handle more than one caller in parallel (one call per calling system)

// It should be possible to do trie-building in parallel as long as we lock critical 
// sections of local_add_read with omp pragmas.  Critical sections being where an
// element is updated. (The primitive has to be read-modify-write, not write.)

// Receiving process should:
// 1) receive command code.
// 2) if more parameters to follow, receive them *from caller*
// 2) If no more parameters to follow, execute local command
// 3) Send results to caller
// 4) If no results, send 'done' acknowlegement to caller

// Note in next iteration, assuming it is safe, we can drop the acknowlegements
// and at some point we hope to use OMP to allow parallel dispatch loops.

// Data:
#define TAG_DATA 1
#define TAG_ACK 2

// Low-level data access commands
#define TAG_SEND_RAW_MEM 3
#define TAG_READ_READ 4
#define TAG_WRITE_READ 5
#define TAG_EXIT_PROGRAM 6

// Application-level Remote Procedure Calls
#define TAG_ADD_READ 7
#define TAG_GET_NEXT_FREE_EDGE 8
#define TAG_OUTPUT_DUPINFO 9
#define TAG_OUTPUT_READ 10
#define TAG_WALK_AND_PRINT_TRIE_INTERNAL 13
#define TAG_DUMP_TRIE 14
/*>*/

/*< Create a virtual array using the memory of as many nodes as are available.
   Use one process per node and give it all available memory.

   We assume that memory rather than CPU is the bottleneck with our algorithm - we
   don't need high concurrency - one processor would be enough, if there were
   sufficient RAM available.  Until the day when we routinely have onboard memories
   measured in Terabytes, this will have to do as a slow approximation.

   (Need to compare to speed using vSMP system built out of 8 nodes)

   This paper may have some relevance: http://www.cs.berkeley.edu/~bonachea/upc/mpi2.html

   If we use higher-level calls (eg add_string_to_trie) then in most cases we'll have one
   cross-machine call per hardware node or less, rather than one for every letter in the READ. */

// initially these were #define'd but now we work them out on the fly...
// Chunk is a power of two.  It is a count of trie cells, not a size of RAM
// There are CHUNKSIZE trie cells per instance of this code.  Multiply by
// mpisize for total size across all nodes.  (MAX_SIZE)
static long long CHUNKBITS, CHUNKSIZE, CHUNKMASK;
/*>*//*>*/

/*< Procedures */

/*< Support procedures. */

static void shut_down_other_nodes(void) // what the user calls
{
  int target_rank;
  long int value = 0L;
  MPI_Status status;

  for (target_rank = 1; target_rank < mpisize; target_rank++) {
    if (target_rank != mpirank) { // not me.  Or 0, unfortunately.
      MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_EXIT_PROGRAM, // The command.
	   MPI_COMM_WORLD);/* always use this */

      // Get acknowlegement that it has been written (to synchronise)
      MPI_Recv(&value,/* message buffer - not used */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   /*rank,*/     MPI_ANY_SOURCE,/* receive from any sender */
	   /* TAG_ACK,*/ MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
    }
  }
}

// Send raw memory via MPI_Send.
// Using this assumes a homogeneous cluster (same byte sex everywhere).  Not a problem for us.
// If used on structs rather than arrays, remember to handle machine-dependent padding issues.

static void send_bytes(int dest, void *memp, long bytes) { // just a type-checked macro
  MPI_Send(memp, bytes, MPI_BYTE, dest, TAG_SEND_RAW_MEM, MPI_COMM_WORLD);
}
/*>*/

/*< remote_*: these procedures call a later rank.  These are mostly internal housekeeping procedures, not application calls */
static void remote_walk_and_print_trie_internal(int target_rank, char *s, EDGE edge, int len)
{
  MPI_Status status;
  long value = 0L; // generic up-front parameter
  long stringlength;

  stringlength = strlen(s)+1;
  // COMMAND CODE
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_WALK_AND_PRINT_TRIE_INTERNAL,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // SEND PARAMETERS - later optimise by sending one param via 'value' above
  MPI_Send(&stringlength,/* message buffer - address part */
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  send_bytes(target_rank, s, stringlength);

  MPI_Send(&edge,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&len,
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

static void remote_dump_trie(int target_rank, char *filename)
{
  MPI_Status status;
  long value = 0L; // generic up-front parameter
  long stringlength;

  stringlength = strlen(filename)+1;
  // COMMAND CODE
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_DUMP_TRIE,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // SEND PARAMETERS - later optimise by sending one param via 'value' above
  MPI_Send(&stringlength,/* message buffer - address part */
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  send_bytes(target_rank, filename, stringlength);

  // Get acknowlegement that it has been written (to synchronise)
  MPI_Recv(&value,/* message buffer - used for len function result */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
}





// what setread calls internally to ask another node to assign to a trie cell
static void remote_setread(int target_rank, INDEX index, CELL value)
{
  long int dummy;
  MPI_Status status;
  //static int debug = 0;
  // write this data in the part of the virtual array held on another node

  // The code is essentially serial - we are only using the other nodes
  // as a cache for memory

  // This is just a sketch from some example pages on the net - NOT the actual code we need:

  // some caveats about sending multiple messages:
  //    http://stackoverflow.com/questions/2024214/mpi-buffered-send-receive-order

  assert(target_rank > mpirank);

  MPI_Send(&dummy,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_WRITE_READ,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&index,/* message buffer - address part */
	   1,/* one data item */
	   MPI_LONG_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // BEWARE!!!!  On the system I'm using, sizeof(long) == sizeof(long long).
  // So the options to select word size are not being properly tested.

  // I do have an older 32-bit machine which supports "long long" in
  // software, and I would expect it to complain about a few of these
  // assertions.

  // What we *should* be doing is use:
  // (sizeof(INDEX) == sizeof(long long) ? MPI_LONG_LONG : MPI_LONG)
  // in the MPI arguments.

  assert(sizeof(value.edge[0]) == sizeof(long));
  MPI_Send(&value.edge,/* message buffer - data part */
	   sizeof(value.edge)/sizeof(value.edge[0]),/* number of data items */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // Get acknowlegement that it has been written (to synchronise)
  // We are being very careful here to avoid parallelism
  // Next iteration once all the code is working reliably will be
  // to fire & forget (ie remove the call below)
  MPI_Recv(&value,/* message buffer - not used */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
}

// what getread calls internally to ask another node to supply us with a trie cell
static void remote_getread(int target_rank, INDEX index, CELL *valuep)
{
  MPI_Status status;
  long dummy;

  // read data from the part of the virtual array held on another node

  // The code is essentially serial - we are only using the other nodes
  // as a cache for memory

  // This is just a sketch from some example pages on the net - NOT the actual code we need:

  MPI_Send(&dummy,/* message buffer - not actually used */
	   1,/* one data item */
	   MPI_LONG,/* data item is a long integer */
	   target_rank,/* destination process rank */
	   TAG_READ_READ,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  MPI_Send(&index,/* message buffer - address part */
	   1,/* one data item */
	   MPI_LONG_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  assert(sizeof(valuep->edge[0]) == sizeof(long));
  MPI_Recv(valuep,/* message buffer */
	   sizeof(valuep->edge)/sizeof(valuep->edge[0]),/* number of data items */
	   MPI_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */
}
/*>*/

/*< (internal support for accessing virtual array) */
// setread and getread are unlikely to be used except possibly in debugging.  They allow us
// to access individual elements of the virtual array regardless of which rank the data is
// stored on.  Not to be used for high numbers of accesses as this is a relatively
// inefficient mechanism (compared to add_read which takes great advantage of locality)

// These do come in useful when prototyping a new procedure but once the logic has been
// worked out with an implementation in Node 0, it is almost always more efficient to
// use our 'migrating process' trick to run the logic local to whichever node hosts the
// relevant data.

static void setread(INDEX index, CELL value) // what the user calls
{
  int target_rank;

  target_rank = index >> CHUNKBITS;
  if (target_rank == mpirank) {
    int i;
    for (i = 0; i < 5; i++) trie_cell[index & CHUNKMASK].edge[i] = value.edge[i];
  } else {
    if (target_rank < mpirank) {
      fprintf(stderr, "PROGRAM BUG: Node %d requested access to trie_cell[%lld] on node %d - assert that we never feed backwards...\n", mpirank, index, target_rank);
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (target_rank >= mpisize) {
      fprintf(stderr, "ERROR: array bounds exceeded!  Requested access to trie_cell[%lld]\n", index);
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    remote_setread(target_rank, index, value);
  }
}

static void getread(INDEX index, CELL *valuep) // what the user calls
{
  int target_rank;

  target_rank = index >> CHUNKBITS;
  if (target_rank == mpirank) {
    int i;
    for (i = 0; i < 5; i++) valuep->edge[i] = trie_cell[index & CHUNKMASK].edge[i];
  } else {
    if (target_rank < mpirank) {
      fprintf(stderr, "PROGRAM BUG: Node %d requested access to trie_cell[%lld] on node %d - assert that we never feed backwards...\n", mpirank, index, target_rank);
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if (target_rank >= mpisize) {
      fprintf(stderr, "ERROR: array bounds exceeded!  Requested access to trie_cell[%lld]\n", index);
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    remote_getread(target_rank, index, valuep);
  }
}
/*>*/

//------------------------------------------------------------------------

/*< local_*: these procedures all execute locally to the caller.  This
    is where the actual work is done.  However the users do not call
    these directly (except possibly in a few cases where some low-level
    optimisation is a necessity). <b><tt>local_add_read</tt></b> is the
    primary work-horse procedure that implement our algorithm. */

static int add_read(char *s, EDGE edge, long read_number, int len);
static INDEX get_next_free_edge(void);

// If we allow multi-threading, this procedure is where all the critical
// sections will be.  Any updates to cells must be done carefully as
// protected read/modify/writes.
static int local_add_read(char *s, EDGE edge, long read_number, int len)
{
  int c;

  assert((edge >> CHUNKBITS) == mpirank);
  assert((*s != '\0') && (*s != '\n') && (*s != '\r'));

//for (;;) {  // on a single node, can optimise tail-recusion by making this a loop...
    c = *s++;

    letters++; // for diag only - NOT SHARED - WON'T WORK IN DISTRIBUTED IMPLEMENTATION!
    freq[c]++; // letters and freq may need to be long longs

    if (c == 'A') c = _A_;
    else if (c == 'C') c = _C_;
    else if (c == 'G') c = _G_;
    else if (c == 'T') c = _T_;
    else c = _N_; // some other char

    if ((*s == '\0') || (*s == '\n') || (*s == '\r')) {
      if (trie_cell[edge&CHUNKMASK].edge[c]&ENDS_WORD) {
        long original_read = trie_cell[edge&CHUNKMASK].edge[c]&EDGE_MASK;

        // Take note of 100% overlaps (dups) as we find them, and don't enter
	// the second or later copies in the trie.
        fprintf(duplicates, "%ld:0 %ld\n", original_read, read_number);
        if (ferror(duplicates)) {
          fprintf(stderr, "\n\n************* add_read() (duplicates) failed, %s\n", strerror(errno));
          shut_down_other_nodes();
          MPI_Finalize();
          exit(EXIT_FAILURE);
	}
        dups++; // for diag only - NOT SHARED - WON'T WORK IN DISTRIBUTED IMPLEMENTATION!
      } else {
        assert((trie_cell[edge&CHUNKMASK].edge[c]&EDGE_MASK) == 0);
        trie_cell[edge&CHUNKMASK].edge[c] = (ENDS_WORD | read_number);
      }

      return len+1;
    }

    // WE NO LONGER SUPPORT STRINGS OF DIFFERING LENGTHS

    if (trie_cell[edge&CHUNKMASK].edge[c] == 0LL) { // we need to extend

      INDEX new_edge = get_next_free_edge(); // Might return an edge from another compute node! // TO DO: return an *EMPTY* edge

      setread(new_edge, empty); // slight overhead if local. optimise later.

      trie_cell[edge&CHUNKMASK].edge[c] = new_edge;
      if (new_edge >= MAX_SIZE) {
        fprintf(stderr, "Ran out of free edges after %d reads (last_used_edge = %lld, MAX_SIZE = %lld)\n", seq, new_edge, MAX_SIZE);
        shut_down_other_nodes();
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    } else { // we need to merge
      // ... with existing trie_cell[edge&CHUNKMASK].edge[c]
    }

    // This add_read call is actually the smart bit of the algorithm because it hands over to
    // another node only when it has to.  An alternative choice would have been to fetch/store
    // each individual trie cell, but the overhead of that was horrendous.

    // CALL NEXT PROCESSOR IF 'edge' IS TOO HIGH.
    // (use of tmp is just so we can look at it with diags if need be)

    return add_read(/* modified */ s, trie_cell[edge&CHUNKMASK].edge[c] /* may be on another node */, read_number, len+1);

    //edge = trie_cell[edge].edge[c]&EDGE_MASK;  // Local version can be optimised by removal of tail recursion...
//}                                              // and converted into a loop.
}

/*>*/

/*< remote_*: these procedures call a later rank.  These are user-application calls. */
static INDEX remote_get_next_free_edge(target_rank)
{
  MPI_Status status;
  long value = 0L; // generic up-front parameter
  INDEX free_edge;

  assert(target_rank != mpirank);

  // COMMAND CODE
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_GET_NEXT_FREE_EDGE,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // NO PARAMS

  assert(sizeof(INDEX) == sizeof(long long));
  MPI_Recv(&free_edge,/* message buffer - used for len function result */
	   1,/* one data item */
	   MPI_LONG_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */

  // Get acknowlegement that call is complete (to synchronise)
  MPI_Recv(&value,/* message buffer - unused */
	   1,/* one data item */
	   MPI_LONG,/* response is a long int for now (temp) */
	   MPI_ANY_SOURCE,/* receive from any sender */
	   MPI_ANY_TAG,/* receive any type of message */
	   MPI_COMM_WORLD,/* always use this */
	   &status);/* info about received message */

  return free_edge;
}

static int remote_add_read(long target_rank, char *s, long edge, long read_number, int len) // pass to another node
{
  MPI_Status status;
  long value = 0L; // generic up-front parameter
  long stringlength;

  assert(target_rank != mpirank);

  stringlength = strlen(s)+1;

  // COMMAND CODE
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_ADD_READ,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // SEND PARAMETERS - later optimise by sending one param via 'value' above
  MPI_Send(&stringlength,/* message buffer - address part */
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  send_bytes(target_rank, s, stringlength);

  assert(sizeof(EDGE) == sizeof(long));
  MPI_Send(&edge,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  assert(sizeof(read_number) == sizeof(long));
  MPI_Send(&read_number,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  assert(sizeof(len) == sizeof(int));
  MPI_Send(&len,
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

  // If len were not needed we could fire & forget, by removing the Recv above...
  return (int)value;
}

static void remote_output_read(long target_rank, char *s, EDGE readindex) // pass to another node
{
  MPI_Status status;
  long value = 0L; // generic up-front parameter
  long stringlength;

  assert(target_rank != mpirank);

  stringlength = strlen(s)+1;

  // COMMAND CODE
  MPI_Send(&value,/* message buffer - data part and trigger to perform operation */
	   1,/* one data item */
	   MPI_LONG,/* data value is a long integer */
	   target_rank,/* destination process rank */
	   TAG_OUTPUT_READ,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  // SEND PARAMETERS - later optimise by sending one param via 'value' above
  MPI_Send(&stringlength,
	   1,/* one data item */
	   MPI_LONG,/* index is a long long integer */
	   target_rank,/* destination process rank */
	   TAG_DATA,/* user chosen message tag */
	   MPI_COMM_WORLD);/* always use this */

  send_bytes(target_rank, s, stringlength);

  MPI_Send(&readindex,
	   1,/* one data item */
	   MPI_LONG_LONG,/* data value is a long integer */
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
/*>*/

/*< accept_*: these are part of the dispatcher on the receiving node.
   They accept the RPC from the earlier node, save its parameters into
   local variables, call the local procedure to perform the action,
   and finally they return any results or success codes back to the caller  */

// NOTE: In later releases there is great scope for optimisation here, both
// from parallelism & pipelining, and from taking a short-cut in returning
// results when there is more than one node in a chain of calls. (MPI hacks)

static void accept_get_next_free_edge(int caller)
{
  INDEX new_edge;
  new_edge = get_next_free_edge(); // get the result.  Will either be local or punted downstream.
  MPI_Send(&new_edge, 1, MPI_LONG_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege and return result
}

static void accept_add_read(int myrank, long value, MPI_Status status) // receive request from RPC mechanism
{
  char *s;
  EDGE edge;
  long read_number, stringlength;
  int len;
  int caller;

  // RECEIVE PARAMETERS (possible optional first parameter passed in as 'value')
  MPI_Recv(&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  caller = status.MPI_SOURCE;

  s = malloc(stringlength);
  if (s == NULL) fprintf(stderr, "CRAP!  Failed to malloc... why\?\?\?\n");
  MPI_Recv(s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  assert(sizeof(EDGE) == sizeof(long));
  MPI_Recv(&edge, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  assert(sizeof(read_number) == sizeof(long));
  MPI_Recv(&read_number, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  assert(sizeof(len) == sizeof(int));
  MPI_Recv(&len, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  {
    long long int target_rank = (long long int)edge >> CHUNKBITS;
    assert(target_rank == mpirank);
  }

  value = add_read(s, edge, read_number, len);
  free(s); s = NULL;

  MPI_Send(&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege and return result
}

static void walk_and_print_trie_internal(char *s, EDGE edge, int len);
static void dump_trie(char *filename);
static void accept_walk_and_print_trie_internal(int myrank, long value, MPI_Status status) // receive request from RPC mechanism
{
  char s[MAX_LINE];
  EDGE edge;
  long stringlength;
  int len;
  int caller;

  MPI_Recv(&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  caller = status.MPI_SOURCE;
  MPI_Recv(s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&edge, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&len, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  walk_and_print_trie_internal(s, edge, len);
  MPI_Send(&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege and return result
}

static void accept_dump_trie(int myrank, long value, MPI_Status status) // receive request from RPC mechanism
{
  char filename[MAX_LINE];
  long stringlength;
  int caller;

  MPI_Recv(&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  caller = status.MPI_SOURCE;
  MPI_Recv(filename, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  dump_trie(filename);
  MPI_Send(&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege and return result
}

static void output_read(char *s, EDGE readindex);
static void accept_output_read(int myrank, long value, MPI_Status status) // receive request from RPC mechanism
{
  char *s;
  EDGE readindex;
  long stringlength;
  int caller;

  // RECEIVE PARAMETERS (possible optional first parameter passed in as 'value')
  MPI_Recv(&stringlength, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  caller = status.MPI_SOURCE;

  s = malloc(stringlength); // TO DO: may be better off stack...
  if (s == NULL) fprintf(stderr, "CRAP!  Failed to malloc... why\?\?\?\n");
  MPI_Recv(s, stringlength, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&readindex, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  output_read(s, readindex);

  free(s); s = NULL;

  MPI_Send(&value, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege and return result
}
/*>*/

// ------------------------------------------------------------------------------

/*< High-level remote procedure calls.  These look like normal procedures
    which could be called by a user program.  Depending on where the data
    is stored that they are accesing, the procedure will either be called
    locally, or the RPC stub will punt the call to the appropriate node
    where the data can be found. */


static void output_read(char *s, EDGE readindex)
{
  time_t curtime;
  if (mpirank == mpisize-1) {
    static int printed = 0;
    fprintf(sorted_and_unique_reads, "%s %12lld\n", s, readindex);
    if (ferror(sorted_and_unique_reads)) {
      fprintf(stderr, "\n\n************* output_read() failed, %s\n", strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    printed++;
    if ((printed % 1000000) == 0) {
      time(&curtime); fprintf(stderr, "%d unique and sorted reads written back at %s", printed, ctime(&curtime));
    }
  } else {
    remote_output_read(mpisize-1, s, readindex);
  }
}

static INDEX get_next_free_edge(void)
{
  if (((last_used_edge+(INDEX)1) >> CHUNKBITS) != mpirank) { // this node is all out of edges...
    INDEX edge;
    static int next_guy = 1;
    static int init = FALSE;
    // Ask the next guy.  If he's full too, he'll ask the guy after that.
    if (!init) {
      // optimise so that we only ever make one call, rather than a chain of calls
      next_guy = mpirank+1; init = TRUE;
      if (mpirank == mpisize) {
        fprintf(stderr, "ERROR: not enough RAM for this input file (%d * %lld cells used).  Try resubmtting with some more processors.\n", mpisize, CHUNKSIZE);
      }
    }
    edge = remote_get_next_free_edge(next_guy);
    next_guy = (edge >> CHUNKBITS); // we remember the last guy in the chain
    return edge;
  } else return ++last_used_edge; // local
}

static int add_read(char *s, EDGE edge, long read_number, int len) // what the user calls...
{
  int len2;
  long long int target_rank = (long long)edge >> CHUNKBITS;

  if (len == 0) { // initial call, not a recursive call with a partial-string
    assert(edge == ROOT_CELL);
    if (read_number > EDGE_MASK) { // too large a read number to store in an EDGE.
      fprintf(stderr, "maketrie: too many READs! (%ld)  Limit is %lld\n", read_number, EDGE_MASK);
      assert(read_number <= EDGE_MASK);
    }
    seq++;
  }

  if (target_rank == mpirank) {
    // the caller could optimise by calling local_add_read directly if it knew it
    // was still within the local chunk.  However the overhead is low (especially
    // if we use the extension "static inline") and this is *much* safer.
    len2 = local_add_read(s, edge, read_number, len); // can optimise more knowing that trie_cell[ROOT_CELL] is always rank 0
    if (len == 0) length[len2]++;  // hacky, and ought to be removed some day.
    return len2;
  } else {
    return remote_add_read((long)target_rank, s, edge, read_number, len);
  }

}

// ------------------------------------------------------------------------------

// For diagnostics only.  Don't call this for mega-large tries!

static void walk_and_print_trie_internal(char *s, EDGE edge, int len)
{
  int i;
  int target_rank = edge >> CHUNKBITS;

  if (target_rank != mpirank) {
    remote_walk_and_print_trie_internal(target_rank, s, edge, len);
    return;
  }

  s[len+1] = '\0';
  for (i = 0; i < 5; i++) {
    s[len] = trt[i];
    if (trie_cell[edge & CHUNKMASK].edge[i]&ENDS_WORD) {
      // the final cell points to the sequence number from original 'read_number'
      // when it was entered and duplicate removal performed.  We'll need this
      // when reading the reads back in for overlap detection.
      output_read(s, trie_cell[edge & CHUNKMASK].edge[i]&EDGE_MASK);
    } else if (trie_cell[edge & CHUNKMASK].edge[i]) { // not final letter, and this letter is present with more to follow
      walk_and_print_trie_internal(s, trie_cell[edge & CHUNKMASK].edge[i], len+1);
    }
  }

}

static void dump_trie(char *filename)
{
  time_t curtime;
  int rc;
  fprintf(stderr, "maketrie[%d]: ", mpirank);
  if (mpirank == 0) {
    trie_file = fopen(filename, "w");
    fprintf(stderr, "Writing");
  } else {
    trie_file = fopen(filename, "a");
    fprintf(stderr, "Appending");
  }
  time(&curtime); fprintf(stderr, " to dumped trie %s at %s\n", filename, ctime(&curtime));
  if (trie_file == NULL) {
    fprintf(stderr, "maketrie[%d]: Cannot save trie to %s - %s\n", mpirank, filename, strerror(errno));
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  if (last_used_edge == (mpirank * CHUNKSIZE - 1)) return /* mpirank-1 */;  // this processor node has not been used
                                                   // TO DO: return this processor's rank back up to rank 0
  fwrite(trie_cell, last_used_edge+1 - (mpirank * CHUNKSIZE), sizeof(CELL), trie_file);
  if (ferror(trie_file)) {
    fprintf(stderr, "maketrie[%d]: Error saving trie to %s - %s\n", mpirank, filename, strerror(errno));
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  time(&curtime); fprintf(stderr, "Written at %s\n", ctime(&curtime));
  if (trie_file) {
    rc = fclose(trie_file); trie_file = NULL;
  if (rc == EOF) {
    fprintf(stderr, "maketrie[%d]: Error saving trie to %s - %s\n", mpirank, filename, strerror(errno));
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  }

  if (last_used_edge == ((mpirank+1) * CHUNKSIZE - 1)) {
    fprintf(stderr, "Not done.  Asking next rank to continue...\n");
    // This node is full so ask the next node to append its data to the trie file as well
    if (mpirank != mpisize-1) /* return */ remote_dump_trie(mpirank+1, filename);
  }
  /* return mpirank; */
}

static void walk_and_print_trie(void) {
  char s[MAX_LINE];
  time_t curtime;

  time(&curtime); fprintf(stderr, "Printing sorted reads at %s", ctime(&curtime));
  walk_and_print_trie_internal(s, ROOT_CELL, 0);
  if (mpirank == mpisize-1) {
    if (sorted_and_unique_reads) {
    int rc = fclose(/* output */sorted_and_unique_reads); sorted_and_unique_reads = NULL;
    if (rc == EOF) {
      fprintf(stderr, "maketrie[%d]: Error closing sorted output - %s\n", mpirank, strerror(errno));
    }
    }
  }
  time(&curtime); fprintf(stderr, "Printing sorted reads complete at %s", ctime(&curtime));
}
/*>*//*>*/

int main(int argc, char **argv)
{
  /*< (local declarations) */
  time_t curtime;
  char fname[1024];
  char line[MAX_LINE];
  int i, c, rc, lineno = 1, number_of_lengths = 0;
  long read_number = 0;
  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  /*>*/

  time(&curtime); fprintf(stderr, "Program started at %s", ctime(&curtime));
  /*< Set up MPI environment.  We also intend to use OMP later, although currently the use of it is trivial. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  fprintf(stderr, "I am rank %d of world size %d\n", mpirank, mpisize);
  MPI_Get_processor_name(processor_name, &namelen);
  if (processor_name && strchr(processor_name, '.')) *strchr(processor_name, '.') = '\0';
  /*>*/

  /*< Command-line parameter handling and file opening */
  if ((mpirank == 0) && (argc > 2)) {
    fprintf(stderr, "warning: extra parameter %s ignored...\n", argv[2]);
  }

  if (argc >= 2) { // command-line handling could be cleaner.

    if (mpirank == 0) {
      read_file = fopen(argv[1], "r");
      if (read_file == NULL) {
        fprintf(stderr, "maketrie: cannot open input \"%s\" - %s\n", argv[1], strerror(errno));
        shut_down_other_nodes();
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
      fprintf(stderr, "Input: %s\n", argv[1]);
    }

    assert(strlen(argv[1]) + 10 < 1024); // Sorry... I'm in a hurry...

    /*< Open all files. */
    // TO DO: use the (mpisize-1) trick to output the duplicates to a single file rather than a separate file per rank.
    // There may be a minimal extra cost in doing this but it will simplify file handling a lot.  May even
    // output the overlaps to the same file later.
    sprintf(fname, "%s-dups-%05d", argv[1], mpirank); // are 10,000 nodes enough for you?
    duplicates = fopen(fname, "w");
    if (duplicates == NULL) {
      fprintf(stderr, "maketrie on rank %d: cannot open output \"%s\" - %s\n", mpirank, fname, strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    fprintf(stderr, "Output: %s\n", fname);

    if (mpirank == mpisize-1) {
    sprintf(fname, "%s-sorted", argv[1]);
    sorted_and_unique_reads = fopen(fname, "w");
    if (sorted_and_unique_reads == NULL) {
      fprintf(stderr, "maketrie[%d]: cannot create sorted output \"%s\" - %s\n", mpirank, fname, strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    fprintf(stderr, "Sorted READ Output: %s\n", fname);

    sprintf(fname, "%s-rejects", argv[1]);
    rejects = fopen(fname, "w");
    if (rejects == NULL) {
      fprintf(stderr, "maketrie[%d]: cannot create reject file \"%s\" - %s\n", mpirank, fname, strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    fprintf(stderr, "Reject Output: %s\n", fname);
    }
    // Might be easier to gather and print stats as we output.  However there's currently no 
    // real need to output the sorted reads and I may remove this later.
    // Also should consider dumping the trie as a binary structure for later fast reloading...
    // (and to allow diagnostic utilities to print subtrees, e.g. to explain overlaps in contigs)

    if (mpirank == 0) {
    sprintf(fname, "%s-index", argv[1]);
    read_index = fopen(fname, "w");
    if (read_index == NULL) {
      fprintf(stderr, "maketrie: cannot create index \"%s\" - %s\n", fname, strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    fprintf(stderr, "READ Index Output: %s\n", fname);
    }
    /*>*/
  } else {
    if (mpirank == 0) fprintf(stderr, "syntax: maketrie input.fastq\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  /*>*/

  /*< Some OS-dependent (Linux) code to let us grab as much memory as possible in a single chunk. */

  // DO NOT throw this away and attempt to build the trie using memory claimed from malloc.
  // Sequential memory allocation is a high-level optimisation and critical to the speed of the program

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
      fprintf(stderr, "Discovered %lld processors per node\n", PROCESSORS_PER_NODE);
    }
    if (PROCESSORS_PER_NODE == 0ULL) PROCESSORS_PER_NODE = CORES_PER_NODE;

    // At least at TACC, omp_get_max_threads() defaults to 1 unless you
    // explicitly set OMP_MAX_THREADS in the job script
    // which we do only for the case where we assign one task per node.
    // This assumption probably won't hold elsewhere and also we don't
    // yet know the number of cores per node... can't seem to find it
    // from any standard MPI or OMP call... currently hard coded...
    // (omp_get_num_procs doesn't seem to do it...)

    TASKS_PER_NODE = PROCESSORS_PER_NODE/(long long)omp_get_max_threads();
    fprintf(stderr, "Node %s, Rank %d, and running %lld ranks on this node.   <-------------------------------\n",
           processor_name, mpirank, TASKS_PER_NODE);

    meminfo = fopen("/proc/meminfo", "r");
    if (meminfo) for (;;) {
      int count;
      char *s = fgets(line, 1023 /*MAX_LINE*/, meminfo);
      if (s == NULL) break;
      //fprintf(stderr, "checking %s", line);
      count = sscanf(line, "MemTotal:     %lld kB", &memsize);
      if (count == 1) {
        //fprintf(stderr, "got memsize %lld kB\n", memsize);
        memsize *= 1024ULL; // re-scale to bytes
        memsize /= TASKS_PER_NODE; // availability per processor
        memsize /= (long long)sizeof(CELL); // convert to cell count

        CHUNKBITS = 1ULL;
        for (;;) {
          if ((1ULL << CHUNKBITS) >= memsize) break;
          CHUNKBITS++;
	}
        CHUNKBITS -= 1ULL;
        fprintf(stderr, "rounding down memsize to %lldM cells per core (%lld bits), ie %lldM cells per node\n",
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
    fprintf(stderr, "Node %d: trying calloc of %lld cells of %d bytes each.\n", mpirank, CHUNKSIZE, (int)sizeof(CELL));
    trie_cell = calloc(CHUNKSIZE, sizeof(CELL)); // trie_cell[0..CHUNKSIZE]
    if (trie_cell == NULL) {
      CHUNKBITS -= 1ULL;
      CHUNKSIZE = (1ULL<<CHUNKBITS);
      CHUNKMASK = (CHUNKSIZE-1ULL);
    } else break; // malloc successful
  }

#ifdef MULTINODE_DEBUG100
  // we want to test on 4 nodes because of a bug that seemed to trigger when add_read
  // had to access the third node of 4.  use with 100 letter test genome (fake100)
  free(trie_cell);
  CHUNKBITS = 8;
  CHUNKSIZE = (1ULL<<CHUNKBITS);
  CHUNKMASK = (CHUNKSIZE-1ULL);
  trie_cell = calloc(CHUNKSIZE, sizeof(CELL)); // trie_cell[0..CHUNKSIZE]
#endif

#ifdef MULTINODE_DEBUG1K
  // we want to test on 4 nodes because of a bug that seemed to trigger when add_read
  // had to access the third node of 4.  use with 100 letter test genome (fake100)
  free(trie_cell);
  CHUNKBITS = 9;
  CHUNKSIZE = (1ULL<<CHUNKBITS);
  CHUNKMASK = (CHUNKSIZE-1ULL);
  trie_cell = calloc(CHUNKSIZE, sizeof(CELL)); // trie_cell[0..CHUNKSIZE]
#endif

  if (trie_cell == NULL) {
    fprintf(stderr,
            "rpctest: rank %d unable to allocate array of %lld longs\n", mpirank, CHUNKSIZE);
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

  MAX_SIZE = CHUNKSIZE*mpisize; // no of elements in total virtual read array
  fprintf(stderr, "setting MAX_SIZE to %lld (%lld * %d)\n",  MAX_SIZE, CHUNKSIZE, mpisize);

  /*>*/

  /*< Housekeeping. */
  // Originally I just did these on Node 0 but realise now that they're really global.

  // these are used for information only.  They need work as they are not YET shared, and
  // some of the data is coming from ranks 1..n

  for (i = 0; i < 256; i++) freq[i] = 0;
  for (i = 0; i < MAX_LINE; i++) length[i] = 0;
  for (i = 0; i < 5; i++) empty.edge[i] = 0LL; // Need to do away with 'empty' and in fact rely on calloc so that get_next_free_cell does no init

  // INITIALISE ROOT_CELL.  This is why last_used_edge is initialised to ROOT_CELL and not 0.
  // this could be replaced by ROOT_CELL = get_next_free_cell() ????
  for (i = 0; i < 5; i++) trie_cell[ROOT_CELL].edge[i] = 0LL;
#ifdef TWONODE_DEBUG
    last_used_edge = CHUNKSIZE-100ULL; // debugging trick to force us to go over an edge boundary and use another node to store some of the data...
#endif
  fprintf(stderr, "last_used_edge: %lld,  CHUNKSIZE: %lld,  MAX_SIZE: %lld\n", last_used_edge, CHUNKSIZE, MAX_SIZE);
  /*>*/

  if (mpirank == 0) {
    // --------------------------- MAIN BODY OF CODE ON PRIMARY PROCESSOR ----------------------------------
    fprintf(stderr,
            "\nCombined system is using %lldM trie edges distributed across %d ranks\n\n",
            (long long)mpisize * (CHUNKSIZE >> 24ULL), mpisize
           );


    /*< PHASE ONE: BUILD A DATA STRUCTURE USING ALL READS */

  // The input handling pretty much assumes input is well-formed, for speed of reading.
  // (although could be faster with some I/O tricks.)  Use of a program like fastqc should
  // enforce clean input?

  // Read data in fastq format:
  for (;;) {
    int len;
    off_t read_start;
    char *s;

    // Build an index to map read_number to the place in the file where the text for that read can be found.
    // size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);

    read_start = ftello(read_file); if (read_index) fwrite(&read_start, sizeof(off_t), 1, read_index); // indexed by read_number

    // @SRR022868.923/1
    s = fgets(line, MAX_LINE, read_file);
    if (s == NULL) break;
    // ATCTACTTACTGGAAGTTTAATTTGAGTAAATTGTTATCCAGTCATTCGTTAGAACTCCTTATAGTACTTATACCNNNNNNNNNNNNNNNNNNNNNNNNNN
    lineno++; fgets(line, MAX_LINE, read_file);
    s = strchr(line, '\n'); if (s) *s = '\0'; // trim trailing newline
    len = add_read(line, ROOT_CELL, read_number++, 0);

    // we could, if we wanted to, have a second ROOT_NODE (eg ROOT_NODE_REVERSED) which
    // we could use to create a similar trie, but first reversing the input line.  This
    // would let us find overlaps in both directions which would simplify finding the
    // start of a contig.

    // +
    lineno++; fgets(line, MAX_LINE, read_file);
    if (line[0] != '+') {
      fprintf(stderr, "Input data format error in READ file line %d\n", lineno);
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    // II*II?IIIIIIIIII8;III2I@I;I4E+2'I?'F>>$.II.&0@+)I:,6(,2#+*#+#*.'('''+$54%$*!!!!!!!!!!!!!!!!!!!!!!!!!!
    lineno++; fgets(line, MAX_LINE, read_file);
    lineno++;
    if ((read_number % 1000000) == 0) {
      time_t curtime;
      time(&curtime); fprintf(stderr, "%ld READs loaded at %s", read_number, ctime(&curtime));
    }
    if (read_number == 0x7FFFFFFF) {
      fprintf(stderr, "maketrie: an assumption was wrong.  We have an input file with more than %d READs.  Code fix needed.\n", 0x7FFFFFFF);
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE); // need to make the number of READs a long long.
    }
  }
  if (duplicates) {
  rc = fclose(/* output */ duplicates); duplicates = NULL;
  if (rc == EOF) {
    fprintf(stderr, "maketrie: error closing %s-duplicates - %s\n", argv[1], strerror(errno));
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  }

  if (read_file) {
  rc = fclose(/* input */ read_file); read_file = NULL;
  if (rc == EOF) {
    fprintf(stderr, "maketrie: error closing %s - %s\n", argv[1], strerror(errno));
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  }

  if (read_index) {
  rc = fclose(/* output */ read_index); read_index = NULL;
  if (rc == EOF) {
    fprintf(stderr, "maketrie: error closing %s-index - %s\n", argv[1], strerror(errno));
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  }

  fprintf(stderr, "\nread trie built using %lld nodes (%0.0f%% of capacity)\n",
          last_used_edge, 100.0 * last_used_edge / MAX_SIZE); /* TO DO: THIS IS WRONG - last_used_edge only applies to rank 0, not last used rank. */
  fprintf(stderr, "\nTotal of %d reads indexed and sorted, including %d (%0.0f%%) duplicates (dup count is temporarily inaccurate when using multiple nodes)\n",
          seq, dups, dups * 100.0 / seq);
  fprintf(stderr, "\nFrequencies:\n");
  for (c = 0; c < 256; c++) if (freq[c]) fprintf(stderr, "   %c  %ld\n", c, freq[c]);
  for (i = 0; i < MAX_LINE; i++) {
    if (length[i]) {
      number_of_lengths++; read_length = i;
    }
  }
  fprintf(stderr, "\n");
  if (number_of_lengths == 0) {
    fprintf(stderr, "Error: No READs found!  Bad input file?\n");
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else if (number_of_lengths != 1) {
    fprintf(stderr, "Error: this code does not handle READs of differing lengths\n");
    fprintf(stderr, "\nWe found READs of lengths:\n");
    for (i = 0; i < MAX_LINE; i++) {
      if (length[i]) {
        fprintf(stderr, "     %d  (%ld)", i, length[i]);
      }
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "\nPlease clean the data first with a program like 'fastqc'.\n\n");
    shut_down_other_nodes();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    fprintf(stderr, "READ length: %d\n", read_length);
  }
  fprintf(stderr, "\n");
  fflush(stderr);
  /*>*/

    walk_and_print_trie();
    sprintf(fname, "%s-edges", argv[1]);
    dump_trie(fname);

    // having built our key data structure, now we find the overlaps, using an algorithm
    // for parallel string searching roughly based on
    // http://www.gtoal.com/wordgames/bigspell/multiscan.c.html

    if (rejects) {
    rc = fclose(/* output */ rejects); rejects = NULL;
    if (rc == EOF) {
      fprintf(stderr, "maketrie: cannot close reject file \"%s\" - %s\n", argv[1], strerror(errno));
      shut_down_other_nodes();
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    }

    //fprintf(stderr, "Asking other nodes to shut down...\n");
    time(&curtime); fprintf(stderr, "Program complete at %s", ctime(&curtime));
    shut_down_other_nodes();

    // or we could go on from here while the trie is in memory to use the trie to
    // stitch the overlapping fragments into a complete genome. :-)

    free(trie_cell); trie_cell = NULL;


  } else {
    /*< ALL OTHER PROCESSORS RUN AN ACCEPT/DISPATCH LOOP FOR REMOTE PROCEDURE CALLS... */
    // (currently single-threaded.  An MPI guru could write this better.)

    // Also it is VERY IMPORTANT to note that Node 0 does *not* run a copy of
    // this dispatcher.  (because it is single-threaded too.)  Would be nice if
    // that option were available, eg to allow callbacks.

    // One optimisation trick planned for the far future...  if A calls B and B calls C,
    // then C could pass the result back directly to A, allowing B to return to its
    // displatch loop and accept more commands.  Which will be fine if it's always A
    // that initiates certain types of call... - just need to be careful with
    // specifying 'MPI_ANY_SOURCE' rather than 'caller'.

    INDEX index;
    CELL read_value;
    long longvalue;
    MPI_Status status;
    //    int debug = 0, rdebug = 0, wdebug = 0;

    for (i = 0; i < 5; i++) empty.edge[i] = 0LL; // Has to be initialised in these instances as well...
    last_used_edge = mpirank * CHUNKSIZE - 1; // IMPORTANT TO INITIALISE TO *virtual* INDEX

    for (;;) {
      int caller;

      // Accept an RPC request.  Get the command code.
      MPI_Recv(&longvalue, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      caller = status.MPI_SOURCE;

      /*
       * Check the tag of the received message and dispatch accordingly
       */
      if (status.MPI_TAG == TAG_READ_READ) {          // RETURN AN ELEMENT
        MPI_Recv(&index, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        getread(index, &read_value);
        MPI_Send(&read_value.edge, sizeof(read_value.edge)/sizeof(read_value.edge[0]), MPI_LONG, caller, 0, MPI_COMM_WORLD); // return the value

      } else if (status.MPI_TAG == TAG_WRITE_READ) {  // WRITE TO AN ELEMENT
        MPI_Recv(&index, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        assert(sizeof(read_value.edge[0]) == sizeof(long));
        MPI_Recv(&read_value.edge, sizeof(read_value.edge)/sizeof(read_value.edge[0]), MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        setread(index, read_value);
        MPI_Send(&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege

      } else if (status.MPI_TAG == TAG_ADD_READ) {  // CALL add_read()
	accept_add_read(mpirank, longvalue, status);

      } else if (status.MPI_TAG == TAG_GET_NEXT_FREE_EDGE) {  // CALL get_next_free_edge()
        accept_get_next_free_edge(caller);
        MPI_Send(&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege

      } else if (status.MPI_TAG == TAG_OUTPUT_READ) {  // CALL output_read()
	accept_output_read(mpirank, longvalue, status);

      } else if (status.MPI_TAG == TAG_WALK_AND_PRINT_TRIE_INTERNAL) {  // CALL walk_and_print_trie_internal();
	accept_walk_and_print_trie_internal(mpirank, longvalue, status);

      } else if (status.MPI_TAG == TAG_DUMP_TRIE) {  // CALL dump_trie_internal
        accept_dump_trie(mpirank, longvalue, status);

      } else if (status.MPI_TAG == TAG_EXIT_PROGRAM) {  // Kill clients
        fprintf(stderr, "Node %d asked to exit\n", mpirank);
        MPI_Send(&longvalue, 1, MPI_LONG, caller, 0, MPI_COMM_WORLD); // Acknowlege
        fprintf(stderr, "Node %d exit acknowleged\n", mpirank);
        break;

      } else {
        // UNKNOWN - CODING ERROR?
      }
    }
    fprintf(stderr, "Node %d exiting cleanly.  local base = %lld,  last_used_edge = %lld,  local maximum = %lld\n",
                     mpirank, mpirank * CHUNKSIZE, last_used_edge, (mpirank+1) * CHUNKSIZE);
    if (duplicates) {
    rc = fclose(/* output */ duplicates); duplicates = NULL;
    if (rc == EOF) {
      fprintf(stderr, "maketrie: error closing %s-duplicates - %s\n", argv[1], strerror(errno));
    }
    }

    if (trie_cell != NULL) free(trie_cell); trie_cell = NULL;
    /*>*/

  }

  MPI_Finalize();

  exit(EXIT_SUCCESS);
  return EXIT_FAILURE;
}
