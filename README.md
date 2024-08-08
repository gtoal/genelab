# genelab
A fast de novo Overlap-Layout-Consensus DNA assembler in C. Single machine or cluster scalable.

# The short README
Download the .c files and the Makefile.<br/>
Install an MPI suite if not done already: <tt>sudo apt install libmpich-dev</tt> may be sufficient.<br/>
Type <tt>make</tt><br/>
Experiment with your .fastq files of shotgun DNA reads. Commands are <a href="https://github.com/gtoal/genelab/blob/main/README.md#the-software">below</a>.<br/>
<hr>

# The LONG README
## My <em>de novo</em> DNA assembler
This is a <em>de novo</em> "Overlap-Layout-Consensus" DNA assembler which I wrote some years ago, but didn't put online at the
time.  It had some significant advantages over what was available in 2013 (mostly De Bruijn graph assemblers) - for example,
we were actually able to assemble our Bacillus thuringiensis sequence on one single Windows PC with this software, even though we had
access to a supercomputer.

I've not kept up with this field, having had the luxury of being able to retire early,
soon after this project completed, so it may be that its advantages are now a standard part of
current software.  I really don't know.  But I've decided to upload it to github just
in case it's still useful to anyone.  By the way, I don't think I ever gave it a cutesy name - I found
my old source code in a directory named "genelab", so I guess I'll just call it that when I upload
it to github.

Here's what it did and how it works: first, shotgun sequencing breaks your DNA down into
short segments called &lsquo;reads&rsquo; - I forget the exact number of pairs - somewhere between 70
and 80 as far as I remember with the system we used.  Just to make talking about it easier
I'll use 70 as the example length in this description.  These reads are also called 'k-mers'
and in our case, we had some data that came from a system where k=37 and others where k=70.

And before I start, I'ld like to take two quick linguistic asides for folks who like myself were not already
steeped in the terminology of Bioinformatics (or Computer Science, if you <em>are</em> a bioinformaticist)...
Firstly: the sequence output by the sequencing machine
is called a &lsquo;read&rsquo; - that's pronounced the same as &lsquo;reed&rsquo; as in a saxophone or clarinet (&lsquo;ɹiːd&rsquo;
in IPA notation), not like &lsquo;red&rsquo; as in &lsquo;I read a good book last week&rsquo;. &lsquo;read&rsquo; the
noun and &lsquo;read&rsquo; the present tense verb are both pronounced the same.  Past-tense
&lsquo;read&rsquo; (red, or &lsquo;rɛd&rsquo; in IPA notation) is not likely to be used much in this document,
but you will have to read carefully at times to distinguish the noun from the verb!;
and secondly: I refer to a data structure used in Computer Science as a &lsquo;trie&rsquo;, with the
word coming from the middle of the expression &lsquo;data re<em>trie</em>val structure&rsquo;!  Yet for some
reason, this is properly pronounced as &lsquo;try&rsquo; - but no-one will fault you if you prefer
to hear &lsquo;tree&rsquo; in your head as you read it, especially since it is actually a form of a
tree and indeed I will probably refer to it as such below except in a few cases where it
may be instructive to draw a distinction.

I believe nowadays the &lsquo;short&rsquo; reads can be a lot longer - it doesn't matter to this
software - the max size is arbitrary.  It doesn't need to be as short as 70 in order
to run efficiently.

At the time I wrote this the majority of DNA assembly code used De Bruijn graphs to
join reads into longer contiguous sequences.  This struck me as a very strange and
inefficient algorithm+data structure to use - whereas to me it was blatantly obvious
that the data structure to use should be a trie.  Of course it was only blatantly
obvious to me because I had been using that data structure for years in the fields
of word games (Scrabble, Boggle etc) and spelling checking and correction, and the
functionality required in the DNA field was basically the same as in spelling correction.

So I built a &lsquo;trie&rsquo; from these sequences - a very simple data structure where the root
node of the data structure has 4 pointers - one for each base letter, and that node
points to the node for the second base in a 70-base read.  So it is a 4-way branching
tree which can be 70 items deep.

The code to build this tree is deliberately simple!  Following Tony Hoare's maxim from
the 60's, &ldquo;There are two methods in software design.  One is to make the program so simple,
there are obviously no errors.  The other is to make it so complicated, there are no
obvious errors.&rdquo;, I designed the software to be so simple that it should not have mistakes
and should be maintainable by any programmer.

During the project development I did discover one or two less well-known open source DNA
assemblers which did indeed use a tree data structure rather than De Bruijn graphs, and
we did try their code before completing our own, but even the best of those programs
had serious problems with running out of memory and were too slow to be useful.

So with the trie data structure, all the heavy computation was done up front in the
conversion from millions of flat 70-character reads into a 4-way 70-deep tree. (By the way
the final node at the end of a sequence of 70 bases would be a count of the number of times
that specific 70-character read was found - information which comes in useful later)

So the basic algorithm can be executed quite easily on a desktop PC, and if your organism
is quite small and your PC has a reasonably large RAM, you can do the entire job on a
single PC.

However if you are assembling say a human genome, one Windows desktop doesn't quite cut it.

BUT... the tree-building algorithm is what we computer folks call &lsquo;embarrasingly parallel&rsquo;!

Its simple to describe: Take 4 PC's and send all the reads whose first pair is &lsquo;C&rsquo; to one of
those PCs, all the reads starting with &lsquo;G&rsquo; to the next one, and so on.  Then very trivially
(with no runtime overhead other than copying the quarter-trees back to a single machine) you
can join the 4 sub-trees together into a single tree which is structured <em>identically</em> to
what you would have built on a larger system in one run.

And I just used a factor of 4 as an easy way to explain the parallel speed-up.  In fact by
taking more base pairs at a time, you could farm the building out to any arbitrary number
of computers you have available, eg reads which start with &lsquo;CC&rsquo; go to computer #1, with &lsquo;CG&rsquo;
to computer #2 and so on, for a distribution among 16 systems. Or 64. Or 256. Etc.

I should explain though that the primary bottleneck on the individual systems is reading
and writing the data (disk or network limited) and the amount of RAM available to hold
the sub-tree being built.  It is <em>not</em> the CPU power, so although each system may have
several (4 or 8 or more) CPU cores, we really only need one.  Running two cores on one
machine doesn't make a noticable difference - having twice as much RAM does.

So... when I ported the desktop version to our cluster in order to handle larger genomes,
I had a bit of a revelation in the middle of the developing the code which caused me to
structure the parallelism a completely different way!  I don't expect any other DNA
assemblers are likely to work this way.  Here's the story:

The 'traditional' way to handle splitting this large problem into smaller chunks
would be to run the same algorithm on all the machines, and combine the results
at the end.  This is quite reasonable but involves multiple machines all reading
the big input data file and just keeping the items that are destined for that
specific machine.  This has a high I/O overhead from all these systems reading
the same data off a shared disk.

The thing about this data structure is that the critical path is the amount
of RAM available.  Doubling the RAM on one computer gives a huge improvement
compared to doubling the CPU power available.  My first improvement over the
traditional way to break down the problem was to have the extra systems in the
cluster merely act as second-level RAM storage for the primary machine - effectively
an explicit virtual memory that was faster to store in and retrieve from than
disk.  Each system would be given a virtual base address for its share of
the combined memory, and it would be those disjoint virtual addresses that
would be stored in the tree nodes rather than real RAM addresses. (Another way
to look at a pointer in a node would be as a combination of a machine ID number
and a RAM offset within that machine's memory.  When saved to disk the RAM
offset would then become a file offset relative to the start of that particular
machine's 'chunk' in the combined tree file.)

This worked but it was while implementing this that I found the perfect solution:
A single read is added to the tree one letter at a time, as you step through
the read one letter at a time, you add to the tree one node at a time.  In the
previous scheme, looking at a node and adding to a node was a remote call to
another machine (once the available memory in the primary machine was full),
but that was when I realised... the initial computer doesn't have to work out
how to build the subsequent nodes once its memory is full - it can simply
take the remaining unprocessed part of the read and pass it on to the next
free machine!  This turns out to be beautifully elegant, because if you cut the
tree building off at a fixed tree depth, all nodes below that depth will be
hanging off the same node, so no information needs to be passed between computers
in order to assemble the final tree from nodes scattered all over the place -
all you need to do is concatenate the entire memory block from each system into
one file in the correct order.  Fast and elegant.

So, to summarize, since the description above was a bit rambling, ...

The 'root' of the tree handles every read sequentially.  It stores the first
&lt;N&gt; characters of each read in its own memory, and when it hits the &lt;N&gt;+1<sup>th</sup>
character in that read, it passes the entire remaining string on to another
processor.  That other processor runs identical software, except instead of the
root of tree that it builds representing the full 70-letter sequences, it represents
a smaller tree for smaller string lengths, and when the trees from all the systems
are amalgamated, that smaller tree fits in to the master tree at exactly the place
where the node that caused the switch to another processor happened.

The recursive call to insert the remainder of a read in the memory of the
other system could just have easily run on the master system, with the data
being written to the remote system via a remote procedure call, but by letting the
remote system actually perform the tree insertion itself directly, the communication
cost between systems is slashed to almost nothing!

In this way we take advantage of all the extra memory of the multiple systems
in a supercomputer cluster.  True, we don't take much advantage of CPU parallelism -
basically only one core in the whole cluster is doing useful work at any one time -
but that doesn't matter.  The problem is RAM-limited, not CPU-limited.

(That said, we initially did a little opportunistic CPU parallelism (using the
"OMP" compiler options) where we could, in maketrie and find-overlap, but it never
gave one of those 'order of magnitude' wins that the overall system design achieved,
and risked making the code less obvious to understand, so it was only done when
the parallelism did not add to the complexity of the code and actually I think it
may have subsequently been removed. It's the sort of optimisation you put back in only
once development is complete and users start asking for a little more speed)

This algorithm makes it possible to sequence huge (human-sized) genomes extremely
quickly, depending on how many cluster machines and consequently how much RAM you
have available to you.


The code in this repository implements the multiprocessing above by using the MPI primitives
common to supercomputer clusters - but really the communications part of the final version
of the algorithm is so simple it could be implemented using regular individual Linux computers,
with the data being passed on to child systems using RPC calls or even something as simple as an ssh command-line!

Note that having one root system, which is the one which sequentially works through the file
of reads one by one, means that the disk I/O is minimised and you don't need one of those fancy distributed
filing systems that clusters usually have.  If your lab has Linux PCs with 64Gb ram, you
could handle a human genome with a few of those PCs in an acceptable time.  By the way
it's relatively easy to turn a few standard Linux systems into an MPI cluster.  The code
here could be used with, say, 32 Pi 4 systems with 8Gb each to create a DNA assembly
engine with 256Gb of RAM.  And if that isn't enough, you can also break the problem
down sequentially - build half the tree in one run, then the other half in a second
run, and again just trivially join the trees together.  You would build the trees for
every read that started with C or G in the first run, and for every read that started
with A or T in the second run, and combining the two runs consists simply of creating
a single node for C,G,A,and T strings which points to the subtrees starting with the
second letter of each string.

There is an implementation detail I need to mention here which is <em>crucial</em>.
The memory space allocated for each node is taken from a simple flat array of nodes.
The links between nodes are simple integer indexes into this array rather than
memory addresses.  It's the simplest data structure you can imagine for storing a
tree, and by taking the nodes from a linear array, you can save the entire tree
by doing a single <em>write</em> operation on the whole array!  You must *not*
use 'malloc' to allocate these nodes as then writing them out to file from being
scattered all over memory becomes a tremendously expensive operation.  The use of
integer indices rather than actual memory pointers removes the need to relocate
the pointers when you write the tree out to disk.  This again is a big speed win.

And now, once the tree data structure is built, the fun begins!

The basic algorithm is simply to take one 70-char read, and find at what offset any other
70-char read best aligns to it - with the overlap <em>usually</em> being exact, but also allowing
up to a specified number of mis-matches to allow for errors in the reads.  Reads where the
overlap is too short (I think we ended up using a minimum overlap of 14 bases, determined by
experience rather than any theory) were not used in the assembly of longer contiguous chains.

In terms of algorithmic complexity, matching N strings to N strings - for millions of strings -
sounds like a monumentally large &lsquo;Order (N<sup>2</sup>)&rsquo; computation, and with the wrong data structure
it would be, but with a Trie, comparing one string against many is an <em>extremely</em> cheap computation - the same cost as
one single string compare.  Admittedly we do have to do that string compare roughly 70 times, with
the comparion string being slid aross the known strings one position at a time.
Done naively, the determination of all overlaps is
N<sup>2</sup> where N is the number of reads.  But with the trie data structure, the algorithmic complexity
is linear in N and proportional to L where L is the length of a read. (70, in our example, although because we reject overlaps
of less than 14 characters, really L=56).

It may not be obvious how N<sup>2</sup> complexity turns into N complexity... the answer is because the N<sup>2</sup>
is really N<sup>2</sup> in space <em>and</em> time, not just time.  And by using the trie data structure
we are effectively using N space and N time, ... so it's still N<sup>2</sup> in space and time, but
we only care about the N in time part.

The result of all this is that the hard work is done up front in the 'make trie' phase,
and assembling contiguous areas done on demand with a constant time algorithm in the 'find overlaps' code,
with the result that contig creation is <em>fast</em>.  <em>Really fast</em>.  And you can
do other valuable things with the now pre-computed data structure as well.  Say you are
looking for a specific sequence in the soup of DNA reads.  You don't have to build up all
the contigous sequences first in order to look for a match - you can simply take a short part of
your target sequence that you want to match against (say from the middle) and use that to build
sequences on the fly using the data
in the tree to <em>immediately</em> determine if the target sequence exists within your data.

The larger sequences which a DNA assembler creates from the reads are called &lsquo;contigs&rsquo;,
short for &lsquo;contiguous sequences.&rsquo;  These can be thousands of pairs long, <em>but</em> they are
only possible to determine when nearby sequences overlap unambiguously.
<p>You need to be aware that there is an
intrinsic problem that no software can solve, when either
A: there is an identical 70-character sequence at two or more places within the genome.
If this happens there is no way to know which sequence that follows one of these duplicated areas goes with which
sequence that precedes it; the answer can be had from an external database of similar organisms, but if this is
the first time such an organism has been sequenced, there's really no way to tell from that data.
Or, B: there is effectively random data which appears between contigs that is not consistent enough
to determine that two given contigs straddle this gap.  Sorry, I'm a computer programmer,
not a biologist, so I don't know (or rather I was told, but I don't remember) the proper
words to describe this issue.  But I'm sure you biologists will know what I'm referring to.

So... this software alone cannot guarantee to reconstruct an entire genome, but since the length of the raw sequences
that this code can process is arbitrary, you can augment your short reads with longer reads from a
different system, and this code can stitch those together with the short reads, straddling the areas on each end
of some 70-char sequence which occurred twice, by using a longer sequence that doesn't repeat,
containing that 70-letter sequence as part of it.

At the time I was working on this, sequences of say 300 were available using a different technology,
but they had a far greater error rate (where an error is a mis-identified letter) than the short
reads we usually used.


So within the limitations of the read length of your technology, this DNA assembly code does do an amazing
job of building those reads into the minimum number of large contigs.  These can be viewed with many well-developed
programs on a PC.  As I said, it was many years ago when I worked on this - I vaguely
recollect the names AMOS, BLAST, and <a href="https://ics.hutton.ac.uk/tablet/">Tablet</a>? - I'm sure other programs have come along since then
to replace them.

We wrote this for one particular project, and my collaborator - the
biologist who actually understands this stuff and explained as much of it as he could
to allow me to write the software - did not feel that the software was worth sharing at
the time due to one factor: it was written by a programmer for a programmer to use for
one specific project (<a href="https://pubmed.ncbi.nlm.nih.gov/26781916/"><em>Cry-like
genes, in an uncommon gene configuration, produce a crystal that localizes within the
exosporium when expressed in an acrystalliferous strain of Bacillus thuringiensis</em></a>)
&mdash; and so it did not have a GUI and wasn't especially accessible by biologists, lacking all of
the minor tweaks that biologists expect - it was just a basic assembly engine that did
one thing, but did it quickly and well, but was driven by a programmer working from the Linux command-line, hand-in-hand
with a biologist telling him what he needed.  He believed (probably rightly) that to be useful, software like this had to be usable
by biologists without any programmers on hand to run it for them.  My own feelings on that matter
are that either the biologists learn to code, or some other programmer with more interest
in user interfaces can take this code and make it idiot-proof.  Assuming the rest of the
world hasn't caught up with its performance yet, it'll be worth the effort.  If the current
equivalent software now outperforms this, well, look on it as an interesting historical
oddity :-)

So, in case the algorithms and data are still worth exploring some 10 years later, I've
decided to just release the code as it was at the end of our project, and let other programmers
experiment with it.  If it is still as advantageous over equivalent software as it was
in 2013, you're entirely free to build on it and add any bells and whistles that would
be necessary to gain acceptance by biologists.  I'm pretty sure the underlying framework
is sound.  When you read the code, don't groan about how obvious and simple some of it
is - it's meant to be that way.

I'm not planning to provide any support for this code, but don't let that hold you back
from asking questions if you really have to.  I'll answer those that I can.  I may have
a few other related utilities on an old drive somewhere.  The code here is primarily just
the basic engine.

As far as I remember, the software understands .afg and .fastq formats. Maybe others. It was a long time ago.
All those formats are quite simple and I'm sure there are conversion tools for all of them nowadays.
And having just had a quick scan through the source code, I've remembered that the trie-builder
also creates an index on the fly that allows you to present any 70-letter read as a string
and it will find the location of that read in the initial .fastq text file.  It must have
been useful for something but I'm damned if I can remember any more!  Although that index
can be quite large (up to 24Gb with the limitations of linux in 2013) looking up a read
(either an exact match or a match allowing for a small number of errors) was very fast,
because the index was memory mapped into the Linux address space and the disk was only
read from as each node in a 70-level deep tree was accessed.  Due to the order of how
these nodes were stored, there were maybe only 2 or three disk accesses per read.  The
database wasn't structured <em>explicitly</em> as indexed-sequential, but the creation
of the index was done in such a way that it behaved as if it were :-)

Having read over the source code again a few other details are coming back to me: I
mentioned that this was just the basic overlap/contig-finding code and didn't have
a lot of bells &amp; whistles... I remember now that we used another assembler package
called AMOS for those - it *did* have all the bells and whistles that my bioinformaticist
colleague wanted, so we would sometimes use parts of that suite in conjunction with
ours - I remember now that I changed some of my internal file formats to use the
ones compatible with AMOS, such as the descriptions of where overlaps occurred.  At the
time when I wrote this, our code outperformed AMOS by a huge margin, but AMOS was relatively
new code in 2013.  I can easily imagine that it might outperform us 10 years later.

A small detail: I said we used a 4-way branching tree to store the data.  In truth we
use a 5-way branching tree.  There's support for a 'N' entry in the k-mers, where the
correct value could not be determined.  We don't use it yet but I allowed for it in
the data to future-proof the software.  True, we do waste 25% of RAM that way, but
25% isn't enough to worry about. So our nodes consist of 5 64-bit integers instead of 4.

## The software  
Here are the programs we were using when our project was complete and we mothballed
the sources.  There is some earlier code and a few ancilliary utilities in the <a href="https://gtoal.com/genelab/old/old/">old</a>
subdirectory on my web server. The links below take you to prettified html versions of the code - the
raw C code is downloadable from the file listing at the top of this page (currently http://gtoal.com/genelab/ )
<ul>
<li><a href="http://gtoal.com/genelab/.html/maketrie.c.html">maketrie</a>: This builds the trie data structure and and index back into the file of raw k-mer reads.  By pre-building these data structures, all subsequent lookups (eg for contig building) are fast.<br/>
<a href="https://gtoal.com/genelab/.html/maketrie-stampede.c.html">This version</a> of Maketrie includes a lot more comments - it is an html file structured using "folds" so you can drill down to specific areas you are interested in.  Functionally it should be the same as the stripped version. (There were so many more comments than code that paradoxically it was getting harder to maintain the C source because of the comments!)
<tt>syntax: maketrie input.fastq</tt></li>
<li><a href="http://gtoal.com/genelab/.html/findoverlaps.c.html">findoverlaps</a>: this takes a file of k-mer reads which have been converted into trie form and locates all
the overlaps between all the k-mers and outputs the results in a file of overlap information.  It can be compiled
to use a simple internal format or to write out the same format that the AMOS suite accepts.  Maketrie and findoverlaps can be combined into a single program - we just kept them separate at first to speed up the software development cycle when only working on one of the two components.<br/><tt>syntax: findoverlaps input.fastq</tt></li>
<li><a href="http://gtoal.com/genelab/.html/glocate.c.html">glocate</a>: this builds a single contig which contains the k-mer passed as a parameter on the command line. To extract a complete contig, you would run this forward using a k-mer that you are trying to locate, and then run it again backwards using the reverse-complement of the k-mer, and join the left and right extensions together to get the complete contig.  Also you would use this for de novo genome assembly by effectively taking random k-mers as a starting point in order to find each independent contig. (In practice the search would pick subsequent k-mers from the remaining ones which had not yet been assigned to a contig, rather than picking completely at random)<br/>Example: <tt>glocate ~gtoal/genelab/data/40kreads-schliesky.fastq AAACCAGCAGATCCAGCACCAACGACGACGACATCAGTCTCAGCATAAGTGATCATATCCGTCATGTACCTTCTCGTCATCTCACGGGACACGATCGATTC</tt></li>
<li><a href="http://gtoal.com/genelab/.html/locate_read.c.html">locate_read</a>: given a single read as a command-line parameter, this locates that read (and identical copies of it) in the original fastq data file.<br/>Example: <tt>locate_read ~gtoal/genelab/data/40kreads-schliesky.fastq AAACCAGCAGATCCAGCACCAACGACGACGACATCAGTCTCAGCATAAGTGATCATATCCGTCATGTACCTTCTCGTCATCTCACGGGACACGATCGATTC</tt></li>
<li><a href="http://gtoal.com/genelab/.html/nearmatch.c.html">nearmatch</a>: given a k-mer parameter on the command-line, it outputs all the reads which match it, allowing for a small number of mis-matched letters.  Using the 'projectname-edges' trie file and the 'projectname-index' index back into the file of raw k-mer reads, this lookup is effectively instantaneous and does not require a large RAM to work.<br/><tt>nearmatch ~gtoal/genelab/data/40kreads-schliesky.fastq AAACCAGCAGATCCAGCACCAACGACGACGACATCAGTCTCAGCATAAGTGATCATATCCGTCATGTACCTTCTCGTCATCTCACGGGACACGATCGATTC</tt></li>
<li><a href="http://gtoal.com/genelab/.html/makeafg.c.html">makeafg</a>: Pass an entire contig on the command-line via a fastq file, and the file of reads that it was extracted from, and this generates a .afg representing that contig and containing all the reads which match 100% (or with a small allowed number of errors) so that it can be viewed
in 'tablet' etc. with the individual k-mers aligned to the contig.<br/><tt>makeafg 256seq.fastq 256seq.fastq-NEWCONTIG</tt><br/>You may also need makeafg.sh</li>
<li><a href="http://gtoal.com/genelab/.html/rcomp.c.html">rcomp</a>: this was a trivial little command-line utility to reverse and complement a k-mer.  It was easier
to keep the interfaces of programs like 'glocate' simple and just invoke them twice, once with each
direction of a k-mer as a parameter, than to have a bunch of options in every program to add the reverse-complement
form implicitly.<br/><tt>syntax: rcomp GATTACA</tt></li>
</ul>

There are some obvious programs missing, such as extending glocate to find all contigs and do a complete genome assembly.  For our project where the entire genome consisted of very few contigs and we knew what we were looking for, we never had a need for that software and with academic deadlines we didn't have an incentive to write code we wouldn't be using, but all the basic code to create that is in place and it would be a relatively small exercise to complete the suite. (which is one of the reasons I'm putting this on github in the hope that someone will do just that.)

One final word: all the code above <em>will</em> run on a single computer.  It does not need a cluster.  However that computer does need the mpi software to be installed.  It turns out that <a href="https://mpitutorial.com/tutorials/installing-mpich2/">installing MPI for a single machine isn't that hard</a> - in fact I just typed <tt>sudo apt install libmpich-dev</tt> and I was able to compile maketrie with <tt>make</tt> immediately (modulo a few compiler warnings that are new since 2013).  Unfortunately we converted our code to run in the MPI environment quite early on and I don't think we still have any version of the trie builder or overlap locator that doesn't refer to the MPI library.  Note there are some driving scripts for clusters under the 'old' subdirectory, and there's a custom version of 'maketrie' for our cluster in the sources above.<br/>
PS I would love to see this running on a cluster of $75 Raspberry Pi 4 machines with 8Gb each.  The 'big' systems we used were 32Gb Intel rack-mounted cards, but they cost a lot more than 4 * $75!


Regards,


Graham
