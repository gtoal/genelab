# You will probably need to issue:
#    module load  mpi/openmpi
# outside of the makefile before running this.

all: maketrie findoverlaps glocate nearmatch
	@echo "# All up to date.  Now try: qsub fake10m-maketrie.job; qsub fake10m-findoverlaps.job "

print:
	ctohtml maketrie.c > maketrie.c.html
	ctohtml findoverlaps.c > findoverlaps.c.html
	ctohtml glocate.c > glocate.c.html
	ctohtml nearmatch.c > nearmatch.c.html
	ctohtml locate_read.c > locate_read.c.html
	ctohtml makeafg.c > makeafg.c.html
	ctohtml rcomp.c > rcomp.c.html
	ctohtml maketrie-stampede.c > maketrie-stampede.c.html

glocate: glocate.c
	cc -o glocate glocate.c
	cp glocate ~/bin/

nearmatch: nearmatch.c
	cc -o nearmatch nearmatch.c
	cp nearmatch ~/bin/

maketrie: maketrie.c
	mpicc -Wall -fopenmp -g -o maketrie maketrie.c
#	cc -O3 -fopenmp -o maketrie -Wall maketrie.c -I/usr/mpi/gcc/openmpi-1.4.3/include -L/usr/mpi/gcc/openmpi-1.4.3/lib64 -lmpi

findoverlaps: findoverlaps.c
	mpicc -Wall -fopenmp -g -o findoverlaps findoverlaps.c
#	cc -O3 -fopenmp -o findoverlaps -Wall findoverlaps.c -I/usr/mpi/gcc/openmpi-1.4.3/include -L/usr/mpi/gcc/openmpi-1.4.3/lib64 -lmpi
