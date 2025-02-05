SHELL = /bin/sh

OS=Linux

# libraries ...
Linux_LIB = -qmkl -qopenmp -lpthread -lfftw3_omp -lfftw3 -lm
Linux_DYLIB = -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_thread -lmkl_core -lm

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64
Linux_CFLAGS = $(Linux_OPTFLAGS) -O3 -qopenmp
#Linux_CFLAGS = -g -qopenmp -O0
#Linux_CFLAGS = -g -O0 # uncomment this line to debug/valgrind

MAINNAM = Filter

# compiler ...
Linux_CC = cc

OBJECTS = \
	main.o read.o init.o save.o projectors.o norm.o rand.o interpolate.o strain.o \
	energy.o coeff.o hamiltonian.o ortho.o Hmat.o write.o nerror.o periocic.o\
	vector.o dipole.o angular.o pseudopotential.o aux.o phamiltonian.o ghaltonian.o

OMP = $(OBJECTS) filter_omp.o

FAST = $(OBJECTS) filter_fast.o

MPI = $(OBJECTS) filter_mpi.o

PERIODIC = $(OBJECTS) filter_omp.o kspace.o

CUBEOBJECTS = \
	makecube.o write.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o vector.o strain.o

# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c  

# Compile the main executable Filter.x and link libraries
$(MAINNAM): $(PERIODIC)
	$(Linux_CC) -o $(MAINNAM)_periodic.x $(PERIODIC) $(LIB)

mpi: $(MPI)
	$(Linux_CC) -o $(MAINNAM)_mpi.x $(MPI) $(LIB)

fast: $(FAST)
	$(Linux_CC) -o $(MAINNAM)_fast.x $(FAST) $(LIB)

omp: $(OMP)
	$(Linux_CC) -o $(MAINNAM).x $(OMP) $(LIB)

# If desired, run make cube to build the makecube.x executable
cube: $(CUBEOBJECTS)
	$(Linux_CC) -o makecube.x $(CUBEOBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
