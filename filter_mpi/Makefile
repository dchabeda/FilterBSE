SHELL = /bin/sh

OS=Linux

# libraries ...
Linux_LIB = -qmkl -qopenmp -lpthread -lfftw3_omp -lfftw3 -lm
Linux_DYLIB = -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_thread -lmkl_core -lm

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64
Linux_CFLAGS = $(Linux_OPTFLAGS) -O3 -qopenmp -Wall
#Linux_CFLAGS = -g -qopenmp -O0 -Wall
#Linux_CFLAGS = -g -Wno-language-extension-token -O0 -Wall # uncomment this line to debug/valgrind

MAINNAM = Filter

# compiler ...
Linux_CC = cc

OBJECTS = \
	main.o mod_init.o mod_mem.o mod_pot.o mod_filter.o mod_ortho.o mod_diag.o mod_output.o\
	read.o init.o save.o projectors.o norm.o rand.o interpolate.o strain.o coeff.o \
	energy.o hamiltonian.o ortho.o Hmat.o write.o nerror.o integral.o phamiltonian.o \
	vector.o aux.o init_periodic.o

FAST = $(OBJECTS) filter_fast.o

CUBEOBJECTS = \
	makecube.o write.o size.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o vector.o strain.o

# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c  

# Compile the main executable Filter.x and link libraries
$(MAINNAM): $(FAST)
	$(Linux_CC) -o $(MAINNAM)_mpi.x $(OMP) $(LIB)


# If desired, run make cube to build the makecube.x executable
cube: $(CUBEOBJECTS)
	$(Linux_CC) -o makecube.x $(CUBEOBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
