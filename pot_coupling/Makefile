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

MAINNAM = pot_coupling

# compiler ...
Linux_CC = cc

OBJECTS = \
        read_bse.o read_filter.o init_bse.o init_filter.o save.o projectors.o interpolate.o \
        pot_coupling.o write.o nerror.o strain.o basis.o vector.o pot_mat_elems.c hamiltonian.o


# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c  

# Compile the main executable Filter.x and link libraries
$(MAINNAM): $(OBJECTS) 
	$(Linux_CC) -o $(MAINNAM).x $(OBJECTS) $(LIB)


clean:
	/bin/rm *.o *.x
