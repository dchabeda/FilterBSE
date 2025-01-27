SHELL = /bin/sh

OS=Linux

# libraries ...
Linux_LIB = -qmkl -qopenmp -lpthread -lfftw3 -lfftw3_threads -lfftw3_omp
Linux_DYLIB = -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_thread -lmkl_core -lm

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64
Linux_CFLAGS = $(Linux_OPTFLAGS) -O3 -qopenmp
#Linux_CFLAGS = -g -qopenmp -O0
#Linux_CFLAGS = -g -O0 # uncomment this line to debug/valgrind

MAINNAM = BSE

# compiler ...
Linux_CC = cc

# obhect files
OBJECTS = \
	main.o write.o save.o read.o basis.o init.o dipole.o coulomb.o hartree.o angular.o bethe-salpeter.o diag.o optical.o

# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c 


# Compile the main executable BSE.x and link libraries
$(MAINNAM): $(OBJECTS) 
	$(Linux_CC) -o $(MAINNAM).x $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
