SHELL = /bin/sh

OS=Linux

# libraries ...
Linux_LIB = -qmkl -qopenmp -lpthread -lfftw3 -lfftw3_threads -lfftw3_omp
Linux_DYLIB = -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_thread -lmkl_core -lm

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
Linux_OPTFLAGS = -O3 -DMKL_ILP64 -DFFTW_FFT -mavx2
Linux_CFLAGS = $(Linux_OPTFLAGS) -qopenmp
#Linux_CFLAGS = -g -O0 # uncomment this line to debug/valgrind

MAINNAM = bse_real

# compiler ...
Linux_CC = cc

# obhect files
OBJECTS = \
	main.o aux.o read.o write.o init.o basis.o save.o \
	mod_init.o spin.o mod_dipole.o dipole.o angular.o ang_mom_operators.o \
	mod_pot.o coulomb.o hartree.o mod_kernel.o \
	mod_bse.o bethe-salpeter.o diag.o optical.o
# mod_pot.o mod_dipole.o mod_kernel.o write.o save.o read.o basis.o init.o dipole.o coulomb.o\
#	hartree.o angular.o bethe-salpeter.o diag.o optical.o aux.o

# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c 


# Compile the main executable BSE.x and link libraries
$(MAINNAM): $(OBJECTS) 
	$(Linux_CC) -o $(MAINNAM).x $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
