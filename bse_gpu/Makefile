SHELL = /bin/sh

OS=Linux

# libraries ...
Linux_LIB = -lpthread -lfftw3 -lfftw3_threads -lfftw3_omp -L$(CRAY_LIBSCI_PREFIX_DIR)/lib -lsci_nvidia
Linux_DYLIB = -lm

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64 -DFFTW_FFT
#Linux_CFLAGS = $(Linux_OPTFLAGS) -O3 -msse4 -g -Minfo=accel -gpu=debug -traceback -lineinfo -mp=gpu -gpu=cc80 -fast -lnvToolsExt
#Linux_CFLAGS = -g -qopenmp -O0
Linux_CFLAGS = -g -O0 # uncomment this line to debug/valgrind

MAINNAM = BSE

# compiler ...
Linux_CC = cc

# obhect files
OBJECTS = \
	main.o mod_init.o spin.o write.o save.o read.o basis.o init.o dipole.o coulomb.o\
	hartree.o angular.o bethe-salpeter.o diag.o optical.o aux.o

# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c 


# Compile the main executable BSE.x and link libraries
$(MAINNAM): $(OBJECTS) 
	$(Linux_CC) -o $(MAINNAM).x $(OBJECTS) $(LIB) -lnvToolsExt

clean:
	/bin/rm *.o *.x
