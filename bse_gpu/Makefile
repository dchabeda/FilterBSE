SHELL = /bin/sh

OS=Linux

# libraries ...
FFTW_PREFIX = $(CRAY_FFTW_PREFIX)/x86_64
Linux_LIB = -lpthread -L$(FFTW_PREFIX)/lib/ -lfftw3 -lfftw3_threads -lfftw3_omp -L$(CRAY_LIBSCI_PREFIX_DIR)/lib -lsci_nvidia
Linux_DYLIB = -lm
FFTW_INCLUDE = -I$(CRAY_FFTW_PREFIX)/x86_64/include/

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
FFT_FLAGS = -I$(CRAY_FFTW_PREFIX)/x86_64/include

Linux_OPTFLAGS = -DMKL_ILP64 -DFFTW_FFT
# GPU offload flags. These MUST appear on BOTH the compile and the link step:
# NVHPC uses relocatable device code (RDC) by default, so the final link has to
# run the device linker (nvlink). Omitting them at link time produces a binary
# whose CUDA module fails to load at runtime ("Could not run target region").
Linux_GPUFLAGS = -mp=gpu -gpu=cc80
Linux_CFLAGS = $(Linux_OPTFLAGS) -O3 -g -Minfo=all $(Linux_GPUFLAGS) -fast $(FFT_FLAGS)
#Linux_CFLAGS = $(Linux_OPTFLAGS) -O3 -march=native -Mvect=simd -g -Minfo=all -gpu=debug -traceback -lineinfo -mp=gpu -gpu=cc80 -fast -target-accel=nvidia80 -lnvToolsExt
#Linux_CFLAGS = -g -qopenmp -O0
#Linux_CFLAGS = -g -O0 # uncomment this line to debug/valgrind

MAINNAM = BSE

# compiler ...
Linux_CC = cc

# obhect files
OBJECTS = \
	main.o mod_init.o spin.o mod_pot.o mod_dipole.o mod_kernel.o write.o save.o read.o basis.o init.o dipole.o coulomb.o\
	hartree.o angular.o bethe-salpeter.o diag.o optical.o aux.o

# compilation ...
# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c 


# Compile the main executable BSE.x and link libraries
$(MAINNAM): $(OBJECTS)
	$(Linux_CC) $(Linux_GPUFLAGS) -o $(MAINNAM).x $(OBJECTS) $(LIB) $(FFTW_INCLUDE)

clean:
	/bin/rm *.o *.x
