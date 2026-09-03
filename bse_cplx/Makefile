SHELL = /bin/sh

# ============================================================================
#  bse_cplx : NVHPC + libsci build WITH optional GPU offload of the Coulomb
#             kernel (coulomb.c).  Native `double complex` typing is preserved
#             everywhere; the GPU path is selected at RUNTIME (input.par
#             `gpuAccel = 1`, auto-gated by device availability) and can be
#             forced off with `OMP_TARGET_OFFLOAD=DISABLED`.
#
#  Toolchain (source in the same shell, no pipes -- a pipe spawns a subshell
#  and the module swap is lost):
#       module load PrgEnv-nvidia cray-fftw cray-libsci
#
#  This produces bse_cplx_gpu.x.  The original Intel/MKL CPU-only build is
#  preserved verbatim in Makefile_mkl (produces bse_cplx.x); that paradigm is
#  also what bse_real uses.  The two builds share .o names, so run
#  `make clean` when switching toolchains.
#
#  GPU flags -mp=gpu -gpu=cc80 MUST appear on BOTH the compile and the link
#  step: NVHPC uses relocatable device code (RDC) by default, so the final
#  link has to run the device linker (nvlink).  Omitting them at link builds a
#  binary whose every `omp target` region aborts at runtime with
#  "Could not run target region on device 0".
# ============================================================================

OS = Linux

# libraries ...
FFTW_PREFIX  = $(CRAY_FFTW_PREFIX)/x86_64
FFTW_INCLUDE = -I$(FFTW_PREFIX)/include/
Linux_LIB    = -lpthread -L$(FFTW_PREFIX)/lib/ -lfftw3 -lfftw3_threads -lfftw3_omp -L$(CRAY_LIBSCI_PREFIX_DIR)/lib -lsci_nvidia
Linux_DYLIB  = -lm

LIB = $(Linux_DYLIB) $(Linux_LIB)

# flags ...
# USE_LIBSCI      -> fd.h pulls <lapacke.h> (cray-libsci) instead of <mkl.h>,
#                    and diag.c uses omp_set_num_threads instead of mkl_set_num_threads.
# USE_GPU_OFFLOAD -> compiles the OpenMP-target offload path in coulomb.c.
Linux_OPTFLAGS = -DFFTW_FFT -DUSE_LIBSCI -DUSE_GPU_OFFLOAD
Linux_GPUFLAGS = -mp=gpu -gpu=cc80
Linux_CFLAGS   = $(Linux_OPTFLAGS) -O3 -fast $(Linux_GPUFLAGS) -Minfo=mp $(FFTW_INCLUDE)

MAINNAM = bse_cplx_gpu

# compiler ...
Linux_CC = cc

# object files
OBJECTS = \
	main.o aux.o read.o write.o init.o basis.o save.o \
	mod_init.o spin.o mod_dipole.o dipole.o angular.o ang_mom_operators.o \
	mod_pot.o coulomb.o hartree.o mod_kernel.o \
	mod_bse.o bethe-salpeter.o diag.o optical.o

# compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c

# link the executable and libraries (GPU flags required here too -- RDC)
$(MAINNAM): $(OBJECTS)
	$(Linux_CC) $(Linux_GPUFLAGS) -o $(MAINNAM).x $(OBJECTS) $(LIB)

clean:
	/bin/rm -f *.o *.x
