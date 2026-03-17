SHELL = /bin/sh

OS = Linux
# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
MKLROOT = /opt/intel/compilers_and_libraries/linux/mkl

INTELLINK = ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a          \
            -Wl,--start-group                                          \
                ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a           \
                ${MKLROOT}/lib/intel64/libmkl_intel_thread.a          \
                ${MKLROOT}/lib/intel64/libmkl_core.a                  \
                ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a  \
            -Wl,--end-group                                            \
            -liomp5 -lpthread -lm -ldl

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
Linux_DYLIB =
Linux_LIB   = ${INTELLINK} -L${HOME}/usr/lib -lfftw3

LIB = $(Linux_DYLIB) $(Linux_LIB)

# ------------------------------------------------------------------------------
# Flags
# ------------------------------------------------------------------------------
Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3
Linux_CFLAGS   = $(Linux_OPTFLAGS) -qopenmp -static
#Linux_CFLAGS  = -g -qopenmp -O0
#Linux_CFLAGS  = -g -O0  # uncomment this line to debug/valgrind

# ------------------------------------------------------------------------------
# Compiler
# ------------------------------------------------------------------------------
Linux_CC = icx

# ------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------
MAINNAM = bse_cplx

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
