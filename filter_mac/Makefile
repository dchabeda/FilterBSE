# Makefile Configuration

# Shell and OS
SHELL = /bin/sh
OS = Linux

# Compiler
ARCHFLAGS = -arch arm64
CC = clang $(ARCHFLAGS)

# Optimization Flags
OPTFLAGS = -O3 -march=native
CFLAGS = $(OPTFLAGS) -Xclang -fopenmp

# Include and Lib Directories (customize these paths)
# With library paths and flags
HOMEBREW_CELLAR  = /opt/homebrew/Cellar

# OpenMP
OPENMP_INCLUDE   = $(HOMEBREW_CELLAR)/libomp/20.1.1/include/
OPENMP_LIB       = $(HOMEBREW_CELLAR)/libomp/20.1.1/lib/
OPENMP_CFLAGS    = -I$(OPENMP_INCLUDE)
OPENMP_LDFLAGS   = -L$(OPENMP_LIB) -lomp

# OpenMPI
OPENMPI_INCLUDE  = $(HOMEBREW_CELLAR)/open-mpi/5.0.7/include/
OPENMPI_LIB      = $(HOMEBREW_CELLAR)/open-mpi/5.0.7/lib/
OPENMPI_CFLAGS   = -I$(OPENMPI_INCLUDE)
OPENMPI_LDFLAGS  = -L$(OPENMPI_LIB) -lmpi

# LAPACK (Linear Algebra subroutines)
LAPACK_INCLUDE   = $(HOMEBREW_CELLAR)/lapack/3.12.1/include/
LAPACK_LIB       = $(HOMEBREW_CELLAR)/lapack/3.12.1/lib/
LAPACK_CFLAGS    = -I$(LAPACK_INCLUDE)
LAPACK_LDFLAGS   = -L$(LAPACK_LIB) -llapack

# FFTW (Fast Fourier Transform)
FFTW_INCLUDE     = $(HOMEBREW_CELLAR)/fftw/3.3.10_2/include/
FFTW_LIB         = $(HOMEBREW_CELLAR)/fftw/3.3.10_2/lib/
FFTW_CFLAGS      = -I$(FFTW_INCLUDE)
FFTW_LDFLAGS     = -L$(FFTW_LIB) -lfftw3 -lfftw3_omp -lfftw3_threads 

# Other Linux Libraries
LINUX_LIBS = -lpthread -lm

# Combine all flags and libraries
ALL_CFLAGS = $(CFLAGS) $(OPENMP_CFLAGS) $(OPENMPI_CFLAGS) $(LAPACK_CFLAGS) $(FFTW_CFLAGS)
ALL_LDFLAGS = $(OPENMP_LDFLAGS) $(OPENMPI_LDFLAGS) $(LAPACK_LDFLAGS) $(FFTW_LDFLAGS)
ALL_LIBS = $(ALL_LDFLAGS) $(LINUX_LIBS)

# Main target name
MAINNAM = filter.cplx

# Object files
OBJECTS = \
	main.o \
	mod_init.o mod_mem.o mod_pot.o \
	mod_filter.o mod_ortho.o mod_diag.o mod_sigma.o \
	mod_output.o \
	read.o write.o init.o save.o nerror.o aux.o \
	projectors.o norm.o rand.o interpolate.o strain.o coeff.o vector.o \
	energy.o hamiltonian.o phamiltonian.o \
	filter.o ortho.o Hmat.o  \
	init_periodic.o integral.o


CUBEOBJECTS = \
	makecube.o write.o size.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o vector.o strain.o

# Compilation rules
.SUFFIXES: .c .o

.c.o:
	$(CC) $(ALL_CFLAGS) -c $*.c

$(MAINNAM).x: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(ALL_LIBS)

# Phony targets
.PHONY: clean

clean:
	/bin/rm -f *.o *.x
