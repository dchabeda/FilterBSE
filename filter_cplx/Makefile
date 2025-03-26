# Makefile Configuration

# Shell and OS
SHELL = /bin/sh
OS = Linux

# Compiler
Linux_CC = clang $(ARCHFLAGS)
ARCHFLAGS = -arch arm64

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
LAPACK_LDFLAGS   = -L$(LAPACK_LIB) -llapacke

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
MAINNAM = Filter

# Object files
OBJECTS = \
	main.o mod_init.o mod_mem.o mod_pot.o mod_filter.o mod_ortho.o mod_diag.o mod_output.o \
	read.o init.o save.o projectors.o norm.o rand.o interpolate.o strain.o coeff.o \
	energy.o hamiltonian.o ortho.o Hmat.o write.o nerror.o integral.o phamiltonian.o \
	vector.o aux.o filter.o


CUBEOBJECTS = \
	makecube.o write.o size.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o vector.o strain.o

# Compilation rules
.SUFFIXES: .c .o

.c.o:
	$(CC) $(ALL_CFLAGS) -c $*.c

# Compilation targets
$(MAINNAM)_omp.x: $(OMP)

mpi: $(MPI)
	$(Linux_CC) -o $(MAINNAM)_mpi.x $(MPI) $(Linux_LIB)

fast: $(FAST)
	$(Linux_CC) -o $(MAINNAM)_fast.x $(FAST) $(Linux_LIB)

cube: $(CUBEOBJECTS)
	$(Linux_CC) -o makecube.x $(CUBEOBJECTS) $(Linux_LIB)

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

# Phony targets
.PHONY: clean $(MAINNAM) mpi fast cube

clean:
	/bin/rm -f *.o *.x