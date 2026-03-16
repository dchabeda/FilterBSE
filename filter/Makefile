SHELL = /bin/sh
OS    = Linux

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
MAINNAM = Filter

OBJECTS = \
    main.o          \
    read.o          \
    init.o          \
    save.o          \
    projectors.o    \
    norm.o          \
    rand.o          \
    interpolate.o   \
    strain.o        \
    energy.o        \
    coeff.o         \
    hamiltonian.o   \
    filter.o        \
    ortho.o         \
    Hmat.o          \
    write.o         \
    nerror.o        \
    vector.o        \
    dipole.o        \
    angular.o       \
    pseudopotential.o \
    aux.o

CUBEOBJECTS = \
    makecube.o      \
    write.o         \
    size.o          \
    init.o          \
    nerror.o        \
    read.o          \
    save.o          \
    interpolate.o   \
    rand.o          \
    norm.o          \
    projectors.o    \
    vector.o        \
    strain.o

GENBSEOBJECTS = \
    gen_bse_psi.o   \
    write.o         \
    size.o          \
    init.o          \
    nerror.o        \
    read.o          \
    save.o          \
    interpolate.o   \
    rand.o          \
    norm.o          \
    projectors.o    \
    vector.o        \
    strain.o

# ------------------------------------------------------------------------------
# Rules
# ------------------------------------------------------------------------------

# Compilation rule for all .c files
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c

# Compile the main executable Filter.x and link libraries
$(MAINNAM): $(OBJECTS)
	$(Linux_CC) -o $(MAINNAM).x $(Linux_CFLAGS) $(OBJECTS) $(LIB)

# If desired, run make cube to build the makecube.x executable
cube: $(CUBEOBJECTS)
	$(Linux_CC) -o makecube.x $(CUBEOBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
