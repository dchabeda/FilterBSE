SHELL = /bin/sh

OS=Linux

# location of mkl .a files on perlmutter: /opt/intel/oneapi/mkl/latest/lib/intel64/
MKLROOT = /opt/intel/oneapi/mkl/latest/lib/intel64

# libraries ...
INTEL_LINK = -Wl,--start-group $(MKLROOT)/libmkl_intel_ilp64.a $(MKLROOT)/libmkl_lapack95_ilp64.a $(MKLROOT)/libmkl_blas95_ilp64.a $(MKLROOT)/libmkl_intel_thread.a $(MKLROOT)/libmkl_core.a -Wl,--end-group

Linux_LIB = $(INTEL_LINK) -qopenmp -lpthread -lfftw3
Linux_DYLIB = -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_thread -lmkl_core -lm
Linux_CLIB = -lm

LIB = $(Linux_DYLIB) $(Linux_LIB) $(Linux_CLIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64 -O3
Linux_CFLAGS = $(Linux_OPTFLAGS) -g -O3 -fopenmp #-I/usr/local/include
#Linux_CFLAGS = -g -O0
Linux_FFLAGS = $(Linux_OPTFLAGS)

MAINNAM = Filter

# compiler ...
Linux_CC = cc
Linux_FF = ftn
Linux_LD = cc

OBJECTS = \
    main.o read.o init.o save.o projectors.o norm.o rand.o interpolate.o strain.o \
    energy.o coeff.o hamiltonian.o filter.o ortho.o Hmat.o write.o nerror.o \
    vector.o dipole.o angular.o pseudopotential.o

TESTOBJECTS = \
    testmain.o angular.o norm.o nerror.o write.o init.o size.o rand.o \
    interpolate.o projectors.o read.o

CUBEOBJECTS = \
    makecube.o write.o size.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o vector.o strain.o

POTMATOBJECTS = \
    pot_mat_elems.o write.o size.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o \
    hamiltonian.o strain.o vector.o

# compilation ...

all: $(MAINNAM)

.f.o:
	$(Linux_FF) $(Linux_FFLAGS) -c $*.f
.c.o:
	$(Linux_CC) -DOS_$(OS) $(Linux_CFLAGS) -c $*.c  

$(MAINNAM): $(OBJECTS) 
	$(Linux_LD) -o $(MAINNAM).x $(OBJECTS) $(LIB)

test: $(TESTOBJECTS)
	$(Linux_LD) -o test.x $(TESTOBJECTS) $(LIB)

cube: $(CUBEOBJECTS)
	$(Linux_LD) -o makecube.x $(CUBEOBJECTS) $(LIB)

potmat: $(POTMATOBJECTS)
	$(Linux_LD) -o pot_mat_elems.x $(POTMATOBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
