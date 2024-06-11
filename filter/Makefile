SHELL = /bin/sh

MKLROOT = /opt/intel/compilers_and_libraries/linux/mkl
INTELLINK =  ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

Linux_LIB =  ${INTELLINK} -L${HOME}/usr/lib -lfftw3
Linux_MYLIB = 
Linux_CLIB = -lm
Linux_FCLIB = 	

CYGWIN_NT-6.3_LIB = -llapack -lblas -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-6.3_MYLIB = 
CYGWIN_NT-6.3_CLIB =  -lm
CYGWIN_NT-6.3_FCLIB = 	

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
#Linux_OPTFLAGS = -fast -qopenmp -static
Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3  -qopenmp -static
#Linux_OPTFLAGS = -O0 -g
Linux_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW -I${HOME}/usr/include -I /opt/rocks/include -I"${MKLROOT}/include" -Wall -pedantic 
Linux_FFLAGS = $(${OS}_OPTFLAGS)

CYGWIN_NT-6.3_OPTFLAGS = -O3 -ffast-math 
CYGWIN_NT-6.3_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
CYGWIN_NT-6.3_FFLAGS = $(${OS}_OPTFLAGS) 

MAINNAM = Filter

#compiler ...
Linux_CC = icx
Linux_FF = ifx
Linux_LD = icx
CYGWIN_NT-6.3_CC = gcc
CYGWIN_NT-6.3_FF = gfortran
CYGWIN_NT-6.3_LD = gcc

OBJECTS = \
	main.o read.o init.o save.o projectors.o norm.o rand.o interpolate.o strain.o \
	energy.o coeff.o hamiltonian.o filter.o ortho.o Hmat.o write.o nerror.o \
        vector.o dipole.o angular.o pseudopotential.o

TESTOBJECTS = \
	      testmain.o angular.o norm.o nerror.o write.o init.o size.o rand.o\
	      interpolate.o projectors.o read.o

CUBEOBJECTS = \
		makecube.o write.o size.o init.o nerror.o read.o save.o interpolate.o rand.o norm.o projectors.o
# compilation ...

.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_${OS} $(${OS}_CFLAGS) -c  $*.c  

all: $(MAINNAM)

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)


test: $(TESTOBJECTS)
	$(${OS}_LD) -o test.x $(${OS}_CFLAGS) $(TESTOBJECTS) $(LIB)

cube: $(CUBEOBJECTS)
	$(${OS}_LD) -o makecube.x $(${OS}_CFLAGS) $(CUBEOBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
