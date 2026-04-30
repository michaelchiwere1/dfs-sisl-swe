
FC = gfortran
OPT = -O3
FFLAGS= -fopenmp -frecord-marker=4 $(OPT)


#FC = ifort
#OPT = -O3 -xCORE-AVX512 -mcmodel=medium
#FFLAGS= -qopenmp -assume byterecl -std $(OPT)


INCLUDE=

OBJS =   \
dcosqb.o \
dcosqf.o \
dcosqi.o \
dcost.o  \
dcosti.o \
dcsqb1.o \
dcsqf1.o \
dradb2.o \
dradb3.o \
dradb4.o \
dradb5.o \
dradbg.o \
dradf2.o \
dradf3.o \
dradf4.o \
dradf5.o \
dradfg.o \
drfftb.o \
drfftf.o \
drffti.o \
drftb1.o \
drftf1.o \
drfti1.o \
dsinqb.o \
dsinqf.o \
dsinqi.o \
dsint.o  \
dsinti.o \

all: libbihar.a

libbihar.a: $(OBJS)
	ar cr libbihar.a $(OBJS)

clean:
	rm -rf *.o *.mod libbihar.a

.SUFFIXES: 
.SUFFIXES:      .f .o
.f.o:	
	$(FC) $(FFLAGS) $(INCLUDE) -c $*.f -o $*.o

