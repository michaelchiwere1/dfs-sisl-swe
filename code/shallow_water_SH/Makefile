
FC = gfortran
#OPT = -g -Wall -O -fbounds-check -Wuninitialized -fbacktrace
OPT =  -O3 -fmax-stack-var-size=32640
FFLAGS= -fopenmp -frecord-marker=4 $(OPT)
LFLAGS= -fopenmp  -lfftw3 -lfftw3_omp -lgomp -lstdc++  $(HOME)/finufft/lib/libfinufft.so -lblas -llapack 


#FC = ifort
#OPT = -O3 -xCORE-AVX512 -mcmodel=medium
#FFLAGS= -qopenmp -assume byterecl -std $(OPT)
#LFLAGS= -qopenmp


INCLUDE = -I /usr/include/ -I  $(HOME)/finufft/include -L /usr/bin

FRAME= \
e_time.o \
prm_phconst.o \
monit_spectrum.o \
monit_norm.o \
monit.o \
grads.o \
xderiv.o \
fft.o    \
legendre.o  \
uv2rotdiv.o \
rotdiv2uv.o \
finufft_mod.o

all: sw_sh adv_sh test_laplacian_helmholtz test_truncate

sw_sh : $(FRAME) sw_sh.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a 
	$(FC) $(LFLAGS) $(FRAME) sw_sh.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a ~/finufft/lib/libfinufft.so $(LFLAGS) -o sw_sh

#sw_eul_sh : $(FRAME) sw_eul_sh.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a
#	$(FC) $(LFLAGS) $(FRAME) adv_sh.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a -o sw_eul_sh

adv_sh : $(FRAME) adv_sh.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a
	$(FC) $(LFLAGS) $(FRAME) adv_sh.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a -o adv_sh

test_laplacian_helmholtz : $(FRAME) test_laplacian_helmholtz.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a
	$(FC) $(LFLAGS) $(FRAME) test_laplacian_helmholtz.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a -o test_laplacian_helmholtz

test_truncate : $(FRAME) test_truncate.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a
	$(FC) $(LFLAGS) $(FRAME) test_truncate.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a -o test_truncate

clean:
	rm -rf *.o *.mod sw_sh adv_sh test_truncate test_laplacian_helmholtz *~

.SUFFIXES: 
finufft_mod.o:
	 $(FC) $(FFLAGS) $(INCLUDE) -c   $(HOME)/finufft/include/finufft_mod.f90 -o finufft_mod.o
.SUFFIXES:      .F90 .o
.F90.o:	
	$(FC) $(FFLAGS) $(INCLUDE) -c $*.F90 -o $*.o
