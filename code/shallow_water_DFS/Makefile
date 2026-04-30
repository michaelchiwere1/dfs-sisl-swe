FC := gfortran
OPT := -O3
MODDIR := .mod
DEPDIR := .deps

#FFLAGS := -fopenmp -frecord-marker=4 $(OPT) -J$(MODDIR) -I$(MODDIR) -cpp
#OPT = -g -O0 -fcheck=all -fbacktrace -Wall
FFLAGS = -fopenmp -frecord-marker=4 $(OPT) -J$(MODDIR) -I$(MODDIR) -cpp

CPPFLAGS := -I/usr/include -I$(HOME)/finufft/include
LDLIBS := -lfftw3 -lfftw3_omp -lgomp -lstdc++ $(HOME)/finufft/lib-static/libfinufft.a -lstdc++

FRAME := \
prm_phconst.o \
com_dfs.o \
e_time.o \
legendre.o \
monit_spectrum.o \
monit_norm.o \
monit.o \
grads.o \
gmean_to_zero_dfs.o \
tridiagonal.o \
tridiagonal_cmplx.o \
pentadiagonal.o \
pentadiagonal_4.o \
pentadiagonal_cmplx.o \
ninediagonal_4.o \
helmholtz_dfs_old.o \
helmholtz_dfs.o \
helmholtz_dfs_cmplx_old.o \
helmholtz_dfs_cmplx.o \
hdiff4_dfs.o \
xderiv.o \
yderiv_dfs_old.o \
uv2divrot_dfs_old.o \
uv2divrot_dfs.o \
divrot2uv_dfs_old.o \
divrot2uv_dfs.o \
fft_y.o \
fft.o \
finufft_mod.o \
comput2DCoeffs.o \
newdfscoeffs.o



SW_DFS := $(FRAME) sw_dfs.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a
ADV_DFS := $(FRAME) adv_dfs.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a
TEST := $(FRAME) test.o ../bihar/libbihar.a ../ispack-3.0.1/libispack3.a

PROGRAMS := sw_dfs adv_dfs test

.PHONY: all clean
all: $(PROGRAMS)

sw_dfs: $(SW_DFS)
	$(FC) -fopenmp -o $@ $^ $(LDLIBS)

adv_dfs: $(ADV_DFS)
	$(FC) -fopenmp -o $@ $^ $(LDLIBS)

test: $(TEST)
	$(FC) -fopenmp -o $@ $^ $(LDLIBS)

$(MODDIR) $(DEPDIR):
	mkdir -p $@

%.o: %.F90 | $(MODDIR) $(DEPDIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

%.o: %.f90 | $(MODDIR) $(DEPDIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

%.o: %.f | $(MODDIR) $(DEPDIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

finufft_mod.o: $(HOME)/finufft/include/finufft_mod.f90 | $(MODDIR) $(DEPDIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@


# Fortran module dependency ordering for make -j

legendre.o: e_time.o
monit_spectrum.o: prm_phconst.o legendre.o

gmean_to_zero_dfs.o: com_dfs.o
xderiv.o: com_dfs.o
yderiv_dfs_old.o: com_dfs.o

fft.o: e_time.o
fft_y.o: com_dfs.o prm_phconst.o tridiagonal.o pentadiagonal.o e_time.o

uv2divrot_dfs_old.o: prm_phconst.o yderiv_dfs_old.o com_dfs.o \
                     pentadiagonal_4.o ninediagonal_4.o pentadiagonal.o \
                     tridiagonal.o e_time.o

uv2divrot_dfs.o: prm_phconst.o com_dfs.o pentadiagonal_4.o ninediagonal_4.o \
                 pentadiagonal.o tridiagonal.o uv2divrot_dfs_old.o e_time.o

divrot2uv_dfs_old.o: yderiv_dfs_old.o
divrot2uv_dfs.o: prm_phconst.o com_dfs.o pentadiagonal.o tridiagonal.o \
                 gmean_to_zero_dfs.o divrot2uv_dfs_old.o e_time.o

helmholtz_dfs_old.o: prm_phconst.o com_dfs.o pentadiagonal.o tridiagonal.o
helmholtz_dfs.o: prm_phconst.o com_dfs.o pentadiagonal.o tridiagonal.o \
                 helmholtz_dfs_old.o

helmholtz_dfs_cmplx_old.o: prm_phconst.o com_dfs.o pentadiagonal_cmplx.o \
                           tridiagonal_cmplx.o
helmholtz_dfs_cmplx.o: prm_phconst.o com_dfs.o pentadiagonal_cmplx.o \
                       tridiagonal_cmplx.o helmholtz_dfs_cmplx_old.o

hdiff4_dfs.o: com_dfs.o helmholtz_dfs_cmplx.o

newdfscoeffs.o: fft.o fft_y.o uv2divrot_dfs.o divrot2uv_dfs.o

sw_dfs.o: prm_phconst.o fft.o fft_y.o xderiv.o com_dfs.o uv2divrot_dfs.o \
          divrot2uv_dfs.o helmholtz_dfs.o hdiff4_dfs.o monit.o \
          monit_spectrum.o monit_norm.o grads.o e_time.o finufft_mod.o \
          comput2DCoeffs.o newdfscoeffs.o


adv_dfs.o: prm_phconst.o fft.o fft_y.o xderiv.o com_dfs.o uv2divrot_dfs.o \
           divrot2uv_dfs.o helmholtz_dfs.o hdiff4_dfs.o monit.o \
           monit_spectrum.o monit_norm.o grads.o e_time.o

test.o: com_dfs.o prm_phconst.o fft.o fft_y.o divrot2uv_dfs.o \
        uv2divrot_dfs.o xderiv.o helmholtz_dfs.o hdiff4_dfs.o \
        gmean_to_zero_dfs.o monit.o grads.o


clean:
	rm -rf a.out $(PROGRAMS) *.o *.mod $(MODDIR) $(DEPDIR) *~
