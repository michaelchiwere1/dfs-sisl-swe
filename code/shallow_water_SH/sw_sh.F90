program sw_sh

  !! Shallow water model using spherical harmonics

  use prm_phconst, only : PI, PI2, ER, OMG, GRAV
  use legendre, only : legendre__ini, legendre__g2w_uv, legendre__g2w,  &
   &                   legendre__w2g_uv, legendre__w2g_dy, legendre__w2g
  use fft, only : fft_run_towave, fft_run_togrid
  use uv2rotdiv, only : uv2rotdiv_run
  use rotdiv2uv, only : rotdiv2uv_run
  use xderiv, only : xderiv_run
  use monit, only : monit__ini, monit__output
  use monit_spectrum, only : monit_spectrum__ini, monit_spectrum__output
  use monit_norm, only : monit_norm__ini, monit_norm__output
  use grads, only : grads__outxy
  use omp_lib
  use e_time, only : e_time__start, e_time__end, e_time__output
  use finufft_mod
!
  implicit none
  include 'fftw3.f03'
!  integer,parameter :: JCN_INITIAL = 1  !! Williamson test case 1 (Advection)
  integer,parameter :: JCN_INITIAL = 2  !! Williamson test case 2
! integer,parameter :: JCN_INITIAL = 5  !! Williamson test case 5
 ! integer,parameter :: JCN_INITIAL = 6  !! Williamson test case 6
 ! integer,parameter :: JCN_INITIAL = 10 !! Galewsky 2004 test case
!  integer,parameter :: JCN_INITIAL = 11 !! Galewsky 2004 like test case

 ! integer,parameter :: JCN_PROG_H = 0 !! h*=h-hs is prognostic variable
  integer,parameter :: JCN_PROG_H = 1 !! h is prognostic variable
  
  integer,parameter :: JCN_HDIFF = 0  !! No diffusion
!  integer,parameter :: JCN_HDIFF = 1  !! 2nd order diff. for Galewsky test case
 ! integer,parameter :: JCN_HDIFF = 2  !! 4th order hyper diffusion
!
 ! integer,parameter :: JCN_DEPARTURE = 0  !! Use extrapolated winds
  integer,parameter :: JCN_DEPARTURE = 2  !! Use predicted winds

!  integer,parameter :: JCN_MONIT = 0  !! No monitor output
  integer,parameter :: JCN_MONIT = 1  !! Write monitor output

!  integer,parameter :: JCN_MONIT_SPECTRUM = 0  !! No kinetic energy spectrum output
  integer,parameter :: JCN_MONIT_SPECTRUM = 1   !! Write kinetic energy spectrum output
  
  integer,parameter :: INTHR_MONIT = 24   !! Interval of monitor output (hour)

 ! integer,parameter :: IMAX = 16
 ! integer,parameter :: JMAX = 8 
 ! integer,parameter :: NMAX = 6
 ! real(8),parameter :: TIMESTEP = 3600.0d0   !! 1 hour
!
 ! integer,parameter :: IMAX = 128            !! About 300km resolution
 ! integer,parameter :: JMAX = 64             !! J is from North to South
 ! integer,parameter :: NMAX = 62             !! Linear grid
 ! integer,parameter :: NMAX = 42             !! Quadric grid
 ! real(8),parameter :: TIMESTEP = 1200.0d0  !! 1 hour
!
  integer,parameter :: IMAX = 160            !! 250km resolution
  integer,parameter :: JMAX = 80             !! J is from North to South
  integer,parameter :: NMAX = 78             !! Linear grid
!  integer,parameter :: NMAX = 53             !! Quadric grid
 ! real(8),parameter :: TIMESTEP = 3600.0d0/128.0d0   !! 1 hour
!
!  integer,parameter :: IMAX = 320            !! 120km resolution
 ! integer,parameter :: JMAX = 160            !! J is from North to South
 ! integer,parameter :: NMAX = 158            !! Linear grid
!  integer,parameter :: NMAX = 106            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 1800.0d0   !! 30 min.
!
 ! integer,parameter :: IMAX = 640            !! 60km resolution
!  integer,parameter :: JMAX = 320            !! J is from North to South
!  integer,parameter :: NMAX = 318            !! Linear grid
!  integer,parameter :: NMAX = 213            !! Quadric grid
  real(8),parameter :: TIMESTEP = 60.0d0   !! 20 min.
!
 ! integer,parameter :: IMAX = 1920           !! 20km resolution
!  integer,parameter :: JMAX = 960            !! J is from North to South
!  integer,parameter :: NMAX = 958            !! Linear grid
!  integer,parameter :: NMAX = 639            !! Quadric grid
 ! real(8),parameter :: TIMESTEP = 300.0d0    !! 10 min.
!
!  integer,parameter :: IMAX = 3840           !! 10km resolution
!  integer,parameter :: JMAX = 1920           !! J is from North to South
!  integer,parameter :: NMAX = 1918           !! Linear grid
!  integer,parameter :: NMAX = 1279           !! Quadric grid
!  integer,parameter :: NMAX = 959            !! Cubic grid
!  real(8),parameter :: TIMESTEP = 360.0d0    !! 6 min.
!
!  integer,parameter :: IMAX = 7680           !! 5km resolution
!  integer,parameter :: JMAX = 3840           !! J is from North to South
!  integer,parameter :: NMAX = 3838           !! Linear grid
!  integer,parameter :: NMAX = 2559           !! Quadric grid
!  integer,parameter :: NMAX = 1919           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 225.0d0    !! 225 sec.
!
!  integer,parameter :: IMAX = 15360          !! 2.6km resolution
!  integer,parameter :: JMAX = 7680           !! J is from North to South
!  integer,parameter :: NMAX = 7678           !! Linear grid
!  integer,parameter :: NMAX = 5119           !! Quadric grid
!  integer,parameter :: NMAX = 3839           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 144.0d0    !! 144 sec.
!
!  integer,parameter :: IMAX = 20480          !! 2.0km resolution
!  integer,parameter :: JMAX = 10240          !! J is from North to South
!  integer,parameter :: NMAX = 10238          !! Linear grid
!  integer,parameter :: NMAX = 6826           !! Quadric grid
!  integer,parameter :: NMAX = 5119           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 120.0d0    !! 120 sec.
!
!  integer,parameter :: IMAX = 30720          !! 1.3km resolution
!  integer,parameter :: JMAX = 15360          !! J is from North to South
!  integer,parameter :: NMAX = 15358          !! Linear grid   
!  integer,parameter :: NMAX = 10239          !! Quadric grid
!  integer,parameter :: NMAX = 7679           !! Cubic grid
  !real(8),parameter :: TIMESTEP = 90.0d0     !! 90 sec.
!
 !  integer, parameter :: IT = 0                      !! LG interpolation
   integer, parameter :: IT = 1                       !! DFS interpolation 

!  integer,parameter :: IVECLEN = 14  !! For test
  integer,parameter :: IVECLEN = min(IMAX,512)     !! For scalar machine
  integer,parameter :: NDIV = (IMAX-1)/IVECLEN + 1 !! For scalar machine
  
!  integer,parameter :: NDIV = 1                    !! For vector machine
!  integer,parameter :: IVECLEN = (IMAX-1)/NDIV + 1 !! For vector machine

  integer,parameter :: MMAX = NMAX      !! MMAX <= NMAX

  integer,parameter :: NNUM = NMAX+1
  integer,parameter :: MNUM = MMAX+1
  integer,save :: MNWAV
  integer,save :: MNWAV_UV
!
  integer,parameter :: IMAX2=IMAX/2 !! IMAX must be even.
!
  integer,parameter :: MGN_X = 15
  integer,parameter :: MGN_Y = 15
!
  real(8),parameter :: DLON = PI2/IMAX
!
  integer,save :: MNSTART(0:MMAX+1)
  integer,save :: MNSTART_UV(0:MMAX+1)
  
  real(8),save :: WGAUSS(7)
  real(8),save :: XGAUSS(7)
!
  real(8),save :: ALON(IMAX)
  real(8),save :: ALON1(IMAX)
  real(8),save :: SINLON(IMAX)
  real(8),save :: COSLON(IMAX)
  real(8),save :: ALAT(JMAX)
  real(8),save :: ALAT_SL(-6:JMAX+7)   !! latitude for semi-Lagrangian
  real(8),save :: YLAT(JMAX)
  real(8),save :: COSLAT(JMAX)
  real(8),save :: COSLAT_INV(JMAX)
  real(8),save :: SINLAT(JMAX)
!
  real(8),save :: ACOSLAT_INV(JMAX)
  real(8),save :: ACOS2LAT_INV(JMAX)
!
  real(8),save :: WEIGHT(JMAX)  !! Latitudinal weight (Half of Gaussian weight)
!
  real(8),save,allocatable :: ERFNN1(:)
  integer,save,allocatable :: NTOTAL(:)
  integer,save,allocatable :: MWAVE(:)
!
! -----------------------------------------------------------------
!
  call e_time__start(1,"shallow water")

  call initialize
!
  call main
!
  call e_time__end(1,"shallow water")
  call e_time__output
!
!******************************************************************
contains
!******************************************************************

subroutine initialize
!
  integer :: i,j
  integer :: m, mn
  real(8) :: an
!
! -----------------------------------------------------------------------------
!
  write(6,*) "JCN_INITIAL        =", JCN_INITIAL
  write(6,*) "JCN_PROG_H         =", JCN_PROG_H
  write(6,*) "JCN_HDIFF          =", JCN_HDIFF
  write(6,*) "JCN_DEPARTURE      =", JCN_DEPARTURE
  write(6,*) "JCN_MONIT          =", JCN_MONIT
  write(6,*) "JCN_MONIT_SPECTRUM =", JCN_MONIT_SPECTRUM
  write(6,*) "INTHR_MONIT        =", INTHR_MONIT
  write(6,*) "IMAX     =", IMAX 
  write(6,*) "JMAX     =", JMAX
  write(6,*) "NMAX     =", NMAX
  write(6,*) "TIMESTEP =", TIMESTEP
!
  call legendre__ini                          &
   &( IMAX, JMAX, NNUM, MNUM,                 &!IN
   &  MNWAV, MNWAV_UV, MNSTART, MNSTART_UV,   &!OUT
   &  WEIGHT, alat, sinlat, coslat, coslat_inv )   !OUT
!
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
  do i=1,IMAX
     ALON(i) = DLON*(i-1)    !## 0 <= alon < 2*PI
     if ( ALON(i) >= PI ) then
        ALON1(i) = ALON(i)-2*PI  !## -PI <= alon1 < PI
     else
        ALON1(i) = ALON(i)       !## -PI <= alon1 < PI
     end if
     SINLON(i) = sin(ALON(i))
     COSLON(i) = cos(ALON(i))
  end do
 !$OMP END DO
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
     ACOSLAT_INV(j) = COSLAT_INV(j)/ER
     ACOS2LAT_INV(j) = COSLAT_INV(j)*COSLAT_INV(j)/ER
     ALAT_SL(j) = ALAT(j)
     YLAT(j) = ALAT(j)*180.0d0/PI
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
  ALAT_SL(0)      = PI/2.0d0
  ALAT_SL(JMAX+1) = -PI/2.0d0
  do j=1,6
    ALAT_SL(-j) = ALAT_SL(0) + ( ALAT_SL(0)-ALAT_SL(j) )
    ALAT_SL(JMAX+1+j) = -ALAT_SL(-j)
  end do

  call output_weight_lat( IMAX, JMAX, YLAT, WEIGHT )

!  write(6,*) 'alat=',alat
!  write(6,*) 'lat='
!  write(6,'(5F10.3)') -ylat
!  write(6,'(5F10.3)') -alat*180.0d0/PI
!xx  write(6,*) 'alat_sl=',alat_sl
!xx  write(6,*) 'sinlat=',sinlat
!xx  write(6,*) 'acos(sinlat)=',acos(sinlat)
!xx  write(6,*) 'coslat=',coslat
!xx  stop 234
!
!xx  NMAX=MNUM-1
  
  allocate( NTOTAL(MNWAV) )
  allocate( MWAVE(MNWAV) )
  allocate( ERFNN1(MNWAV) )

!$OMP PARALLEL default(SHARED), private(m,mn,an)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MMAX
    do mn=MNSTART(m),MNSTART(m+1)-1
      NTOTAL(mn) = m + mn - MNSTART(m)
      MWAVE(mn) = m
      an = NTOTAL(mn)
      ERFNN1(mn) = an*(an+1.0d0)/(ER*ER)  !## n(n+1)/(a*a)
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
  XGAUSS(1) = -0.94910791234275852453d0
  XGAUSS(2) = -0.74153118559939443986d0
  XGAUSS(3) = -0.40584515137739716691d0
  XGAUSS(4) =  0.0d0
  XGAUSS(5) =  0.40584515137739716691d0
  XGAUSS(6) =  0.74153118559939443986d0
  XGAUSS(7) =  0.94910791234275852453d0

  WGAUSS(1) = 0.1294849661688696932d0
  WGAUSS(2) = 0.27970539148927666793d0
  WGAUSS(3) = 0.38183005050511894494d0
  WGAUSS(4) = 0.41795918367346938776d0
  WGAUSS(5) = 0.38183005050511894494d0
  WGAUSS(6) = 0.27970539148927666793d0
  WGAUSS(7) = 0.1294849661688696932d0

!  do j=1,JMAX/2+1
!    hosei(j) = PI - 2.0d0*atan2(1.0d0+cos(DLON), sin(DLON)*SINLAT(j))
!    hosei(JMAX+1-j) = -hosei(j)
!  end do
!  hosei(0) = PI - 2.0d0*atan2(1.0d0+cos(DLON), sin(DLON))
!  hosei(JMAX+1) = -hosei(0)
!
!xx  write(6,*) 'hosei=',hosei
!xx  stop 234

!  write(6,*) "ddddddddd"

end subroutine initialize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine output_weight_lat( imax, jmax, ylat, weight )
  integer,intent(in) :: imax
  integer,intent(in) :: jmax
  real(8),intent(in) :: ylat(jmax)
  real(8),intent(in) :: weight(jmax)

  real(4) :: data_4byte(jmax)
  character(len=20) :: filename   !! Output file name
  character(len=30) :: fdata
  character(len=30) :: fctl
  integer :: j,iunit

  filename = 'weight_lat'
  fdata = trim(filename)//'.dr'
  fctl  = trim(filename)//'.ctl'

    iunit=51
    open(iunit,file=fctl,form='formatted',access='sequential')
!    write(iunit,'(a)') 'OPTIONS big_endian'
    write(iunit,'(a)') 'DSET ^'//fdata
    write(iunit,'(a)') 'TITLE '//filename
    write(iunit,'(a)') 'UNDEF -9.99E33'
    write(iunit,'(a,f14.10)') 'XDEF 1  LINEAR  0.0 ',360.0d0/IMAX    
    write(iunit,'(a,i6,a)') 'YDEF ',jmax,'  LEVELS'
    write(iunit,'( (5x,5(f14.10,1x)) )') -ylat
    write(iunit,'(a)') 'ZDEF 1     LEVELS  1000.0'
    write(iunit,'(a)') 'TDEF 1  LINEAR  00Z01JAN2000 1HR'
    write(iunit,'(a)') 'VARS  1'
    write(iunit,'(a)') 'wgt   0  0   latitudinal weight'
    write(iunit,'(a)') 'ENDVARS'
    close(iunit)
    !
    data_4byte(:) = weight(:)
    !
    open( iunit, file=fdata, form='unformatted', access='direct', &
     &    recl=4*jmax )
    write(iunit,rec=1) data_4byte
    close(iunit)

end subroutine output_weight_lat


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine main
!
  real(8) :: qphi(2,mnwav)               !! GRAV*height (geopotential height)
  real(8) :: qrot(2,mnwav)               !! Vorticity
  real(8) :: qdiv(2,mnwav)               !! Divergence
!
  real(8),save :: phi0(IMAX,JMAX)        !! GRAV*height
  real(8),save :: phi_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y) !! With halo region 
!
  real(8),save :: phis (IMAX,JMAX)       !! Surface geopotential height
  real(8),save :: phisx(IMAX,JMAX)
  real(8),save :: phisy(IMAX,JMAX)
!
  real(8),save :: u0(IMAX,JMAX)          !! Zonal wind
  real(8),save :: v0(IMAX,JMAX)          !! Meridional wind
  real(8),save :: u_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)  !! With halo region
  real(8),save :: v_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),save :: u_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),save :: v_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
!
  real(8),save :: um_adv(IMAX,JMAX)      !! Zonal wind for advection (past)
  real(8),save :: vm_adv(IMAX,JMAX)      !! Zonal wind for advection (past)
!
  real(8),save :: u0_adv(IMAX,JMAX)      !! Zonal wind for advection (present)
  real(8),save :: v0_adv(IMAX,JMAX)      !! Zonal wind for advection (present)
!
  real(8),save :: dudtm(IMAX,JMAX)       !! d(u)/dt (past)
  real(8),save :: dvdtm(IMAX,JMAX)       !! d(v)/dt (past)
  real(8),save :: dphidtm(IMAX,JMAX)     !! d(phi)/dt (past)

  real(8),save :: dudtm_li(IMAX,JMAX)    !! d(u)/dt linear term (past)
  real(8),save :: dvdtm_li(IMAX,JMAX)    !! d(v)/dt linear term (past)
  real(8),save :: dphidtm_li(IMAX,JMAX)  !! d(u)/dt linear term (past)
!
  real(8),save :: dudt0(IMAX,JMAX)       !! d(u)/dt (present)
  real(8),save :: dvdt0(IMAX,JMAX)       !! d(v)/dt (present)
  real(8),save :: dphidt0(IMAX,JMAX)     !! d(phi)/dt (present)

  real(8),save :: dudt0_li(IMAX,JMAX)    !! d(u)/dt linear term (present)
  real(8),save :: dvdt0_li(IMAX,JMAX)    !! d(v)/dt linear term (present)
  real(8),save :: dphidt0_li(IMAX,JMAX)  !! d(phi)/dt linear term (present)
!
  integer,save :: ii(IMAX,JMAX)          !! Grid position near departure point for semi-Lag.
  integer,save :: jj(IMAX,JMAX)          !! Grid position near departure point for semi-Lag.
  real(8),save :: xi(IMAX,JMAX)          !! Interpolation in longitudinal direction
  real(8),save :: yj(IMAX,JMAX)          !! Latitude of departure point
  real(8),save :: cosdtheta(IMAX,JMAX)   !! Rotation of the wind for semi-Lag.
  real(8),save :: sindtheta(IMAX,JMAX)   !! Rotation of the wind for semi-Lag.
!
  real(8),save :: f_8byte(IMAX,JMAX)
!
  real(8),save :: phi_zonal(JMAX)
  real(8),save :: u_zonal(JMAX)
  real(8),save :: work(JMAX)
  
  real(8),save :: work1(IMAX)
  real(8),save :: work2(IMAX)
  real(8),save :: work3(IMAX)
  real(8),save :: work4(IMAX)
!
  real(8) :: phibar
  real(8) :: dta, dtb, beta, betauv
  real(8) :: dt, dt_prev
  real(8) :: alon_c, alat_c, rr, r, phis0, hh0, uu0
  real(8) :: aa,bb,cc,uu
  real(8) :: thetac, ramdac, ramda, alpha0, cosa, sina
  real(8) :: alat0, alat1, alat2, umax, ss, ww, xx, en, f, hhat
  real(8) :: ww1, ww2
  real(8) :: gmean, gmean_phis, gmean_mass_0, gmean_energy_0, gmean_vorticity_0, gmean_divergence_0
!
  integer :: iterx, kt_end, ntimes_hour, iadv_only
  integer :: kind_phi, kind_phis, kind_uv, nstep, nstep_start, intv_monit
  integer :: istep, i, j, nn, kk, mn, idiv, i1, i2, ntmax
  integer :: idepart_prev
  integer :: istep_monit = 0
!
  real(8),parameter :: omg_6 = 7.848d-6
  real(8),parameter :: ak_6  = omg_6
  real(8),parameter :: h0_6  = 8.0d3
  real(8),parameter :: ir_6  = 4
!
  ! DFS interpolation declarations
  real(8), save, allocatable :: lld(:, :)
  real(8), save, allocatable :: ttd(:, :)
  real(8), allocatable:: udhalo(:,:),vdhalo(:,:),hdhalo(:,:),wdhalo(:,:)
  real(8), allocatable :: thd(:), lbd(:), th(:)
  complex*16, allocatable :: coeffs2d(:,:), coeffs1d(:)
  real(8), allocatable :: tmph(:),tmpu(:),tmpv(:),tmpw(:)
  integer, allocatable :: ipiv(:)
  complex*16, allocatable :: fu(:), L(:,:), cn(:,:)
  integer :: tmpindex, iret
  integer(8) :: N1, N2
  real(8) :: vij, uij, vv, uv
  integer :: nthreads
!
  integer :: ier, iflag, ntrans, dim, ttype,k,info
  integer(8) :: M, plan, planF
  integer(8) :: nmodes(3)
  real(8), pointer :: dummy => null()
  real(8) :: tol 
  complex*16, parameter :: img = dcmplx(0,1)
! 
! create the options struct
  type(finufft_opts), pointer :: opts => null()
!
!
  call e_time__start(2,"main")

  if ( JCN_INITIAL == 1 ) then
    iadv_only = 1   !! Advection
    kind_phi  = 1
    kind_phis = 1
    kind_uv   = 2
    kt_end    = 24*12  !! 12 days
  else if ( jcn_initial == 2 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 2
    kind_phis = 1
    kind_uv   = 2
    kt_end    = 24*5*24  !! 5 days
  else if ( jcn_initial == 5 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 5
    kind_phis = 2
    kind_uv   = 5
    kt_end    = 24*15  !! 15 days
  else if ( jcn_initial == 6 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 6
    kind_phis = 1
    kind_uv   = 6
    kt_end    = 24*14  !! 14 days
  else if ( jcn_initial == 10 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 10
    kind_phis = 1
    kind_uv   = 10
    kt_end    = 24*6*5   !! 6 days
  else if ( jcn_initial == 11 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 11
    kind_phis = 1
    kind_uv   = 11
    kt_end    = 24*6   !! 6 days
  else
    write(6,*) "Error: jcn_initial=",jcn_initial," is not supported."
    stop 999
  end if
!
  ntimes_hour = int( (3600.0d0+1.0d-10)/TIMESTEP )
  if ( abs( 3600.0d0 - TIMESTEP*ntimes_hour ) > 1.0d-10 ) then
    write(6,*) "Error: 3600.0/TIMESTEP should be integer."
    write(6,*) "       TIMESTEP =", TIMESTEP
    stop 999
  end if
!
  nstep_start = 1    !# nstep_start>=1
  nstep = kt_end*ntimes_hour + nstep_start
!
  intv_monit = INTHR_MONIT*ntimes_hour
!
  beta   = 1.00d0    !## 1.0 <= beta <= 1.2? 1.4?
  betauv = 1.00d0    !## 1.0 <= beta <= 1.2? 1.4?
!
! =====================================================================
!
! ####   Initial values   #####
!
  if ( kind_uv == 2 ) then
!    alpha0=0.0d0             !! Around the Equator
!    alpha0=0.05d0            !! (minor shift)
    alpha0=PI/2-0.05d0       !! (minor shift)
!    alpha0=PI/2              !! Over the Poles
  else
    alpha0 = 0.0d0
  end if

! ------------------------------
  if(kind_phi==0) then
! ------------------------------
    phibar = 1.0d5
   !$OMP PARALLEL default(SHARED), private(i,j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        phi0(i,j)=phibar
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
! ------------------------------
  else if(kind_phi==1) then
! ------------------------------
    ramdac=3.0d0*PI/2
!    ramdac=2.5*PI/2
    thetac=0.0d0
!    thetac=0.5*PI/2
    rr=ER/3
    phibar = 1100.0d0*grav
   !$OMP PARALLEL default(SHARED), private(i,j,ramda,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ramda=ALON(i)
        r = ER*acos( sin(thetac)*SINLAT(j) + cos(thetac)*COSLAT(j)*cos(ramda-ramdac) )
        if ( r < rr ) then
          phi0(i,j) = 1000.0d0*grav/2*( 1.0d0 + cos(PI*r/rr) )
        else
          phi0(i,j) = 0.0d0
        end if
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

! ------------------------------
  else if(kind_phi==2) then
! ------------------------------
    uu0=2.0d0*PI*ER/(12.0d0*24*60*60)  !! 2*PI*ER/(12days) = about 40m/s
    phibar = 3.2d4
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
       do i=1,IMAX
          phi0(i,j) = 2.94d4 - (ER*OMG*uu0 + 0.5d0*uu0*uu0)                              &
           &            *((-COSLON(i)*COSLAT(j)*sin(alpha0) + SINLAT(j)*cos(alpha0))**2.0d0 )
       end do
    end do
 !$OMP END DO
!$OMP END PARALLEL

! ------------------------------
  else if (kind_phi==5) then
! ------------------------------
    phibar = 6500.0d0*GRAV
    hh0 = 5960.0d0
    uu0 = 20.0d0
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        phi0(i,j) = GRAV*hh0 - ( ER*omg*uu0 + uu0*uu0/2 )*SINLAT(j)**2
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  else if (kind_phi==6) then 
! ------------------------------
    phibar = 11000.0d0*GRAV
!$OMP PARALLEL default(SHARED), private(i,j,aa,bb,cc)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
      aa = omg_6*0.5d0*( 2.0d0*omg + omg_6 )*COSLAT(j)**2   &
       &   + 0.25d0*ak_6**2*COSLAT(j)**(2*ir_6)             &
       &     *( (ir_6+1)*COSLAT(j)**2 + (2*ir_6**2-ir_6-2)  &
       &        - 2.0d0*ir_6**2/COSLAT(j)**2          )
      bb = 2.0d0*( omg + omg_6 )*ak_6/((ir_6+1)*(ir_6+2))*COSLAT(j)**ir_6 &
       &   *( (ir_6**2+2*ir_6+2) - (ir_6+1)**2*COSLAT(j)**2 )
      cc = 0.25d0*ak_6**2*COSLAT(j)**(2*ir_6)*( (ir_6+1)*COSLAT(j)**2 - (ir_6+2) )
      do i=1,IMAX
        phi0(i,j) = GRAV*h0_6 + ER**2*( aa + bb*cos(ALON(i)*ir_6)  &
         &                             + cc*cos(ALON(i)*ir_6*2) )
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  else if(kind_phi==10 .or. kind_phi==11) then
! ------------------------------
    umax = 80.0d0
    alat0 = PI/7.0d0
    alat1 = PI/2.0d0 - alat0
    en = exp( -4.0d0/(alat1-alat0)**2 )
    phi_zonal(JMAX) = 0.0d0  !! Near South Pole
!$OMP PARALLEL default(SHARED), private(i,j,aa,bb,nn,xx,uu,f,ss)
 !$OMP DO schedule(STATIC)
      do j=JMAX-1,1,-1
        aa = ( ALAT(j) + ALAT(j+1) )*0.5d0
        bb = ( ALAT(j) - ALAT(j+1) )*0.5d0   !! ALAT(j+1) < ALAT(j)
        work(j) = 0.0d0
        do nn=1,7
           xx = aa + bb*XGAUSS(nn)
           if ( alat0 < xx .and. xx < alat1 ) then
              uu = umax/en*exp( 1.0d0/((xx-alat0)*(xx-alat1)) )
           else
              uu = 0.0d0
           end if
           f  = 2.0d0*OMG*sin(xx)
           work(j) = work(j) + ER*uu*( f + tan(xx)*uu/ER )*WGAUSS(nn)*0.5d0
        end do
        work(j) = work(j)*( ALAT(j) - ALAT(j+1) )
     end do   
 !$OMP END DO
!$OMP END PARALLEL
     do j=JMAX-1,1,-1
        phi_zonal(j) = phi_zonal(j+1) - work(j) !! Integral
     end do
     ss = sum( phi_zonal(1:JMAX)*WEIGHT(1:JMAX) )/sum( WEIGHT(1:JMAX) )
!$OMP PARALLEL default(SHARED), private(j)
 !$OMP DO schedule(STATIC)
     do j=1,JMAX
        phi_zonal(j) = phi_zonal(j) - ss + GRAV*10000.0d0 !! global mean = 10000m
     end do
 !$OMP END DO
!$OMP END PARALLEL
     
     alat2 = PI/4.0d0
     aa = 1.0d0/3.0d0
     bb = 1.0d0/15.0d0
     hhat = 120.0d0
     if ( kind_phi == 11 ) then
       phibar = 22000.0d0*GRAV
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
       do j=1,JMAX
         do i=1,IMAX
           phi0(i,j) = phi_zonal(j) + phi_zonal(JMAX+1-j) &
            &  + GRAV*hhat*COSLAT(j)*exp(-(ALON1(i)/aa)**2)*exp(-((alat2-ALAT(j))/bb)**2) &
            &  + GRAV*hhat*COSLAT(j)*exp(-(ALON1(i)/aa)**2)*exp(-((-alat2-ALAT(j))/bb)**2)
         end do
       end do
 !$OMP END DO
!$OMP END PARALLEL
     else
       phibar = 11000.0d0*GRAV
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
       do j=1,JMAX
         do i=1,IMAX
           phi0(i,j) = phi_zonal(j) &
            &  + GRAV*hhat*COSLAT(j)*exp(-(ALON1(i)/aa)**2)*exp(-((alat2-ALAT(j))/bb)**2)
         end do
       end do
 !$OMP END DO
!$OMP END PARALLEL
     end if
! ------------------------------
  else
! ------------------------------
!$OMP PARALLEL default(SHARED), private(j)
 !$OMP DO schedule(STATIC)
     do j=1,JMAX
       phi0(:,j) = phibar
     end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  end if
! ------------------------------
!
!
! ------------------------------
  if(kind_phis==1) then
! ------------------------------
!$OMP PARALLEL default(SHARED), private(j)
 !$OMP DO schedule(STATIC)
     do j=1,JMAX
       phis(:,j) = 0.0d0
     end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  else
! ------------------------------
    alon_c = PI*3/2
    alat_c = PI/6
    rr     = PI/9
    phis0  = 2000.0d0*GRAV
!$OMP PARALLEL default(SHARED), private(i,j,r)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        r = min( sqrt( (ALON(i)-alon_c)**2 + (ALAT(j)-alat_c)**2 ), rr )
        phis(i,j) = phis0*(1.0d0-r/rr)
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  end if
! ------------------------------
!
! ------------------------------
  if (kind_uv.eq.0) then
! ------------------------------
    !$OMP PARALLEL default(SHARED), private(i,j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
           u0(i,j)=0.0d0
           v0(i,j)=0.0d0
        end do
     end do
    !$OMP END DO
    !$OMP END PARALLEL
! ------------------------------
  else if (kind_uv.eq.2) then
! ------------------------------
    uu0=2.0d0*PI*ER/(12.0d0*24*60*60)  !! 2*PI*ER/(12days) = about 40m/s
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        u0(i,j) = uu0*( COSLAT(j)*cos(alpha0) + SINLAT(j)*COSLON(i)*sin(alpha0) )
        v0(i,j) = -uu0*SINLON(i)*sin(alpha0)
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  else if (kind_uv.eq.5) then
! ------------------------------
    uu0 = 20.0d0
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        u0(i,j)=uu0*COSLAT(j)
        v0(i,j)=0.0d0
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  else if (kind_uv.eq.6) then
! ------------------------------
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        u0(i,j)=ER*( omg_6*COSLAT(j) + ak_6*COSLAT(j)**(ir_6-1)            &
         &                               *(ir_6*SINLAT(j)**2-COSLAT(j)**2)   &
         &                               *cos(ALON(i)*ir_6)               )
        v0(i,j)=-ER*ak_6*ir_6*COSLAT(j)**(ir_6-1)*SINLAT(j)*sin(ir_6*ALON(i))
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
! ------------------------------
  else if (kind_uv.eq.10 .or. kind_uv == 11) then
! ------------------------------
     umax = 80.0d0
     alat0 = PI/7.0d0
     alat1 = PI/2.0d0 - alat0
     en = exp( -4.0d0/(alat1-alat0)**2 )
!     write(6,*) "en=",en
     do j=1,JMAX
       if ( alat0 < ALAT(j) .and. ALAT(j) < alat1 ) then
         u_zonal(j) = umax/en*exp( 1.0d0/((ALAT(j)-alat0)*(ALAT(j)-alat1)) )
       else
         u_zonal(j) = 0.0d0
       end if
     end do
     if ( kind_uv == 11 ) then
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
       do j=1,JMAX
         do i=1,IMAX
           u0(i,j) = u_zonal(j) + u_zonal(JMAX+1-j)
           v0(i,j) = 0.0d0
         end do
       end do
 !$OMP END DO
!$OMP END PARALLEL
     else
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
       do j=1,JMAX
         do i=1,IMAX
           u0(i,j) = u_zonal(j)
           v0(i,j) = 0.0d0
         end do
       end do
 !$OMP END DO
!$OMP END PARALLEL
     end if
! ------------------------------
  else
! ------------------------------
    write(6,*) "Error: kind_uv=",kind_uv," is not supported."
    stop 999
! ------------------------------
  end if
! ------------------------------

  cosa = cos(alpha0)
  sina = sin(alpha0)
!

  if ( iadv_only == 1 ) then
   !$OMP PARALLEL default(SHARED), private(i,j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        u0_adv(i,j) = u0(i,j)
        v0_adv(i,j) = v0(i,j)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    call grid_to_wave_to_grid &
     &( phi0,       &!INOUT
     &  qphi  )      !WORK

  else
    if ( JCN_PROG_H == 0 ) then
     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          phi0(i,j) = phi0(i,j) - phis(i,j)
        end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    end if

    call fft_run_towave &
     &( IMAX, JMAX,     &!IN
     &  phis        )    !OUT

    call legendre__g2w &
     &( phis,          &!IN
     &  qphi  )         !OUT

    call legendre__w2g &
     &( qphi,          &!IN
     &  phis    )       !OUT
    call legendre__w2g_dy &
     &( qphi,        &!IN
     &  phisy      )  !OUT
    call xderiv_run        &
     &( IMAX, JMAX, MNUM, &!IN
     &  phis,              &!IN
     &  phisx    )          !OUT

    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  phis      )      !INOUT
    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  phisx       )    !INOUT
    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  phisy         )  !INOUT
   
    call grid_to_wave      &
     &( u0, v0, phi0,      &!IN
     &  qrot, qdiv, qphi )  !OUT

    call tendency                           &
     &( qrot, qdiv, qphi,                   &!IN
     &  phis, phisx, phisy,                 &!IN
     &  phibar,                             &!IN
     &  u0, v0, phi0,                       &!OUT
     &  um_adv, vm_adv,                     &!OUT
     &  dudtm, dvdtm, dphidtm,              &!OUT
     &  dudtm_li, dvdtm_li, dphidtm_li )     !OUT
  end if

  if ( JCN_MONIT == 1 ) then
    if ( JCN_PROG_H == 0 ) then
     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          f_8byte(i,j) = ( phi0(i,j) + phis(i,j) )/GRAV
        end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    else
     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          f_8byte(i,j) = phi0(i,j)/GRAV
        end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    end if

    ntmax = kt_end/INTHR_MONIT + 1
    call monit_norm__ini( "norm", INTHR_MONIT, ntmax, IMAX, JMAX, WEIGHT, f_8byte )
    call monit__ini( INTHR_MONIT, ntmax, IMAX, JMAX, YLAT ) !IN
    if ( JCN_MONIT_SPECTRUM == 1 ) then
      call monit_spectrum__ini( "spectrum", INTHR_MONIT, ntmax, MNWAV, NMAX, ER )
    end if
  end if
!
! ===================================================================

  dt = -999.0d0
!
  if(IT .eq. 1) then
    N1 = JMAX * 2
    N2 = IMAX
    allocate(thd(JMAX*IMAX),lbd(JMAX*IMAX),th(N1),tmpv(IMAX*JMAX))
    allocate(udhalo(JMAX,IMAX),vdhalo(JMAX, IMAX),hdhalo(JMAX,IMAX),wdhalo(JMAX, IMAX))
    allocate(coeffs2d(N1,N2),coeffs1d(N1*N2),tmph(IMAX*JMAX),tmpu(IMAX*JMAX),tmpw(IMAX*JMAX))
    allocate(ipiv(N1),L(N1,N2),cn(N1,N2),fu(IMAX*JMAX))
    allocate(lld(IMAX,JMAX), ttd(IMAX,JMAX))
!
!    set mandotary parameters FINUFFT guru interface
    iflag = 1
    tol = 1d-14
    ntrans = 1
    ttype = 2
    dim = 2
    !allocate(nmodes(3))
    nmodes(1) = N1
    nmodes(2) = N2
    M = JMAX*IMAX
    nthreads=omp_get_max_threads()
  
  !
  ! FFTW3 shifting coefficients
    call shiftfftw(JMAX, IMAX, cn)
  !
  ! coeffients FourierS
    th(1:JMAX) =  ALAT !(JMAX:1:-1)
    th(JMAX+1:2*JMAX) = alat - PI
    do j=1, 2*JMAX
      do k =-JMAX, JMAX-1
        i = k + JMAX + 1
        L(j,i) = cdexp(img*k*th(j))
      end do
    end do
    call zgetrf(N1, N1, L, N1, ipiv, info)
  end if
  
!
  do istep = 1, nstep

    call e_time__start(3,"main_cp_past")
    
    dt_prev = dt

    if ( istep == 1 ) then
      dt = TIMESTEP/2**(nstep_start-1)
      istep_monit = 0
    else if ( istep <= nstep_start ) then
      dt = TIMESTEP/2**(nstep_start-istep+1)
      istep_monit = -1
    else if ( istep == nstep_start + 1 ) then
      dt = TIMESTEP
      istep_monit = 1
    else
      dt = TIMESTEP
      istep_monit = istep_monit + 1

     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          um_adv(i,j)  = u0_adv(i,j)
          vm_adv(i,j)  = v0_adv(i,j)
          dudtm(i,j)   = dudt0(i,j)
          dvdtm(i,j)   = dvdt0(i,j)
          dphidtm(i,j) = dphidt0(i,j)
          dudtm_li(i,j)   = dudt0_li(i,j)
          dvdtm_li(i,j)   = dvdt0_li(i,j)
          dphidtm_li(i,j) = dphidt0_li(i,j)
        end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    end if
!
    dta = dt*0.5d0
    dtb = dt*0.5d0
!
    call e_time__end(3,"main_cp_past")
    
!    write(6,*) 'n       = ',n
!    write(6,*) 'istep_monit = ',istep_monit
!
!   ------------------------------------------------------
!
    if ( iadv_only == 1 ) then
      call grid_to_wave_to_grid &
       &( phi0,       &!INOUT
       &  qphi  )      !WORK

    else
      call tendency                           &
       &( qrot, qdiv, qphi,                   &!IN
       &  phis, phisx, phisy,                 &!IN
       &  phibar,                             &!IN
       &  u0, v0, phi0,                       &!OUT
       &  u0_adv, v0_adv,                     &!OUT
       &  dudt0, dvdt0, dphidt0,              &!OUT
       &  dudt0_li, dvdt0_li, dphidt0_li )     !OUT
    end if
!
!   -------------------------------------------------------
!
    if ( JCN_MONIT == 1 .and. istep_monit >= 0 .and. mod( istep_monit, intv_monit ) == 0 ) then
!
      call e_time__start(15,"output")
!
      nn = istep_monit/intv_monit + 1
      write(6,*) "nn = ", nn
!
      if ( JCN_PROG_H == 0 ) then
        !$OMP PARALLEL default(SHARED), private(i,j)
        !$OMP DO schedule(STATIC)
         do j=1,JMAX
            do i=1,IMAX
               f_8byte(i,j) = ( phi0(i,j) + phis(i,j) )/GRAV
            end do
         end do
        !$OMP END DO
        !$OMP END PARALLEL
      else
        !$OMP PARALLEL default(SHARED), private(i,j)
        !$OMP DO schedule(STATIC)
         do j=1,JMAX
            do i=1,IMAX
               f_8byte(i,j) = phi0(i,j)/GRAV
            end do
         end do
        !$OMP END DO
        !$OMP END PARALLEL
      end if
!
      if ( JCN_PROG_H == 0 ) then
        !$OMP PARALLEL default(SHARED), private(i,j)
        !$OMP DO schedule(STATIC)
         do j=1,JMAX
            do i=1,IMAX
               f_8byte(i,j) = ( phi0(i,j) + phis(i,j) )/GRAV
            end do
         end do
        !$OMP END DO
        !$OMP END PARALLEL
      else
        !$OMP PARALLEL default(SHARED), private(i,j)
        !$OMP DO schedule(STATIC)
         do j=1,JMAX
            do i=1,IMAX
               f_8byte(i,j) = phi0(i,j)/GRAV
            end do
         end do
        !$OMP END DO
        !$OMP END PARALLEL
      end if
      call monit__output( 1, nn, f_8byte )
      call monit__output( 2, nn, u0 )
      call monit__output( 3, nn, v0 )
      call monit_norm__output( nn, f_8byte )
      
      call cal_gmean &
       &( f_8byte,   &
       &  gmean)      !! Global mean of h
      if ( nn == 1 ) then
         call cal_gmean &
          &( phis,      &
          &  gmean_phis )
         gmean_mass_0 = gmean - gmean_phis/GRAV
      end if
      gmean = gmean - gmean_phis/GRAV  !! Global mean of mass
      write(6,*) "Global mean of mass = ",gmean, "(", (gmean-gmean_mass_0)/gmean_mass_0, ")"

     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
         do i=1,IMAX
            f_8byte(i,j) = 0.5d0*( ( f_8byte(i,j) - phis(i,j)/GRAV )*( u0(i,j)**2 + v0(i,j)**2 ) &
             &                     + GRAV*( f_8byte(i,j)**2 - (phis(i,j)/GRAV)**2 ) )
         end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL

      call cal_gmean &
       &( f_8byte,   &
       &  gmean)
      if ( nn == 1 ) then
         gmean_energy_0 = gmean
      end if
      write(6,*) "Global mean of energy = ",gmean, "(", (gmean-gmean_energy_0)/gmean_energy_0, ")"

      if ( JCN_INITIAL >= 2 ) then
      
        if ( JCN_MONIT_SPECTRUM == 1 ) then
          call monit_spectrum__output( nn, NTOTAL, MWAVE, MNSTART, qrot, qdiv )
        end if
!
        call legendre__w2g &
         &( qrot,          &!IN
         &  f_8byte    )    !OUT
        call fft_run_togrid      &
         &( IMAX, JMAX,          &!IN
         &  f_8byte  )            !INOUT
        call monit__output( 4, nn, f_8byte )

        call cal_gmean &
         &( f_8byte,   &
         &  gmean)
        if ( nn == 1 ) then
           gmean_vorticity_0 = gmean
        end if
        write(6,*) "Global mean of vorticity = ",gmean, "(", (gmean-gmean_vorticity_0)/gmean_vorticity_0, ")"
!
        call legendre__w2g &
         &( qdiv,          &!IN
         &  f_8byte    )    !OUT
        call fft_run_togrid      &
         &( IMAX, JMAX,          &!IN
         &  f_8byte  )            !INOUT
        call monit__output( 5, nn, f_8byte )

!        call cal_gmean &
!         &( f_8byte,   &
!         &  gmean)
!        if ( nn == 1 ) then
!           gmean_divergence_0 = gmean
!        end if
!        write(6,*) "Global mean of divergence = ",gmean, "(", (gmean-gmean_divergence_0)/gmean_divergence_0, ")"
      end if

      call e_time__end(15,"output")
!
    end if


    call e_time__start(4,"main_integ1")

    if ( iadv_only == 1 ) then

     !$OMP PARALLEL default(SHARED), private(j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          u_adv_halo(i,j) = u0_adv(i,j)
          v_adv_halo(i,j) = v0_adv(i,j)
        end do
        do i=1,IMAX
          u_halo(i,j) = u0(i,j)
          v_halo(i,j) = v0(i,j)
          phi_halo(i,j) = phi0(i,j)
        end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL

    else
    
      if ( JCN_DEPARTURE == 0 ) then
       !$OMP PARALLEL default(SHARED), private(i,j)
       !$OMP DO schedule(STATIC)
        do j=1,JMAX
          do i=1,IMAX
            u_adv_halo(i,j) = -um_adv(i,j) + u0_adv(i,j)*2.0d0
            v_adv_halo(i,j) = -vm_adv(i,j) + v0_adv(i,j)*2.0d0
          end do
          do i=1,IMAX          
            u_halo(i,j) = u0(i,j) + ( -dudtm(i,j) + dudt0(i,j)*2.0d0 )*dtb + ( dudtm_li(i,j) - dudt0_li(i,j) )*betauv*dtb &
             &            + 2*omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
            v_halo(i,j) = v0(i,j) + ( -dvdtm(i,j) + dvdt0(i,j)*2.0d0 )*dtb + ( dvdtm_li(i,j) - dvdt0_li(i,j) )*betauv*dtb &
             &            - 2*omg*ER*sina*SINLON(i)
            phi_halo(i,j) = phi0(i,j) + ( -dphidtm(i,j) + dphidt0(i,j)*2.0d0 )*dtb + ( dphidtm_li(i,j) - dphidt0_li(i,j) )*beta*dtb
          end do
        end do
       !$OMP END DO
       !$OMP END PARALLEL

      else
       !$OMP PARALLEL default(SHARED), private(i,j,ww1,ww2)
       !$OMP DO schedule(STATIC)
        do j=1,JMAX
          do i=1,IMAX
            ww1 = ( -dudtm(i,j) + dudt0(i,j)*2.0d0 )*dtb + 2.0d0*omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
            u_adv_halo(i,j) = u0(i,j) + ww1*0.5d0
            u_halo(i,j) = u0(i,j) + ( dudtm_li(i,j) - dudt0_li(i,j) )*betauv*dtb + ww1
            ww2 = ( -dvdtm(i,j) + dvdt0(i,j)*2.0d0 )*dtb - 2.0d0*omg*ER*sina*SINLON(i)
            v_adv_halo(i,j) = v0(i,j) + ww2*0.5d0
            v_halo(i,j) = v0(i,j) + ( dvdtm_li(i,j) - dvdt0_li(i,j) )*betauv*dtb + ww2
            phi_halo(i,j) = phi0(i,j) + ( -dphidtm(i,j) + dphidt0(i,j)*2.0d0 )*dtb + ( dphidtm_li(i,j) - dphidt0_li(i,j) )*beta*dtb
          end do
        end do
       !$OMP END DO
       !$OMP END PARALLEL

      end if
    end if
    if (IT .eq. 0) then
    call set_halo_uv &
     &( u_adv_halo )  !INOUT
    call set_halo_uv &
     &( v_adv_halo )  !INOUT
     
    call set_halo_uv &
     &( u_halo )  !OUT
    call set_halo_uv &
     &( v_halo )  !OUT

    call set_halo  &
     &( phi_halo )  !OUT
!
    call e_time__end(4,"main_integ1")
    
    call e_time__start(5,"main_semilag")
     
    if ( dt_prev == dt ) then
      iterx = 3         !! Number of iterations
      idepart_prev = 1  !! Use departure point data of previous step
    else
      iterx = 4         !! Number of iterations
      idepart_prev = 0  !! Not use departure point data of previous step
    end if
    
   !$OMP PARALLEL default(SHARED), private( j, idiv, i1, i2, work1, work2, work3, work4 )
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do idiv = 1, NDIV
       !! i is from i1 to i2.
        i1 = 1 + (idiv-1)*IVECLEN
        i2 = min( IMAX, idiv*IVECLEN )
        
        call departure                       &
         &( j, i1, i2, iadv_only,            &!IN
         &  iterx, idepart_prev,             &!IN
         &  dt, dta,                         &!IN
         &  cosa, sina,                      &!IN
         &  u0_adv(1,j), v0_adv(1,j),        &!IN
         &  dudt0(1,j), dvdt0(1,j),          &!IN
         &  u_adv_halo, v_adv_halo,          &!IN
         &  work1, work2, work3, work4,      &!INOUT
         &  ii(1,j), jj(1,j),                &!INOUT
         &  xi(1,j), yj(1,j),                &!INOUT
         &  ttd(1,j), lld(1,j),                &!INOUT
         &  cosdtheta(1,j), sindtheta(1,j) )  !INOUT

        call lag3             &
         &( i1, i2,           &
         &  ii(1,j), jj(1,j), &!IN
         &  xi(1,j), yj(1,j), &!IN
         &  phi_halo,         &!IN
         &  phi0(1,j)   )      !INOUT
         
        call lag5_uv                        &
         &( i1, i2,                         &!IN
         &  ii(1,j), jj(1,j),               &!IN
         &  xi(1,j), yj(1,j),               &!IN
         &  cosdtheta(1,j), sindtheta(1,j), &!IN
         &  u_halo, v_halo,                 &!IN
         &  u0(1,j), v0(1,j)      )          !INOUT
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
  else if(IT .eq. 1) then
    !$OMP PARALLEL DO private(i,j)
    do j = 1, JMAX
      do i = 1, IMAX
        udhalo(j,i) = -SINLON(i)*u_halo(i,j) -  COSLON(i)*SINLAT(j)*v_halo(i,j)
        vdhalo(j,i) = COSLON(i)*u_halo(i,j) - SINLON(i)*SINLAT(j)*v_halo(i,j)
        wdhalo(j,i) = COSLAT(j)*v_halo(i,j)
        hdhalo(j,i) = phi_halo(i,j) 
      end  do
    end do
  !$OMP END PARALLEL DO
  !  
    call set_halo_uv &
     &( u_adv_halo )  !INOUT
    call set_halo_uv &
     &( v_adv_halo )  !INOUT
     
    call set_halo_uv &
     &( u_halo )  !OUT
    call set_halo_uv &
     &( v_halo )  !OUT
    call set_halo  &
     &( phi_halo )  !OUT
!
    call e_time__end(4,"main_integ1")
    
    call e_time__start(5,"main_semilag")
     
    if ( dt_prev == dt ) then
      iterx = 3         !! Number of iterations
      idepart_prev = 1  !! Use departure point data of previous step
    else
      iterx = 4         !! Number of iterations
      idepart_prev = 0  !! Not use departure point data of previous step
    end if
    
   !$OMP PARALLEL default(SHARED), private( j, idiv, i1, i2, work1, work2, work3, work4 )
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do idiv = 1, NDIV
       !! i is from i1 to i2.
        i1 = 1 + (idiv-1)*IVECLEN
        i2 = min( IMAX, idiv*IVECLEN )
        
        call departure                       &
         &( j, i1, i2, iadv_only,            &!IN
         &  iterx, idepart_prev,             &!IN
         &  dt, dta,                         &!IN
         &  cosa, sina,                      &!IN
         &  u0_adv(1,j), v0_adv(1,j),        &!IN
         &  dudt0(1,j), dvdt0(1,j),          &!IN
         &  u_adv_halo, v_adv_halo,          &!IN
         &  work1, work2, work3, work4,      &!INOUT
         &  ii(1,j), jj(1,j),                &!INOUT
         &  xi(1,j), yj(1,j),                &!INOUT
         &  ttd(1,j), lld(1,j),                &!INOUT
         &  cosdtheta(1,j), sindtheta(1,j) )  !INOUT

        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
!
  !  thd(1:M) = reshape(ttd, [IMAX*JMAX])
  !  lbd(1:M) = reshape(lld, [IMAX*JMAX])
    !$OMP PARALLEL default(SHARED), private(i,j,tmpindex)
    !$OMP DO schedule(STATIC)
    do j=1, JMAX
      do i=1, IMAX
        tmpindex = (j-1)*IMAX + i
        thd(tmpindex) = ttd(i,j)
        lbd(tmpindex) = lld(i,j)
      end do
    end do
  !$OMP END DO
   !$OMP END PARALLEL
!
!   interpolate u
    ! Set FINUFFT Plan
    call finufft_makeplan(ttype,dim,nmodes,iflag,ntrans,&
    & tol,plan,opts,ier)

    call finufft_setpts(plan,M,thd,lbd,dummy,dummy,&
      &     dummy,dummy,dummy,ier)
    call doubleUp(udhalo, coeffs2d)
    call lgcoeffs(coeffs2d, cn, nthreads, ipiv, L)
    coeffs1d = reshape(coeffs2d, [N1*N2])
    
    call finufft_execute(plan,fu,coeffs1d,ier)
    tmpu = dble(fu)
   
    
!
! interpolate v
    call doubleUp(vdhalo, coeffs2d)
    call lgcoeffs(coeffs2d, cn, nthreads, ipiv, L)
    coeffs1d = reshape(coeffs2d, [N1*N2])
    call finufft_execute(plan,fu,coeffs1d,ier)
    tmpv = dble(fu)
!
! interpolate w
    call doubleUp(wdhalo, coeffs2d)
    call lgcoeffs(coeffs2d, cn, nthreads, ipiv, L)
    coeffs1d = reshape(coeffs2d, [N1*N2])
    call finufft_execute(plan,fu,coeffs1d,ier)
    tmpw = dble(fu)
!
! interpolate h
    call doubleUp(hdhalo, coeffs2d)
    call lgcoeffs(coeffs2d, cn, nthreads, ipiv, L)
    coeffs1d = reshape(coeffs2d, [N1*N2])
    call finufft_execute(plan,fu,coeffs1d,ier)
    tmph = dble(fu)
    call finufft_destroy(plan, ier)
!
! move the quantities back to the grid
    !$OMP PARALLEL DO private(i,j,tmpindex, uij, vij)
    do j=1, JMAX
      do i=1, IMAX
        tmpindex = (j-1)*IMAX + i
        uij = COS(lbd(tmpindex))*tmpv(tmpindex)-SIN(lbd(tmpindex))*tmpu(tmpindex)
        vij = COS(thd(tmpindex))*tmpw(tmpindex)-SIN(lbd(tmpindex))*SIN(thd(tmpindex))*tmpv(tmpindex)& 
 &         -COS(lbd(tmpindex))*sin(thd(tmpindex))*tmpu(tmpindex)
        u0(i,j) = cosdtheta(i,j)*uij - sindtheta(i,j)*vij
        v0(i,j) = sindtheta(i,j)*uij  + cosdtheta(i,j)*vij
      
        phi0(i,j) = tmph(tmpindex)
      end do
    end do
  !$OMP END PARALLEL DO
  end if
    call e_time__end(5,"main_semilag")
   
!
!   ===========================================================
    call e_time__start(6,"main_integ2")
    !
        if ( iadv_only /= 1 ) then
    !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
          do j=1,JMAX
            do i=1,IMAX
              u0(i,j) = u0(i,j) + dudt0(i,j)*dta - betauv*dudt0_li(i,j)*dta &
               &       - 2*omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
              v0(i,j) = v0(i,j) + dvdt0(i,j)*dta - betauv*dvdt0_li(i,j)*dta &
               &       + 2*omg*ER*sina*SINLON(i)
              phi0(i,j) = phi0(i,j) + dphidt0(i,j)*dta - beta*dphidt0_li(i,j)*dta
            end do
          end do
     !$OMP END DO
    !$OMP END PARALLEL
    !
          call grid_to_wave      &
           &( u0, v0, phi0,      &!IN
           &  qrot, qdiv, qphi )  !OUT
    !  
    !$OMP PARALLEL default(SHARED), private(mn,kk)
     !$OMP DO schedule(STATIC)
          do mn = 1,mnwav
            do kk=1,2
              qphi(kk,mn) = ( qphi(kk,mn) - beta*dta*phibar*qdiv(kk,mn) )   &
               &            /( 1.0d0 + ERFNN1(mn)*beta*betauv*phibar*dta**2 )
              qdiv(kk,mn) = betauv*dta*ERFNN1(mn)*qphi(kk,mn) + qdiv(kk,mn)
            end do
          end do
     !$OMP END DO
    !$OMP END PARALLEL
        end if
      
        call e_time__end(6,"main_integ2")
    
        call e_time__start(7,"main_hdiff")
    
        if ( JCN_HDIFF == 1) then
          call hdiff2            &
           &( dt,                &
           &  qrot, qdiv, qphi  )
        else if ( JCN_HDIFF == 2 ) then
          call hdiff4            &
           &( dt,                &
           &  qrot, qdiv, qphi  )
        end if
    
        call e_time__end(7,"main_hdiff")
    !
    !   ------------------------------------------------------
    !
    !
    !    write(6,*) 'f=',f
    !    write(6,*) 'df=',df
    !
    !    stop
    !
      end do
      
      call e_time__end(2,"main")
    !
  if(IT .eq. 1) then
    deallocate(thd,lbd,lld,ttd,tmph,tmpu,tmpv,tmpw)
    deallocate(coeffs1d,coeffs2d,ipiv,L,cn,fu)
    deallocate(udhalo,vdhalo,wdhalo,hdhalo)
  end if
end subroutine main


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine tendency                     &
 &( qrot, qdiv, qphi,                   &!IN
 &  phis, phisx, phisy,                 &!IN
 &  phibar,                             &!IN
 &  u0, v0, phi0,                       &!IN
 &  u0_adv, v0_adv,                     &!OUT
 &  dudt0, dvdt0, dphidt0,              &!OUT
 &  dudt0_li, dvdt0_li, dphidt0_li )     !OUT
!
  real(8),intent(in) :: qrot(2,mnwav)
  real(8),intent(in) :: qdiv(2,mnwav)
  real(8),intent(in) :: qphi(2,mnwav)
  real(8),intent(in) :: phis (IMAX,JMAX)
  real(8),intent(in) :: phisx(IMAX,JMAX)
  real(8),intent(in) :: phisy(IMAX,JMAX)
  real(8),intent(in) :: phibar
  real(8),intent(out) :: u0(IMAX,JMAX)
  real(8),intent(out) :: v0(IMAX,JMAX)
  real(8),intent(out) :: phi0(IMAX,JMAX)
  real(8),intent(out) :: u0_adv(IMAX,JMAX)
  real(8),intent(out) :: v0_adv(IMAX,JMAX)
  real(8),intent(out) :: dudt0(IMAX,JMAX)
  real(8),intent(out) :: dvdt0(IMAX,JMAX)
  real(8),intent(out) :: dphidt0(IMAX,JMAX)
  real(8),intent(out) :: dudt0_li(IMAX,JMAX)
  real(8),intent(out) :: dvdt0_li(IMAX,JMAX)
  real(8),intent(out) :: dphidt0_li(IMAX,JMAX)
!
  real(8) :: div(IMAX,JMAX)
  real(8) :: phi0x(IMAX,JMAX)
  real(8) :: phi0y(IMAX,JMAX)
  
  real(8) :: qu(2,mnwav_uv)
  real(8) :: qv(2,mnwav_uv)
!
  integer :: i,j
!
! ===========================================================

  call e_time__start(11,"tendency")
!
  call rotdiv2uv_run     &!##  qrot,qdiv -> qu, qv
   &( mnwav, MNUM,   &!IN
   &  mnstart, mnstart_uv, &!IN
   &  ER,             &!IN
   &  qrot, qdiv,        &!IN
   &  qu, qv          )   !OUT
!
  call legendre__w2g_uv   &!##  qu (Q) -> um (W)
   &( qu,      &!IN
   &  u0     )  !OUT
!
  call legendre__w2g_uv &!##  qu (Q) -> um (W)
   &( qv,        &!IN
   &  v0     )    !OUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  u0           )   !INOUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  v0           )   !INOUT
!
  call legendre__w2g &
   &( qdiv,          &!IN
   &  div        )    !OUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  div          )   !INOUT
!
  call legendre__w2g &
   &( qphi,          &!IN
   &  phi0        )   !OUT
!
  call xderiv_run        &
   &( IMAX, JMAX, MNUM, &!IN
   &  phi0,              &!IN
   &  phi0x      )        !OUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  phi0        )    !INOUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  phi0x      )    !INOUT
   
  call legendre__w2g_dy &
   &( qphi,             &!IN
   &  phi0y       )      !OUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  phi0y       )    !INOUT
!
  if ( JCN_PROG_H == 0 ) then
   !$OMP PARALLEL default(SHARED), private(i,j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
           dudt0(i,j) = -(phi0x(i,j) + phisx(i,j))*ACOSLAT_INV(j)
           dvdt0(i,j) = -(phi0y(i,j) + phisy(i,j))*ACOSLAT_INV(j)
           dudt0_li(i,j) = -phi0x(i,j)*COSLAT_INV(j)/ER
           dvdt0_li(i,j) = -phi0y(i,j)*COSLAT_INV(j)/ER
           dphidt0(i,j)  = -phi0(i,j)*div(i,j)
           dphidt0_li(i,j) = -phibar*div(i,j)
           u0(i,j) = u0(i,j)*COSLAT_INV(j)
           v0(i,j) = v0(i,j)*COSLAT_INV(j)
           u0_adv(i,j) = u0(i,j)
           v0_adv(i,j) = v0(i,j)
        end do
     end do
    !$OMP END DO
   !$OMP END PARALLEL

  else
   !$OMP PARALLEL default(SHARED), private(i,j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
           dudt0(i,j) = -phi0x(i,j)*ACOSLAT_INV(j)
           dvdt0(i,j) = -phi0y(i,j)*ACOSLAT_INV(j)
           dudt0_li(i,j) = dudt0(i,j)
           dvdt0_li(i,j) = dvdt0(i,j)
           dphidt0(i,j) = - ( phi0(i,j)-phis(i,j) )*div(i,j)            &
            &             + ( u0(i,j)*phisx(i,j) + v0(i,j)*phisy(i,j) ) &
            &               *ACOS2LAT_INV(j)
           dphidt0_li(i,j) = -phibar*div(i,j)
           u0(i,j) = u0(i,j)*COSLAT_INV(j)
           v0(i,j) = v0(i,j)*COSLAT_INV(j)
           u0_adv(i,j) = u0(i,j)
           v0_adv(i,j) = v0(i,j)
        end do
     end do
    !$OMP END DO
   !$OMP END PARALLEL

  end if
!
!
!xx     if ( istep_monit >= 1 ) then
!xx
!xx        f_4byte(1:IMAX,1:JMAX) = dudtm(1:IMAX,JMAX:1:-1)
!xx        write(12,rec=istep_monit*3-2) f_4byte
!xx
!xx        f_4byte(1:IMAX,1:JMAX) = dvdtm(1:IMAX,JMAX:1:-1)
!xx        write(12,rec=istep_monit*3-1) f_4byte
!xx
!xx        f_4byte(1:IMAX,1:JMAX) = dphidtm(1:IMAX,JMAX:1:-1)
!xx        write(12,rec=istep_monit*3  ) f_4byte
!xx
!xx      end if

  call e_time__end(11,"tendency")
!
end subroutine tendency


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine grid_to_wave &
 &( u, v, phi,         &
 &  qrot, qdiv, qphi  )
!
  real(8),intent(inout) :: u(IMAX,JMAX)
  real(8),intent(inout) :: v(IMAX,JMAX)
  real(8),intent(inout) :: phi(IMAX,JMAX)
  real(8),intent(out) :: qrot(2,mnwav)
  real(8),intent(out) :: qdiv(2,mnwav)
  real(8),intent(out) :: qphi(2,mnwav)
  
  real(8) :: qu(2,mnwav_uv)
  real(8) :: qv(2,mnwav_uv)
  integer :: i,j
!
! ==========================================================

  call e_time__start(12,"grid_to_wave")
!
!$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    do i=1,IMAX
     !###   1/(a*cos(lat)) * u   ###
      u(i,j) = ACOSLAT_INV(j)*u(i,j)
     !###   1/(a*cos(lat)) * v   ###
      v(i,j) = ACOSLAT_INV(j)*v(i,j)
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
  call fft_run_towave &
   &( IMAX, JMAX,     &!IN
   &  u       )        !OUT
!
  call fft_run_towave &
   &( IMAX, JMAX,     &!IN
   &  v       )        !OUT
!
  call legendre__g2w_uv &
   &( u,     &!IN
   &  qu   )  !OUT
!
  call legendre__g2w_uv   &
   &( v,     &!IN
   &  qv   )  !OUT
!
  call uv2rotdiv_run       &
   &( mnwav, MNUM,        &!IN
   &  mnstart, mnstart_uv, &!IN
   &  qu, qv,              &!IN
   &  qrot, qdiv   )        !OUT 1/(a*cos2(lat))*(U,V) -> rot, div
!
! ----------------------------------------------------
!
  call fft_run_towave &
   &( IMAX, JMAX,     &!IN
   &  phi         )    !OUT
!
  call legendre__g2w &
   &( phi,           &!IN
   &  qphi  )         !OUT

  call e_time__end(12,"grid_to_wave")
!
end subroutine grid_to_wave


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine grid_to_wave_to_grid &
 &( phi,         &
 &  qphi  )
!
  real(8),intent(inout) :: phi(IMAX,JMAX)
  real(8),intent(out) :: qphi(2,mnwav)

  call e_time__start(12,"grid_to_wave")
!
  call fft_run_towave &
   &( IMAX, JMAX,     &!IN
   &  phi         )    !OUT
  call legendre__g2w &
   &( phi,           &!IN
   &  qphi  )         !OUT
!
  call legendre__w2g &
   &( qphi,          &!IN
   &  phi   )         !OUT
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  phi         )    !OUT

  call e_time__end(12,"grid_to_wave")
!
end subroutine grid_to_wave_to_grid


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine hdiff2  &
 &( dt,            &
 &  qrot, qdiv, qphi  )
!
  real(8),intent(in) :: dt
  real(8),intent(inout) :: qrot(2,mnwav)
  real(8),intent(inout) :: qdiv(2,mnwav)
  real(8),intent(inout) :: qphi(2,mnwav)
!
  real(8) :: nyu
  integer :: l, k
!
! ==========================================================
!
  nyu = 1.0d5
!
!$OMP PARALLEL default(SHARED), private(l,k)
 !$OMP DO schedule(STATIC)
  do l=2,mnwav
    do k=1,2
      qrot(k,l) = qrot(k,l)/( 1.0d0 + nyu*dt*ERFNN1(l) )
      qdiv(k,l) = qdiv(k,l)/( 1.0d0 + nyu*dt*ERFNN1(l) )
      qphi(k,l) = qphi(k,l)/( 1.0d0 + nyu*dt*ERFNN1(l) )
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
end subroutine hdiff2



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine hdiff4  &
 &( dt,            &
 &  qrot, qdiv, qphi  )
!
  real(8),intent(in) :: dt
  real(8),intent(inout) :: qrot(2,mnwav)
  real(8),intent(inout) :: qdiv(2,mnwav)
  real(8),intent(inout) :: qphi(2,mnwav)
!
  real(8) :: tau, coef
  integer :: l, k
!
! ==========================================================
!
 !! 4th order hyper diffusion
  tau  = 7.2D0*(107.0D0/(NMAX+1))**2
  coef = dt*ER**4/(tau*3600.0D0*(real(NMAX,kind=8)*(NMAX+1))**2)
!
!$OMP PARALLEL default(SHARED), private(l,k)
 !$OMP DO schedule(STATIC)
  do l=2,mnwav
    do k=1,2
      qrot(k,l) = qrot(k,l)/( 1.0d0 + coef*ERFNN1(l)**2 )
      qdiv(k,l) = qdiv(k,l)/( 1.0d0 + coef*ERFNN1(l)**2 )
      qphi(k,l) = qphi(k,l)/( 1.0d0 + coef*ERFNN1(l)**2 )
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
end subroutine hdiff4



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine cal_gmean &
 &( data,           &
 &  gmean)
!
  real(8),intent(in) :: data(IMAX,JMAX)
  real(8),intent(out) :: gmean       !! Global mean
!
  real(8) :: work(JMAX)
  real(8),save :: denominator = -1.0d0
!
  real(8) :: ww
  integer :: i,j
!
! ====================================================================
!
  if ( denominator < 0.0d0 ) then
    denominator = 0.0d0
    do j=1,JMAX
      denominator = denominator + WEIGHT(j)
    end do
    denominator = denominator*IMAX
  end if
  
!$OMP PARALLEL default(SHARED), private(i,j,ww)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    ww=0.0d0
    do i=1,IMAX
      ww = ww + data(i,j)
    end do
    work(j) = ww
  end do
 !$OMP END DO
!$OMP END PARALLEL
  
  gmean = 0.0d0
  do j=1,JMAX
    gmean = gmean + work(j)*WEIGHT(j)
  end do
!
  gmean = gmean/denominator
!
end subroutine cal_gmean



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine departure      &
  &( j, i1, i2, iadv_only, &
  &  iterx, idepart_prev,  &!IN
  &  dt, dta,              &!IN
  &  cosa, sina,           &!IN
  &  u0_adv, v0_adv,       &!IN
  &  dudt0, dvdt0,         &!IN
  &  u_adv_halo, v_adv_halo,   &!IN
  &  utmp, vtmp, uwork, vwork, &!IN
  &  ii, jj,               &!OUT
  &  xi, yj,               &!OUT
  &  tt, ll,               &!OUT
  &  cosdtheta, sindtheta ) !OUT
 !
   integer,intent(in) :: j
   integer,intent(in) :: i1
   integer,intent(in) :: i2
   integer,intent(in) :: iadv_only
   integer,intent(in) :: iterx
   integer,intent(in) :: idepart_prev
   real(8),intent(in) :: dt
   real(8),intent(in) :: dta
   real(8),intent(in) :: cosa
   real(8),intent(in) :: sina
   real(8),intent(in) :: u0_adv(IMAX)
   real(8),intent(in) :: v0_adv(IMAX)
   real(8),intent(in) :: dudt0(IMAX)
   real(8),intent(in) :: dvdt0(IMAX)
   real(8),intent(in) :: u_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
   real(8),intent(in) :: v_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
   real(8),intent(inout) :: utmp(IMAX)
   real(8),intent(inout) :: vtmp(IMAX)
   real(8),intent(inout) :: uwork(IMAX)
   real(8),intent(inout) :: vwork(IMAX)
   integer,intent(inout) :: ii(IMAX)
   integer,intent(inout) :: jj(IMAX)
   real(8),intent(inout) :: xi(IMAX)
   real(8),intent(inout) :: yj(IMAX)
   real(8),intent(inout) :: cosdtheta(IMAX)
   real(8),intent(inout) :: sindtheta(IMAX)
   real(8),intent(inout) :: tt(IMAX)
   real(8),intent(inout) :: ll(IMAX)
 !
   real(8),parameter :: ddd = 0.5d0  !## 0.5 <= delta <= 1.0,  delta=2.0/3.0 is best?
   real(8),parameter :: ddd2 = 1.0d0 - ddd
 
   real(8) :: dt_ER
   real(8) :: uu, vv, rdis2, sindis_rdis, cosdis, sinlatd, coslatd_cosdlond, coslatd_sindlond
   real(8) :: coslatd, ww1, ww2, ww3, ww4, cosdlond, sindlon_d, alond, alatd, yyy
 
 !
   integer :: iter
   integer :: i
 !
 ! ===============================================================
   
   if ( iadv_only /= 1 .and. JCN_DEPARTURE /= 0 ) then
     do i=i1,i2
       uwork(i) = 0.5d0*dudt0(i)*dta - omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
       vwork(i) = 0.5d0*dvdt0(i)*dta + omg*ER*sina*SINLON(i)
     end do
   end if
 !
   iter_loop : do iter = 1, iterx
 !
     if ( idepart_prev == 0 .and. iter == 1 ) then
       !! Not use departure point data of previous step
       do i=i1,i2
          utmp(i) = u0_adv(i)
          vtmp(i) = v0_adv(i)
       end do
 
     else
       !! Use departure point data of previous step
       call lag5_uv                &
        &( i1, i2,                 &!IN
        &  ii, jj,                 &!IN
        &  xi, yj,                 &!IN
        &  cosdtheta, sindtheta,   &!IN
        &  u_adv_halo, v_adv_halo, &!IN
        &  utmp, vtmp      )        !INOUT
 
       if ( iadv_only == 1 ) then
           do i=i1,i2
             utmp(i) = ( utmp(i) + u0_adv(i)      )*0.5d0
             vtmp(i) = ( vtmp(i) + v0_adv(i)      )*0.5d0
           end do
       else
         if ( JCN_DEPARTURE == 0 ) then
           do i=i1,i2
             utmp(i) = ( utmp(i) + u0_adv(i) )*0.5d0
             vtmp(i) = ( vtmp(i) + v0_adv(i) )*0.5d0
           end do
         else
           do i=i1,i2
             utmp(i) = utmp(i) + uwork(i)
             vtmp(i) = vtmp(i) + vwork(i)
           end do
         end if
       end if
     end if
 
     dt_er = dt/ER
     ww4 = 1.0d0/(ALAT_SL(1)-ALAT_SL(JMAX))*(JMAX-1)
 
     do i=i1,i2
         uu = utmp(i)*dt_er
         vv = vtmp(i)*dt_er
 
         rdis2 = uu**2 + vv**2
         sindis_rdis = 1.0d0 - rdis2/6.0d0   !! sin(rdis)/rdis
         cosdis      = 1.0d0 - rdis2*0.5d0
         sinlatd     = SINLAT(j)*cosdis - vv*COSLAT(j)*sindis_rdis
         coslatd_cosdlond = COSLAT(j)*cosdis + vv*SINLAT(j)*sindis_rdis
         coslatd_sindlond = -uu*sindis_rdis
         coslatd     = sqrt( coslatd_cosdlond**2 + coslatd_sindlond**2 )
         alatd       = atan2(sinlatd,coslatd)  
         
         ww1 = sign(0.5d0,coslatd+1.0d-50) - sign(0.5d0,coslatd-1.0d-50)
         ww2 = 1.0d0/( ww1 + coslatd )
         cosdlond = (ww1+coslatd_cosdlond)*ww2
         sindlon_d =      coslatd_sindlond *ww2
         alond     = ALON(i) + atan2(sindlon_d,cosdlond) + PI*2.0d0
         
         ww3 = 1.0d0/( 1.0d0 + cosdis )
         cosdtheta(i) = ww3*( COSLAT(j)*coslatd + (1.0d0+SINLAT(j)*sinlatd)*cosdlond )
         sindtheta(i) = ww3*( (SINLAT(j)+sinlatd)*sindlon_d )
 
         xi(i) = alond/DLON + 1.0d0
         ii(i) = int( xi(i) )
         xi(i) = xi(i) - ii(i)
         ii(i) = mod(ii(i)-1,IMAX) + 1
 !           
         yj(i) = alatd
         yyy = (ALAT_SL(1)-alatd)*ww4 + 1.0d0
         jj(i) = int( yyy )
         !! ALAT_SL(jj(i,j)) >= alatd >= ALAT_SL(jj(i,j)+1) should be satisfied.
         jj(i) = jj(i) + int( sign(0.501D0,ALAT_SL(jj(i))-alatd) + sign(0.501D0,ALAT_SL(jj(i)+1)-alatd) )

        
         if(iter == iterx) then
           tt(i) = alatd
           ll(i) = alond-2.0d0*PI
         end if
     end do
 !
   end do iter_loop
 
 end subroutine departure
 



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine lag3 &
 &( i1, i2,     &
 &  ii, jj,     &!IN
 &  xi, yj,     &!IN
 &  phi_halo,   &!IN
 &  phi   )      !INOUT
!
  integer,intent(in) :: i1
  integer,intent(in) :: i2
  integer,intent(in) :: ii(IMAX)
  integer,intent(in) :: jj(IMAX)
  real(8),intent(in) :: xi(IMAX)
  real(8),intent(in) :: yj(IMAX)
  real(8),intent(in) :: phi_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(inout) :: phi(IMAX)
!
  real(8) :: a1,a2,a3,a4,xx,yy
  real(8) :: phi1,phi2,phi3,phi4
!
  integer :: i, ni, nj
!
! ===========================================================================

!xx  call e_time__start(6,"lag3")
!
!xx!$OMP PARALLEL default(SHARED), private(i,j,xx,yy,ni,nj,a1,a2,a3,a4,phi1,phi2,phi3,phi4)
!xx !$OMP DO schedule(STATIC)
!xx  do j = 1, JMAX
    do i = i1,i2
!
      xx = xi(i)
      yy = yj(i)
      ni = ii(i)
      nj = jj(i)
!
      a1 = (xx  )*(xx-1)*(xx-2)/((-1)*(-2)*(-3))
      a2 = (xx+1)*(xx-1)*(xx-2)/(( 1)*(-1)*(-2))
      a3 = (xx+1)*(xx  )*(xx-2)/(( 2)*( 1)*(-1))
      a4 = (xx+1)*(xx  )*(xx-1)/(( 3)*( 2)*( 1))
!
      phi1 =  a1*phi_halo(ni-1,nj-1) + a2*phi_halo(ni  ,nj-1) &
       &    + a3*phi_halo(ni+1,nj-1) + a4*phi_halo(ni+2,nj-1)
!
      phi2 =  a1*phi_halo(ni-1,nj  ) + a2*phi_halo(ni  ,nj  ) &
       &    + a3*phi_halo(ni+1,nj  ) + a4*phi_halo(ni+2,nj  )
!
      phi3 =  a1*phi_halo(ni-1,nj+1) + a2*phi_halo(ni  ,nj+1) &
       &    + a3*phi_halo(ni+1,nj+1) + a4*phi_halo(ni+2,nj+1)
!
      phi4 =  a1*phi_halo(ni-1,nj+2) + a2*phi_halo(ni  ,nj+2) &
       &    + a3*phi_halo(ni+1,nj+2) + a4*phi_halo(ni+2,nj+2)
!
      a1 = (yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+1))*(yy-ALAT_SL(nj+2))        &
       &   / ( (ALAT_SL(nj-1)-ALAT_SL(nj  ))*(ALAT_SL(nj-1)-ALAT_SL(nj+1)) &
       &      *(ALAT_SL(nj-1)-ALAT_SL(nj+2)) )
      a2 = (yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj+1))*(yy-ALAT_SL(nj+2))        &
       &   / ( (ALAT_SL(nj)-ALAT_SL(nj-1))*(ALAT_SL(nj)-ALAT_SL(nj+1))     &
       &      *(ALAT_SL(nj)-ALAT_SL(nj+2))  )
      a3 = (yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+2))        &
       &   / ( (ALAT_SL(nj+1)-ALAT_SL(nj-1))*(ALAT_SL(nj+1)-ALAT_SL(nj  )) &
       &      *(ALAT_SL(nj+1)-ALAT_SL(nj+2)) )
      a4 = (yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+1))        &
       &   / ( (ALAT_SL(nj+2)-ALAT_SL(nj-1))*(ALAT_SL(nj+2)-ALAT_SL(nj  )) &
       &      *(ALAT_SL(nj+2)-ALAT_SL(nj+1)) )
!
      phi(i) =  a1*phi1 + a2*phi2 + a3*phi3 + a4*phi4
!
!      phi(i,j) = max( phi(i,j), min( phi_halo(ni,nj  ),phi_halo(ni+1,nj  ),   &
!       &                             phi_halo(ni,nj+1),phi_halo(ni+1,nj+1) ) )
!
!      phi(i,j) = min( phi(i,j), max( phi_halo(ni,nj  ),phi_halo(ni+1,nj  ),   &
!       &                             phi_halo(ni,nj+1),phi_halo(ni+1,nj+1) ) )
!
    end do
!xx  end do
!xx !$OMP END DO
!xx!$OMP END PARALLEL

!xx  call e_time__end(6,"lag3")
!
end subroutine lag3


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine lag3_uv        &
 &( i1, i2,               &
 &  ii, jj,               &!IN
 &  xi, yj,               &!IN
 &  cosdtheta, sindtheta, &!IN
 &  u_halo, v_halo,       &!IN
 &  u, v      )            !INOUT
!
  integer,intent(in) :: i1
  integer,intent(in) :: i2
  integer,intent(in) :: ii(IMAX)
  integer,intent(in) :: jj(IMAX)
  real(8),intent(in) :: xi(IMAX)
  real(8),intent(in) :: yj(IMAX)
  real(8),intent(in) :: cosdtheta(IMAX)
  real(8),intent(in) :: sindtheta(IMAX)
  real(8),intent(in) :: u_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(in) :: v_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(inout) :: u(IMAX)
  real(8),intent(inout) :: v(IMAX)
!
!
  integer :: i,ni,nj
  real(8) :: a1,a2,a3,a4
  real(8) :: u1,u2,u3,u4,v1,v2,v3,v4,uu,vv
  real(8) :: xx,yy
!
! ===========================================================================

!xx  call e_time__start(8,"lag3_uv")
!
! ------------------------------------------------------------------------
!
!
!xx!$OMP PARALLEL default(SHARED), private(i,j,ni,nj,a1,a2,a3,a4,u1,u2,u3,u4,v1,v2,v3,v4) &
!xx!$OMP  & private(uu1,uu2,uu3,uu4,uu5,uu6,uu,vv1,vv2,vv3,vv4,vv,xx,yy )
!xx !$OMP DO schedule(STATIC)
!xx  do j = 1, JMAX
    do i = i1,i2
!
      xx = xi(i)
      yy = yj(i)
      ni = ii(i)
      nj = jj(i)
!
      a1 = (yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+1))*(yy-ALAT_SL(nj+2))        &
       &   / ( (ALAT_SL(nj-1)-ALAT_SL(nj  ))*(ALAT_SL(nj-1)-ALAT_SL(nj+1)) &
       &      *(ALAT_SL(nj-1)-ALAT_SL(nj+2)) )
      a2 = (yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj+1))*(yy-ALAT_SL(nj+2))        &
       &   / ( (ALAT_SL(nj  )-ALAT_SL(nj-1))*(ALAT_SL(nj  )-ALAT_SL(nj+1)) &
       &      *(ALAT_SL(nj  )-ALAT_SL(nj+2)) )
      a3 = (yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+2))        &
       &   / ( (ALAT_SL(nj+1)-ALAT_SL(nj-1))*(ALAT_SL(nj+1)-ALAT_SL(nj  )) &
       &      *(ALAT_SL(nj+1)-ALAT_SL(nj+2)) )
      a4 = ( yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+1))       &
       &   / ( (ALAT_SL(nj+2)-ALAT_SL(nj-1))*(ALAT_SL(nj+2)-ALAT_SL(nj  )) &
       &      *(ALAT_SL(nj+2)-ALAT_SL(nj+1)) )
!
      u1 =   a1*u_halo(ni-1,nj-1) + a2*u_halo(ni-1,nj  )  &
       &   + a3*u_halo(ni-1,nj+1) + a4*u_halo(ni-1,nj+2)
!
      u2 =   a1*u_halo(ni  ,nj-1) + a2*u_halo(ni  ,nj  )  &
       &   + a3*u_halo(ni  ,nj+1) + a4*u_halo(ni  ,nj+2)
!
      u3 =   a1*u_halo(ni+1,nj-1) + a2*u_halo(ni+1,nj  )  &
       &   + a3*u_halo(ni+1,nj+1) + a4*u_halo(ni+1,nj+2)
!
      u4 =   a1*u_halo(ni+2,nj-1) + a2*u_halo(ni+2,nj  )  &
       &   + a3*u_halo(ni+2,nj+1) + a4*u_halo(ni+2,nj+2)
!
!
      v1 =   a1*v_halo(ni-1,nj-1) + a2*v_halo(ni-1,nj  )  &
       &   + a3*v_halo(ni-1,nj+1) + a4*v_halo(ni-1,nj+2)
!
      v2 =   a1*v_halo(ni  ,nj-1) + a2*v_halo(ni  ,nj  )  &
       &   + a3*v_halo(ni  ,nj+1) + a4*v_halo(ni  ,nj+2)
!
      v3 =   a1*v_halo(ni+1,nj-1) + a2*v_halo(ni+1,nj  )  &
       &   + a3*v_halo(ni+1,nj+1) + a4*v_halo(ni+1,nj+2)
!
      v4 =   a1*v_halo(ni+2,nj-1) + a2*v_halo(ni+2,nj  )  &
       &   + a3*v_halo(ni+2,nj+1) + a4*v_halo(ni+2,nj+2)
!
      a1 = (xx  )*(xx-1)*(xx-2)/((-1)*(-2)*(-3))
      a2 = (xx+1)*(xx-1)*(xx-2)/(( 1)*(-1)*(-2))
      a3 = (xx+1)*(xx  )*(xx-2)/(( 2)*( 1)*(-1))
      a4 = (xx+1)*(xx  )*(xx-1)/(( 3)*( 2)*( 1))
!
!      hosei_d = hosei(nj)*(1.0d0-yy) + hosei(nj+1)*yy
!      hosei1 = hosei_d*(-1.0d0-xx)
!      hosei2 = hosei_d*(      -xx)
!      hosei3 = hosei_d*( 1.0d0-xx)
!      hosei4 = hosei_d*( 2.0d0-xx)
!      cos_hosei1 = 1.0d0-0.5d0*hosei1**2
!      cos_hosei2 = 1.0d0-0.5d0*hosei2**2
!      cos_hosei3 = 1.0d0-0.5d0*hosei3**2
!      cos_hosei4 = 1.0d0-0.5d0*hosei4**2
!
!      uu1 = u1*cos_hosei1 - v1*hosei1
!      uu2 = u2*cos_hosei2 - v2*hosei2
!      uu3 = u3*cos_hosei3 - v3*hosei3
!      uu4 = u4*cos_hosei4 - v4*hosei4
!
!      vv1 = v1*cos_hosei1 + u1*hosei1
!      vv2 = v2*cos_hosei2 + u2*hosei2
!      vv3 = v3*cos_hosei3 + u3*hosei3
!      vv4 = v4*cos_hosei4 + u4*hosei4
!
!      uu =  a1*uu1 + a2*uu2 + a3*uu3 + a4*uu4
!      vv =  a1*vv1 + a2*vv2 + a3*vv3 + a4*vv4
!
      uu =  a1*u1 + a2*u2 + a3*u3 + a4*u4
      vv =  a1*v1 + a2*v2 + a3*v3 + a4*v4
!
      u(i) = cosdtheta(i)*uu - sindtheta(i)*vv
      v(i) = sindtheta(i)*uu + cosdtheta(i)*vv
!
    end do
!xx  end do
!xx !$OMP END DO
!xx!$OMP END PARALLEL

!xx  call e_time__end(8,"lag3_uv")
!
end subroutine lag3_uv



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine lag5_uv        &
 &( i1, i2,            &!IN
 &  ii, jj,               &!IN
 &  xi, yj,               &!IN
 &  cosdtheta, sindtheta, &!IN
 &  u_halo, v_halo,       &!IN
 &  u, v      )            !INOUT
!
  integer,intent(in) :: i1
  integer,intent(in) :: i2
  integer,intent(in) :: ii(IMAX)
  integer,intent(in) :: jj(IMAX)
  real(8),intent(in) :: xi(IMAX)
  real(8),intent(in) :: yj(IMAX)
  real(8),intent(in) :: cosdtheta(IMAX)
  real(8),intent(in) :: sindtheta(IMAX)
  real(8),intent(in) :: u_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(in) :: v_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(inout) :: u(IMAX)
  real(8),intent(inout) :: v(IMAX)
!
  integer :: i,ni,nj
  real(8) :: a1,a2,a3,a4,a5,a6
  real(8) :: u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6,uu,vv
  real(8) :: xx,yy,yyy
!
! ===========================================================================

!xx  call e_time__start(9,"lag5_uv")
!
! ------------------------------------------------------------------------
!
!xx!$OMP PARALLEL default(SHARED), private(i,j,ni,nj,a1,a2,a3,a4,a5,a6,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6) &
!xx!$OMP  & private(uu1,uu2,uu3,uu4,uu5,uu6,uu,vv1,vv2,vv3,vv4,vv5,vv6,vv,xx,yy,yyy )
!xx !$OMP DO schedule(STATIC)
!xx  do j = 1, JMAX
    do i = i1,i2
!
      xx = xi(i)
      yy = yj(i)
      ni = ii(i)
      nj = jj(i)
!
      a1 = ( yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+1))       &
       &   *(yy-ALAT_SL(nj+2))*(yy-ALAT_SL(nj+3))                          &
       &   / ( (ALAT_SL(nj-2)-ALAT_SL(nj-1))*(ALAT_SL(nj-2)-ALAT_SL(nj  )) &
       &      *(ALAT_SL(nj-2)-ALAT_SL(nj+1))*(ALAT_SL(nj-2)-ALAT_SL(nj+2)) &
       &      *(ALAT_SL(nj-2)-ALAT_SL(nj+3)) )
      a2 = ( yy-ALAT_SL(nj-2))*(yy-ALAT_SL(nj  ))*(yy-ALAT_SL(nj+1))       &
       &   *(yy-ALAT_SL(nj+2))*(yy-ALAT_SL(nj+3))                          &
       &   / ( (ALAT_SL(nj-1)-ALAT_SL(nj-2))*(ALAT_SL(nj-1)-ALAT_SL(nj  )) &
       &      *(ALAT_SL(nj-1)-ALAT_SL(nj+1))*(ALAT_SL(nj-1)-ALAT_SL(nj+2)) &
       &      *(ALAT_SL(nj-1)-ALAT_SL(nj+3)) )
      a3 = ( yy-ALAT_SL(nj-2))*(yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj+1))       &
       &   *(yy-ALAT_SL(nj+2))*(yy-ALAT_SL(nj+3))                          &
       &   / ( (ALAT_SL(nj)-ALAT_SL(nj-2))*(ALAT_SL(nj)-ALAT_SL(nj-1))     &
       &      *(ALAT_SL(nj)-ALAT_SL(nj+1))*(ALAT_SL(nj)-ALAT_SL(nj+2))     &
       &      *(ALAT_SL(nj)-ALAT_SL(nj+3)) )
      a4 = ( yy-ALAT_SL(nj-2))*(yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))       &
       &   *(yy-ALAT_SL(nj+2))*(yy-ALAT_SL(nj+3))                          &
       &   / ( (ALAT_SL(nj+1)-ALAT_SL(nj-2))*(ALAT_SL(nj+1)-ALAT_SL(nj-1)) &
       &      *(ALAT_SL(nj+1)-ALAT_SL(nj  ))*(ALAT_SL(nj+1)-ALAT_SL(nj+2)) &
       &      *(ALAT_SL(nj+1)-ALAT_SL(nj+3)) )
      a5 = ( yy-ALAT_SL(nj-2))*(yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))       &
       &   *(yy-ALAT_SL(nj+1))*(yy-ALAT_SL(nj+3))                          &
       &   / ( (ALAT_SL(nj+2)-ALAT_SL(nj-2))*(ALAT_SL(nj+2)-ALAT_SL(nj-1)) &
       &      *(ALAT_SL(nj+2)-ALAT_SL(nj  ))*(ALAT_SL(nj+2)-ALAT_SL(nj+1)) &
       &      *(ALAT_SL(nj+2)-ALAT_SL(nj+3)) )
      a6 = ( yy-ALAT_SL(nj-2))*(yy-ALAT_SL(nj-1))*(yy-ALAT_SL(nj  ))       &
       &   *(yy-ALAT_SL(nj+1))*(yy-ALAT_SL(nj+2))                          &
       &   / ( (ALAT_SL(nj+3)-ALAT_SL(nj-2))*(ALAT_SL(nj+3)-ALAT_SL(nj-1)) &
       &      *(ALAT_SL(nj+3)-ALAT_SL(nj  ))*(ALAT_SL(nj+3)-ALAT_SL(nj+1)) &
       &      *(ALAT_SL(nj+3)-ALAT_SL(nj+2)) )
!
      u1 =   a1*u_halo(ni-2,nj-2) + a2*u_halo(ni-2,nj-1)  &
       &   + a3*u_halo(ni-2,nj  ) + a4*u_halo(ni-2,nj+1)  &
       &   + a5*u_halo(ni-2,nj+2) + a6*u_halo(ni-2,nj+3)
!
      u2 =   a1*u_halo(ni-1,nj-2) + a2*u_halo(ni-1,nj-1)  &
       &   + a3*u_halo(ni-1,nj  ) + a4*u_halo(ni-1,nj+1)  &
       &   + a5*u_halo(ni-1,nj+2) + a6*u_halo(ni-1,nj+3)
!
      u3 =   a1*u_halo(ni  ,nj-2) + a2*u_halo(ni  ,nj-1)  &
       &   + a3*u_halo(ni  ,nj  ) + a4*u_halo(ni  ,nj+1)  &
       &   + a5*u_halo(ni  ,nj+2) + a6*u_halo(ni  ,nj+3)
!
      u4 =   a1*u_halo(ni+1,nj-2) + a2*u_halo(ni+1,nj-1)  &
       &   + a3*u_halo(ni+1,nj  ) + a4*u_halo(ni+1,nj+1)  &
       &   + a5*u_halo(ni+1,nj+2) + a6*u_halo(ni+1,nj+3)
!
      u5 =   a1*u_halo(ni+2,nj-2) + a2*u_halo(ni+2,nj-1)  &
       &   + a3*u_halo(ni+2,nj  ) + a4*u_halo(ni+2,nj+1)  &
       &   + a5*u_halo(ni+2,nj+2) + a6*u_halo(ni+2,nj+3)
!
      u6 =   a1*u_halo(ni+3,nj-2) + a2*u_halo(ni+3,nj-1)  &
       &   + a3*u_halo(ni+3,nj  ) + a4*u_halo(ni+3,nj+1)  &
       &   + a5*u_halo(ni+3,nj+2) + a6*u_halo(ni+3,nj+3)
!
      v1 =   a1*v_halo(ni-2,nj-2) + a2*v_halo(ni-2,nj-1)  &
       &   + a3*v_halo(ni-2,nj  ) + a4*v_halo(ni-2,nj+1)  &
       &   + a5*v_halo(ni-2,nj+2) + a6*v_halo(ni-2,nj+3)
!
      v2 =   a1*v_halo(ni-1,nj-2) + a2*v_halo(ni-1,nj-1)  &
       &   + a3*v_halo(ni-1,nj  ) + a4*v_halo(ni-1,nj+1)  &
       &   + a5*v_halo(ni-1,nj+2) + a6*v_halo(ni-1,nj+3)
!
      v3 =   a1*v_halo(ni  ,nj-2) + a2*v_halo(ni  ,nj-1)  &
       &   + a3*v_halo(ni  ,nj  ) + a4*v_halo(ni  ,nj+1)  &
       &   + a5*v_halo(ni  ,nj+2) + a6*v_halo(ni  ,nj+3)
!
      v4 =   a1*v_halo(ni+1,nj-2) + a2*v_halo(ni+1,nj-1)  &
       &   + a3*v_halo(ni+1,nj  ) + a4*v_halo(ni+1,nj+1)  &
       &   + a5*v_halo(ni+1,nj+2) + a6*v_halo(ni+1,nj+3)
!
      v5 =   a1*v_halo(ni+2,nj-2) + a2*v_halo(ni+2,nj-1)  &
       &   + a3*v_halo(ni+2,nj  ) + a4*v_halo(ni+2,nj+1)  &
       &   + a5*v_halo(ni+2,nj+2) + a6*v_halo(ni+2,nj+3)
!
      v6 =   a1*v_halo(ni+3,nj-2) + a2*v_halo(ni+3,nj-1)  &
       &   + a3*v_halo(ni+3,nj  ) + a4*v_halo(ni+3,nj+1)  &
       &   + a5*v_halo(ni+3,nj+2) + a6*v_halo(ni+3,nj+3)
!
      a1 = (xx+1)*(xx  )*(xx-1)*(xx-2)*(xx-3)/((-1)*(-2)*(-3)*(-4)*(-5))
      a2 = (xx+2)*(xx  )*(xx-1)*(xx-2)*(xx-3)/(( 1)*(-1)*(-2)*(-3)*(-4))
      a3 = (xx+2)*(xx+1)*(xx-1)*(xx-2)*(xx-3)/(( 2)*( 1)*(-1)*(-2)*(-3))
      a4 = (xx+2)*(xx+1)*(xx  )*(xx-2)*(xx-3)/(( 3)*( 2)*( 1)*(-1)*(-2))
      a5 = (xx+2)*(xx+1)*(xx  )*(xx-1)*(xx-3)/(( 4)*( 3)*( 2)*( 1)*(-1))
      a6 = (xx+2)*(xx+1)*(xx  )*(xx-1)*(xx-2)/(( 5)*( 4)*( 3)*( 2)*( 1))
!
      yyy = ( yy-ALAT_SL(nj) )/( ALAT_SL(nj+1)-ALAT_SL(nj) )
!xx      yyy=(ALAT_SL(jj(i,j))-yj(i,j))/(ALAT_SL(jj(i,j))-ALAT_SL(jj(i,j)+1))
!
!      hosei_d = hosei(nj)*(1.0d0-yyy) + hosei(nj+1)*yyy
!      hosei1 = hosei_d*(-2.0d0-xx)
!      hosei2 = hosei_d*(-1.0d0-xx)
!      hosei3 = hosei_d*(      -xx)
!      hosei4 = hosei_d*( 1.0d0-xx)
!      hosei5 = hosei_d*( 2.0d0-xx)
!      hosei6 = hosei_d*( 3.0d0-xx)
!      cos_hosei1 = 1.0d0-0.5d0*hosei1**2
!      cos_hosei2 = 1.0d0-0.5d0*hosei2**2
!      cos_hosei3 = 1.0d0-0.5d0*hosei3**2
!      cos_hosei4 = 1.0d0-0.5d0*hosei4**2
!      cos_hosei5 = 1.0d0-0.5d0*hosei5**2
!      cos_hosei6 = 1.0d0-0.5d0*hosei6**2
!
!      uu1 = u1*cos_hosei1 - v1*hosei1
!      uu2 = u2*cos_hosei2 - v2*hosei2
!      uu3 = u3*cos_hosei3 - v3*hosei3
!      uu4 = u4*cos_hosei4 - v4*hosei4
!      uu5 = u5*cos_hosei5 - v5*hosei5
!      uu6 = u6*cos_hosei6 - v6*hosei6
!
!      vv1 = v1*cos_hosei1 + u1*hosei1
!      vv2 = v2*cos_hosei2 + u2*hosei2
!      vv3 = v3*cos_hosei3 + u3*hosei3
!      vv4 = v4*cos_hosei4 + u4*hosei4
!      vv5 = v5*cos_hosei5 + u5*hosei5
!      vv6 = v6*cos_hosei6 + u6*hosei6
!
!      uu =  a1*uu1 + a2*uu2 + a3*uu3  &
!       &  + a4*uu4 + a5*uu5 + a6*uu6
!
!      vv =  a1*vv1 + a2*vv2 + a3*vv3  &
!       &  + a4*vv4 + a5*vv5 + a6*vv6
!
      uu =  a1*u1 + a2*u2 + a3*u3  &
       &  + a4*u4 + a5*u5 + a6*u6
!
      vv =  a1*v1 + a2*v2 + a3*v3  &
       &  + a4*v4 + a5*v5 + a6*v6
!
      u(i) = cosdtheta(i)*uu - sindtheta(i)*vv
      v(i) = sindtheta(i)*uu + cosdtheta(i)*vv
!
    end do
!xx  end do
!xx !$OMP END DO
!xx!$OMP END PARALLEL

!xx  call e_time__end(9,"lag5_uv")
!
end subroutine lag5_uv

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine set_halo   &
 & ( data_halo  )      !INOUT
!
  real(8),intent(inout) :: data_halo(1-MGN_X:IMAX+MGN_X,1-MGN_Y:JMAX+MGN_Y)
!
  logical,save :: first = .true.
  real(8),save :: a1,a2
  real(8) :: zmean11,zmean12,zmean21,zmean22,zmean31,zmean32,zmean41,zmean42
  real(8) :: data_np, data_sp
  integer :: m, i, j
!
! =========================================================================
!
  if ( first ) then
!
    first = .false.
!
    a1 = (ALAT_SL(0)-ALAT_SL(-2))*(ALAT_SL(0)-ALAT_SL(1))*(ALAT_SL(0)-ALAT_SL(2))        &
     &   /( (ALAT_SL(-1)-ALAT_SL(-2))*(ALAT_SL(-1)-ALAT_SL(1))*(ALAT_SL(-1)-ALAT_SL(2)) )
    a2 = (ALAT_SL(0)-ALAT_SL(-1))*(ALAT_SL(0)-ALAT_SL(1))*(ALAT_SL(0)-ALAT_SL(2))        &
     &   /( (ALAT_SL(-2)-ALAT_SL(-1))*(ALAT_SL(-2)-ALAT_SL(1))*(ALAT_SL(-2)-ALAT_SL(2)) )
!
!xx    write(6,*) 'a1,a2,a3=',a1,a2,a3
!xx    stop 344
!
  end if
!
!$OMP PARALLEL default(SHARED), private(m,i,j)
  do m=1,MGN_Y-1
   !$OMP DO schedule(STATIC)
    do i=1,IMAX2
      data_halo(i      ,-m) = data_halo(i+IMAX2,m)
      data_halo(i+IMAX2,-m) = data_halo(i      ,m)
      data_halo(i      ,JMAX+1+m) = data_halo(i+IMAX2,JMAX+1-m)
      data_halo(i+IMAX2,JMAX+1+m) = data_halo(i      ,JMAX+1-m)
    end do
   !$OMP END DO
  end do
! !$OMP DO schedule(STATIC)
!  do i=1,IMAX2
!    data_halo(i      ,0) = a2*(data_halo(i,1)+data_halo(i,-1))  &
!     &                     + a1*(data_halo(i,2)+data_halo(i,-2))
!!    data_halo(i+IMAX2,0) = data_halo(i,0)
!    data_halo(i      ,JMAX+1) = a2*(data_halo(i,JMAX+1+1)+data_halo(i,JMAX+1-1))  &
!     &                          + a1*(data_halo(i,JMAX+1+2)+data_halo(i,JMAX+1-2))
!!    data_halo(i+IMAX2,JMAX+1) = data_halo(i,JMAX+1)
!  end do
! !$OMP END DO
 !$OMP SECTIONS
  !$OMP SECTION
  zmean11 = sum( data_halo(1:IMAX2,1) )/IMAX2
  !$OMP SECTION
  zmean12 = sum( data_halo(IMAX2+1:IMAX,1) )/IMAX2
  !$OMP SECTION
  zmean21 = sum( data_halo(1:IMAX2,2) )/IMAX2
  !$OMP SECTION
  zmean22 = sum( data_halo(IMAX2+1:IMAX,2) )/IMAX2
  !$OMP SECTION
  zmean31 = sum( data_halo(1:IMAX2,JMAX-1) )/IMAX2
  !$OMP SECTION
  zmean32 = sum( data_halo(IMAX2+1:IMAX,JMAX-1) )/IMAX2
  !$OMP SECTION
  zmean41 = sum( data_halo(1:IMAX2,JMAX) )/IMAX2
  !$OMP SECTION
  zmean42 = sum( data_halo(IMAX2+1:IMAX,JMAX) )/IMAX2
 !$OMP END SECTIONS
  data_np = a1*(zmean11+zmean12) + a2*(zmean21+zmean22)
  data_sp = a1*(zmean41+zmean42) + a2*(zmean31+zmean32)
 !$OMP DO schedule(STATIC)
  do i=1,IMAX
    data_halo(i,0)      = data_np
    data_halo(i,JMAX+1) = data_sp
  end do
 !$OMP END DO
 !$OMP DO schedule(STATIC)
  do j=1-MGN_Y,JMAX+MGN_Y
    do m=1,MGN_X
      data_halo(1-m,j)    = data_halo(IMAX+1-m,j)
      data_halo(IMAX+m,j) = data_halo(m,j)
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
end subroutine set_halo


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine set_halo_uv  &
 & ( data_halo )   !INOUT
!
  real(8),intent(out) :: data_halo(1-MGN_X:IMAX+MGN_X,1-MGN_Y:JMAX+MGN_Y)
!
  logical,save :: first = .true.
  real(8),save :: a1,a2
  real(8) :: zonal_cos11, zonal_cos21, zonal_cos31, zonal_cos41
  real(8) :: zonal_sin11, zonal_sin21, zonal_sin31, zonal_sin41
  real(8) :: zonal_cos12, zonal_cos22, zonal_cos32, zonal_cos42
  real(8) :: zonal_sin12, zonal_sin22, zonal_sin32, zonal_sin42
  real(8) :: zonal_np_cos, zonal_np_sin, zonal_sp_cos, zonal_sp_sin
  integer :: m, i, j
!
! =========================================================================
!
  if ( first ) then
!
    first = .false.

    a1 = (ALAT_SL(0)-ALAT_SL(-2))*(ALAT_SL(0)-ALAT_SL(1))*(ALAT_SL(0)-ALAT_SL(2))        &
     &   /( (ALAT_SL(-1)-ALAT_SL(-2))*(ALAT_SL(-1)-ALAT_SL(1))*(ALAT_SL(-1)-ALAT_SL(2)) )
    a2 = (ALAT_SL(0)-ALAT_SL(-1))*(ALAT_SL(0)-ALAT_SL(1))*(ALAT_SL(0)-ALAT_SL(2))        &
     &   /( (ALAT_SL(-2)-ALAT_SL(-1))*(ALAT_SL(-2)-ALAT_SL(1))*(ALAT_SL(-2)-ALAT_SL(2)) )
!
!xx    write(6,*) 'a1,a2,a3=',a1,a2,a3
!xx    stop 344
!
  end if
!
!$OMP PARALLEL default(SHARED), private(m,i,j)
  do m=1,MGN_Y-1
   !$OMP DO schedule(STATIC)
    do i=1,IMAX2
      data_halo(i      ,-m) = -data_halo(i+IMAX2,m)
      data_halo(i+IMAX2,-m) = -data_halo(i      ,m)
      data_halo(i      ,JMAX+1+m) = -data_halo(i+IMAX2,JMAX+1-m)
      data_halo(i+IMAX2,JMAX+1+m) = -data_halo(i      ,JMAX+1-m)
    end do
   !$OMP END DO
  end do
 !$OMP DO schedule(STATIC)
  do i=1,IMAX2
    data_halo(i      ,0) = a1*(data_halo(i,1)+data_halo(i,-1))  &
     &                     + a2*(data_halo(i,2)+data_halo(i,-2))
    data_halo(i+IMAX2,0) = - data_halo(i,0)
    data_halo(i      ,JMAX+1) = a1*(data_halo(i,JMAX+1+1)+data_halo(i,JMAX+1-1))  &
     &                          + a2*(data_halo(i,JMAX+1+2)+data_halo(i,JMAX+1-2))
    data_halo(i+IMAX2,JMAX+1) = - data_halo(i,JMAX+1)
  end do
 !$OMP END DO

 !$OMP SECTIONS
  !$OMP SECTION
  zonal_cos11 = sum( data_halo(1:IMAX2,1)*COSLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos12 = sum( data_halo(IMAX2+1:IMAX,1)*COSLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin11 = sum( data_halo(1:IMAX2,1)*SINLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin12 = sum( data_halo(IMAX2+1:IMAX,1)*SINLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos21 = sum( data_halo(1:IMAX2,2)*COSLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos22 = sum( data_halo(IMAX2+1:IMAX,2)*COSLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin21 = sum( data_halo(1:IMAX2,2)*SINLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin22 = sum( data_halo(IMAX2+1:IMAX,2)*SINLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos31 = sum( data_halo(1:IMAX2,JMAX-1)*COSLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos32 = sum( data_halo(IMAX2+1:IMAX,JMAX-1)*COSLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin31 = sum( data_halo(1:IMAX2,JMAX-1)*SINLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin32 = sum( data_halo(IMAX2+1:IMAX,JMAX-1)*SINLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos41 = sum( data_halo(1:IMAX2,JMAX)*COSLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_cos42 = sum( data_halo(IMAX2+1:IMAX,JMAX)*COSLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin41 = sum( data_halo(1:IMAX2,JMAX)*SINLON(1:IMAX2) )*2.0d0/IMAX2
  !$OMP SECTION
  zonal_sin42 = sum( data_halo(IMAX2+1:IMAX,JMAX)*SINLON(IMAX2+1:IMAX) )*2.0d0/IMAX2
 !$OMP END SECTIONS

  zonal_np_cos = a1*(zonal_cos11+zonal_cos12) + a2*(zonal_cos21+zonal_cos22)
  zonal_np_sin = a1*(zonal_sin11+zonal_sin12) + a2*(zonal_sin21+zonal_sin22)
  zonal_sp_cos = a1*(zonal_cos41+zonal_cos42) + a2*(zonal_cos31+zonal_cos32)
  zonal_sp_sin = a1*(zonal_sin41+zonal_sin42) + a2*(zonal_sin31+zonal_sin32)

 !$OMP DO schedule(STATIC)
  do i=1,IMAX
     data_halo(i,0) = zonal_np_cos*COSLON(i) + zonal_np_sin*SINLON(i)
     data_halo(i,JMAX+1) = zonal_sp_cos*COSLON(i) + zonal_sp_sin*SINLON(i)
  end do
 !$OMP END DO

 !$OMP DO schedule(STATIC)
  do j=1-MGN_Y,JMAX+MGN_Y
    do m=1,MGN_X
      data_halo(1-m,j)    = data_halo(IMAX+1-m,j)
      data_halo(IMAX+m,j) = data_halo(m,j)
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
end subroutine set_halo_uv


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! subroutine for computing the Fourier coefficients
subroutine lgcoeffs(F, cn, nthreads, ipiv, A)
  complex*16, intent(inout) :: F(:,:)
  complex*16, intent(in) :: cn(:,:), A(:,:)
  integer, intent(in) :: ipiv(:), nthreads

  integer :: m, n, k, iret, info
  complex*16, allocatable :: tmp(:,:), tmp1d(:)
  integer(8) :: plan
!    
  n = size(F,1)
  m = size(F,2)
  allocate(tmp1d(m))
  allocate(tmp(n,m))
  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(nthreads)
  call dfftw_plan_dft_1d(plan, m, tmp1d, tmp1d, FFTW_FORWARD, FFTW_ESTIMATE)
  do k = 1,n
      tmp1d = F(k,:)
      call dfftw_execute_dft(plan, tmp1d, tmp1d)
      tmp(k,:) = (/tmp1d(m/2+1:m),tmp1d(1:m/2)/) / m
  end do
  deallocate(tmp1d)
  call dfftw_destroy_plan(plan)
!   solve the linear system using the blas
  call zgetrs('N', n, m, A, n, ipiv, tmp, n, info)
!  
  if (info .ne. 0) then
      print*,'No Succes s'
  end if
!
  F = tmp * cn
  !k = 1
  !F(1:k,:) = 0.0d0; F(:, 1:k) = 0.0d0;
  !F(n-k+1:n,:) = 0.0d0; Tcoeff(:,m-k+1:m) = 0.0d0;
  deallocate(tmp)
end subroutine lgcoeffs
!
! Subroutine for doubling up the function
subroutine doubleUp(F, dF)
  real(8), intent(in) :: F(:,:)
  complex*16, intent(inout) :: dF(:,:)
!
  integer :: m, n
  n = size(F, 1)
  m = size(F, 2)
 
  dF(1:n,1:m/2) = F(:,m/2+1:m)
  dF(1:n,m/2+1:m) = F(:,1:m/2)
  dF(n+1:2*n,:) = F(n:1:-1,:)
  
end subroutine doubleUp
!
! subroutine for initializing the shifts
 subroutine shiftfftw(n, m, cn)
        integer, intent(in) :: m, n
        complex*16, intent(inout) :: cn(:,:)
        complex*16, allocatable :: ck(:)
        integer :: k,j
!
        allocate(ck(m))
        !$OMP PARALLEL DO private(j,k)
        do k=-m/2,m/2-1
            j = k + m/2 + 1
            ck(j) = (-1d0) ** (k)
        end do
        !$OMP END PARALLEL DO        
        cn = spread(ck,1,2*n)
        deallocate(ck)
  end subroutine shiftfftw

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
end program sw_sh
