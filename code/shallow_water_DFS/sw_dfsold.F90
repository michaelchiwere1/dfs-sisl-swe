program sw_dfs

  !! Semi-implicit semi-Lagrangian shallow water model using double Fourier series

  use prm_phconst, only : PI, PI2, ER, OMG, GRAV
  use fft, only : fft__ini, fft__towave, fft__togrid
  use fft_y, only : fft_y__ini, fft_y__g2w, fft_y__g2w_uv, fft_y__w2g, fft_y__w2g_dy, fft_y__w2g_uv
  use xderiv, only : xderiv_run
  use com_dfs, only : com_dfs__ini
  use uv2divrot_dfs, only : uv2divrot_dfs__run, uv2divrot_dfs__uv2chipsi, uv2divrot_dfs__laplacian
  use divrot2uv_dfs, only : divrot2uv_dfs__run
  use helmholtz_dfs, only : helmholtz_dfs__run, helmholtz_dfs__like, helmholtz_dfs__hdiff2
  use hdiff4_dfs, only : hdiff4_dfs__run
  use monit, only : monit__ini, monit__output
  use monit_spectrum, only : monit_spectrum__ini, monit_spectrum__output
  use monit_norm, only : monit_norm__ini, monit_norm__output
  use grads, only : grads__outxy
  use e_time, only : e_time__start, e_time__end, e_time__output
  use vals
!
  implicit none
!
!  integer,parameter :: JCN_DFS = 0 !! Old method and Cheong's basis functions
!  integer,parameter :: JCN_DFS = 1 !! New method to calculate expantion coefficients
  integer,parameter :: JCN_DFS = 2 !! New method to calculate expantion coefficients and New basis functions

!  integer,parameter :: JCN_GRID = 0   !! colat=PI*(j+0.5)/JMAX0, j=0,...,JMAX0-1
  integer,parameter :: JCN_GRID = 1   !! colat=PI*j/JMAX0, j=0,...,JMAX0   (JMAX=JMAX0+1)
!  integer,parameter :: JCN_GRID = -1  !! colat=PI*j/JMAX0, j=1,...,JMAX0-1 (JMAX=JMAX0-1)
 
 ! integer,parameter :: JCN_INITIAL = 1  !! Williamson test case 1 (Advection)
  integer,parameter :: JCN_INITIAL = 2  !! Williamson test case 2
 ! integer,parameter :: JCN_INITIAL = 5  !! Williamson test case 5
 ! integer,parameter :: JCN_INITIAL = 6  !! Williamson test case 6
 ! integer,parameter :: JCN_INITIAL = 10 !! Galewsky 2004 test case
 ! integer,parameter :: JCN_INITIAL = 11 !! Galewsky 2004 like test case

 ! integer,parameter :: JCN_PROG_H = 0 !! h*=h-hs is prognostic variable
  integer,parameter :: JCN_PROG_H = 1 !! h is prognostic variable

  integer,parameter :: JCN_HDIFF = 0  !! No diffusion
 ! integer,parameter :: JCN_HDIFF = 1  !! 2nd order diff. for Galewsky test case
!  integer,parameter :: JCN_HDIFF = 2  !! 4th order hyper diffusion
  
 ! integer,parameter :: JCN_ZONALFILTER = 0 !! No zonal filter
  integer,parameter :: JCN_ZONALFILTER = 1 !! Use zonal Fourier filter
  integer,parameter :: M0 = 20             !! For zonal Fourier filter

 ! integer,parameter :: JCN_DEPARTURE = 0  !! Use extrapolated winds
  integer,parameter :: JCN_DEPARTURE = 1  !! Use predicted winds

!  integer,parameter :: JCN_MONIT = 0  !! No monitor output
  integer,parameter :: JCN_MONIT = 1  !! Write monitor output

!  integer,parameter :: JCN_MONIT_SPECTRUM = 0  !! No kinetic energy spectrum output
  integer,parameter :: JCN_MONIT_SPECTRUM = 1   !! Write kinetic energy spectrum output.
                                                !! This needs Legendre transform.
  
  integer,parameter :: INTHR_MONIT = 24  !! Interval of monitor output (hour)

 ! integer, parameter :: IT = 0    !! Lagrange spatial interpolation
  integer, parameter ::  IT = 1    !! dfs spatial interpolation

  integer, parameter :: VELOCITYD = 3   !! 3D velocity field
 ! integer, parameter :: VELOCITYD = 2   !! 2D velocity field

 ! integer,parameter :: IMAX = 16*32
!  integer,parameter :: JMAX = 8 * 32 + JCN_GRID
 ! integer,parameter :: NMAX = 7
 ! integer, parameter :: NMAX = JMAX - 1
 ! real(8),save :: TIMESTEP = 3600.0d0         !! 1 hour
!
 ! integer,parameter :: IMAX = 128            !! About 300km resolution
 ! integer,parameter :: JMAX = 64+JCN_GRID
 ! integer,parameter :: NMAX = 63             !! Linear grid
!  integer,parameter :: NMAX = 42            !! Quadric grid
  real(8),save :: TIMESTEP = 3600.0d0 / 1d0       !! 1 hour
!
  integer,parameter :: IMAX = 160            !! 250km resolution
  integer,parameter :: JMAX = 80+JCN_GRID
  integer,parameter :: NMAX = 78             !! Linear grid
!  integer,parameter :: NMAX = 53             !! Quadric grid
 ! real(8),save :: TIMESTEP = 3600.0d0        !! 1 hour
!
 ! integer,parameter :: IMAX = 320            !! 120km resolution
 ! integer,parameter :: JMAX = 160+JCN_GRID
!  integer,parameter :: NMAX = 159            !! Linear grid
!  integer,parameter :: NMAX = 106            !! Quadric grid
 ! real(8),parameter :: TIMESTEP = 1800.0d0   !! 30 min.
!
  !integer,parameter :: IMAX = 640            !! 60km resolution
 ! integer,parameter :: JMAX = 320+JCN_GRID
 ! integer,parameter :: NMAX = 319            !! Linear grid
!  integer,parameter :: NMAX = 213            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 1200.0d0   !! 20 min.
!
!  integer,parameter :: IMAX = 1920           !! 20km resolution
!  integer,parameter :: JMAX = 960+JCN_GRID   !! J is from North to South
!  integer,parameter :: NMAX = 959            !! Linear grid
!  integer,parameter :: NMAX = 639            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 600.0d0    !! 10 min.
!
!  integer,parameter :: IMAX = 3840           !! 10km resolution
!  integer,parameter :: JMAX = 1920+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 1919           !! Linear grid
!  integer,parameter :: NMAX = 1279           !! Quadric grid
!  integer,parameter :: NMAX = 959            !! Cubic grid
!  real(8),parameter :: TIMESTEP = 360.0d0    !! 6 min.
!
!  integer,parameter :: IMAX = 7680           !! 5km resolution
!  integer,parameter :: JMAX = 3840+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 3839           !! Linear grid
!  integer,parameter :: NMAX = 1919           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 225.0d0    !! 225 sec.
!
!  integer,parameter :: IMAX = 15360          !! 2.6km resolution
!  integer,parameter :: JMAX = 7680+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 7679           !! Linear grid
!  integer,parameter :: NMAX = 5119           !! Quadric grid
!  integer,parameter :: NMAX = 3839           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 144.0d0    !! 144 sec.
!
!  integer,parameter :: IMAX = 20480          !! 2.0km resolution
!  integer,parameter :: JMAX = 10240+JCN_GRID !! J is from North to South
!  integer,parameter :: NMAX = 10239          !! Linear grid
!  integer,parameter :: NMAX = 6826           !! Quadric grid
!  integer,parameter :: NMAX = 5119           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 120.0d0    !! 120 sec.
!
!  integer,parameter :: IMAX = 30720          !! 1.3km resolution
!  integer,parameter :: JMAX = 15360+JCN_GRID !! J is from North to South
!  integer,parameter :: NMAX = 15359          !! Linear grid   
!  integer,parameter :: NMAX = 10239          !! Quadric grid
!  integer,parameter :: NMAX = 7679           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 90.0d0     !! 90 sec.
!
!  integer,parameter :: IVECLEN = 14  !! For test
  integer,parameter :: IVECLEN = min(IMAX,512)     !! For scalar machine
  integer,parameter :: NDIV = (IMAX-1)/IVECLEN + 1 !! For scalar machine
  
!  integer,parameter :: NDIV = 1                    !! For vector machine
!  integer,parameter :: IVECLEN = (IMAX-1)/NDIV + 1 !! For vector machine

  integer,parameter :: MMAX = NMAX
  integer,parameter :: NNUM = NMAX + 1
!
  integer,parameter :: IMAX2=IMAX/2 !! IMAX must be even.
!
  integer,parameter :: MGN_X = 5
  integer,parameter :: MGN_Y = 5
!
  real(8),parameter :: DLON = PI2/IMAX
  real(8),parameter :: DLAT = PI/(JMAX-JCN_GRID)
!
  real(8),save :: ALON(IMAX)
  real(8),save :: ALON1(IMAX)
  real(8),save :: SINLON(IMAX)
  real(8),save :: COSLON(IMAX)
  real(8),save :: ALAT(JMAX)
  real(8),save :: YLAT(JMAX)
  real(8),save :: COSLAT(JMAX)
  real(8),save :: COSLAT_INV(JMAX)
  real(8),save :: SINLAT(JMAX)
!
  real(8),save :: ACOS2LAT_INV(JMAX)
  real(8),save :: ACOSLAT_INV(JMAX)
  real(8),save :: WEIGHT(JMAX)
  
  real(8),save :: WGAUSS(7)
  real(8),save :: XGAUSS(7)
!
! ==================================================================
!
  call e_time__start(1,"shallow water")

  call initialize
!
  call main
!
  call e_time__end(1,"shallow water")
  call e_time__output
!
!
!******************************************************************
contains
!******************************************************************


subroutine initialize
!
  integer :: i,j
!
! -----------------------------------------------------------------------------
!
  write(6,*) "JCN_DFS            =", JCN_DFS
  write(6,*) "JCN_GRID           =", JCN_GRID
  write(6,*) "JCN_INITIAL        =", JCN_INITIAL
  write(6,*) "JCN_PROG_H         =", JCN_PROG_H
  write(6,*) "JCN_HDIFF          =", JCN_HDIFF
  write(6,*) "JCN_ZONALFILTER    =", JCN_ZONALFILTER
  write(6,*) "M0                 =", M0
  write(6,*) "JCN_DEPARTURE      =", JCN_DEPARTURE
  write(6,*) "JCN_MONIT          =", JCN_MONIT
  write(6,*) "JCN_MONIT_SPECTRUM =", JCN_MONIT_SPECTRUM
  write(6,*) "INTHR_MONIT        =", INTHR_MONIT
  write(6,*) "IMAX     =", IMAX 
  write(6,*) "JMAX     =", JMAX
  write(6,*) "NMAX     =", NMAX
  write(6,*) "TIMESTEP =", TIMESTEP
!  
  call com_dfs__ini( JCN_DFS, JCN_GRID, IMAX, JMAX, MMAX, NMAX )
  call fft_y__ini( IMAX, JMAX, MMAX, alat, weight ) !IN
  
!$OMP PARALLEL default(SHARED), private(j,i)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    COSLAT(j) = cos(ALAT(j))
    SINLAT(j) = sin(ALAT(j))
    COSLAT_INV(j) = 1.0d0/COSLAT(j)
    YLAT(j) = ALAT(j)*180.0d0/PI
    ACOS2LAT_INV(j) = COSLAT_INV(j)**2/ER
    ACOSLAT_INV(j) = COSLAT_INV(j)/ER
  end do
 !$OMP END DO

 !$OMP DO schedule(STATIC)
  do i=1,IMAX
    ALON(i) = DLON*(i-1)     !##   0 <= ALON  < 2*PI
    if ( ALON(i) >= PI ) then
       ALON1(i) = ALON(i)-2*PI  !## -PI <= ALON1 < PI
    else
       ALON1(i) = ALON(i)       !## -PI <= ALON1 < PI
    end if
    SINLON(i) = sin(ALON(i))
    COSLON(i) = cos(ALON(i))
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call output_weight_lat( IMAX, JMAX, YLAT, WEIGHT )

!  write(6,*) 'YLAT(:)=',YLAT(:)
!  write(6,*) 'WEIGHT(:)=',WEIGHT(:)

  call fft__ini( JCN_ZONALFILTER, M0, IMAX, JMAX, MMAX, coslat )

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
    write(iunit,'(a,i6,a,f14.10,a,f14.10)') 'YDEF ',JMAX,'  LINEAR  ',ylat(JMAX),' ',ylat(JMAX-1)-ylat(JMAX)
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
  !! Expansion coefficients
  !!  Real part: [cosine omponent]/2, Imaginary part: -[sine component]/2 for m >= 1
  complex(8),save :: qphi(NNUM,0:MMAX)  !! GRAV*height (geopotential height)
  complex(8),save :: qrot(NNUM,0:MMAX)  !! Vorticity
  complex(8),save :: qdiv(NNUM,0:MMAX)  !! Divergence
!
  real(8),save :: phi0(IMAX,JMAX)       !! GRAV*height
  
  real(8),save :: phi_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y) !! With halo region 
!
  real(8),save :: phis (IMAX,JMAX)      !! Surface geopotential height
  real(8),save :: phisx(IMAX,JMAX)
  real(8),save :: phisy(IMAX,JMAX)
!
  real(8),save :: u0(IMAX,JMAX)         !! Zonal wind
  real(8),save :: v0(IMAX,JMAX)         !! Meridional wind
!
  real(8),save :: u_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y) !! With halo region 
  real(8),save :: v_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  
  real(8),save :: u_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),save :: v_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
!
  real(8),save :: um_adv(IMAX,JMAX)     !! Zonal wind for advection (past)
  real(8),save :: vm_adv(IMAX,JMAX)     !! Zonal wind for advection (past)
!
  real(8),save :: u0_adv(IMAX,JMAX)     !! Zonal wind for advection (present)
  real(8),save :: v0_adv(IMAX,JMAX)     !! Zonal wind for advection (present)
!
  real(8),save :: dudtm(IMAX,JMAX)      !! d(u)/dt (past)
  real(8),save :: dvdtm(IMAX,JMAX)      !! d(v)/dt (past)
  real(8),save :: dphidtm(IMAX,JMAX)    !! d(phi)/dt (past)
  
  real(8),save :: dudtm_li(IMAX,JMAX)   !! d(u)/dt linear term (past)
  real(8),save :: dvdtm_li(IMAX,JMAX)   !! d(v)/dt linear term (past)
  real(8),save :: dphidtm_li(IMAX,JMAX) !! d(u)/dt linear term (past)
!
  real(8),save :: dudt0(IMAX,JMAX)      !! d(u)/dt (present)
  real(8),save :: dvdt0(IMAX,JMAX)      !! d(v)/dt (present)
  real(8),save :: dphidt0(IMAX,JMAX)    !! d(phi)/dt (present)
  
  real(8),save :: dudt0_li(IMAX,JMAX)   !! d(u)/dt linear term (present)
  real(8),save :: dvdt0_li(IMAX,JMAX)   !! d(v)/dt linear term (present)
  real(8),save :: dphidt0_li(IMAX,JMAX) !! d(phi)/dt linear term (present)
!
  integer,save :: ii(IMAX,JMAX)    !! Grid position near departure point for semi-Lag.
  integer,save :: jj(IMAX,JMAX)    !! Grid position near departure point for semi-Lag.
  real(8),save :: xi(IMAX,JMAX)    !! Interpolation in longitudinal direction
  real(8),save :: yj(IMAX,JMAX)    !! Interpolation in latitudinal direction
  real(8),save :: cosdtheta(IMAX,JMAX) !! Rotation of the wind for semi-Lag.
  real(8),save :: sindtheta(IMAX,JMAX) !! Rotation of the wind for semi-Lag.
!
  real(8),save :: f_8byte(IMAX,JMAX)
  real(8),save :: phi_zonal(JMAX)
  real(8),save :: u_zonal(JMAX)
  
  real(8) :: workj(JMAX)
  real(8) :: work1(IMAX)
  real(8) :: work2(IMAX)
  real(8) :: work3(IMAX)
  real(8) :: work4(IMAX)

! Initialize for dfs Interpolation
  real(8), save :: ll(IMAX, JMAX)
  real(8), save :: tt(IMAX, JMAX)
  real(8), allocatable :: ud_halo(:,:), vd_halo(:,:), hd_halo(:,:), wd_halo(:,:)
  real(8), save :: th_d(IMAX*JMAX)  ! departure 1d
  real(8), save :: lb_d(IMAX*JMAX)
  complex(8), allocatable :: coeffs2du(:,:), coeffs1du(:)  !! dfs coefficiensts
  complex(8), allocatable :: coeffs2dv(:,:), coeffs1dv(:)
  complex(8), allocatable :: coeffs2dh(:,:), coeffs1dh(:)
  complex(8), allocatable :: coeffs2dw(:,:), coeffs1dw(:)
  real(8), allocatable :: tmph(:), tmpu(:), tmpv(:), tmpw(:)
  integer :: tmpindex
  integer(8) :: N1, N2
  real(8) :: vij, uij

  real(8) :: tgrid, tvals, tfin, tcopv,tcoeff, tstart, tfinal
  integer, parameter :: nthreads = 8  ! Number of threads for fftw3
  real(8) :: tmpcosd
!
  real(8) :: phibar
  real(8) :: dta, dtb, beta, betauv
  real(8) :: dt, dt_prev
  real(8) :: alon_c, alat_c, rr, r, phis0, hh0, uu0
  real(8) :: aa, bb, cc, ww1, ww2
  real(8) :: thetac, ramdac, ramda, alpha0, cosa, sina
  real(8) :: alat0, alat1, alat2, umax, ss, ww, xx, uu, en, f, hhat
  real(8) :: gmean, gmean_phis, gmean_mass_0, gmean_energy_0, gmean_vorticity_0, gmean_divergence_0
!
  integer :: iterx, kt_end, ntimes_hour, iadv_only
  integer :: kind_phi, kind_phis, kind_uv, nstep, nstep_start, intv_monit
  integer :: istep, i, j, nn, idiv, i1, i2, idepart_prevstep, ntmax
  integer :: istep_monit = 0
!
  real(8),parameter :: omg_6 = 7.848d-6
  real(8),parameter :: ak_6  = omg_6
  real(8),parameter :: hh0_6  = 8.0d3
  real(8),parameter :: ir_6  = 4

  call e_time__start(2,"main")

  if ( JCN_INITIAL == 1 ) then
    iadv_only = 1   !! Advection
    kind_phi  = 1
    kind_phis = 1
    kind_uv   = 2
    kt_end    = 24*12  !! 12 days
  else if ( JCN_INITIAL == 2 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 2
    kind_phis = 1
    kind_uv   = 2
    kt_end    = 24*5  !! 5 days
  else if ( JCN_INITIAL == 5 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 5
    kind_phis = 2
    kind_uv   = 5
    kt_end    = 24*15  !! 15 days
  else if ( JCN_INITIAL == 6 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 6
    kind_phis = 1
    kind_uv   = 6
    kt_end    = 24*14  !! 14 days
  else if ( JCN_INITIAL == 10 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 10
    kind_phis = 1
    kind_uv   = 10
    kt_end    = 24*6  !! 6 days
  else if ( JCN_INITIAL == 11 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 11
    kind_phis = 1
    kind_uv   = 11
    kt_end    = 24*6   !! 6 days
  else
    write(6,*) "Error: JCN_INITIAL=",JCN_INITIAL," is not supported."
    stop 999
  end if
!
  ntimes_hour = int( (3600.0d0+1.0d-20)/TIMESTEP )
  if ( abs( 3600.0d0 - TIMESTEP*ntimes_hour ) > 1.0d-20 ) then
    write(6,*) "Error: 3600.0/TIMESTEP should be integer."
    write(6,*) "       TIMESTEP =", TIMESTEP
    stop 999
  end if
!
  nstep_start = 1    !# nstep_start >= 1
  nstep = kt_end*ntimes_hour + nstep_start
!
  intv_monit = INTHR_MONIT*ntimes_hour
!
  beta   = 1.0d0    !! 1.0 <= beta <= 1.2? 1.4?
  betauv = 1.0d0    !! 1.0 <= beta <= 1.2? 1.4?
! =====================================================================
!
!  open( 11, file = 'data.dr',  form = 'UNFORMATTED',     &
!   &        access = 'DIRECT', recl = 4*(IMAX)*(JMAX)  )
!!
!  open( 12, file = 'data2.dr',  form = 'UNFORMATTED',     &
!   &        access = 'DIRECT', recl = 4*(IMAX)*(JMAX)  )

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
         phi0(i,j) = GRAV*hh0_6 + ER**2*( aa + bb*cos(ALON(i)*ir_6)  &
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
        workj(j) = 0.0d0
        do nn=1,7
           xx = aa + bb*XGAUSS(nn)
           if ( alat0 < xx .and. xx < alat1 ) then
              uu = umax/en*exp( 1.0d0/((xx-alat0)*(xx-alat1)) )
           else
              uu = 0.0d0
           end if
           f  = 2.0d0*OMG*sin(xx)
           workj(j) = workj(j) + ER*uu*( f + tan(xx)*uu/ER )*WGAUSS(nn)*0.5d0
        end do
        workj(j) = workj(j)*( ALAT(j) - ALAT(j+1) )
     end do   
    !$OMP END DO
    !$OMP END PARALLEL
     do j=JMAX-1,1,-1
        phi_zonal(j) = phi_zonal(j+1) - workj(j) !! Integral
     end do
     ss = sum( phi_zonal(1:JMAX)*WEIGHT(1:JMAX) )/sum( WEIGHT(1:JMAX) )
    !$OMP PARALLEL default(SHARED), private(j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        phi_zonal(j) = phi_zonal(j) - ss + GRAV*10000.0d0 !! global mean = 10000m
     end do
    !$OMP END DO
    !$OMP END PARALLEL
     
!     write(6,*) "phi_zonal(:) = ",phi_zonal(:)
  
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
!           phi(i,j) = phi_zonal(j)
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
    write(6,*) "Error: kind_phi=",kind_phi," is not supported."
    stop 999
! ------------------------------
  end if
! ------------------------------
!
! ###  surface g.p.h.
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
  else if ( kind_phis == 2 ) then
! -------------------------
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
  else
! ------------------------------
    write(6,*) "Error: kind_phis=",kind_phis," is not supported."
    stop 999
! ------------------------------
  end if
! ------------------------------
   
! ------------------------------
  if (kind_uv == 0) then
! ------------------------------
    !$OMP PARALLEL default(SHARED), private(i,j,ww)
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
  else if (kind_uv == 2) then
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
  else if (kind_uv == 5) then
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
  else if (kind_uv == 6) then
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
  else if (kind_uv == 10 .or. kind_uv == 11) then
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
!
  cosa = cos(alpha0)
  sina = sin(alpha0)

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
     &( phi0,   &!INOUT
     &  qphi  )  !WORK

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
!
    call fft__towave &
     &( phis      )      !INOUT
!
    call fft_y__g2w &
     &( phis,    &!IN
     &  qphi  )   !INOUT
!
    call fft_y__w2g  &
     &( qphi,    &!IN
     &  phis   )  !OUT
    call fft_y__w2g_dy &
     &( qphi,      &!IN
     &  phisy  )    !INOUT
    call xderiv_run              &
     &( phis, phisy, coslat_inv, &!IN
     &  phisx          )          !OUT
!
    call fft__togrid &
     &( phis        )    !INOUT
    call fft__togrid &
     &( phisx       )    !INOUT
    call fft__togrid &
     &( phisy      )     !INOUT
   
!   call grads__outxy( "phisy", 1, IMAX, JMAX, ylat, phisy(1:IMAX,1:JMAX) ) !IN
!
    call grid_to_wave      &
     &( u0, v0, phi0,      &!IN
     &  qrot, qdiv, qphi )  !OUT
    if ( JCN_DFS >= 1 ) then
      call uv2divrot_dfs__laplacian(qrot) !INOUT !! Stream function => vorticity
      call uv2divrot_dfs__laplacian(qdiv) !INOUT !! Velocity potential => divergence
    end if
!
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
    call monit_norm__ini( "norm", INTHR_MONIT, ntmax, IMAX, JMAX, weight, f_8byte )
    call monit__ini( "data", INTHR_MONIT, ntmax, IMAX, JMAX, YLAT ) !IN
    if ( JCN_MONIT_SPECTRUM == 1 ) then
      call monit_spectrum__ini( "spectrum", INTHR_MONIT, ntmax, JCN_DFS, NMAX, MMAX, ER )
    end if
  end if

!   call grads__outxy( "u0", 1, IMAX, JMAX, ylat, u0(1:IMAX,1:JMAX) ) !IN
!   call grads__outxy( "v0", 1, IMAX, JMAX, ylat, v0(1:IMAX,1:JMAX) ) !IN
!   call grads__outxy( "phi0", 1, IMAX, JMAX, ylat, phi0(1:IMAX,1:JMAX) ) !IN
!
!
! =====================================================================
! #####   Start integration   #####
! =====================================================================
!
!
!  call cal_mean &
!   &( phi,      &
!   &  phi_mean)
!  write(6,*) 'phi_mean = ',phi_mean
!
! ===================================================================
!
  dt = -999.0d0
  if(IT == 1) then
    N1 = 2*(JMAX-JCN_GRID)
    N2 = IMAX
    allocate(ud_halo(JMAX,IMAX))
    allocate(vd_halo(JMAX,IMAX))
    allocate(hd_halo(JMAX,IMAX))
    allocate(wd_halo(JMAX,IMAX))
    allocate(coeffs1du(2*(JMAX-JCN_GRID)*IMAX))
    allocate(coeffs1dv(2*(JMAX-JCN_GRID)*IMAX))
    allocate(coeffs1dh(2*(JMAX-JCN_GRID)*IMAX))
    allocate(coeffs1dw(2*(JMAX-JCN_GRID)*IMAX))
    allocate(tmpu(JMAX*IMAX))
    allocate(tmph(JMAX*IMAX))
    allocate(tmpv(JMAX*IMAX))
    allocate(tmpw(JMAX*IMAX))
    tgrid =0d0; tvals=0d0; tfin=0d0; tcopv=0d0; tcoeff=0d0;
    
  end if



  do istep = 1, nstep
!
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
     !$OMP PARALLEL default(SHARED), private(j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        um_adv(:,j)     = u0_adv(:,j) 
        vm_adv(:,j)     = v0_adv(:,j) 
        dudtm(:,j)      = dudt0(:,j)
        dvdtm(:,j)      = dvdt0(:,j)
        dphidtm(:,j)    = dphidt0(:,j)
        dudtm_li(:,j)   = dudt0_li(:,j)
        dvdtm_li(:,j)   = dvdt0_li(:,j)
        dphidtm_li(:,j) = dphidt0_li(:,j)
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    end if
!
    dta = dt*0.5d0
    dtb = dt*0.5d0
!
    call e_time__end(3,"main_cp_past")
!
!   ------------------------------------------------------
!
!
    if ( iadv_only == 1 ) then

      call grid_to_wave_to_grid &
       &( phi0,   &!INOUT
       &  qphi  )  !WORK

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
      call monit__output( 1, nn, f_8byte )      
      call monit__output( 2, nn, u0 )
      call monit__output( 3, nn, v0 )
      call monit_norm__output( nn, f_8byte )
          
      call cal_gmean &
       &( f_8byte,   &
       &  gmean)      !! Global mean of h
      if ( nn == 1 ) then
         call cal_gmean &
          &( phis,   & 
          &  gmean_phis)
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
!
!
      if ( JCN_INITIAL >= 2 ) then
        if ( JCN_MONIT_SPECTRUM == 1 ) then
          call monit_spectrum__output( nn, qrot, qdiv )
        end if
        
        call fft_y__w2g    &
         &( qrot,                &!IN
         &  f_8byte  )            !OUT
        call fft__togrid      &
         &( f_8byte  )            !INOUT
        call monit__output( 4, nn, f_8byte )

        call cal_gmean &
         &( f_8byte,   &
         &  gmean)
        if ( nn == 1 ) then
           gmean_vorticity_0 = gmean
        end if
        write(6,*) "Global mean of vorticity = ",gmean, "(", (gmean-gmean_vorticity_0)/gmean_vorticity_0, ")"
!
        call fft_y__w2g    &
         &( qdiv,                &!IN
         &  f_8byte  )            !OUT
        call fft__togrid      &
         &( f_8byte  )            !INOUT
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
!
!   -------------------------------------------------------
    
    call e_time__start(4,"main_integ1")
!
    if ( iadv_only == 1 ) then
     !$OMP PARALLEL default(SHARED), private(i,j)
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
            u_halo(i,j) = u0(i,j) + (-dudtm(i,j) + dudt0(i,j)*2.0d0)*dtb + ( dudtm_li(i,j) - dudt0_li(i,j) )*betauv*dtb &
             &       + 2.0d0*omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
            v_halo(i,j) = v0(i,j) + (-dvdtm(i,j) + dvdt0(i,j)*2.0d0)*dtb + ( dvdtm_li(i,j) - dvdt0_li(i,j) )*betauv*dtb &
             &       - 2.0d0*omg*ER*sina*SINLON(i)
            phi_halo(i,j) = phi0(i,j) + (-dphidtm(i,j) + dphidt0(i,j)*2.0d0)*dtb + ( dphidtm_li(i,j) - dphidt0_li(i,j) )*beta*dtb
          end do
        end do
       !$OMP END DO
       !$OMP END PARALLEL
      else
       !$OMP PARALLEL default(SHARED), private(i,j,ww1,ww2)
       !$OMP DO schedule(STATIC)
        do j=1,JMAX
          do i=1,IMAX
            ww1 = (-dudtm(i,j) + dudt0(i,j)*2.0d0)*dtb + 2.0d0*omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
            u_adv_halo(i,j) = u0(i,j) + 0.5*ww1
            u_halo(i,j) = u0(i,j) + ( dudtm_li(i,j) - dudt0_li(i,j) )*betauv*dtb + ww1
             
            ww2 = (-dvdtm(i,j) + dvdt0(i,j)*2.0d0)*dtb - 2.0d0*omg*ER*sina*SINLON(i)
            v_adv_halo(i,j) = v0(i,j) + 0.5*ww2
            v_halo(i,j) = v0(i,j) + ( dvdtm_li(i,j) - dvdt0_li(i,j) )*betauv*dtb + ww2
             
            phi_halo(i,j) = phi0(i,j) + (-dphidtm(i,j) + dphidt0(i,j)*2.0d0)*dtb + ( dphidtm_li(i,j) - dphidt0_li(i,j) )*beta*dtb
          end do
        end do
       !$OMP END DO
       !$OMP END PARALLEL
      end if
    end if
!   
    if (IT == 0) then
    call set_halo_uv &
     &( u_adv_halo )  !INOUT
    call set_halo_uv &
     &( v_adv_halo )  !INOUT 

    call set_halo_uv  &
     &( u_halo )  !INOUT
    call set_halo_uv  &
     &( v_halo )  !OUT
 
    call set_halo &
     &( phi_halo ) !INOUT
!
    call e_time__end(4,"main_integ1")
    
    call e_time__start(5,"main_semilag")
!
    if ( dt_prev == dt ) then
      iterx = 3             !! Number of iterations
      idepart_prevstep = 1  !! Use departure point data of previous step
    else
      iterx = 4             !! Number of iterations
      idepart_prevstep = 0  !! Not use departure point data of previous step
    end if

   !$OMP PARALLEL default(SHARED), private( j, idiv, i1, i2, work1, work2, work3, work4 )
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do idiv = 1, NDIV
       !! i is from i1 to i2.
        i1 = 1 + (idiv-1)*IVECLEN
        i2 = min( IMAX, idiv*IVECLEN )
        
        call departure                      &
         &( j, i1, i2, iadv_only,           &!IN
         &  iterx, idepart_prevstep,        &!IN
         &  dt, dta,                        &!IN
         &  cosa, sina,                     &!IN
         &  u0_adv, v0_adv,                 &!IN
         &  dudt0, dvdt0,                   &!IN
         &  u_adv_halo, v_adv_halo,         &!IN
         &  work1, work2, work3, work4,     &!INOUT
         &  ii(1,j), jj(1,j),               &!INOUT
         &  xi(1,j), yj(1,j),               &!INOUT
         &  tt(1,j), ll(1,j),               &!INOUT
         &  cosdtheta(1,j), sindtheta(1,j) ) !INOUT
!
        call lag3               &
         &( i1, i2,             &!IN
         &  ii(1,j), jj(1,j),   &!IN
         &  xi(1,j), yj(1,j),   &!IN
         &  phi_halo,           &!IN
         &  phi0(1,j)   )        !INOUT

        call lag5_uv                        &
         &( i1, i2,                         &!IN
         &  ii(1,j), jj(1,j),               &!IN
         &  xi(1,j), yj(1,j),               &!IN
         &  cosdtheta(1,j), sindtheta(1,j), &!IN
         &  u_halo, v_halo,                 &!IN
         &  u0(1,j), v0(1,j)      )          !OUT
      end do   
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    call e_time__end(5,"main_semilag")
  !
! interpolating using dfs method 
  else if (IT == 1 ) then
    if(VELOCITYD .eq. 3) then
      !$OMP PARALLEL DO private(i,j,tmpcosd)
      do j = 1, JMAX
        do i = 1, IMAX
          ud_halo(j,i) = -SINLON(i)*u_halo(i,j) -  COSLON(i)*SINLAT(j)*v_halo(i,j)
          vd_halo(j,i) = COSLON(i)*u_halo(i,j) - SINLON(i)*SINLAT(j)*v_halo(i,j)
          wd_halo(j,i) = COSLAT(j)*v_halo(i,j)
          hd_halo(j,i) = phi_halo(i,j) 
        end  do
      end do
    !$OMP END PARALLEL DO

    else
      !!$OMP PARALLEL DO private(i,j,tmpcosd)
      do j = 1, JMAX
        tmpcosd = DCOS(ALAT(j))
        do i = 1, IMAX
          ud_halo(j,i) = u_halo(i,j)*tmpcosd
          vd_halo(j,i) = v_halo(i,j)*tmpcosd
          hd_halo(j,i) = phi_halo(i,j) 
        end  do
      end do
   ! !$OMP END PARALLEL DO
    end if

    call set_halo_uv &
     &( u_adv_halo )  !INOUT
    call set_halo_uv &
     &( v_adv_halo )  !INOUT 
     call set_halo_uv  &
     &( u_halo )  !INOUT
    call set_halo_uv  &
     &( v_halo )  !OUT
 
              
  !
      call e_time__end(4,"main_integ1")
      call e_time__start(5,"main_semilag")
  !
      if ( dt_prev == dt ) then
        iterx = 3             !! Number of iterations
        idepart_prevstep = 1  !! Use departure point data of previous step
      else
        iterx = 4             !! Number of iterations
        idepart_prevstep = 0  !! Not use departure point data of previous step
      end if
  
     !$OMP PARALLEL default(SHARED), private( j, idiv, i1, i2, work1, work2, work3, work4 )
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do idiv = 1, NDIV
         !! i is from i1 to i2.
          i1 = 1 + (idiv-1)*IVECLEN
          i2 = min( IMAX, idiv*IVECLEN )
          
          call departure                      &
           &( j, i1, i2, iadv_only,           &!IN
           &  iterx, idepart_prevstep,        &!IN
           &  dt, dta,                        &!IN
           &  cosa, sina,                     &!IN
           &  u0_adv, v0_adv,                 &!IN
           &  dudt0, dvdt0,                   &!IN
           &  u_adv_halo, v_adv_halo,         &!IN
           &  work1, work2, work3, work4,     &!INOUT
           &  ii(1,j), jj(1,j),               &!INOUT
           &  xi(1,j), yj(1,j),               &!INOUT
           &  tt(1,j), ll(1,j),               &!INOUT
           &  cosdtheta(1,j), sindtheta(1,j) ) !INOUT

         
          end do   
      end do
     !$OMP END DO
     !$OMP END PARALLEL

  ! 
      tstart = omp_get_wtime()
      !$OMP PARALLEL DO private(i,j,tmpindex)
      do j=1, JMAX
        do i=1, IMAX
          tmpindex = (j-1)*IMAX + i
          th_d(tmpindex) = tt(i,j)
          lb_d(tmpindex) = mod(ll(i,j),2.0d0*pi)
                 end do
      end do
    !$OMP END PARALLEL DO
      tfinal = omp_get_wtime()
      tgrid = tgrid + (tfinal - tstart)
  !    print*, lb_d
  !
  ! compute the coefficients
      tstart = omp_get_wtime()
      call vals2coeffs(ud_halo, coeffs2du, IMAX,JMAX, JCN_GRID, nthreads)
      call vals2coeffs(vd_halo, coeffs2dv, IMAX, JMAX, JCN_GRID, nthreads)
      call vals2coeffs(wd_halo, coeffs2dw, IMAX, JMAX, JCN_GRID, nthreads)
      call vals2coeffs(hd_halo, coeffs2dh, IMAX, JMAX, JCN_GRID, nthreads)
      tfinal = omp_get_wtime()
      tvals = tvals + (tfinal - tstart)

   
! 
      tstart = omp_get_wtime()
      coeffs1dh = reshape(coeffs2dh, [2*(JMAX-JCN_GRID)*IMAX])
      coeffs1du = reshape(coeffs2du, [2*(JMAX-JCN_GRID)*IMAX])
      coeffs1dv = reshape(coeffs2dv, [2*(JMAX-JCN_GRID)*IMAX])
      coeffs1dw = reshape(coeffs2dw, [2*(JMAX-JCN_GRID)*IMAX])
      tfinal = omp_get_wtime()
      tcoeff = tcoeff + (tfinal - tstart)

      tstart = omp_get_wtime()
      call dfsInterp(coeffs1du,tmpu, N1, N2, th_d, lb_d)
      call dfsInterp(coeffs1dv, tmpv, N1, N2, th_d, lb_d)
      call dfsInterp(coeffs1dh, tmph, N1, N2, th_d, lb_d)
      call dfsInterp(coeffs1dw, tmpw, N1, N2, th_d, lb_d)
      tfinal = omp_get_wtime()
      tfin = tfin + (tfinal - tstart)
     
! transform the values to grid
      tstart = omp_get_wtime()
      if( VELOCITYD .eq. 3 ) then
        !$OMP PARALLEL DO private(i,j,tmpindex, uij, vij)
      do j=1, JMAX
        do i=1, IMAX
          tmpindex = (j-1)*IMAX + i
          uij = COS(lb_d(tmpindex))*tmpv(tmpindex)-SIN(lb_d(tmpindex))*tmpu(tmpindex)
          vij = COS(th_d(tmpindex))*tmpw(tmpindex)-SIN(lb_d(tmpindex))*SIN(th_d(tmpindex))*tmpv(tmpindex)& 
   &         -COS(lb_d(tmpindex))*sin(th_d(tmpindex))*tmpu(tmpindex)
          u0(i,j) = cosdtheta(i,j)*uij - sindtheta(i,j)*vij
          v0(i,j) = sindtheta(i,j)*uij  + cosdtheta(i,j)*vij
        
          phi0(i,j) = tmph(tmpindex)
          print*, 'here'
        end do
      end do
    !$OMP END PARALLEL DO
      tfinal = omp_get_wtime()
      tcopv = tcopv + (tfinal - tstart)
      else
      !$OMP PARALLEL DO private(i,j,tmpindex, tmpcosd)
      do j=1, JMAX
        do i=1, IMAX
          tmpindex = (j-1)*IMAX + i
          tmpcosd  = Dcos(th_d(tmpindex))
          u0(i,j) = (cosdtheta(i,j)*tmpu(tmpindex) - sindtheta(i,j)*tmpv(tmpindex))/tmpcosd
          v0(i,j) = (sindtheta(i,j)*tmpu(tmpindex)  + cosdtheta(i,j)*tmpv(tmpindex))/tmpcosd
          phi0(i,j) = tmph(tmpindex)
        end do
      end do
    !$OMP END PARALLEL DO
    end if
 !
      call e_time__end(5,"main_semilag")
    else

  end if
!
!   ===========================================================
!
    call e_time__start(6,"main_integ2")
    
    if ( iadv_only /= 1 ) then
     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          u0(i,j) = u0(i,j) + dudt0(i,j)*dta - betauv*dudt0_li(i,j)*dta             &
           &       - 2.0d0*omg*ER*( cosa * COSLAT(j) + sina * SINLAT(j) * COSLON(i) )
          v0(i,j) = v0(i,j) + dvdt0(i,j)*dta - betauv*dvdt0_li(i,j)*dta             &
           &       + 2.0d0*omg*ER*sina*SINLON(i)
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
      call semi_implicit               &
         &( dta, beta, betauv, phibar, &!IN
         &  qrot, qdiv, qphi )          !INOUT

      call hdiffusion( dt, qrot, qdiv, qphi )
!
    end if
  
    call e_time__end(6,"main_integ2")
!
  end do

  call e_time__end(2,"main")
  if(IT == 1) then
    
    deallocate(ud_halo)
    deallocate(vd_halo)
    deallocate(hd_halo)
    deallocate(wd_halo)
    deallocate(coeffs1du)
    deallocate(coeffs1dv)
    deallocate(coeffs1dh)
    deallocate(coeffs1dw)
    deallocate(tmpu)
    deallocate(tmph)
    deallocate(tmpv)
    deallocate(tmpw)
  end if
!
  PRINT*, "TIME GRID", tgrid
  PRINT*, "TIME  TVALS", tvals
  PRINT*, "TIME TFIN", tfin
  PRINT*, "TIME COPYV", tcopv
  PRINT*, "TIME COPYCOEFFS",tcoeff
  PRINT*, "NSTEP", nstep
end subroutine main



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine tendency                     &
 &( qrot, qdiv, qphi,                   &!IN
 &  phis, phisx, phisy,                 &!IN
 &  phibar,                             &!IN
 &  u0, v0, phi0,                       &!OUT
 &  u0_adv, v0_adv,                     &!OUT
 &  dudt0, dvdt0, dphidt0,              &!OUT
 &  dudt0_li, dvdt0_li, dphidt0_li   )   !OUT
!
  complex(8),intent(inout) :: qrot(NNUM,0:MMAX)
  complex(8),intent(inout) :: qdiv(NNUM,0:MMAX)
  complex(8),intent(inout) :: qphi(NNUM,0:MMAX)
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
  complex(8) :: qu(NNUM,0:MMAX)
  complex(8) :: qv(NNUM,0:MMAX)
!
!xx  real(8) :: rot(IMAX,JMAX)
  real(8) :: div(IMAX,JMAX)
!xxx  real(8) :: umx(IMAX,JMAX)
  real(8) :: phi0x(IMAX,JMAX)
  real(8) :: phi0y(IMAX,JMAX)
!
  integer :: i,j
!xxx  integer :: n

  call e_time__start(11,"tendency")
!
! ===========================================================
!

  call divrot2uv_dfs__run     &!##  qrot,qdiv -> qu, qv
   &( qdiv, qrot,        &!IN
   &  qu, qv        )     !OUT

  call fft_y__w2g_uv  &
   &( qu,      &!IN
   &  u0     )  !OUT
  call fft_y__w2g_uv  &
   &( qv,      &!IN
   &  v0     )  !OUT

  call fft__togrid &
   &( u0         )     !INOUT
  call fft__togrid &
   &( v0         )     !INOUT
!
  call fft_y__w2g &
   &( qdiv,       &!IN
   &  div    )     !OUT
  call fft__togrid &
   &( div       )      !INOUT
!
  call fft_y__w2g &
   &( qphi,       &!IN
   &  phi0    )    !OUT
  call fft_y__w2g_dy &
   &( qphi,       &!IN
   &  phi0y    )   !OUT
  call xderiv_run              &
   &( phi0, phi0y, coslat_inv, &!IN
   &  phi0x          )          !OUT
!
  call fft__togrid &
   &( phi0       )     !INOUT
  call fft__togrid &
   &( phi0x      )     !INOUT
  call fft__togrid &
   &( phi0y      )     !INOUT
!
  if ( JCN_PROG_H == 0 ) then

   !$OMP PARALLEL default(SHARED), private(i,j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
           dudt0(i,j) = -(phi0x(i,j) + phisx(i,j))/ER
           dvdt0(i,j) = -(phi0y(i,j) + phisy(i,j))/ER
           dudt0_li(i,j) = -phi0x(i,j)/ER
           dvdt0_li(i,j) = -phi0y(i,j)/ER

           dphidt0(i,j)  = -phi0(i,j)*div(i,j)
           dphidt0_li(i,j) = -phibar*div(i,j)

           u0_adv(i,j) = u0(i,j)
           v0_adv(i,j) = v0(i,j)
        end do
     end do
    !$OMP END DO
   !$OMP END PARALLEL

  else
!
   !$OMP PARALLEL default(SHARED), private(i,j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
           dudt0(i,j) = -phi0x(i,j)/ER
           dvdt0(i,j) = -phi0y(i,j)/ER
           dudt0_li(i,j) = dudt0(i,j)
           dvdt0_li(i,j) = dvdt0(i,j)
           dphidt0(i,j) = - ( phi0(i,j) - phis(i,j) )*div(i,j)            &
            &             + ( u0(i,j)*phisx(i,j) + v0(i,j)*phisy(i,j) )/ER
           dphidt0_li(i,j) = -phibar*div(i,j)
           u0_adv(i,j) = u0(i,j)
           v0_adv(i,j) = v0(i,j)
        end do
     end do
    !$OMP END DO
   !$OMP END PARALLEL

  end if

  call e_time__end(11,"tendency")
!
end subroutine tendency


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine grid_to_wave &
 &( u, v, phi,          &
 &  qrot, qdiv, qphi  )
!
  real(8),intent(inout) :: u(IMAX,JMAX)
  real(8),intent(inout) :: v(IMAX,JMAX)
  real(8),intent(inout) :: phi(IMAX,JMAX)
  complex(8),intent(out) :: qrot(NNUM,0:MMAX)
  complex(8),intent(out) :: qdiv(NNUM,0:MMAX)
  complex(8),intent(out) :: qphi(NNUM,0:MMAX)
!
  complex(8) :: qu(NNUM,0:MMAX)
  complex(8) :: qv(NNUM,0:MMAX)
!
!xxx  real(4) :: f_4byte(IMAX,JMAX)

  call e_time__start(12,"grid_to_wave")
!
! ==========================================================
!
!xx  do j=1,JMAX
!xx    do i=1,IMAX
!xx     !###   1/(a*cos2(lat)) * U   ###
!xx      u(i,j) = ACOS2LAT_INV(j)*u(i,j)
!xx     !###   1/(a*cos2(lat)) * V   ###
!xx      v(i,j) = ACOS2LAT_INV(j)*v(i,j)
!xx    end do
!xx  end do
   
!  write(6,*)
!  write(6,*) "0:u(IMAX/4,1:JMAX)=",u(IMAX/4,1:JMAX)
!  write(6,*) "0:v(IMAX/4,1:JMAX)=",v(IMAX/4,1:JMAX)
!
  call fft__towave &
   &( phi       )      !INOUT
  call fft_y__g2w &
   &( phi,        &!IN
   &  qphi    )    !INOUT
!
  call fft__towave &
   &( u          )     !INOUT
  call fft__towave &
   &( v          )     !INOUT

  call fft_y__g2w_uv     &
   &( u,       &!IN
   &  qu    )   !INOUT
  call fft_y__g2w_uv     &
   &( v,       &!IN
   &  qv    )   !INOUT
  
  if ( JCN_DFS >= 1 ) then
    call uv2divrot_dfs__uv2chipsi &
     &( qu, qv,        &!IN
     &  qdiv, qrot  )   !OUT  U, V -> velocity potential, stream function
  else
    call uv2divrot_dfs__run &
     &( qu, qv,        &!IN
     &  qdiv, qrot  )   !OUT  U, V -> divergence, vorticity
  end if
!

!  write(6,*) "qrot(:,1)=",qrot(:,1)
!  write(6,*) "qdiv(:,1)=",qdiv(:,1)
  

!
! ----------------------------------------------------
!
   
!  write(6,*) "qphi(:,1)=",qphi(:,1)

  call e_time__end(12,"grid_to_wave")
!
end subroutine grid_to_wave



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine grid_to_wave_to_grid &
 &( phi,          &
 &  qphi  )
!
  real(8),intent(inout) :: phi(IMAX,JMAX)
  complex(8),intent(out) :: qphi(NNUM,0:MMAX)

  call e_time__start(12,"grid_to_wave")
!
! ==========================================================
!
  call fft__towave &
   &( phi     )    !INOUT
  call fft_y__g2w &
   &( phi,        &!IN
   &  qphi    )    !INOUT

  call fft_y__w2g  &
   &( qphi,    &!IN
   &  phi    )  !OUT
  call fft__togrid &
   &( phi    )  !INOUT

!
! ----------------------------------------------------
!
   
!  write(6,*) "qphi(:,1)=",qphi(:,1)

  call e_time__end(12,"grid_to_wave")
!
end subroutine grid_to_wave_to_grid


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine semi_implicit       &
 &( dt, beta, betauv, phibar,  &!IN
 &  qrot, qdiv, qphi )          !INOUT
!
  real(8),intent(in) :: dt
  real(8),intent(in) :: beta
  real(8),intent(in) :: betauv
  real(8),intent(in) :: phibar
  complex(8),intent(inout) :: qrot(NNUM,0:MMAX)
  complex(8),intent(inout) :: qdiv(NNUM,0:MMAX)
  complex(8),intent(inout) :: qphi(NNUM,0:MMAX)
!
  complex(8) :: qlapphi(NNUM,0:MMAX)
!
  real(8) :: coef
!
  integer :: m,j

  call e_time__start(13,"semi_implicit")
  
  if ( JCN_DFS >= 1 ) then
  
    call uv2divrot_dfs__laplacian(qrot)  !! Stream fuinction => Vorticity

   !$OMP PARALLEL default(SHARED), private(m,j)
   !$OMP DO schedule(STATIC)
    do m = 0,MMAX
      do j = 1,NNUM
        qdiv(j,m) = qdiv(j,m) - betauv*dt*qphi(j,m)  !! Velocity potential
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    coef = beta*betauv*phibar*dt**2

    call helmholtz_dfs__like &
     &( dt, coef,     &!IN
     &  qdiv    )      !INOUT   !! Velocity potential => Divergence

   !$OMP PARALLEL default(SHARED), private(m,j)
   !$OMP DO schedule(STATIC)
    do m = 0,MMAX
      do j = 1,NNUM
        qphi(j,m) = qphi(j,m) - beta*dt*phibar*qdiv(j,m)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
  
  else
    
   !$OMP PARALLEL default(SHARED), private(m,j)
   !$OMP DO schedule(STATIC)
    do m = 0,MMAX
      do j = 1,NNUM
        qlapphi(j,m) = qphi(j,m)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
  
    call uv2divrot_dfs__laplacian(qlapphi)        !INOUT

   !$OMP PARALLEL default(SHARED), private(m,j)
   !$OMP DO schedule(STATIC)
    do m = 0,MMAX
      do j = 1,NNUM
        qdiv(j,m) = qdiv(j,m) - betauv*dt*qlapphi(j,m)  !! Divergence
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    coef = beta*betauv*phibar*dt**2

    call helmholtz_dfs__run &!harmonic spectral filter
     &( dt, coef,     &!IN
     &  qdiv    )  !INOUT

   !$OMP PARALLEL default(SHARED), private(m,j)
   !$OMP DO schedule(STATIC)
    do m = 0,MMAX
      do j = 1,NNUM
        qphi(j,m) = qphi(j,m) - beta*dt*phibar*qdiv(j,m)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
  end if
!
  call e_time__end(13,"semi_implicit")
!
end subroutine semi_implicit



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine hdiffusion( dt, qrot, qdiv, qphi )
  real(8),intent(in) :: dt
  complex(8),intent(inout) :: qrot(NNUM,0:MMAX)
  complex(8),intent(inout) :: qdiv(NNUM,0:MMAX)
  complex(8),intent(inout) :: qphi(NNUM,0:MMAX)

  real(8) :: coef, tau, factor
  call e_time__start(14,"hdiffusion")

  if ( JCN_HDIFF == 1 ) then
    !         if ( kind_phi == 8 ) then
                !! 2nd order diffusion for Galewsky (2004) test case
                coef = 1.0d5*dt
                call helmholtz_dfs__hdiff2 &
                 &( dt, coef,     &!IN
                 &  qdiv    )  !INOUT
                call helmholtz_dfs__hdiff2 &
                 &( dt, coef,     &!IN
                 &  qrot    )  !INOUT
                call helmholtz_dfs__hdiff2 &
                 &( dt, coef,     &!IN
                 &  qphi    )  !INOUT
    !         end if
          else if ( JCN_HDIFF == 2 ) then
             !! 4th order hyper diffusion
             tau  = 7.2D0*(107.0D0/(NMAX+1))**2
             factor = 1.0d0
             coef = factor*dt*ER**4/(tau*3600.0D0*(real(NMAX,kind=8)*(NMAX+1))**2)
             call hdiff4_dfs__run &
              &( dt, coef,     &!IN
              &  qdiv    )  !INOUT
             call hdiff4_dfs__run &
              &( dt, coef,     &!IN
              &  qrot    )  !INOUT
             call hdiff4_dfs__run &
              &( dt, coef,     &!IN
              &  qphi    )  !INOUT
          end if
    
      call e_time__end(14,"hdiffusion")
    
 
    
end subroutine hdiffusion

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine cal_gmean &
 &( data,           &
 &  gmean)
!
  real(8),intent(in) :: data(IMAX,JMAX)
  real(8),intent(out) :: gmean      !! Global mean
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
!
    denominator = 0.0d0
    do j=1,JMAX
      denominator = denominator + WEIGHT(j)
    end do
    denominator = denominator*IMAX
!
!xx    write(6,*) 'COSLAT(:)/WEIGHT(:)=',COSLAT(:)/WEIGHT(:)
!
  end if
!
!
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



subroutine departure          &
 &( j, i1, i2, iadv_only,     &!IN
 &  iterx, idepart_prevstep,  &!IN
 &  dt, dta,                  &!IN
 &  cosa, sina,               &!IN
 &  u0_adv, v0_adv,           &!IN
 &  dudt0, dvdt0,             &!IN
 &  u_adv_halo, v_adv_halo,   &!IN
 &  utmp, vtmp, uwork, vwork, &!OUT
 &  ii, jj,                   &!OUT
 &  xi, yj,                   &!OUT
 &  tt, ll,                   &!OUT
 &  cosdtheta, sindtheta )     !OUT
!
  integer,intent(in) :: j
  integer,intent(in) :: i1
  integer,intent(in) :: i2
  integer,intent(in) :: iadv_only
  integer,intent(in) :: iterx
  integer,intent(in) :: idepart_prevstep
  real(8),intent(in) :: dt
  real(8),intent(in) :: dta
  real(8),intent(in) :: cosa
  real(8),intent(in) :: sina
  real(8),intent(in) :: u0_adv(IMAX,JMAX)
  real(8),intent(in) :: v0_adv(IMAX,JMAX)
  real(8),intent(in) :: dudt0(IMAX,JMAX)
  real(8),intent(in) :: dvdt0(IMAX,JMAX)
  real(8),intent(in) :: u_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(in) :: v_adv_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
!
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
  real(8), intent(inout) :: tt(IMAX)
  real(8), intent(inout) :: ll(IMAX)
!
  real(8) :: dt_ER
  real(8) :: rdis2, uu, vv, cosdis, alond, alatd
  real(8) :: sindis_rdis, sinlatd, coslatd_cosdlond, coslatd_sindlond
  real(8) :: coslatd, ww1, ww2, ww3, cosdlond, sindlond
! call dfsInterp(coeffs1dh, tmph, N1, N2, th_d, lb_d)
  integer :: iter
  integer :: i
  
!xx  call e_time__start(6,"departure")
!
! ===============================================================
!


!xx!$OMP PARALLEL default(SHARED), private( i, j, rdis2, uu, vv, cosdis, alond, alatd ) &
!xx!$OMP  & private( sindis_rdis, sinlatd, coslatd_cosdlond, coslatd_sindlond ) &
!xx!$OMP  & private( coslatd, ww1, ww2, ww3, cosdlond, sindlond ) &
!xx!$OMP  & private( iter, dt_ER, utmp, vtmp )
!xx!$OMP DO schedule(STATIC)
!xx  do j=1,JMAX
  
    if ( iadv_only /= 1 .and. JCN_DEPARTURE /= 0 ) then
      do i=i1,i2
        uwork(i) = 0.5d0*dudt0(i,j)*dta    &
         &         - omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
        vwork(i) = 0.5d0*dvdt0(i,j)*dta    &
         &         + omg*ER*sina*SINLON(i)
      end do
    end if

    iter_loop : do iter = 1, iterx
!
      if ( idepart_prevstep == 0 .and. iter == 1 ) then
        !! Not use departure point data of previous step
        do i=i1,i2
          utmp(i) = u0_adv(i,j)
          vtmp(i) = v0_adv(i,j)
        end do
      else
        !! Use departure point data of previous step
        call lag3_uv                &
         &( i1, i2,                 &!IN
         &  ii, jj,                 &!IN
         &  xi, yj,                 &!IN
         &  cosdtheta, sindtheta,   &!IN
         &  u_adv_halo, v_adv_halo, &!IN
         &  utmp, vtmp      )      !INOUT
!
        if ( iadv_only == 1 ) then
          do i=i1,i2
            utmp(i) = ( utmp(i) + u0_adv(i,j) )*0.5d0
            vtmp(i) = ( vtmp(i) + v0_adv(i,j) )*0.5d0
          end do
!
        else

          if ( JCN_DEPARTURE == 0 ) then
            do i=i1,i2
              utmp(i) = ( utmp(i) + u0_adv(i,j) )*0.5d0
              vtmp(i) = ( vtmp(i) + v0_adv(i,j) )*0.5d0
            end do

          else
            do i=i1,i2
              utmp(i) = utmp(i) + 0.5d0*dudt0(i,j)*dta    &
               &          - omg*ER*( cosa*COSLAT(j) + sina*SINLAT(j)*COSLON(i) )
              vtmp(i) = vtmp(i) + 0.5d0*dvdt0(i,j)*dta    &
               &          + omg*ER*sina*SINLON(i)
            end do

          end if
        end if
      end if

      dt_er = dt/ER
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
        !tt(i) = atan2(sinlatd,coslatd) 

        ww1 = sign(0.5d0,coslatd+1.0d-50) - sign(0.5d0,coslatd-1.0d-50)
        ww2 = 1.0d0/( ww1 + coslatd )
        cosdlond = ( ww1 + coslatd_cosdlond )*ww2
        sindlond = coslatd_sindlond *ww2
        alond    = ALON(i) + atan2(sindlond,cosdlond) + PI*2.0d0
        !ll(i) =  ALON(i) + atan2(sindlond,cosdlond) + PI*2.0d0
        ww3 = 1.0d0/( 1.0d0 + cosdis )
        cosdtheta(i) = ww3*( COSLAT(j)*coslatd + (1.0d0+SINLAT(j)*sinlatd)*cosdlond )
        sindtheta(i) = ww3*( (SINLAT(j)+sinlatd)*sindlond )

        xi(i) = alond/DLON + 1.0d0
        !xi(i) = ll(i)/DLON + 1.0d0
        ii(i) = int( xi(i) )
        xi(i) = xi(i) - ii(i)
        ii(i) = mod(ii(i)-1,IMAX) + 1

        yj(i) = -alatd/DLAT + (JMAX+1)*0.5d0
       ! yj(i) = -tt(i)/DLAT + (JMAX+1)*0.5d0
        jj(i) = int( yj(i) )
        yj(i) = yj(i) - jj(i)

        if(iter == iterx) then
          tt(i) = alatd
          ll(i) = alond
        end if
        
      end do
      

    end do iter_loop
  
!xx  end do
!xx!$OMP END DO
!xx!$OMP END PARALLEL
  
!xx  call e_time__end(6,"departure")
!
!
end subroutine departure



!***************************************************************************


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
!
  real(8),intent(in) :: phi_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(inout) :: phi(IMAX)
!
  real(8) :: xx1,xx2,xx3,xx4,xxf,xxg
  real(8) :: yy1,yy2,yy3,yy4,yyf,yyg
  real(8) :: phi1,phi2,phi3,phi4
!
  integer :: i,ni,nj
!
!  call e_time__start(7,"lag3")
! ===========================================================================
!
!xx!$OMP PARALLEL default(SHARED), private(i,j,ni,nj,xx1,xx2,xx3,xx4,xxf,xxg,yy1,yy2,yy3,yy4,yyf,yyg,phi1,phi2,phi3,phi4)
!xx !$OMP DO schedule(STATIC)
!xx  do j = 1, JMAX
    do i = i1,i2
 
      ni = ii(i)
      nj = jj(i)   
    
      xx1 = 2.0d0-xi(i)
      xx2 = 1.0d0-xi(i)
      xx3 = xi(i)
      xx4 = 1.0d0+xi(i)
      xxf = xx1*xx4*0.5d0
      xxg = -xx2*xx3/6.0d0
      xx1 = xx1*xxg
      xx2 = xx2*xxf
      xx3 = xx3*xxf
      xx4 = xx4*xxg
!
      yy1 = 2.0d0-yj(i)
      yy2 = 1.0d0-yj(i)
      yy3 = yj(i)
      yy4 = 1.0d0+yj(i)
      yyf = yy1*yy4*0.5d0
      yyg = -yy2*yy3/6.0d0
      yy1 = yy1*yyg
      yy2 = yy2*yyf
      yy3 = yy3*yyf
      yy4 = yy4*yyg

      phi1 = xx1*phi_halo(ni-1,nj-1) + xx2*phi_halo(ni  ,nj-1)  &
       &     + xx3*phi_halo(ni+1,nj-1) + xx4*phi_halo(ni+2,nj-1)
      phi2 = xx1*phi_halo(ni-1,nj  ) + xx2*phi_halo(ni  ,nj  )  &
       &     + xx3*phi_halo(ni+1,nj  ) + xx4*phi_halo(ni+2,nj  )
      phi3 = xx1*phi_halo(ni-1,nj+1) + xx2*phi_halo(ni  ,nj+1)  &
       &     + xx3*phi_halo(ni+1,nj+1) + xx4*phi_halo(ni+2,nj+1)
      phi4 = xx1*phi_halo(ni-1,nj+2) + xx2*phi_halo(ni  ,nj+2)  &
       &     + xx3*phi_halo(ni+1,nj+2) + xx4*phi_halo(ni+2,nj+2)
    
      phi(i) = yy1*phi1 + yy2*phi2 + yy3*phi3 + yy4*phi4
    end do
!xx  end do
!xx !$OMP END DO
!xx!$OMP END PARALLEL
!
!  call e_time__end(7,"lag3")
!
end subroutine lag3


!***************************************************************************


subroutine lag3_uv        &
 &( i1, i2,               &!IN
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
!
  real(8),intent(in) :: u_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(in) :: v_halo(1-MGN_X:IMAX+MGN_X, 1-MGN_Y:JMAX+MGN_Y)
  real(8),intent(inout) :: u(IMAX)
  real(8),intent(inout) :: v(IMAX)
!
  integer :: i,ni,nj
  real(8) :: u1,u2,u3,u4,v1,v2,v3,v4
  real(8) :: uu,vv
  real(8) :: xx1,xx2,xx3,xx4,xxf,xxg,yy1,yy2,yy3,yy4,yyf,yyg
!
!xx  call e_time__start(9,"lag3_uv")
! ===========================================================================
!
!xx!$OMP PARALLEL default(SHARED), private(i,j,ni,nj,u1,u2,u3,u4,v1,v2,v3,v4,uu,vv) &
!xx!$OMP  & private(xx1,xx2,xx3,xx4,xxf,xxg,yy1,yy2,yy3,yy4,yyf,yyg)
!xx !$OMP DO schedule(STATIC)
!xx  do j = 1, JMAX
    do i = i1, i2
!
      ni = ii(i)
      nj = jj(i)
      
      xx1 = 2.0d0-xi(i)
      xx2 = 1.0d0-xi(i)
      xx3 = xi(i)
      xx4 = 1.0d0+xi(i)
      xxf = xx1*xx4*0.5d0
      xxg = -xx2*xx3/6.0d0
      xx1 = xx1*xxg
      xx2 = xx2*xxf
      xx3 = xx3*xxf
      xx4 = xx4*xxg
!
      yy1 = 2.0d0-yj(i)
      yy2 = 1.0d0-yj(i)
      yy3 = yj(i)
      yy4 = 1.0d0+yj(i)
      yyf = yy1*yy4*0.5d0
      yyg = -yy2*yy3/6.0d0
      yy1 = yy1*yyg
      yy2 = yy2*yyf
      yy3 = yy3*yyf
      yy4 = yy4*yyg

      u1 =  xx1*u_halo(ni-1,nj-1) + xx2*u_halo(ni  ,nj-1)  &
       &  + xx3*u_halo(ni+1,nj-1) + xx4*u_halo(ni+2,nj-1)
      u2 =  xx1*u_halo(ni-1,nj  ) + xx2*u_halo(ni  ,nj  )  &
       &  + xx3*u_halo(ni+1,nj  ) + xx4*u_halo(ni+2,nj  )
      u3 =  xx1*u_halo(ni-1,nj+1) + xx2*u_halo(ni  ,nj+1)  &
       &  + xx3*u_halo(ni+1,nj+1) + xx4*u_halo(ni+2,nj+1)
      u4 =  xx1*u_halo(ni-1,nj+2) + xx2*u_halo(ni  ,nj+2)  &
       &  + xx3*u_halo(ni+1,nj+2) + xx4*u_halo(ni+2,nj+2)

      v1 =  xx1*v_halo(ni-1,nj-1) + xx2*v_halo(ni  ,nj-1)  &
       &  + xx3*v_halo(ni+1,nj-1) + xx4*v_halo(ni+2,nj-1)
      v2 =  xx1*v_halo(ni-1,nj  ) + xx2*v_halo(ni  ,nj  )  &
       &  + xx3*v_halo(ni+1,nj  ) + xx4*v_halo(ni+2,nj  )
      v3 =  xx1*v_halo(ni-1,nj+1) + xx2*v_halo(ni  ,nj+1)  &
       &  + xx3*v_halo(ni+1,nj+1) + xx4*v_halo(ni+2,nj+1)
      v4 =  xx1*v_halo(ni-1,nj+2) + xx2*v_halo(ni  ,nj+2)  &
       &  + xx3*v_halo(ni+1,nj+2) + xx4*v_halo(ni+2,nj+2)
!
      uu = yy1*u1 + yy2*u2 + yy3*u3 + yy4*u4
      vv = yy1*v1 + yy2*v2 + yy3*v3 + yy4*v4
!
      u(i)  = cosdtheta(i)*uu  - sindtheta(i)*vv
      v(i)  = sindtheta(i)*uu  + cosdtheta(i)*vv
!
    end do
!xx  end do
!xx !$OMP END DO
!xx!$OMP END PARALLEL
!
! ------------------------------------------------------------------------
!  
!xx  call e_time__end(9,"lag3_uv")
!
end subroutine lag3_uv


!***************************************************************************


subroutine lag5_uv        &
 &( i1, i2,               &!IN
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
  real(8) :: w1,w2,w3,a1,a2,a3,a4,a5,a6
  real(8) :: u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6
  real(8) :: uu,vv
  real(8) :: xx,yy
!
!xx  call e_time__start(10,"lag5_uv")
! ===========================================================================
!
!xx!$OMP PARALLEL default(SHARED), private(i,j,ni,nj,w1,w2,w3,a1,a2,a3,a4,a5,a6,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6) &
!xx!$OMP  & private(uu,vv,xx,yy)
!
!xx !$OMP DO schedule(STATIC)
!xx  do j = 1, JMAX

    do i = i1,i2
!
      xx = xi(i)
      yy = yj(i)
      ni = ii(i)
      nj = jj(i)

      w1 = (xx+2)*(xx+1)
      w2 = (xx  )*(xx-1)
      w3 = (xx-2)*(xx-3)
      a1 = (xx+1)*w2*w3/((-1)*(-2)*(-3)*(-4)*(-5))
      a2 = (xx+2)*w2*w3/(( 1)*(-1)*(-2)*(-3)*(-4))
      a3 = w1*(xx-1)*w3/(( 2)*( 1)*(-1)*(-2)*(-3))
      a4 = w1*(xx  )*w3/(( 3)*( 2)*( 1)*(-1)*(-2))
      a5 = w1*w2*(xx-3)/(( 4)*( 3)*( 2)*( 1)*(-1))
      a6 = w1*w2*(xx-2)/(( 5)*( 4)*( 3)*( 2)*( 1))

!      a1 = (xx+1)*(xx  )*(xx-1)*(xx-2)*(xx-3)/((-1)*(-2)*(-3)*(-4)*(-5))
!      a2 = (xx+2)*(xx  )*(xx-1)*(xx-2)*(xx-3)/(( 1)*(-1)*(-2)*(-3)*(-4))
!      a3 = (xx+2)*(xx+1)*(xx-1)*(xx-2)*(xx-3)/(( 2)*( 1)*(-1)*(-2)*(-3))
!      a4 = (xx+2)*(xx+1)*(xx  )*(xx-2)*(xx-3)/(( 3)*( 2)*( 1)*(-1)*(-2))
!      a5 = (xx+2)*(xx+1)*(xx  )*(xx-1)*(xx-3)/(( 4)*( 3)*( 2)*( 1)*(-1))
!      a6 = (xx+2)*(xx+1)*(xx  )*(xx-1)*(xx-2)/(( 5)*( 4)*( 3)*( 2)*( 1))
!
      u1 =   a1*u_halo(ni-2,nj-2) + a2*u_halo(ni-1,nj-2)  &
       &   + a3*u_halo(ni  ,nj-2) + a4*u_halo(ni+1,nj-2)  &
       &   + a5*u_halo(ni+2,nj-2) + a6*u_halo(ni+3,nj-2)
!
      u2 =   a1*u_halo(ni-2,nj-1) + a2*u_halo(ni-1,nj-1)  &
       &   + a3*u_halo(ni  ,nj-1) + a4*u_halo(ni+1,nj-1)  &
       &   + a5*u_halo(ni+2,nj-1) + a6*u_halo(ni+3,nj-1)
!
      u3 =   a1*u_halo(ni-2,nj  ) + a2*u_halo(ni-1,nj  )  &
       &   + a3*u_halo(ni  ,nj  ) + a4*u_halo(ni+1,nj  )  &
       &   + a5*u_halo(ni+2,nj  ) + a6*u_halo(ni+3,nj  )
!
      u4 =   a1*u_halo(ni-2,nj+1) + a2*u_halo(ni-1,nj+1)  &
       &   + a3*u_halo(ni  ,nj+1) + a4*u_halo(ni+1,nj+1)  &
       &   + a5*u_halo(ni+2,nj+1) + a6*u_halo(ni+3,nj+1)
!
      u5 =   a1*u_halo(ni-2,nj+2) + a2*u_halo(ni-1,nj+2)  &
       &   + a3*u_halo(ni  ,nj+2) + a4*u_halo(ni+1,nj+2)  &
       &   + a5*u_halo(ni+2,nj+2) + a6*u_halo(ni+3,nj+2)
!
      u6 =   a1*u_halo(ni-2,nj+3) + a2*u_halo(ni-1,nj+3)  &
       &   + a3*u_halo(ni  ,nj+3) + a4*u_halo(ni+1,nj+3)  &
       &   + a5*u_halo(ni+2,nj+3) + a6*u_halo(ni+3,nj+3)
      
      v1 =   a1*v_halo(ni-2,nj-2) + a2*v_halo(ni-1,nj-2)  &
       &   + a3*v_halo(ni  ,nj-2) + a4*v_halo(ni+1,nj-2)  &
       &   + a5*v_halo(ni+2,nj-2) + a6*v_halo(ni+3,nj-2)
!
      v2 =   a1*v_halo(ni-2,nj-1) + a2*v_halo(ni-1,nj-1)  &
       &   + a3*v_halo(ni  ,nj-1) + a4*v_halo(ni+1,nj-1)  &
       &   + a5*v_halo(ni+2,nj-1) + a6*v_halo(ni+3,nj-1)
!
      v3 =   a1*v_halo(ni-2,nj  ) + a2*v_halo(ni-1,nj  )  &
       &   + a3*v_halo(ni  ,nj  ) + a4*v_halo(ni+1,nj  )  &
       &   + a5*v_halo(ni+2,nj  ) + a6*v_halo(ni+3,nj  )
!
      v4 =   a1*v_halo(ni-2,nj+1) + a2*v_halo(ni-1,nj+1)  &
       &   + a3*v_halo(ni  ,nj+1) + a4*v_halo(ni+1,nj+1)  &
       &   + a5*v_halo(ni+2,nj+1) + a6*v_halo(ni+3,nj+1)
!
      v5 =   a1*v_halo(ni-2,nj+2) + a2*v_halo(ni-1,nj+2)  &
       &   + a3*v_halo(ni  ,nj+2) + a4*v_halo(ni+1,nj+2)  &
       &   + a5*v_halo(ni+2,nj+2) + a6*v_halo(ni+3,nj+2)
!
      v6 =   a1*v_halo(ni-2,nj+3) + a2*v_halo(ni-1,nj+3)  &
       &   + a3*v_halo(ni  ,nj+3) + a4*v_halo(ni+1,nj+3)  &
       &   + a5*v_halo(ni+2,nj+3) + a6*v_halo(ni+3,nj+3)
!
      w1 = (yy+2)*(yy+1)
      w2 = (yy  )*(yy-1)
      w3 = (yy-2)*(yy-3)
      a1 = (yy+1)*w2*w3/((-1)*(-2)*(-3)*(-4)*(-5))
      a2 = (yy+2)*w2*w3/(( 1)*(-1)*(-2)*(-3)*(-4))
      a3 = w1*(yy-1)*w3/(( 2)*( 1)*(-1)*(-2)*(-3))
      a4 = w1*(yy  )*w3/(( 3)*( 2)*( 1)*(-1)*(-2))
      a5 = w1*w2*(yy-3)/(( 4)*( 3)*( 2)*( 1)*(-1))
      a6 = w1*w2*(yy-2)/(( 5)*( 4)*( 3)*( 2)*( 1))
      
!      a1 = (yy+1)*(yy  )*(yy-1)*(yy-2)*(yy-3)/((-1)*(-2)*(-3)*(-4)*(-5))
!      a2 = (yy+2)*(yy  )*(yy-1)*(yy-2)*(yy-3)/(( 1)*(-1)*(-2)*(-3)*(-4))
!      a3 = (yy+2)*(yy+1)*(yy-1)*(yy-2)*(yy-3)/(( 2)*( 1)*(-1)*(-2)*(-3))
!      a4 = (yy+2)*(yy+1)*(yy  )*(yy-2)*(yy-3)/(( 3)*( 2)*( 1)*(-1)*(-2))
!      a5 = (yy+2)*(yy+1)*(yy  )*(yy-1)*(yy-3)/(( 4)*( 3)*( 2)*( 1)*(-1))
!      a6 = (yy+2)*(yy+1)*(yy  )*(yy-1)*(yy-2)/(( 5)*( 4)*( 3)*( 2)*( 1))
        
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

!xx  call e_time__end(10,"lag5_uv")
!
end subroutine lag5_uv


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine set_halo &
 & ( data_halo )  !OUT
!
  real(8),intent(inout) :: data_halo(1-MGN_X:IMAX+MGN_X,1-MGN_Y:JMAX+MGN_Y)
!
  real(8) :: zmean11, zmean12, zmean21, zmean22, zmean31, zmean32, zmean41, zmean42
  real(8) :: data_np, data_sp
  integer :: m, i, j
!
!$OMP PARALLEL default(SHARED), private(m,i,j)

  if ( JCN_GRID == 1 ) then
   !$OMP DO schedule(STATIC)
    do m=1,MGN_Y
      do i=1,IMAX2
        data_halo(i      ,1-m) = data_halo(i+IMAX2,1+m)
        data_halo(i+IMAX2,1-m) = data_halo(i      ,1+m)
        data_halo(i      ,JMAX+m) = data_halo(i+IMAX2,JMAX-m)
        data_halo(i+IMAX2,JMAX+m) = data_halo(i      ,JMAX-m)
      end do
    end do
   !$OMP END DO
   
  else if ( JCN_GRID == -1 ) then
   !$OMP DO schedule(STATIC)
    do m=2,MGN_Y
      do i=1,IMAX2
        data_halo(i      ,1-m) = data_halo(i+IMAX2,m-1)
        data_halo(i+IMAX2,1-m) = data_halo(i      ,m-1)
        data_halo(i      ,JMAX+m) = data_halo(i+IMAX2,JMAX+2-m)
        data_halo(i+IMAX2,JMAX+m) = data_halo(i      ,JMAX+2-m)
      end do
    end do
   !$OMP END DO
    
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

    data_np = ( 4.0d0*(zmean11+zmean12) - (zmean21+zmean22) )/6.0d0
    data_sp = ( 4.0d0*(zmean41+zmean42) - (zmean31+zmean32) )/6.0d0
   !$OMP DO schedule(STATIC)
    do i=1,IMAX
      data_halo(i,0) = data_np
      data_halo(i,JMAX+1) = data_sp
    end do
   !$OMP END DO
  else
   !$OMP DO schedule(STATIC)
    do m=1,MGN_Y
      do i=1,IMAX2
        data_halo(i      ,1-m) = data_halo(i+IMAX2,m)
        data_halo(i+IMAX2,1-m) = data_halo(i      ,m)
        data_halo(i      ,JMAX+m) = data_halo(i+IMAX2,JMAX+1-m)
        data_halo(i+IMAX2,JMAX+m) = data_halo(i      ,JMAX+1-m)
      end do
    end do
   !$OMP END DO
  end if

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


!***************************************************************************


subroutine set_halo_uv  &
 & ( data_halo  )   !OUT
!
  real(8),intent(inout) :: data_halo(1-MGN_X:IMAX+MGN_X,1-MGN_Y:JMAX+MGN_Y)
!
  real(8) :: zonal_cos11, zonal_cos21, zonal_cos31, zonal_cos41
  real(8) :: zonal_cos12, zonal_cos22, zonal_cos32, zonal_cos42
  real(8) :: zonal_sin11, zonal_sin21, zonal_sin31, zonal_sin41
  real(8) :: zonal_sin12, zonal_sin22, zonal_sin32, zonal_sin42
  real(8) :: zonal_np_cos, zonal_np_sin, zonal_sp_cos, zonal_sp_sin
  integer :: m, i, j
!
 !$OMP PARALLEL default(SHARED), private(m,i,j)

  if ( JCN_GRID == 1 ) then
   !$OMP DO schedule(STATIC)
    do m=1,MGN_Y
      do i=1,IMAX2
        data_halo(i      ,1-m) = -data_halo(i+IMAX2,1+m)
        data_halo(i+IMAX2,1-m) = -data_halo(i      ,1+m)
        data_halo(i      ,JMAX+m) = -data_halo(i+IMAX2,JMAX-m)
        data_halo(i+IMAX2,JMAX+m) = -data_halo(i      ,JMAX-m)
      end do
    end do
   !$OMP END DO

  else if ( JCN_GRID == -1 ) then

   !$OMP DO schedule(STATIC)
    do m=2,MGN_Y
      do i=1,IMAX2
        data_halo(i      ,1-m) = -data_halo(i+IMAX2,m-1)
        data_halo(i+IMAX2,1-m) = -data_halo(i      ,m-1)
        data_halo(i      ,JMAX+m) = -data_halo(i+IMAX2,JMAX+2-m)
        data_halo(i+IMAX2,JMAX+m) = -data_halo(i      ,JMAX+2-m)
      end do
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

    zonal_np_cos = ( 4.0d0*(zonal_cos11+zonal_cos12) - (zonal_cos21+zonal_cos22) )/6.0d0
    zonal_np_sin = ( 4.0d0*(zonal_sin11+zonal_sin12) - (zonal_sin21+zonal_sin22) )/6.0d0
    zonal_sp_cos = ( 4.0d0*(zonal_cos41+zonal_cos42) - (zonal_cos31+zonal_cos32) )/6.0d0
    zonal_sp_sin = ( 4.0d0*(zonal_sin41+zonal_sin42) - (zonal_sin31+zonal_sin32) )/6.0d0

   !$OMP DO schedule(STATIC)
    do i=1,IMAX
      data_halo(i,0) = zonal_np_cos*COSLON(i) + zonal_np_sin*SINLON(i)
      data_halo(i,JMAX+1) = zonal_sp_cos*COSLON(i) + zonal_sp_sin*SINLON(i)
    end do
   !$OMP END DO
    
    
!    write(6,*)
!    write(6,*) "data(1,-2:2)=",data(1,-2:2)
!    write(6,*) "data(IMAX/4,-2:2)=",data(IMAX/4,-2:2)
!    write(6,*)
!    write(6,*) "data(1,JMAX-1:JMAX+3)=",data(1,JMAX-1:JMAX+3)
!    write(6,*) "data(IMAX/4,JMAX-1:JMAX+3)=",data(IMAX/4,JMAX-1:JMAX+3)
!    write(6,*)
    
  else

   !$OMP DO schedule(STATIC)
    do m=1,MGN_Y
      do i=1,IMAX2
        data_halo(i      ,1-m) = -data_halo(i+IMAX2,m)
        data_halo(i+IMAX2,1-m) = -data_halo(i      ,m)
        data_halo(i      ,JMAX+m) = -data_halo(i+IMAX2,JMAX+1-m)
        data_halo(i+IMAX2,JMAX+m) = -data_halo(i      ,JMAX+1-m)
      end do
    end do
   !$OMP END DO

  end if
  
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


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


end program sw_dfs