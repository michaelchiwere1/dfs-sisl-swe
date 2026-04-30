program adv_dfs

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
!
  implicit none

!  integer,parameter :: JCN_DFS = 0 !! Old method and Cheong's basis functions
!  integer,parameter :: JCN_DFS = 1 !! New method to calculate expantion coefficients
  integer,parameter :: JCN_DFS = 2 !! New method to calculate expantion coefficients and New basis functions

  integer,parameter :: JCN_GRID = 0   !! colat=PI*(j+0.5)/JMAX0, j=0,...,JMAX0-1
!  integer,parameter :: JCN_GRID = 1   !! colat=PI*j/JMAX0, j=0,...,JMAX0   (JMAX=JMAX0+1)
!  integer,parameter :: JCN_GRID = -1  !! colat=PI*j/JMAX0, j=1,...,JMAX0-1 (JMAX=JMAX0-1)
  
  integer,parameter :: JCN_INITIAL = 1  !! Williamson test case 1 (Advection)
!  integer,parameter :: JCN_INITIAL = 2  !! Williamson test case 2
!  integer,parameter :: JCN_INITIAL = 5  !! Williamson test case 5
!  integer,parameter :: JCN_INITIAL = 6  !! Williamson test case 6
!  integer,parameter :: JCN_INITIAL = 10 !! Galewsky 2004 test case

  integer,parameter :: JCN_HDIFF = 0  !! No diffusion
!  integer,parameter :: JCN_HDIFF = 1  !! 2nd order diff. for Galewsky test case
!  integer,parameter :: JCN_HDIFF = 2  !! 4th order hyper diffusion
  
!  integer,parameter :: JCN_ZONALFILTER = 0 !! No zonal filter
  integer,parameter :: JCN_ZONALFILTER = 1 !! Use zonal Fourier filter
  integer,parameter :: M0 = 1              !! For zonal Fourier filter

!  integer,parameter :: JCN_MONIT = 0  !! No monitor output
  integer,parameter :: JCN_MONIT = 1  !! Write monitor output

  integer,parameter :: INTHR_MONIT = 24   !! Interval of monitor output (hour)

!  integer,parameter :: IMAX = 16
!  integer,parameter :: JMAX = 8 + JCN_GRID
!  integer,parameter :: NMAX = 7
!  real(8),save :: TIMESTEP = 3600.0d0         !! 1 hour
!
  integer,parameter :: IMAX = 128            !! About 300km resolution
  integer,parameter :: JMAX = 64+JCN_GRID
  integer,parameter :: NMAX = 42             !! Quadric grid
  real(8),save :: TIMESTEP = 1800.0d0        !! 30 min.
!
!  integer,parameter :: IMAX = 160            !! 250km resolution
!  integer,parameter :: JMAX = 80+JCN_GRID
!  integer,parameter :: NMAX = 53             !! Quadric grid
!  real(8),save :: TIMESTEP = 1800.0d0        !! 30 min.
!
!  integer,parameter :: IMAX = 320            !! 120km resolution
!  integer,parameter :: JMAX = 160+JCN_GRID
!  integer,parameter :: NMAX = 106            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 900.0d0    !! 15 min.
!
!  integer,parameter :: IMAX = 640            !! 60km resolution
!  integer,parameter :: JMAX = 320+JCN_GRID
!  integer,parameter :: NMAX = 213            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 450.0d0    !! 7.5 min.
!
!  integer,parameter :: IMAX = 1920           !! 20km resolution
!  integer,parameter :: JMAX = 960+JCN_GRID   !! J is from North to South
!  integer,parameter :: NMAX = 639            !! Quadric grid             
!  integer,parameter :: NMAX = 959            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 150.0d0    !! 2.5 min.
!  real(8),parameter :: TIMESTEP = 200.0d0    !! 3.33 min.
!  real(8),parameter :: TIMESTEP = 240.0d0    !! 4.0 min.
!  real(8),parameter :: TIMESTEP = 300.0d0    !! 5.0 min.
!
!  integer,parameter :: IMAX = 3840           !! 10km resolution
!  integer,parameter :: JMAX = 1920+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 1279           !! Quadric grid
!  integer,parameter :: NMAX = 959            !! Cubic grid
!  real(8),parameter :: TIMESTEP = 80.0d0     !! 80 sec.
!
!  integer,parameter :: IMAX = 7680           !! 5km resolution
!  integer,parameter :: JMAX = 3840+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 2559           !! Quadric grid
!  integer,parameter :: NMAX = 1919           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 40.0d0     !! 40 sec.
!
!  integer,parameter :: IMAX = 15360          !! 2.6km resolution
!  integer,parameter :: JMAX = 7680+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 5119           !! Quadric grid
!  integer,parameter :: NMAX = 3839           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 20.0d0    !!  20 sec.
!
!  integer,parameter :: IMAX = 20480          !! 2.0km resolution
!  integer,parameter :: JMAX = 10240+JCN_GRID !! J is from North to South
!  integer,parameter :: NMAX = 6826           !! Quadric grid
!  integer,parameter :: NMAX = 5119           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 15.0d0     !! 15 sec.
!
!  integer,parameter :: IMAX = 30720          !! 1.3km resolution
!  integer,parameter :: JMAX = 15360+JCN_GRID !! J is from North to South
!  integer,parameter :: NMAX = 10239          !! Quadric grid
!  integer,parameter :: NMAX = 7679           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 10.0d0     !! 10 sec.
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
  
  real(8),save :: eps_save = 0.05d0  !## 0.0 <= eps_save <= 0.1?
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
  write(6,*) "JCN_HDIFF          =", JCN_HDIFF
  write(6,*) "JCN_ZONALFILTER    =", JCN_ZONALFILTER
  write(6,*) "M0                 =", M0
  write(6,*) "JCN_MONIT          =", JCN_MONIT
  write(6,*) "INTHR_MONIT        =", INTHR_MONIT
  write(6,*) "IMAX     =", IMAX 
  write(6,*) "JMAX     =", JMAX
  write(6,*) "NMAX     =", NMAX
  write(6,*) "TIMESTEP =", TIMESTEP
  
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
  !! Expansion coefficients
  !!  Real part: [cosine omponent]/2, Imaginary part: -[sine component]/2 for m >= 1
  complex(8),save :: qphi(NNUM,0:MMAX)  !! GRAV*height (geopotential height)
!
  real(8),save :: phim(IMAX,JMAX)       !! GRAV*height
  real(8),save :: phi0(IMAX,JMAX)       !! GRAV*height
  real(8),save :: phi(IMAX,JMAX)       !! GRAV*height
!
  real(8),save :: u0_adv(IMAX,JMAX)     !! Zonal wind for advection (present)
  real(8),save :: v0_adv(IMAX,JMAX)     !! Zonal wind for advection (present)
!
  real(8),save :: dphidt0(IMAX,JMAX)    !! d(phi)/dt (present)
!
  real(8),save :: f_8byte(IMAX,JMAX)
  real(8),save :: phi_zonal(JMAX)
  
  real(8) :: workj(JMAX)
  real(8) :: work1(IMAX)
  real(8) :: work2(IMAX)
  real(8) :: work3(IMAX)
  real(8) :: work4(IMAX)
!
  real(8) :: phibar
  real(8) :: beta, betauv
  real(8) :: dt, dt_prev
  real(8) :: alon_c, alat_c, rr, r, phis0, hh0, uu0
  real(8) :: aa, bb, cc, ww1, ww2
  real(8) :: thetac, ramdac, ramda, alpha0, cosa, sina
  real(8) :: alat0, alat1, alat2, umax, ss, ww, xx, uu, en, f, hhat
  real(8) :: eps1
  real(8) :: gmean, gmean_h_0, gmean_energy_0, gmean_vorticity_0
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
    kind_phis = 0
    kind_uv   = 2
    kt_end    = 24*12  !! 12 days
  else if ( JCN_INITIAL == 2 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 2
    kind_phis = 1
    kind_uv   = 2
    kt_end    = 24*15  !! 15 days
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
    kt_end    = 24*6   !! 6 days
  else
    write(6,*) "Error: JCN_INITIAL=",JCN_INITIAL," is not supported."
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
    phibar = 1100.0d0*GRAV
   !$OMP PARALLEL default(SHARED), private(i,j,ramda,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ramda=ALON(i)
        r = ER*acos( sin(thetac)*SINLAT(j) + cos(thetac)*COSLAT(j)*cos(ramda-ramdac) )
        if ( r < rr ) then
          phi0(i,j) = 1000.0d0*GRAV/2*( 1.0d0 + cos(PI*r/rr) )
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
  else
! ------------------------------
    write(6,*) "Error: kind_phi=",kind_phi," is not supported."
    stop 999
! ------------------------------
  end if
! ------------------------------
!
   
! ------------------------------
  if (kind_uv.eq.0) then
! ------------------------------
    !$OMP PARALLEL default(SHARED), private(i,j,ww)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
           u0_adv(i,j)=0.0d0
           v0_adv(i,j)=0.0d0
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
        u0_adv(i,j) = uu0*( COSLAT(j)*cos(alpha0) + SINLAT(j)*COSLON(i)*sin(alpha0) )
        v0_adv(i,j) = -uu0*SINLON(i)*sin(alpha0)
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
        u0_adv(i,j)=uu0*COSLAT(j)
        v0_adv(i,j)=0.0d0
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
! ------------------------------
  else
! ------------------------------
    write(6,*) "Error: kind_uv=",kind_uv," is not supported."
    stop 999
! ------------------------------
  end if
! ------------------------------
!
    call grid_to_wave      &
     &( phi0,      &!IN
     &  qphi )  !OUT
!
    call tendency          &  
     &( qphi,              &!IN
     &  phibar,            &!IN
     &  u0_adv, v0_adv,    &!OUT
     &  phi0, dphidt0 )     !OUT
!
 !$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    do i=1,IMAX
      phim(i,j) = phi0(i,j)
    end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL

  if ( JCN_MONIT == 1 ) then
   !$OMP PARALLEL default(SHARED), private(i,j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        f_8byte(i,j) = phi0(i,j)/GRAV
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    
    ntmax = kt_end/INTHR_MONIT + 1
    call monit_norm__ini( "norm", INTHR_MONIT, ntmax, IMAX, JMAX, weight, f_8byte )
    call monit__ini( "data", INTHR_MONIT, ntmax, IMAX, JMAX, YLAT ) !IN
  end if

!   call grads__outxy( "u0", 2, IMAX, JMAX, ylat, u(1:IMAX,1:JMAX) ) !IN
!   call grads__outxy( "v0", 2, IMAX, JMAX, ylat, v(1:IMAX,1:JMAX) ) !IN
!   call grads__outxy( "phi0", 2, IMAX, JMAX, ylat, phi(1:IMAX,1:JMAX) ) !IN
!
!
! =====================================================================
! #####   Start integration   #####
! =====================================================================
!
!
  dt = -999.0d0

  do istep = 1, nstep
!
    call e_time__start(3,"main_cp_past")

    dt_prev = dt

    if ( istep == 1 ) then
      dt = TIMESTEP/2**nstep_start
      istep_monit = 0
      eps1 = 0.0d0
    else if ( istep <= nstep_start ) then
      dt = TIMESTEP/2**(nstep_start-istep+1)
      istep_monit = -1
      eps1 = 0.0d0
    else if ( istep == nstep_start + 1 ) then
      dt = TIMESTEP
      istep_monit = 1
      eps1 = 0.0d0
    else
      dt = TIMESTEP
      istep_monit = istep_monit + 1
      eps1 = eps_save
      
     !$OMP PARALLEL default(SHARED), private(j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        phim(:,j) = eps1*phim(:,j) + (1.0d0-eps1*2.0d0)*phi0(:,j)
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    end if
!
    call e_time__end(3,"main_cp_past")
!
!   ------------------------------------------------------
!

      call tendency          &
       &( qphi,              &!IN
       &  phibar,            &!IN
       &  u0_adv, v0_adv,    &!OUT
       &  phi0, dphidt0   )   !OUT
!
!   -------------------------------------------------------
!
    if ( JCN_MONIT == 1 .and. istep_monit >= 0 .and. mod( istep_monit, intv_monit ) == 0 ) then

      call e_time__start(15,"output")
!
      nn = istep_monit/intv_monit + 1

     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
         do i=1,IMAX
            f_8byte(i,j) = phi0(i,j)/GRAV
         end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL
     
      call monit__output( 1, nn, f_8byte )
      call monit__output( 2, nn, u0_adv )
      call monit__output( 3, nn, v0_adv )
      call monit_norm__output( nn, f_8byte )
      
      call cal_gmean &
       &( f_8byte,   &
       &  gmean)
      if ( nn == 1 ) then
         gmean_h_0 = gmean
      end if
      write(6,*) "Global mean of height = ",gmean, "(", (gmean-gmean_h_0)/gmean_h_0, ")"
!
      call e_time__end(15,"output")
!
    end if
!
!   -------------------------------------------------------
    
    call e_time__start(4,"main_integ1")

   !$OMP PARALLEL default(SHARED), private(i,j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        phim(i,j) = phim(i,j) + eps1*phi0(i,j)
        phi(i,j) = phim(i,j) + dphidt0(i,j)*2.0d0*dt
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL    
    
      call grid_to_wave      &
       &( phi,      &!IN
       &  qphi )      !OUT            
!
    call e_time__end(4,"main_integ1")
!
!   ===========================================================
!
    call e_time__start(6,"main_hdiff")

    call hdiffusion( dt, qphi )
  
    call e_time__end(6,"main_hdiff")
!
  end do

  call e_time__end(2,"main")
!
end subroutine main



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine tendency                     &
 &( qphi,                &!IN
 &  phibar,              &!IN
 &  u0_adv, v0_adv,      &!IN
 &  phi0, dphidt0      )  !OUT
!
  complex(8),intent(inout) :: qphi(NNUM,0:MMAX)
  real(8),intent(in) :: phibar
  real(8),intent(in) :: u0_adv(IMAX,JMAX)
  real(8),intent(in) :: v0_adv(IMAX,JMAX)
  real(8),intent(out) :: phi0(IMAX,JMAX)
  real(8),intent(out) :: dphidt0(IMAX,JMAX)
!
  real(8) :: phi0x(IMAX,JMAX)
  real(8) :: phi0y(IMAX,JMAX)
!
  integer :: i,j

  call e_time__start(11,"tendency")
!
! ===========================================================
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

 !$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    do i=1,IMAX
      dphidt0(i,j) = ( -u0_adv(i,j)*phi0x(i,j) -v0_adv(i,j)*phi0y(i,j) )/ER
    end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL

  call e_time__end(11,"tendency")
!
end subroutine tendency


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine grid_to_wave &
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
   &( phi       )      !INOUT
  call fft_y__g2w &
   &( phi,        &!IN
   &  qphi    )    !INOUT

  call e_time__end(12,"grid_to_wave")
!
end subroutine grid_to_wave


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine hdiffusion( dt, qphi )
  real(8),intent(in) :: dt
  complex(8),intent(inout) :: qphi(NNUM,0:MMAX)

  real(8) :: coef, tau, factor

  call e_time__start(14,"hdiffusion")

      if ( JCN_HDIFF == 1 ) then
!         if ( kind_phi == 8 ) then
            !! 2nd order diffusion for Galewsky (2004) test case
            coef = 1.0d5*dt
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


end program adv_dfs
