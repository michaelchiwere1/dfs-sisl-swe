program sw_eul_sh

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
  use e_time, only : e_time__start, e_time__end, e_time__output
!
  implicit none
  
!  integer,parameter :: JCN_INITIAL = 1  !! Williamson test case 1 (Advection)
  integer,parameter :: JCN_INITIAL = 2  !! Williamson test case 2
!  integer,parameter :: JCN_INITIAL = 5  !! Williamson test case 5
!  integer,parameter :: JCN_INITIAL = 6  !! Williamson test case 6
!  integer,parameter :: JCN_INITIAL = 10 !! Galewsky 2004 test case
  
   integer,parameter :: JCN_HDIFF = 0  !! No diffusion
!  integer,parameter :: JCN_HDIFF = 1  !! 2nd order diff. for Galewsky test case
!  integer,parameter :: JCN_HDIFF = 2  !! 4th order hyper diffusion

!  integer,parameter :: JCN_MONIT = 0  !! No monitor output
  integer,parameter :: JCN_MONIT = 1  !! Write monitor output
  
  integer,parameter :: INTHR_MONIT = 24  !! Interval of monitor output (hour)

!  integer,parameter :: IMAX = 16
!  integer,parameter :: JMAX = 8
!  integer,parameter :: NMAX = 5
!  real(8),parameter :: TIMESTEP = 3600.0d0   !! 1 hour
!
  integer,parameter :: IMAX = 128            !! About 300km resolution
  integer,parameter :: JMAX = 64             !! J is from North to South
  integer,parameter :: NMAX = 42            !! Quadric grid
  real(8),parameter :: TIMESTEP = 1800.0d0   !! 1 hour
!
!  integer,parameter :: IMAX = 160            !! 250km resolution
!  integer,parameter :: JMAX = 80             !! J is from North to South
!  integer,parameter :: NMAX = 53             !! Quadric grid
!  real(8),parameter :: TIMESTEP = 1800.0d0   !! 1 hour
!
!  integer,parameter :: IMAX = 320            !! 120km resolution
!  integer,parameter :: JMAX = 160            !! J is from North to South
!  integer,parameter :: NMAX = 106            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 900.0d0   !! 30 min.
!
!  integer,parameter :: IMAX = 640            !! 60km resolution
!  integer,parameter :: JMAX = 320            !! J is from North to South
!  integer,parameter :: NMAX = 213            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 450.0d0   !! 15 min.
!
!  integer,parameter :: IMAX = 1920           !! 20km resolution
!  integer,parameter :: JMAX = 960            !! J is from North to South
!  integer,parameter :: NMAX = 639            !! Quadric grid
!  real(8),parameter :: TIMESTEP = 150.0d0    !! 2.5 min.
!  real(8),parameter :: TIMESTEP = 240.0d0    !! 4 min.
!  real(8),parameter :: TIMESTEP = 3600.0d0/14.0d0    !! 4.285714 min.
!  real(8),parameter :: TIMESTEP = 3600.0d0/13.0d0    !! 4.615385 min.
!  real(8),parameter :: TIMESTEP = 300.0d0    !! 5 min.
!
!  integer,parameter :: IMAX = 3840           !! 10km resolution
!  integer,parameter :: JMAX = 1920           !! J is from North to South
!  integer,parameter :: NMAX = 1279           !! Quadric grid
!  integer,parameter :: NMAX = 959            !! Cubic grid
!  real(8),parameter :: TIMESTEP = 80.0d0    !! 80 sec.
!
!  integer,parameter :: IMAX = 7680           !! 5km resolution
!  integer,parameter :: JMAX = 3840           !! J is from North to South
!  integer,parameter :: NMAX = 2559           !! Quadric grid
!  integer,parameter :: NMAX = 1919           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 40.0d0    !! 40 sec.
!
!  integer,parameter :: IMAX = 15360          !! 2.6km resolution
!  integer,parameter :: JMAX = 7680           !! J is from North to South
!  integer,parameter :: NMAX = 7678           !! Linear grid
!  integer,parameter :: NMAX = 5119           !! Quadric grid
!  integer,parameter :: NMAX = 3839           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 20.0d0     !! 20 sec.
!
!  integer,parameter :: IMAX = 20480          !! 2.0km resolution
!  integer,parameter :: JMAX = 10240          !! J is from North to South
!  integer,parameter :: NMAX = 6826           !! Quadric grid
!  integer,parameter :: NMAX = 5119           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 15.0d0     !! 15 sec.
!
!  integer,parameter :: IMAX = 30720          !! 1.3km resolution
!  integer,parameter :: JMAX = 15360          !! J is from North to South
!  integer,parameter :: NMAX = 10239          !! Quadric grid
!  integer,parameter :: NMAX = 7679           !! Cubic grid
!  real(8),parameter :: TIMESTEP = 10.0d0     !! 10 sec.

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
  real(8),save :: GW(JMAX)
!
  real(8),save,allocatable :: ERFNN1(:)
  integer,save,allocatable :: NTOTAL(:)
  
  real(8),save :: eps_save = 0.05d0  !## 0.0 <= eps_save <= 0.1?
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
!
  call legendre__ini                          &
   &( IMAX, JMAX, NNUM, MNUM,                 &!IN
   &  MNWAV, MNWAV_UV, MNSTART, MNSTART_UV,   &!OUT
   &  gw, alat, sinlat, coslat, coslat_inv )   !OUT
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
!
!  write(6,*) 'alat=',alat
  write(6,*) 'lat='
  write(6,'(5F10.3)') -alat*180.0d0/PI
!xx  write(6,*) 'alat_sl=',alat_sl
!xx  write(6,*) 'sinlat=',sinlat
!xx  write(6,*) 'acos(sinlat)=',acos(sinlat)
!xx  write(6,*) 'coslat=',coslat
!xx  stop 234
!
!xx  NMAX=MNUM-1
  
  allocate( NTOTAL(MNWAV) )
  allocate( ERFNN1(MNWAV) )

!$OMP PARALLEL default(SHARED), private(m,mn,an)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MMAX
    do mn=MNSTART(m),MNSTART(m+1)-1
      NTOTAL(mn) = m + mn - MNSTART(m)
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



subroutine main
!
  real(8) :: qphi(2,mnwav)               !! GRAV*height (geopotential height)
!
  real(8),save :: phim(IMAX,JMAX)        !! GRAV*height 
  real(8),save :: phi0(IMAX,JMAX)        !! GRAV*height
  real(8),save :: phi(IMAX,JMAX)         !! GRAV*height
!
  real(8),save :: phis (IMAX,JMAX)       !! Surface GRAV*height 
  real(8),save :: phisx(IMAX,JMAX)
  real(8),save :: phisy(IMAX,JMAX)
!
  real(8),save :: u0_adv(IMAX,JMAX)      !! Zonal wind for advection (present)
  real(8),save :: v0_adv(IMAX,JMAX)      !! Zonal wind for advection (present)
!
  real(8),save :: dphidt0(IMAX,JMAX)     !! d(phi)/dt (present)
  real(8),save :: dudt0(IMAX,JMAX)       !! d(u)/dt (present)
  real(8),save :: dvdt0(IMAX,JMAX)       !! d(v)/dt (present)
!
  real(8),save :: f_8byte(IMAX,JMAX)
!
  real(8),save :: phi_zonal(JMAX)
!
  real(8) :: phibar
  real(8) :: beta, betauv
  real(8) :: dt, dt_prev
  real(8) :: alon_c, alat_c, rr, r, phis0, hh0, uu0
  real(8) :: aa,bb,cc,uu
  real(8) :: thetac, ramdac, ramda, alpha0, cosa, sina
  real(8) :: alat0, alat1, alat2, umax, ss, ww, xx, en, f, hhat
  real(8) :: ww1, ww2
  real(8) :: eps1
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

  call e_time__start(2,"main")

  if ( JCN_INITIAL == 1 ) then
    iadv_only = 1   !! Advection
    kind_phi  = 1
    kind_phis = 0
    kind_uv   = 2
    kt_end    = 24*12  !! 12 days
  else if ( jcn_initial == 2 ) then
    iadv_only = 0   !! Shallow water
    kind_phi  = 2
    kind_phis = 1
    kind_uv   = 2
    kt_end    = 24*15  !! 15 days
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
     ss = sum( phi_zonal(1:JMAX)*GW(1:JMAX) )/sum( GW(1:JMAX) )
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

    call grid_to_wave_advonly &
     &( phi0,  &!IN
     &  qphi )  !OUT

    call tendency_advonly    &
     &( qphi,                &!IN
     &  phibar,              &!IN
     &  u0_adv, v0_adv,      &!IN
     &  phi0, dphidt0       ) !OUT

   !$OMP PARALLEL default(SHARED), private(i,j)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        phim(i,j) = phi0(i,j)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

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

 	 
 	 
! 	 write(6,*) "phi0=",phi0
 	 
 	 
!
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
    call monit_norm__ini( "norm", INTHR_MONIT, ntmax, IMAX, JMAX, gw, f_8byte )
    call monit__ini( INTHR_MONIT, ntmax, IMAX, JMAX, YLAT ) !IN
  end if
!
! ===================================================================

  dt = -999.0d0

  do istep = 1, nstep

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
      
     !$OMP PARALLEL default(SHARED), private(i,j)
     !$OMP DO schedule(STATIC)
      do j=1,JMAX
        do i=1,IMAX
          phim(i,j) = eps1*phim(i,j) + (1.0d0-eps1*2.0d0)*phi0(i,j)
        end do
      end do
     !$OMP END DO
     !$OMP END PARALLEL
    end if
!
    call e_time__end(3,"main_cp_past")
    

    call tendency            &
     &( qphi,                &!IN
     &  phibar,              &!IN
     &  u0_adv, v0_adv,      &!IN
     &  phi0, dphidt0      )  !OUT
!    write(6,*) 'n       = ',n
!    write(6,*) 'istep_monit = ',istep_monit
!
!   -------------------------------------------------------
!
    if ( JCN_MONIT == 1 .and. istep_monit >= 0 .and. mod( istep_monit, intv_monit ) == 0 ) then
!
      call e_time__start(15,"output")
!
      nn = istep_monit/intv_monit + 1
!
        !$OMP PARALLEL default(SHARED), private(i,j)
        !$OMP DO schedule(STATIC)
         do j=1,JMAX
            do i=1,IMAX
               f_8byte(i,j) = phi0(i,j)/GRAV
            end do
         end do
        !$OMP END DO
        !$OMP END PARALLEL
        	
        	write(6,*) "f_8byte=",f_8byte(1,:)
!
      call monit__output( 1, nn, f_8byte )
      call monit__output( 2, nn, u0_adv )
      call monit__output( 3, nn, v0_adv )
      call monit_norm__output( nn, f_8byte )
!
      call e_time__end(15,"output")
!
    end if


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
    
    call grid_to_wave &
     &( phi,       &!IN
     &  qphi     )  !OUT                      
!
    call e_time__end(4,"main_integ1")

    call e_time__start(7,"main_hdiff")

    if ( JCN_HDIFF == 1) then
      call hdiff2            &
       &( dt,                &
       &  qphi  )
    else if ( JCN_HDIFF == 2 ) then
      call hdiff4            &
       &( dt,                &
       &  qphi  )
    end if

    call e_time__end(7,"main_hdiff")
!
  end do
  
  call e_time__end(2,"main")
!
end subroutine main


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine tendency_advonly &
 &( qphi,                &!IN
 &  phibar,              &!IN
 &  u0_adv, v0_adv,      &!IN
 &  phi0, dphidt0      )  !OUT
!
  real(8),intent(in) :: qphi(2,mnwav)
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
!
! ===========================================================

  call e_time__start(11,"tendency")
!
  call legendre__w2g &
   &( qphi,          &!IN
   &  phi0        )   !OUT
!
  call xderiv_run        &
   &( IMAX, JMAX, MNUM,  &!IN
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
 !$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    do i=1,IMAX
      dphidt0(i,j) = ( -u0_adv(i,j)*phi0x(i,j) -v0_adv(i,j)*phi0y(i,j) )*ACOSLAT_INV(j)
    end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL
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
end subroutine tendency_advonly


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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


    call xderiv_run            &
     &( IMAX, JMAX, MNUM, &!IN
     &  u0           &!IN
     &  ux       )    !OUT
!
    call xderiv_run            &
     &( IMAX, JMAX, MNUM, &!IN
     &  v0,         &!IN
     &  vx      )    !OUT
!
    call fft_run_togrid     &
     &( imax, jmax,         &!IN
     &  ux               )   !INOUT
!
    call fft_run_togrid     &
     &( imax, jmax,         &!IN
     &  vx               )   !INOUT
!
  call legendre__w2g &
   &( qdiv,          &!IN
   &  div        )    !OUT
!
  call fft_run_togrid &
   &( IMAX, JMAX,     &!IN
   &  div          )   !INOUT
!
    call legendre__w2g   &
     &( qrot,      &!IN
     &  rot      )  !OUT
!
    call fft_run_togrid &
     &( imax, jmax,     &!IN
     &  rot          )   !INOUT
!

   !$OMP PARALLEL default(SHARED), private(i,j)
    !$OMP DO schedule(STATIC)
     do j=1,JMAX
        do i=1,IMAX
        uy(i,j) = vx(i,j) - earth*rot(i,j)*coslat(j)**2
        vy(i,j) = earth*div(i,j)*coslat(j)**2 - ux(i,j)
     end do
    !$OMP END DO
   !$OMP END PARALLEL
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


 !$OMP PARALLEL default(SHARED), private(i,j)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    do i=1,IMAX
      dphidt0(i,j) = dphidt0(i,j) + ( -u0(i,j)*phi0x(i,j) -v0(i,j)*phi0y(i,j) )*ACOSLAT_INV(j)
      dudt(i,j) = dudt(i,j)                               &
       &          - ( u0(i,j)*ux(i,j) + v0(i,j)*uy(i,j) ) &
       &            *ACOSLAT_INV(j)                       &
       &          + 2.0d0*omg*sinlat(j)*v0(i,j)
      dvdt(i,j) = dvdt(i,j)                                     &
       &          - ( u0(i,j)*vx(i,j) + v0(i,j)*vy(i,j)         &
       &              + sinlat(j)*( u0(i,j)**2 + v0(i,j)**2 ) ) &
       &            *ACOSLAT_INV(j)                             &
       &          - 2.0d0*omg*sinlat(j)*u0(i,j)
    end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL

  call e_time__end(11,"tendency")
!
end subroutine tendency



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


subroutine grid_to_wave_advonly &
 &( phi,         &
 &  qphi   )
!
  real(8),intent(inout) :: phi(IMAX,JMAX)
  real(8),intent(out) :: qphi(2,mnwav)
  
  integer :: i,j
!
! ==========================================================

  call e_time__start(12,"grid_to_wave")
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
end subroutine grid_to_wave_advonly


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



subroutine hdiff2  &
 &( dt,            &
 &  qphi  )
!
  real(8),intent(in) :: dt
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
 &  qphi  )
!
  real(8),intent(in) :: dt
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
      qphi(k,l) = qphi(k,l)/( 1.0d0 + coef*ERFNN1(l)**2 )
    end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL
!
end subroutine hdiff4



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine cal_mean &
 &( data,           &
 &  amean)
!
  real(8),intent(in) :: data(IMAX,JMAX)
  real(8),intent(out) :: amean
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
      denominator = denominator + GW(j)
    end do
    denominator = denominator*IMAX
  end if
  
!$OMP PARALLEL default(SHARED), private(i,j,ww)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    ww=0.0d0
    do i=1,IMAX
      ww = ww + data(i,j)*GW(j)
    end do
    work(j) = ww
  end do
 !$OMP END DO
!$OMP END PARALLEL
  
  amean = 0.0d0
  do j=1,JMAX
    amean = amean + work(j)*GW(j)
  end do
!
  amean = amean/denominator
!
end subroutine cal_mean


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



end program sw_eul_sh
