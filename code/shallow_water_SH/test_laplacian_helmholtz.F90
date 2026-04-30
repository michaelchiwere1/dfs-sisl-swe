program test_laplacian_helmholtz

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

!  integer,parameter :: IMAX = 16
!  integer,parameter :: JMAX = 8 + jcn_grid
!  integer,parameter :: NMAX = 6
!  real(8),parameter :: TIMESTEP = 3600.0d0   !! 1 hour
!
  integer,parameter :: IMAX = 128            !! About 300km resolution
  integer,parameter :: JMAX = 64             !! J is from North to South
!  integer,parameter :: NMAX = 62             !! Linear grid
  integer,parameter :: NMAX = 42            !! Quadric grid
!
!  integer,parameter :: IMAX = 160            !! 250km resolution
!  integer,parameter :: JMAX = 80             !! J is from North to South
!  integer,parameter :: NMAX = 78             !! Linear grid
!  integer,parameter :: NMAX = 53             !! Quadric grid
!
!  integer,parameter :: IMAX = 320            !! 120km resolution
!  integer,parameter :: JMAX = 160            !! J is from North to South
!  integer,parameter :: NMAX = 158            !! Linear grid
!  integer,parameter :: NMAX = 106            !! Quadric grid
!
!  integer,parameter :: IMAX = 640            !! 60km resolution
!  integer,parameter :: JMAX = 320            !! J is from North to South
!  integer,parameter :: NMAX = 318            !! Linear grid
!  integer,parameter :: NMAX = 213            !! Quadric grid
!
!  integer,parameter :: IMAX = 1920           !! 20km resolution
!  integer,parameter :: JMAX = 960            !! J is from North to South
!  integer,parameter :: NMAX = 958            !! Linear grid
!  integer,parameter :: NMAX = 639            !! Quadric grid
!
!  integer,parameter :: IMAX = 3840           !! 10km resolution
!  integer,parameter :: JMAX = 1920           !! J is from North to South
!  integer,parameter :: NMAX = 1918           !! Linear grid
!  integer,parameter :: NMAX = 1279           !! Quadric grid
!  integer,parameter :: NMAX = 959            !! Cubic grid
!
!  integer,parameter :: IMAX = 7680           !! 5km resolution
!  integer,parameter :: JMAX = 3840           !! J is from North to South
!  integer,parameter :: NMAX = 3838           !! Linear grid
!  integer,parameter :: NMAX = 2559           !! Quadric grid
!  integer,parameter :: NMAX = 1919           !! Cubic grid
!
!  integer,parameter :: IMAX = 15360          !! 2.6km resolution
!  integer,parameter :: JMAX = 7680           !! J is from North to South
!  integer,parameter :: NMAX = 7678           !! Linear grid
!  integer,parameter :: NMAX = 5119           !! Quadric grid
!  integer,parameter :: NMAX = 3839           !! Cubic grid
  
!  integer,parameter :: IMAX = 20480          !! 2.0km resolution
!  integer,parameter :: JMAX = 10240          !! J is from North to South
!  integer,parameter :: NMAX = 10238          !! Linear grid
!  integer,parameter :: NMAX = 6826           !! Quadric grid
!  integer,parameter :: NMAX = 5119           !! Cubic grid
!
!  integer,parameter :: IMAX = 30720          !! 1.3km resolution
!  integer,parameter :: JMAX = 15360          !! J is from North to South
!  integer,parameter :: NMAX = 15358          !! Linear grid   
!  integer,parameter :: NMAX = 10239          !! Quadric grid
!  integer,parameter :: NMAX = 7679           !! Cubic grid

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
  real(8),save :: GW(JMAX)
!
  real(8),save,allocatable :: ERFNN1(:)
  integer,save,allocatable :: NTOTAL(:)
!
! -----------------------------------------------------------------
!
  call e_time__start(1,"shallow water")

  call initialize
!
  call main
!
  call e_time__end(1,"shallow water")
!  call e_time__output
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
  write(6,*) "IMAX=",IMAX
  write(6,*) "JMAX=",JMAX
  write(6,*) "NMAX=",NMAX
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
!  write(6,*) 'lat='
!  write(6,'(5F10.3)') -alat*180.0d0/PI
! 
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
  real(8) :: qff(2,mnwav)
  real(8) :: qgg(2,mnwav)
  real(8) :: qhh(2,mnwav)
!
  real(8),save :: ff(IMAX,JMAX)
  real(8),save :: ff2(IMAX,JMAX)
  real(8),save :: lap_ff(IMAX,JMAX)
  real(8),save :: gg(IMAX,JMAX)
  real(8),save :: gg2(IMAX,JMAX)
  real(8),save :: hh(IMAX,JMAX)
  real(8),save :: hh2(IMAX,JMAX)
  real(8),save :: helm_hh(IMAX,JMAX)
  
  real(8),save :: ameanj(JMAX)
  real(8),save :: ameanj2(JMAX)
  
  real(8) :: alonc, alatc, rr, r, eps, ww1, ww2, ww3, ee, ee2
!
  integer :: iterx, kt_end, ntimes_hour, iadv_only
  integer :: kind_phi, kind_phis, kind_uv, nstep, nstep_start, intv_monit
  integer :: istep, i, j, nn, kk, mn, idiv, i1, i2, ntmax
  integer :: idepart_prev
  integer :: istep_monit = 0

  call e_time__start(2,"main")
  
  
    alonc=3.0d0*PI/2                                                          
!    alonc=2.5*PI/2
                                                                                              
!    alatc=0.0d0
    alatc=PI/2-0.05d0
!    alatc=0.5*PI/2

    rr=ER/3
    eps=0.01d0*ER**2
    
   !$OMP PARALLEL default(SHARED), private(i,j,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        r = ER*acos( sin(alatc)*SINLAT(j) + cos(alatc)*COSLAT(j)*cos(ALON(i)-alonc) )
        if ( r < rr ) then
          ww1 = PI*r/rr
          ww2 = PI*ER/rr
          ff(i,j) = 1000.0d0/4*( 1.0d0 + cos(ww1) )**2
          if ( r > 0.0d0 ) then  
            ww3 = -cos(r/ER)/sin(r/ER)*1000.0d0/(2*ER**2)*ww2*(1.0d0+cos(ww1))*sin(ww1)
          else
            ww3 = -1000.0d0*(PI/rr)**2
          end if
          gg(i,j) = ww3 + 1000.0d0/(2*ER**2)*ww2**2*( sin(ww1)**2 - (1.0d0+cos(ww1))*cos(ww1) )
        else
          ff(i,j) = 0.0d0
          gg(i,j) = 0.0d0
        end if
        hh(i,j) = ff(i,j) - eps*gg(i,j)
      end do
      
      do i=1,IMAX
        ff2(i,j) = ff(i,j)
        gg2(i,j) = gg(i,j)
        helm_hh(i,j) = hh(i,j)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL

    call fft_run_towave &
     &( IMAX, JMAX,     &!IN
     &  ff2        )     !OUT
    call legendre__g2w &
     &( ff2,           &!IN
     &  qff    )        !OUT
    call legendre__w2g &
     &( qff,           &!IN
     &  ff2    )        !OUT
    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  ff2      )      !INOUT
     
!$OMP PARALLEL default(SHARED), private(mn,kk)
 !$OMP DO schedule(STATIC)
      do mn = 1,mnwav
        do kk=1,2
          qff(kk,mn) = -ERFNN1(mn)*qff(kk,mn)
        end do
      end do
 !$OMP END DO
!$OMP END PARALLEL
    call legendre__w2g &
     &( qff,           &!IN
     &  lap_ff    )     !OUT
    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  lap_ff    )      !INOUT
     
    call fft_run_towave &
     &( IMAX, JMAX,     &!IN
     &  gg2        )     !OUT
    call legendre__g2w &
     &( gg2,           &!IN
     &  qgg    )        !OUT
    call legendre__w2g &
     &( qgg,           &!IN
     &  gg2    )        !OUT
    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  gg2      )      !INOUT
     
    call fft_run_towave &
     &( IMAX, JMAX,     &!IN
     &  helm_hh    )     !OUT
    call legendre__g2w &
     &( helm_hh,       &!IN
     &  qhh    )        !OUT
!$OMP PARALLEL default(SHARED), private(mn,kk)
 !$OMP DO schedule(STATIC)
    do mn = 1,mnwav
      do kk=1,2
        qhh(kk,mn) = qhh(kk,mn)/( 1.0d0 + eps*ERFNN1(mn) )
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL
    call legendre__w2g &
     &( qhh,           &!IN
     &  helm_hh    )    !OUT
    call fft_run_togrid &
     &( IMAX, JMAX,     &!IN
     &  helm_hh   )      !INOUT                 
     
!    call grads__outxy( "ff", 1, IMAX, JMAX, ylat, ff ) !IN
!    call grads__outxy( "gg", 1, IMAX, JMAX, ylat, gg ) !IN     
!    call grads__outxy( "lap_ff", 1, IMAX, JMAX, ylat, lap_ff ) !IN
!    call grads__outxy( "hh", 1, IMAX, JMAX, ylat, hh ) !IN
!    call grads__outxy( "helm_hh", 1, IMAX, JMAX, ylat, helm_hh ) !IN
        
    ameanj(:) = 0.0d0
    ameanj2(:) = 0.0d0
   !$OMP PARALLEL default(SHARED), private(i,j,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ameanj(j) = ameanj(j) + ( lap_ff(i,j) - gg(i,j) )**2
        ameanj2(j) = ameanj2(j) + gg(i,j)**2
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    ee  = sqrt( sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) ) )
    ee2 = sqrt( sum(GW(:)*ameanj2(:))/( IMAX*sum(GW(:)) ) )
    write(6,*) "RMSE(lap_ff,gg_exact)=",ee,"(",ee/ee2,")"

    ameanj(:) = 0.0d0
    ameanj2(:) = 0.0d0
   !$OMP PARALLEL default(SHARED), private(i,j,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ameanj(j) = ameanj(j) + ( lap_ff(i,j) - gg2(i,j) )**2
        ameanj2(j) = ameanj2(j) + gg2(i,j)**2
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    ee  = sqrt( sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) ) )
    ee2 = sqrt( sum(GW(:)*ameanj2(:))/( IMAX*sum(GW(:)) ) )
    write(6,*) "RMSE(lap_ff,gg_trunc)=",ee,"(",ee/ee2,")"

!    ameanj(:) = 0.0d0
!   !$OMP PARALLEL default(SHARED), private(i,j,r)
!   !$OMP DO schedule(STATIC)
!    do j=1,JMAX
!      do i=1,IMAX
!        ameanj(j) = ameanj(j) + ( gg2(i,j) - gg(i,j) )**2
!      end do
!    end do
!   !$OMP END DO
!   !$OMP END PARALLEL
!    write(6,*) "RMSE(gg2,gg)=",sqrt( sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) ) )
    
    ameanj(:) = 0.0d0
   !$OMP PARALLEL default(SHARED), private(i,j,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ameanj(j) = ameanj(j) + lap_ff(i,j)
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    write(6,*) "mean(lap_ff)=",sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) )

    ameanj(:) = 0.0d0
    ameanj2(:) = 0.0d0
   !$OMP PARALLEL default(SHARED), private(i,j,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ameanj(j) = ameanj(j) + ( helm_hh(i,j) - ff(i,j) )**2
        ameanj2(j) = ameanj2(j) + ff(i,j)**2
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    ee  = sqrt( sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) ) )
    ee2 = sqrt( sum(GW(:)*ameanj2(:))/( IMAX*sum(GW(:)) ) )
    write(6,*) "RMSE(helm_hh,ff_exact)=",ee,"(",ee/ee2,")"

    ameanj(:) = 0.0d0
    ameanj2(:) = 0.0d0
   !$OMP PARALLEL default(SHARED), private(i,j,r)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      do i=1,IMAX
        ameanj(j) = ameanj(j) + ( helm_hh(i,j) - ff2(i,j) )**2
        ameanj2(j) = ameanj2(j) + ff2(i,j)**2
      end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    ee  = sqrt( sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) ) )
    ee2 = sqrt( sum(GW(:)*ameanj2(:))/( IMAX*sum(GW(:)) ) )
    write(6,*) "RMSE(helm_hh,ff_trunc)=",ee,"(",ee/ee2,")"

!    ameanj(:) = 0.0d0
!   !$OMP PARALLEL default(SHARED), private(i,j,r)
!   !$OMP DO schedule(STATIC)
!    do j=1,JMAX
!      do i=1,IMAX
!        ameanj(j) = ameanj(j) + ( ff2(i,j) - ff(i,j) )**2
!      end do
!    end do
!   !$OMP END DO
!   !$OMP END PARALLEL
!    write(6,*) "RMSE(ff2,ff)=",sqrt( sum(GW(:)*ameanj(:))/( IMAX*sum(GW(:)) ) )
     
  call e_time__end(2,"main")
!
end subroutine main


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


end program test_laplacian_helmholtz
