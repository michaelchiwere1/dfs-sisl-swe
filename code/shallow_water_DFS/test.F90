program test
!
  use com_dfs, only : com_dfs__ini
  use prm_phconst, only : ER, PI, PI2, ER, OMG, GRAV
  use fft, only : fft__ini, fft__towave, fft__togrid
  use fft_y, only : fft_y__ini, fft_y__g2w, fft_y__g2w_uv, fft_y__w2g, fft_y__w2g_dy, fft_y__w2g_uv
  use fft_y, only : grid2sin, sin2grid, grid2cos, cos2grid
  use divrot2uv_dfs, only : divrot2uv_dfs__run, divrot2uv_dfs__poisson
  use uv2divrot_dfs, only : uv2divrot_dfs__run, uv2divrot_dfs__laplacian
  use xderiv, only : xderiv_run
  use helmholtz_dfs, only : helmholtz_dfs__run
  use hdiff4_dfs, only : hdiff4_dfs__run
  use gmean_to_zero_dfs, only : gmean_to_zero_dfs__run
  use monit, only : monit__ini, monit__output
  use grads, only : grads__outxy
!
  implicit none

!  integer,parameter :: JCN_DFS = 0 !! Old method and Cheong's basis functions
!  integer,parameter :: JCN_DFS = 1 !! New method to calculate expantion coefficients
  integer,parameter :: JCN_DFS = 2 !! New method to calculate expantion coefficients and New basis functions

!  integer,parameter :: JCN_GRID = 0   !! colat=PI*(j+0.5)/JMAX0, j=0,...,JMAX0-1
!  integer,parameter :: JCN_GRID = 1   !! colat=PI*j/JMAX0, j=0,...,JMAX0   (JMAX=JMAX0+1)
  integer,parameter :: JCN_GRID = -1  !! colat=PI*j/JMAX0, j=1,...,JMAX0-1 (JMAX=JMAX0-1)
  
  integer,parameter :: JCN_ZONALFILTER = 0 !! No zonal filter
!  integer,parameter :: JCN_ZONALFILTER = 1 !! Use zonal Fourier filter
  integer,parameter :: M0 = 20             !! For zonal Fourier filter

!  integer,parameter :: IMAX = 16
!  integer,parameter :: JMAX = 8 + JCN_GRID
!  integer,parameter :: NMAX = 7
!  integer,parameter :: NMAX = 4
!
!  integer,parameter :: IMAX = 128            !! About 300km resolution
!  integer,parameter :: JMAX = 64+JCN_GRID
!  integer,parameter :: NMAX = 63             !! Linear grid
!  integer,parameter :: NMAX = 42            !! Quadric grid
!
!  integer,parameter :: IMAX = 160            !! 250km resolution
!  integer,parameter :: JMAX = 80+JCN_GRID
!  integer,parameter :: NMAX = 79             !! Linear grid
!  integer,parameter :: NMAX = 53             !! Quadric grid
!
!  integer,parameter :: IMAX = 320            !! 120km resolution
!  integer,parameter :: JMAX = 160+JCN_GRID
!  integer,parameter :: NMAX = 159            !! Linear grid
!  integer,parameter :: NMAX = 106            !! Quadric grid
!
!  integer,parameter :: IMAX = 640            !! 60km resolution
!  integer,parameter :: JMAX = 320+JCN_GRID
!  integer,parameter :: NMAX = 319            !! Linear grid
!  integer,parameter :: NMAX = 213            !! Quadric grid
!
  integer,parameter :: IMAX = 1920           !! 20km resolution
  integer,parameter :: JMAX = 960+JCN_GRID   !! J is from North to South
!  integer,parameter :: NMAX = 959            !! Linear grid
  integer,parameter :: NMAX = 639            !! Quadric grid
!
!  integer,parameter :: IMAX = 3840           !! 10km resolution
!  integer,parameter :: JMAX = 1920+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 1919           !! Linear grid
!  integer,parameter :: NMAX = 1279           !! Quadric grid
!  integer,parameter :: NMAX = 959            !! Cubic grid
!
!  integer,parameter :: IMAX = 7680           !! 5km resolution
!  integer,parameter :: JMAX = 3840+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 3839           !! Linear grid
!  integer,parameter :: NMAX = 2559           !! Quadric grid
!  integer,parameter :: NMAX = 1919           !! Cubic grid
!
!  integer,parameter :: IMAX = 15360          !! 2.6km resolution
!  integer,parameter :: JMAX = 7680+JCN_GRID  !! J is from North to South
!  integer,parameter :: NMAX = 7679           !! Linear grid
!  integer,parameter :: NMAX = 5119           !! Quadric grid
!  integer,parameter :: NMAX = 3839           !! Cubic grid
!
!  integer,parameter :: IMAX = 20480          !! 2.0km resolution
!  integer,parameter :: JMAX = 10240+JCN_GRID !! J is from North to South
!  integer,parameter :: NMAX = 10239          !! Linear grid
!  integer,parameter :: NMAX = 6826           !! Quadric grid
!  integer,parameter :: NMAX = 5119           !! Cubic grid
!
!  integer,parameter :: IMAX = 30720          !! 1.3km resolution
!  integer,parameter :: JMAX = 15360+JCN_GRID !! J is from North to South
!  integer,parameter :: NMAX = 15359          !! Linear grid   
!  integer,parameter :: NMAX = 10239          !! Quadric grid
!  integer,parameter :: NMAX = 7679           !! Cubic grid

  integer,parameter :: MMAX = NMAX
  integer,parameter :: NNUM = NMAX + 1
!
  real(8),parameter :: dm = pi2*0.5d0
!  real(8),parameter :: dm = pi2*0.4d0
!
  real(8),parameter :: dlon = pi2/imax
  real(8),parameter :: dlat = pi/(jmax-jcn_grid)
!
  real(8) :: alon(imax)
  real(8) :: sinlon(imax)
  real(8) :: coslon(imax)
  real(8) :: alat(jmax)
  real(8) :: ylat(jmax)
  real(8) :: coslat(jmax)
  real(8) :: coslat_inv(jmax)
  real(8) :: sinlat(jmax)
  
  real(8) :: weight(jmax)
!
! ==================================================================
!
  call initialize
!
  call main
!
!
!******************************************************************
contains
!******************************************************************


subroutine initialize
!
  integer :: i,j
  real(8) :: work(JMAX)
!
! -----------------------------------------------------------------------------
!
  call com_dfs__ini( jcn_dfs, jcn_grid, imax, jmax, mmax, nmax )
  call fft_y__ini( imax, jmax, mmax, alat, weight ) !IN[3],OUT[1]
  
  coslat(:) = cos(alat(:))
  
  call fft__ini( jcn_zonalfilter, M0, imax, jmax, mmax, coslat )
  
  sinlat(:) = sin(alat(:))
  coslat_inv(:) = 1.0d0/coslat(:)
  ylat(:) = alat(:)*180.0d0/PI
!
  do i=1,imax
    alon(i) = dlon*(i-1)    !## 0 <= alon < 2*pi
    sinlon(i) = sin(alon(i))
    coslon(i) = cos(alon(i))
  end do
  
  write(6,*) "ylat(:)=",ylat(:)
!
!
end subroutine initialize




!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine main
!
!
  complex(8) :: qphi(nnum,0:mmax)
  complex(8) :: qphi2(nnum,0:mmax)
  complex(8) :: qu(nnum,0:mmax)
  complex(8) :: qv(nnum,0:mmax)
  complex(8) :: qrot(nnum,0:mmax)
  complex(8) :: qdiv(nnum,0:mmax)
  complex(8) :: qchi(nnum,0:mmax)
!
  real(8) :: phi (imax, jmax)
  real(8) :: phi2 (imax, jmax)
  real(8) :: phiy (imax, jmax)
  real(8) :: phix (imax, jmax)
  real(8) :: rot (imax, jmax)
  real(8) :: div (imax, jmax)
  real(8) :: chi (imax, jmax)
  real(8) :: psi (imax, jmax)
  real(8) :: u (imax, jmax)
  real(8) :: v (imax, jmax)
  real(8) :: work(imax, jmax)
!
  real(4) :: f_4byte(imax,jmax)
!
  real(8),save :: dt_save, alpha_save
!
  real(8) :: alpha, dta, dtb, gamma, eps, eps2
!xx  real(8) :: delta
  real(8) :: dt, cfl
!xxx  real(8) :: du, dv, dphi, div
  real(8) :: velocity, trad, thresh, y, y2
!xxx  real(8) :: dphidt_mean
  real(8) :: alon_c, alat_c, rr, r, phis0, h0, u0, amean
  real(8) :: phi_mean, phi_hosei_mean, rate
  real(8) :: aa,bb,cc, eps_save
  real(8) :: gh0, thetac, ramdac, ramda, alpha0, col
!xxx  real(8) :: uu
!
  integer :: iterx
  integer :: kind_phi, kind_phis, kind_uv, nstep, nstep_start, intv_monit
  integer :: n, i, j, nn, l, h, m, it
!xxx  integer :: kk, mn

  integer :: l_start, l_end, m_start, m_end
!
! =====================================================================
!
!  open( 11, file = 'test.dr',  form = 'UNFORMATTED',     &
!   &        access = 'DIRECT', recl = 4*(imax)*(jmax)  )
   
  call monit__ini( "test", 1, 10, IMAX, JMAX, YLAT ) !IN



  phi(:,:) = 0.0d0
  phi(:,1) = 1.0d0

!  write(6,*) "phi(1,:)=",phi(1,:)
  write(6,*) "mean(phi)=",sum(phi(1,:)*weight(:))

  call fft_y__g2w  &
   &( phi,         &!IN
   &  qphi    )     !INOUT

!  write(6,*) "qphi(:,0)=",qphi(:,0)

  call fft_y__w2g  &
   &( qphi,        &!IN
   &  phi )  !OUT

!  write(6,*) "phi(1,:)=",phi(1,:)
  write(6,*) "mean(phi)=",sum(phi(1,:)*weight(:))
  write(6,*)

 qphi(:,:)= 0.0d0
 qphi(NNUM,0)= 1.0d0

!  write(6,*) "qphi(:,0)=",qphi(:,0)

  call fft_y__w2g  &
   &( qphi,        &!IN
   &  phi )  !OUT

!  write(6,*) "weight(:)=",weight(:)
!  write(6,*) "phi(1,:)=",phi(1,:)
  write(6,*) "mean(phi)=",sum(phi(1,:)*weight(:))
  write(6,*)

!  stop 333

  
 qphi(:,:) = 0.0d0
 
! qphi(1,0)= 6.569924074045746D-006
 qphi(2,0)= 1.149878984631135D-006
! qphi(3,0)= 1.807288516123984D-005
! qphi(4,0)= 1.024822244477884D-006
! qphi(5,0)= 8.184435304486987D-006
 call gmean_to_zero_dfs__run(qphi)

  write(6,*) "qphi(:,0)=",qphi(:,0)

  call fft_y__w2g  &
   &( qphi,        &!IN
   &  phi )  !OUT

!  write(6,*) "weight(:)=",weight(:)
  write(6,*) "phi(1,:)=",phi(1,:)
  write(6,*) "mean(phi)=",sum(phi(1,:)*weight(:))
  write(6,*)

  call fft_y__g2w  &
   &( phi,         &!IN
   &  qphi    )     !INOUT

!  write(6,*) "qphi(:,0)=",qphi(:,0)

  call fft_y__w2g  &
   &( qphi,        &!IN
   &  phi )  !OUT

!  write(6,*) "phi(1,:)=",phi(1,:)
  write(6,*) "mean(phi)=",sum(phi(1,:)*weight(:))
  write(6,*)


!3.282724453977258E-005
!1.934938372716417E-005
!-1.614511230441241E-006
!-6.209535579072673E-006 
!-3.318525782707108E-006
!-6.209535579072672E-006
!-1.614511230441241E-006
!1.934938372716416E-005
!3.282724453977258E-005



  stop 333




  
! ####   Initial values   #####
!
  
  kind_phi = 3

  !
  if ( kind_phi == 2 ) then
    m_start = 3
    m_end = 3
    l_start = 3
    l_end = 3
    call get_spherical_harmonics        &
     &( m_start, m_end, l_start, l_end, &!IN
     &  work, phi )                        !OUT
  else if ( kind_phi == 3 ) then
      m_start = 0
      m_end = 5
      l_start = 0
      l_end = 2
!      m_start = 3
!      m_end = 3
!      l_start = 1
!      l_end = 1
      phi(:,:) = 0.0d0
      do j=1,JMAX
          do m=m_start,m_end
            do l=l_start,l_end
              do i=1,IMAX
                ramda = alon(i)
                col = ( PI/2 - alat(j) )
                if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
!                  phi(i,k,j) = phi(i,k,j) + sin(col)*sin(col)*sin(l*col)*cos(m*ramda)
                  phi(i,j) = phi(i,j) + sin(col)*sin(col)*sin(l*col)*(cos(m*ramda)+sin(m*ramda))
                else if ( m == 0 ) then
!                  phi(i,k,j) = phi(i,k,j) + cos(l*col)*cos(m*ramda)
                  phi(i,j) = phi(i,j) + cos(l*col)*(cos(m*ramda)+sin(m*ramda))
                else if ( mod(m,2) == 1 ) then
!                  phi(i,k,j) = phi(i,k,j) + sin(l*col)*cos(m*ramda)
!                  phi(i,k,j) = phi(i,k,j) + sin(col)*cos(l*col)*cos(m*ramda)
                  phi(i,j) = phi(i,j) + sin(col)*cos(l*col)*(cos(m*ramda)+sin(m*ramda))
                else
!                  phi(i,k,j) = phi(i,k,j) + sin(col)*sin(l*col)*cos(m*ramda)
!                  phi(i,k,j) = phi(i,k,j) + sin(col)*sin(l*col)*(cos(m*ramda)+sin(m*ramda))
!                  phi(i,k,j) = phi(i,k,j) + sin(col)*sin(col)*cos(l*col)*(cos(m*ramda))
                  phi(i,j) = phi(i,j) + sin(col)*sin(col)*cos(l*col)*(cos(m*ramda)+sin(m*ramda))
                end if
              end do
            end do
          end do
      end do

  else if (kind_phi==5) then
    h0 = 5960.0d0
    u0 = 20.0d0
    do j=1,jmax
      do i=1,imax
        phi(i,j) = GRAV*h0 - ( ER*omg*u0 + u0*u0/2 )*sinlat(j)**2
      end do
    end do
  end if


    kind_uv = 1
    
    
    if ( kind_uv == 1 ) then
      !! Williamson's Test Case 1 !!
!      alpha=0.0d0             !! Around the Equator
!      alpha=0.05d0            !! (minor shift)
!      alpha=PI/2-0.05d0       !! (minor shift)
      alpha=PI/2              !! Over the Poles
      velocity=2.0d0*PI*ER/(12*24*60*60)  !! 2*PI*ER/(12days) = about 40m/s
      do j=1,JMAX
          do i=1,IMAX
            ramda = ALON(i)*PI/180.0d0
            u(i,j) = velocity*( COSLAT(j)*cos(alpha) + SINLAT(j)*cos(ramda)*sin(alpha) )
            v(i,j) = -velocity*sin(ramda)*sin(alpha)
          end do
      end do
    else if ( kind_uv == 2 ) then
!      chi(:,:) = 0.0d0
      chi(:,:) = phi(:,:)
!      psi(:,:) = 0.0d0
      psi(:,:) = phi(:,:)


!      write(6,*) "0:chi(1,:)=",chi(1,:)
      
      call fft__towave       &
       &( chi(1:imax,1:jmax)  )    !INOUT

      call fft__towave       &
       &( psi(1:imax,1:jmax)  )    !INOUT

!      write(6,*) "chi(1,:)=",chi(1,:)
!
      call fft_y__g2w  &
       &( chi(1:imax,1:jmax),  &!IN
       &  qdiv    )   !INOUT
!
      call fft_y__g2w  &
       &( psi(1:imax,1:jmax),  &!IN
       &  qrot    )   !INOUT
      !

      m=0
!      write(6,*) "0:qchi(:,m)=",qdiv(:,m)
      
      
      call uv2divrot_dfs__laplacian(qdiv)        !INOUT
      call uv2divrot_dfs__laplacian(qrot)        !INOUT
     

!      write(6,*) "0:qdiv(:,m)=",qdiv(:,m)

      call divrot2uv_dfs__poisson(qdiv,qchi)        !INOUT  

!      write(6,*) "1:qchi(:,m)=",qchi(:,m)
!      qdiv=qchi
!      call uv2divrot_dfs__laplacian(qdiv)        !INOUT
!      write(6,*) "1:qdiv(:,0)=",qdiv(:,0)


      call divrot2uv_dfs__run     &!##  qrot,qdiv -> qu, qv
       &( qdiv, qrot,       &!IN
       &  qu, qv       )     !OUT
       
       
      
!      write(6,*) "0:qv(:,m)=",qv(:,m)
     
      call uv2divrot_dfs__run     &!##  qu, qv -> qrot,qdiv
       &( qu, qv,      &!IN
       &  qdiv, qrot )  !OUT
      
!      write(6,*) "1:qdiv(:,m)=",qdiv(:,m)
     
!      stop 333
     
      call divrot2uv_dfs__run     &!##  qrot,qdiv -> qu, qv
       &( qdiv, qrot,       &!IN
       &  qu, qv       )     !OUT
      
!      write(6,*) "1:qv(:,0)=",qv(:,0)

      call uv2divrot_dfs__run     &!##  qu, qv -> qrot,qdiv
       &( qu, qv,      &!IN
       &  qdiv, qrot )  !OUT
      
!      write(6,*) "2:qdiv(:,0)=",qdiv(:,0)


!      stop 333
      
      call fft_y__w2g_uv  &
       &( qu,                 &!IN
       &  u(1:imax,1:jmax) )  !OUT
      call fft_y__w2g_uv  &
       &( qv,                 &!IN
       &  v(1:imax,1:jmax) )  !OUT
!
      call fft__togrid        &
       &( u(1:imax,1:jmax)   )   !INOUT
      call fft__togrid        &
       &( v(1:imax,1:jmax)   )   !INOUT
!
!     
    end if
!
    phiy(:,:) = 0.0d0
    phi2(:,:) = 0.0d0
!
!      f_4byte(1:imax,1:jmax) = u(1:imax,jmax:1:-1)
!      write(11,rec=1) f_4byte
!
!      f_4byte(1:imax,1:jmax) = v(1:imax,jmax:1:-1)
!      write(11,rec=2) f_4byte
!
!      f_4byte(1:imax,1:jmax) = phi(1:imax,jmax:1:-1)
!      write(11,rec=3) f_4byte
 
      call monit__output( 1, 1, phi )
      call monit__output( 2, 1, u )
      call monit__output( 3, 1, v )
      
      call grads__outxy( "phiy", 1, imax, jmax, ylat, phiy ) !IN  !! zero
      call grads__outxy( "phix", 1, imax, jmax, ylat, phiy ) !IN  !! zero
      call grads__outxy( "phi", 1, imax, jmax, ylat, phi ) !IN
      
 do it=1,2

!  write(6,*) "0:phi(1,1:jmax)=",phi(1,1:jmax)
!  write(6,*) "0:phi(:,2)=",phi(:,2)
      
  call fft__towave       &
   &( phi(1:imax,1:jmax)  )    !OUT

!  write(6,*) "1:phi(1,1:jmax)=",phi(1,1:jmax)
  write(6,*) "1:phi(3,1:jmax)=",phi(3,1:jmax)
!  write(6,*) "1:phi(:,2)=",phi(:,2)
  

  call fft_y__g2w  &
   &( phi(1:imax,1:jmax),  &!IN
   &  qphi    )   !INOUT
      
!      write(6,*) 
!      write(6,*) "1:qphi(:,3)=",qphi(:,3)     
!      call uv2divrot_dfs__laplacian(qphi)        !INOUT
!      qphi2 = qphi
!      write(6,*) "2:qphi2(:,3)=",qphi2(:,3)
!      call poisson2(qphi2,qphi)        !INOUT  
!      write(6,*) "2:qphi(:,3)=",qphi(:,3)
!      write(6,*) 
      
   
!  write(6,*) "qphi(1:jmax,0)=",qphi(1:jmax,0)
!  write(6,*) "qphi(1:jmax,1)=",qphi(1:jmax,1)
!  stop 111
  
!
  call fft__towave       &
   &( u(1:imax,1:jmax)  )    !OUT
  call fft__towave       &
   &( v(1:imax,1:jmax)  )    !OUT
  call fft_y__g2w_uv  &
   &( u(1:imax,1:jmax),  &!IN
   &  qu    )   !INOUT
  call fft_y__g2w_uv  &
   &( v(1:imax,1:jmax),  &!IN
   &  qv    )   !INOUT
!
  call uv2divrot_dfs__run     &!##  qu, qv -> qrot,qdiv
   &( qu, qv,      &!IN
   &  qdiv, qrot )  !OUT
  call divrot2uv_dfs__run     &!##  qrot,qdiv -> qu, qv
   &( qdiv, qrot,       &!IN
   &  qu, qv       )     !OUT
!

   
  write(6,*) "2:qphi(:,0)=",qphi(:,0)
  write(6,*) "2:qphi(:,1)=",qphi(:,1)
  write(6,*) "2:qphi(:,2)=",qphi(:,2)
  write(6,*) "2:qphi(:,3)=",qphi(:,3)
   
  write(6,*) "2:qrot(:,0)=",qrot(:,0)
  write(6,*) "2:qrot(:,1)=",qrot(:,1)
  write(6,*) "2:qrot(:,2)=",qrot(:,2)
  write(6,*) "2:qrot(:,3)=",qrot(:,3)
  write(6,*) "2:qrot(:,4)=",qrot(:,4)

  gamma = 0.01d0*ER*ER
  qphi2(:,:) = qphi(:,:)
  call helmholtz_dfs__run   &!Harmonic spectral filter
   &( 1.0d0, gamma,   &!IN
   &  qphi2     )      !INOUT
  call fft_y__w2g   &
   &( qphi2,                 &!IN
   &  phi2(1:imax,1:jmax) )  !OUT
  call fft__togrid        &
   &( phi2(1:imax,1:jmax) )   !INOUT

  call grads__outxy( "phi2", it, imax, jmax, ylat, phi2 ) !IN

  gamma = 0.001d0*ER*ER*ER*ER
  qphi2(:,:) = qphi(:,:)
  call hdiff4_dfs__run   &!Harmonic spectral filter
   &( 1.0d0, gamma,   &!IN
   &  qphi2     )      !INOUT
  call fft_y__w2g   &
   &( qphi2,                 &!IN
   &  phi2(1:imax,1:jmax) )  !OUT
  call fft__togrid        &
   &( phi2(1:imax,1:jmax) )   !INOUT

  call grads__outxy( "phi4", it, imax, jmax, ylat, phi2 ) !IN
   
  write(6,*) "3:qphi2(:,0)=",qphi2(:,0)
  write(6,*) "3:qphi2(:,1)=",qphi2(:,1)
  write(6,*) "3:qphi2(:,2)=",qphi2(:,2)
  write(6,*) "3:qphi2(:,3)=",qphi2(:,3)
!
  call fft_y__w2g   &
   &( qphi,                 &!IN
   &  phi(1:imax,1:jmax) )  !OUT
   
  write(6,*) "2:phi(:,2)=",phi(:,2)

  call fft_y__w2g_dy   &
   &( qphi,                 &!IN
   &  phiy(1:imax,1:jmax) )  !OUT
   
  call xderiv_run                                          &
   &( phi(1:imax,1:jmax), phiy(1:imax,1:jmax), coslat_inv, &!IN
   &  phix(1:imax,1:jmax)        ) !OUT

  write(6,*) "3:phi(1,1:jmax)=",phi(1,1:jmax)
  write(6,*) "3:phi(3,1:jmax)=",phi(3,1:jmax)
  
   
  call fft__togrid        &
   &( phi(1:imax,1:jmax) )   !INOUT
   
  call fft__togrid        &
   &( phix(1:imax,1:jmax) )   !INOUT
   
  write(6,*) "4:phi(1,1:jmax)=",phi(1,1:jmax)
   
!  write(6,*) "qphi(:,4)=",qphi(:,4)
!  write(6,*) "phi(:,1)=",phi(:,1)
  
!  stop 333
   
   
  call fft__togrid        &
   &( phiy(1:imax,1:jmax) )   !INOUT
   
!  write(6,*) "5:qrot(:,4)=",qrot(:,4)

  call fft_y__w2g   &
   &( qrot,                 &!IN
   &  rot(1:imax,1:jmax) )  !OUT
   
!  write(6,*) "5:rot(7,:)=",rot(9,:)


  write(6,*) "qrot(:,0)=",qrot(:,0)
      write(6,*) "mean(rot)=",sum(rot(1,:)*weight(:))
   
  call fft__togrid        &
   &( rot(1:imax,1:jmax) )   !INOUT
   
  call fft_y__w2g   &
   &( qdiv,                 &!IN
   &  div(1:imax,1:jmax) )  !OUT


      write(6,*) "mean(div)=",sum(div(1,:)*weight(:))

  call fft__togrid        &
   &( div(1:imax,1:jmax) )   !INOUT
   
!
  call fft_y__w2g_uv  &
   &( qu,                 &!IN
   &  u(1:imax,1:jmax) )  !OUT
  call fft_y__w2g_uv  &
   &( qv,                 &!IN
   &  v(1:imax,1:jmax) )  !OUT
  call fft__togrid        &
   &( u(1:imax,1:jmax)   )   !INOUT
  call fft__togrid        &
   &( v(1:imax,1:jmax)   )   !INOUT
!


!
! ---------------------------------------------------
!
!      f_4byte(1:imax,1:jmax) = u(1:imax,jmax:1:-1)
!      write(11,rec=1+it*3) f_4byte
!
!      f_4byte(1:imax,1:jmax) = v(1:imax,jmax:1:-1)
!      write(11,rec=2+it*3) f_4byte
!
!      f_4byte(1:imax,1:jmax) = phi(1:imax,jmax:1:-1)
!      write(11,rec=3+it*3) f_4byte
      
      call monit__output( 1, it+1, phi )
      call monit__output( 2, it+1, u )
      call monit__output( 3, it+1, v )
      call monit__output( 4, it+1, rot )
      call monit__output( 5, it+1, div )
!
 end do
!
!      f_4byte(1:imax,1:jmax) = rot(1:imax,jmax:1:-1)
!      write(11,rec=7) f_4byte
!
!      f_4byte(1:imax,1:jmax) = div(1:imax,jmax:1:-1)
!      write(11,rec=8) f_4byte
!
      call grads__outxy( "phiy", 2, imax, jmax, ylat, phiy ) !IN
!
      call grads__outxy( "phix", 2, imax, jmax, ylat, phix ) !IN
      
      
      call grads__outxy( "phi", 2, imax, jmax, ylat, phi2 ) !IN
      
!
end subroutine main




!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



subroutine get_spherical_harmonics  &
 &( m_start, m_end, l_start, l_end, &!IN
 &  x, d2_x )                        !OUT
!
  integer,intent(in) :: m_start, m_end
  integer,intent(in) :: l_start, l_end
  real(8),intent(out) :: x(imax,jmax)
  real(8),intent(out) :: d2_x(imax,jmax)
!
  real(8) :: pnm1(jmax)
  real(8) :: pnm2(jmax)
  real(8) :: aaa,cc1,cc2,pnm0
  integer :: m, ii, j, n, i, l
!
  d2_x(:,:) = 0.0d0
  x(:,:) = 0.0d0
!
!xx  do m = m_start, m_end
!xx    do j = 1, jmax
!xx      do i = 1, imax
!xx!        aaa = sinlat(j)*coslat(j)**m*cos(m*(alon(i)-dm))
!xx!        d2_x(i,j) = d2_x(i,j) - (m+1)*(m+2)*aaa
!xx!        x(i,j) = x(i,j) + aaa
!xx!
!xx        aaa = coslat(j)**m*cos(m*(alon(i)-dm))
!xx        d2_x(i,j) = d2_x(i,j) - m*(m+1)*aaa
!xx        x(i,j) = x(i,j) + aaa
!xx!
!xx!         d2_x(i,j) = -2.0d0*(1.0d0+3.0d0*cos(2*colat(j)))
!xx!         x(i,j) = 1.0d0/3.0d0+cos(2*colat(j))
!xx!
!xx      end do
!xx    end do
!xx  end do
!
!
  do m = m_start, m_end
!
!    cc1 = sqrt( product( (/(ii*1.0d0,ii=2*m+1, 1,-2)/) )   &
!     &         /product( (/(ii*1.0d0,ii=2*m  , 1,-2)/) ) )
    cc1 = sqrt( product( (/((ii+1)*1.0d0/ii,ii=2*m, 1,-2)/) ) )
!
!    cc2 = sqrt( product( (/(ii*1.0d0,ii=2*m+3, 1,-2)/) )   &
!     &         /product( (/(ii*1.0d0,ii=2*m  , 1,-2)/) ) )
    cc2 = sqrt( product( (/((ii+1)*1.0d0/ii,ii=2*m, 1,-2)/) )*(2*m+3) )
!    cc2 = sqrt( product( (/((ii+3)*1.0d0/ii,ii=2*m, 1,-2)/) ) )
!
    do j = 1, jmax
      pnm1(j) = cc1*coslat(j)**m
      pnm2(j) = cc2*sinlat(j)*coslat(j)**m
    end do
!
    if ( l_start == 0 ) then
      n = m
      do j = 1, jmax
        do i = 1, imax
          aaa = pnm1(j)*cos(m*(alon(i)-dm))
          d2_x(i,j) = d2_x(i,j) - n*(n+1)*aaa
          x(i,j) = x(i,j) + aaa
        end do
      end do
    end if
!
    if ( l_start <= 1 .and. l_end >= 1 ) then
      n = m+1
      do j = 1, jmax
        do i = 1, imax
          aaa = pnm2(j)*cos(m*(alon(i)-dm))
          d2_x(i,j) = d2_x(i,j) - n*(n+1)*aaa
          x(i,j) = x(i,j) + aaa
        end do
      end do
    end if
!
    do l=2,l_start-1
      n = m+l
      do j = 1, jmax
        pnm0    = pnm1(j)
        pnm1(j) = pnm2(j)
        pnm2(j) = - ( sqrt(1.0d0*((n-1)**2-m*m)/(4*(n-1)**2-1))*pnm0   &
         &            - sinlat(j)*pnm1(j)                      )  &
         &          /sqrt(1.0d0*(n*n-m*m)/(4*n*n-1))
      end do   
    end do
!
    do l=max(l_start,2),l_end
      n = m+l
      do j = 1, jmax
        pnm0    = pnm1(j)
        pnm1(j) = pnm2(j)
        pnm2(j) = - ( sqrt(1.0d0*((n-1)**2-m*m)/(4*(n-1)**2-1))*pnm0   &
         &            - sinlat(j)*pnm1(j)                      )  &
         &          /sqrt(1.0d0*(n*n-m*m)/(4*n*n-1))      
        do i = 1, imax
          aaa = pnm2(j)*cos(m*(alon(i)-dm))
          d2_x(i,j) = d2_x(i,j) - n*(n+1)*aaa
          x(i,j) = x(i,j) + aaa
        end do
      end do   
    end do
!
  end do
!
end subroutine get_spherical_harmonics



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


end program test
