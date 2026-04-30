module uv2divrot_dfs_old

  !! Calculate div,rot from U,V

  use prm_phconst,only : ER_INV,ER
  use yderiv_dfs_old, only : yderiv_dfs_old__run2
  use com_dfs, only : jcn_dfs, MMAX, NNUM, NNUMHF1, NNUMHF2,                                      &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS, N_TRUNC_M3_DFS,             &
   &                  N_L2TRUNC1_M0_DFS, N_L2TRUNC2_M0_DFS, N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS, &
   &                  N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS, N_L2TRUNC1_M3_DFS, N_L2TRUNC2_M3_DFS, &
   &                  NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS
  use pentadiagonal_4, only : pentadiagonal_4__allocate, pentadiagonal_4__ini, pentadiagonal_4__solve
  use ninediagonal_4, only : ninediagonal_4__allocate, ninediagonal_4__ini, ninediagonal_4__solve
  use pentadiagonal, only : pentadiagonal__allocate, pentadiagonal__ini, pentadiagonal__solve
  use tridiagonal, only : tridiagonal__allocate, tridiagonal__ini, tridiagonal__solve
  use e_time, only : e_time__start, e_time__end

  implicit none

  private
  public :: uv2divrot_dfs_old__run
  public :: uv2divrot_dfs_old__laplacian

contains


  subroutine uv2divrot_dfs_old__run(qu,qv,qdiv,qrot) !IN,IN,OUT,OUT
    !
    complex(8),intent(IN) :: qu(NNUM,0:MMAX)   !! U
    complex(8),intent(IN) :: qv(NNUM,0:MMAX)   !! V
    complex(8),intent(OUT) :: qdiv(NNUM,0:MMAX) !! Divergence
    complex(8),intent(OUT) :: qrot(NNUM,0:MMAX) !! Rotation
    !
    integer :: m,l
      !
      call yderiv_dfs_old__run2(1,qu,qrot) !IN,IN,OUT   !!cos(lat)@U/@lat
      call yderiv_dfs_old__run2(1,qv,qdiv) !IN,IN,OUT   !!cos(lat)@V/@lat
      !
      !! div = ( cos(lat)@V/@lat + @U/@ramda )*/( er*cos2(lat) )
      !! rot = ( - cos(lat)@U/@lat + @V/@ramda )*/( er*cos2(lat) )
      !
!$OMP PARALLEL default(SHARED), private(m,l)
!$OMP DO schedule(STATIC)
        do m=0,MMAX
          do l=1,NNUM
            qdiv(l,m) = (  qdiv(l,m) + qu(l,m)*m*dcmplx(0.0d0,1.0d0) )*ER_INV
            qrot(l,m) = ( -qrot(l,m) + qv(l,m)*m*dcmplx(0.0d0,1.0d0) )*ER_INV
          end do
        end do
!$OMP END DO
!$OMP END PARALLEL
      !
      call divide_cos2lat_old( qdiv ) !INOUT
      call divide_cos2lat_old( qrot ) !INOUT
    !
  end subroutine uv2divrot_dfs_old__run



  subroutine divide_cos2lat_old( qvar ) !IN,OUT
    !
    complex(8),intent(INOUT) :: qvar(NNUM,0:MMAX)
    !
    real(8),save,allocatable :: a1(:,:)
    real(8),save,allocatable :: b1(:,:)
    real(8),save,allocatable :: c1(:,:)
    !
    real(8),save,allocatable :: a2(:,:)
    real(8),save,allocatable :: b2(:,:)
    real(8),save,allocatable :: c2(:,:)
    !
    logical,save :: first_cos2lat = .true.
    !
    !
    if ( first_cos2lat ) then
      !
      first_cos2lat = .false.
      !
!      call initialize
      !
      allocate( a1(NNUMHF1+1,0:mmax) )
      allocate( b1(NNUMHF1+1,0:mmax) )
      allocate( c1(NNUMHF1+1,0:mmax) )
      !
      allocate( a2(NNUMHF2+1,0:mmax) )
      allocate( b2(NNUMHF2+1,0:mmax) )
      allocate( c2(NNUMHF2+1,0:mmax) )
      !
      call calc_abcd_divide_cos2lat_old  &
       &( a1, b1, c1,    &!OUT
       &  a2, b2, c2  )   !OUT
      !
    end if
    !
    call divide_cos2lat_old_core &!biharmonic spectral filter
     &( a1, b1, c1,       &!IN
     &  a2, b2, c2,       &!IN
     &  qvar       )       !INOUT
    !
  end subroutine divide_cos2lat_old


  subroutine divide_cos2lat_old_core &!biharmonic spectral filter
   &( a1, b1, c1,    &!IN
   &  a2, b2, c2,    &!IN
   &  qvar     )      !INOUT
    !
    real(8),intent(in) :: a1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: b1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: c1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: a2(NNUMHF2+1,0:mmax)
    real(8),intent(in) :: b2(NNUMHF2+1,0:mmax)
    real(8),intent(in) :: c2(NNUMHF2+1,0:mmax)
    complex(8),intent(inout) :: qvar(NNUM,0:MMAX)
    !
!xx#ifndef __SX__
    complex(8) :: g1(NNUMHF1)
    complex(8) :: g2(NNUMHF2)
    complex(8) :: x1(NNUMHF1)
    complex(8) :: x2(NNUMHF2)
!xx#else
!    complex(8) :: g1(0:mmax,NNUMHF)
!    complex(8) :: g2(0:mmax,NNUMHF)
!    complex(8) :: x1(0:mmax,NNUMHF)
!    complex(8) :: x2(0:mmax,NNUMHF)
!xx#endif
    !
    integer :: m,l,l2,nn1,nn2,nn
    complex(8) :: ww,x1aa,x2aa
    !
    ! =============================================================
    !
!xx#ifndef __SX__
!$OMP PARALLEL default(SHARED), private(m,l,l2,ww,x1aa,x2aa,g1,g2,x1,x2,nn1,nn2,nn)
!$OMP DO schedule(STATIC)
      do m=0,mmax
        if ( m == 0 ) then
        else
          do l2=1,NNUMHF1
            l=l2*2-1                 !! l = 1,3,5,...,JMAX-3
            g1(l2) = qvar(l,m)
          end do
          do l2=1,NNUMHF2
            l=l2*2                   !! l = 2,4,6,...,JMAX-2
            g2(l2) = qvar(l,m)
          end do

          !! Truncation
          nn1 = NNUM1_M_DFS(m)
          nn2 = NNUM2_M_DFS(m)

          !! forward substitution
           x1(1) = g1(1)*b1(1,m)
           do l2=2,NNUMHF1
             x1(l2) = ( g1(l2) - a1(l2,m)*x1(l2-1) )*b1(l2,m)
           end do

          !! backward substitution
           do l2=NNUMHF1-1,1,-1
             x1(l2) = x1(l2) - c1(l2,m)*x1(l2+1)
           end do

          !! forward substitution
          x2(1) = g2(1)*b2(1,m)
          do l2=2,NNUMHF2
            x2(l2) = ( g2(l2) - a2(l2,m)*x2(l2-1) )*b2(l2,m)
          end do

          !! backward substitution
          do l2=NNUMHF2-1,1,-1
            x2(l2) = x2(l2) - c2(l2,m)*x2(l2+1)
          end do

          do l2=1,NNUMHF1
            l = l2*2-1
            qvar(l,m)   = x1(l2)    !! l   = 1,3,5,...,JMAX-1
          end do
          do l2=1,NNUMHF2
            l = l2*2-1
            qvar(l+1,m) = x2(l2)    !! l+1 = 2,4,6,...,JMAX
          end do
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine divide_cos2lat_old_core


  subroutine uv2divrot_dfs_old__laplacian( qvar ) !IN,OUT
    !
    complex(8),intent(INOUT) :: qvar(NNUM,0:MMAX)
    !
    real(8),save,allocatable :: a1(:,:)
    real(8),save,allocatable :: b1(:,:)
    real(8),save,allocatable :: c1(:,:)
    !
    real(8),save,allocatable :: a2(:,:)
    real(8),save,allocatable :: b2(:,:)
    real(8),save,allocatable :: c2(:,:)
    !
    real(8),save,allocatable :: aa1(:,:)
    real(8),save,allocatable :: bb1(:,:)
    real(8),save,allocatable :: cc1(:,:)
    !
    real(8),save,allocatable :: aa2(:,:)
    real(8),save,allocatable :: bb2(:,:)
    real(8),save,allocatable :: cc2(:,:)
    !
    logical,save :: first_laplacian = .true.
    !
    if ( first_laplacian ) then
      !
      first_laplacian = .false.
      !
        allocate( a1(NNUMHF1+1,0:MMAX) )
        allocate( b1(NNUMHF1+1,0:MMAX) )
        allocate( c1(NNUMHF1+1,0:MMAX) )
        allocate( a2(NNUMHF2+1,0:MMAX) )
        allocate( b2(NNUMHF2+1,0:MMAX) )
        allocate( c2(NNUMHF2+1,0:MMAX) )
        
        allocate( aa1(NNUMHF1+1,0:MMAX) )
        allocate( bb1(NNUMHF1+1,0:MMAX) )
        allocate( cc1(NNUMHF1+1,0:MMAX) )
        allocate( aa2(NNUMHF2+1,0:MMAX) )
        allocate( bb2(NNUMHF2+1,0:MMAX) )
        allocate( cc2(NNUMHF2+1,0:MMAX) )
      
        call calc_abcd_lap_old &
         &( a1, b1, c1,     &!OUT
         &  a2, b2, c2   )   !OUT
        call calc_abcd_divide_cos2lat_old  &
         &( aa1, bb1, cc1,     &!OUT
         &  aa2, bb2, cc2   )   !OUT
      !
    end if
    !
      call laplacian_core_old &
       &( a1, b1, c1,     &!IN
       &  a2, b2, c2,     &!IN
       &  aa1, bb1, cc1,  &!IN
       &  aa2, bb2, cc2,  &!IN
       &  qvar       )     !INOUT
    !
  end subroutine uv2divrot_dfs_old__laplacian
  

  subroutine laplacian_core_old &!biharmonic spectral filter
   &( a1, b1, c1,       &!IN
   &  a2, b2, c2,       &!IN
   &  aa1, bb1, cc1,    &!IN
   &  aa2, bb2, cc2,    &!IN
   &  qvar     )         !INOUT
    !
    real(8),intent(in) :: a1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: b1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: c1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: a2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: b2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: c2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: aa1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: bb1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: cc1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: aa2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: bb2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: cc2(NNUMHF2+1,0:MMAX)
    complex(8),intent(inout) :: qvar(NNUM,0:MMAX)
    !
!xx#ifndef __SX__
    complex(8) :: g1(NNUMHF1)
    complex(8) :: g2(NNUMHF2)
    complex(8) :: x1(NNUMHF1)
    complex(8) :: x2(NNUMHF2)
!xx#else
!    complex(8) :: g1(0:mmax,NNUMHF)
!    complex(8) :: g2(0:mmax,NNUMHF)
!    complex(8) :: x1(0:mmax,NNUMHF)
!    complex(8) :: x2(0:mmax,NNUMHF)
!xx#endif
    !
    integer :: m,l,l2,nn1,nn2,nn
    complex(8) :: ww,x1aa,x2aa
    !
    ! =============================================================
    !
!xx#ifndef __SX__
!$OMP PARALLEL default(SHARED), private(m,l,l2,nn1,nn2,nn,ww,x1aa,x2aa,g1,g2,x1,x2)
!$OMP DO schedule(STATIC)
      do m=0,mmax
        l2=1
        l=l2*2-1                     !! l = 1
        g1(l2) = qvar(l,m)*b1(l2,m) + qvar(l+2,m)*c1(l2,m)
        l=l+1                        !! l = 2
        g2(l2) = qvar(l,m)*b2(l2,m) + qvar(l+2,m)*c2(l2,m)

        do l2=2,NNUMHF1-1
          l=l2*2-1                 !! l = 3,5,...,JMAX-3
          g1(l2) = qvar(l-2,m)*a1(l2,m)                       &
           &       + qvar(l,m)*b1(l2,m) + qvar(l+2,m)*c1(l2,m)
        end do
        do l2=2,NNUMHF2-1
          l=l2*2                   !! l = 4,6,...,JMAX-2
          g2(l2) = qvar(l-2,m)*a2(l2,m)                       &
           &       + qvar(l,m)*b2(l2,m) + qvar(l+2,m)*c2(l2,m)
        end do

        l2 = NNUMHF1
        l=l2*2-1                     !! l = JMAX-1
        g1(l2) = qvar(l-2,m)*a1(l2,m) + qvar(l,m)*b1(l2,m)

        l2 = NNUMHF2
        l=l2*2                       !! l = JMAX
        g2(l2) = qvar(l-2,m)*a2(l2,m) + qvar(l,m)*b2(l2,m)

        nn1 = NNUMHF1
        nn2 = NNUMHF2

          !! forward substitution
          x1(1) = g1(1)*bb1(1,m)
          do l2=2,NNUMHF1
            x1(l2) = ( g1(l2) - aa1(l2,m)*x1(l2-1) )*bb1(l2,m)
          end do

          !! backward substitution
          do l2=NNUMHF1-1,1,-1
            x1(l2) = x1(l2) - cc1(l2,m)*x1(l2+1)
          end do

          !! forward substitution
          x2(1) = g2(1)*bb2(1,m)
          do l2=2,NNUMHF2
            x2(l2) = ( g2(l2) - aa2(l2,m)*x2(l2-1) )*bb2(l2,m)
          end do

          !! backward substitution
          do l2=NNUMHF2-1,1,-1
            x2(l2) = x2(l2) - cc2(l2,m)*x2(l2+1)
          end do

          do l2=1,NNUMHF1
            l = l2*2-1
            qvar(l,m)   = x1(l2)    !! l   = 1,3,5,...,JMAX-1
          end do
          do l2=1,NNUMHF2
            l = l2*2-1
            qvar(l+1,m) = x2(l2)    !! l+1 = 2,4,6,...,JMAX
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine laplacian_core_old


  subroutine calc_abcd_divide_cos2lat_old  &
   &( a1, b1, c1,        &!OUT
   &  a2, b2, c2   )      !OUT
    !
    real(8),intent(out) :: a1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: b1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: c1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: a2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: b2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: c2(NNUMHF2+1,0:mmax)
    !
    integer :: m,l2,nn
    !
!$OMP PARALLEL default(SHARED), private(m,l2,nn)
!$OMP DO schedule(STATIC)
    do m=0,mmax
      do l2=1,NNUMHF1+1
          a1(l2,m) = - 1.0d0*0.25d0
          b1(l2,m) = + 2.0d0*0.25d0
          c1(l2,m) = - 1.0d0*0.25d0
      end do
      do l2=1,NNUMHF2+1
          a2(l2,m) = - 1.0d0*0.25d0
          b2(l2,m) = + 2.0d0*0.25d0
          c2(l2,m) = - 1.0d0*0.25d0
      end do
      if ( m == 0 ) then
        b2(1,m) = + 1.0d0*0.25d0
        !! c2(m,1) = - 1.0d0*0.25d0
        !!b1(m,1) = + 2.0d0*0.25d0
        !!c1(m,1) = - 1.0d0*0.25d0
        a1(2,m) = - 2.0d0*0.25d0
        !!b1(m,2) = + 2.0d0*0.25d0
        !!c1(m,2) = - 1.0d0*0.25d0
      else
        b1(1,m) = + 3.0d0*0.25d0
      end if

    !
    ! --------------------------------------
    !     LU Decomposition
    ! --------------------------------------
      nn = NNUM1_M_DFS(m)
      b1(1,m) = 1.0d0/b1(1,m)
      c1(1,m) = c1(1,m)*b1(1,m)
      do l2=2,nn-1
        b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
        c1(l2,m) = c1(l2,m)*b1(l2,m)
      end do
      l2=nn
      b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
      c1(l2,m) = 0.0d0
      do l2=nn+1,NNUMHF1+1
        a1(l2,m) = 0.0d0
        b1(l2,m) = 0.0d0
        c1(l2,m) = 0.0d0
      end do

      nn = NNUM2_M_DFS(m)
      b2(1,m) = 1.0d0/b2(1,m)
      c2(1,m) = c2(1,m)*b2(1,m)
      do l2=2,nn-1
        b2(l2,m) = 1.0d0/( b2(l2,m) - a2(l2,m)*c2(l2-1,m) )
        c2(l2,m) = c2(l2,m)*b2(l2,m)
      end do
      l2=nn
      b2(l2,m) = 1.0d0/( b2(l2,m) - a2(l2,m)*c2(l2-1,m) )
      c2(l2,m) = 0.0d0
      do l2=nn+1,NNUMHF2+1
        a2(l2,m) = 0.0d0
        b2(l2,m) = 0.0d0
        c2(l2,m) = 0.0d0
      end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine calc_abcd_divide_cos2lat_old


  subroutine calc_abcd_lap_old    &!
   &( a1, b1, c1,     &!OUT
   &  a2, b2, c2   )   !OUT
    !
    real(8),intent(out) :: a1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: b1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: c1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: a2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: b2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: c2(NNUMHF2+1,0:mmax)
    !
!    real(8) :: a11(0:mmax,NNUMHF)
!    real(8) :: b11(0:mmax,NNUMHF)
!    real(8) :: c11(0:mmax,NNUMHF)
!    real(8) :: a22(0:mmax,NNUMHF)
!    real(8) :: b22(0:mmax,NNUMHF)
!    real(8) :: c22(0:mmax,NNUMHF)
    !
    integer :: m,l2,l
    !
!$OMP PARALLEL default(SHARED), private(m,l2,l)
!$OMP DO schedule(STATIC)
    do m=0,mmax
      if ( m == 0 ) then          !! zonal wavenumber = 0
        do l2=1,NNUMHF1+1
          l=l2*2-2                   !! l = 0,2,4,6,...,JMAX-2,JMAX
          a1(l2,m) = (l-1)*(l-2)*0.25d0/(ER**2)
          b1(l2,m) = -2*l*l*0.25d0/(ER**2)
          c1(l2,m) = (l+1)*(l+2)*0.25d0/(ER**2)
        end do
        do l2=1,NNUMHF2+1
          l=l2*2-1                      !! l = 1,3,5,...,JMAX-1
          a2(l2,m) = (l-1)*(l-2)*0.25d0/(ER**2)
          b2(l2,m) = -2*l*l*0.25d0/(ER**2)
          c2(l2,m) = (l+1)*(l+2)*0.25d0/(ER**2)
        end do
        
      else if ( mod(m,2) == 1 ) then   !! if zonal wavenumber is odd
          do l2=1,NNUMHF1+1
            l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
            a1(l2,m) = (l-1)*(l-2)*0.25d0/(ER**2)
            b1(l2,m) = ( -2*l*l -4*m*m)*0.25d0/(ER**2)
            c1(l2,m) = (l+1)*(l+2)*0.25d0/(ER**2)
            l=l+1                    !! l = 2,4,6,...,JMAX
            a2(l2,m) = (l-1)*(l-2)*0.25d0/(ER**2)
            b2(l2,m) = (-2*l*l -4*m*m)*0.25d0/(ER**2)
            c2(l2,m) = (l+1)*(l+2)*0.25d0/(ER**2)
          end do
          do l2=1,NNUMHF2+1
            l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
            a1(l2,m) = (l-1)*(l-2)*0.25d0/(ER**2)
            b1(l2,m) = ( -2*l*l -4*m*m)*0.25d0/(ER**2)
            c1(l2,m) = (l+1)*(l+2)*0.25d0/(ER**2)
            l=l+1                    !! l = 2,4,6,...,JMAX
            a2(l2,m) = (l-1)*(l-2)*0.25d0/(ER**2)
            b2(l2,m) = (-2*l*l -4*m*m)*0.25d0/(ER**2)
            c2(l2,m) = (l+1)*(l+2)*0.25d0/(ER**2)
          end do
      else                   !! if zonal wavenumber is even, not zero
          do l2=1,NNUMHF1+1
            l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
            a1(l2,m) = l*(l-1)*0.25d0/(ER**2)
            b1(l2,m) = (-2*l*l -4*m*m)*0.25d0/(ER**2)
            c1(l2,m) = l*(l+1)*0.25d0/(ER**2)
          end do
          do l2=1,NNUMHF2+1
            l=l2*2                    !! l = 2,4,6,...,JMAX
            a2(l2,m) = l*(l-1)*0.25d0/(ER**2)
            b2(l2,m) = (-2*l*l -4*m*m)*0.25d0/(ER**2)
            c2(l2,m) = l*(l+1)*0.25d0/(ER**2)
          end do
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL
    !
!    m=1
!    write(6,*) 'calc_abc: m=',m
!    write(6,*) 'a1(m,:)=',a1(m,:)
!    write(6,*) 'b1(m,:)=',b1(m,:)
!    write(6,*) 'c1(m,:)=',c1(m,:)
!    write(6,*) 'd1(m,:)=',d1(m,:)
!    write(6,*) 'a2(m,:)=',a2(m,:)
!    write(6,*) 'b2(m,:)=',b2(m,:)
!    write(6,*) 'c2(m,:)=',c2(m,:)
!    write(6,*) 'd2(m,:)=',d2(m,:)
    !
  end subroutine calc_abcd_lap_old

end module uv2divrot_dfs_old
