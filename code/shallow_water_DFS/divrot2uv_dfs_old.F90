module divrot2uv_dfs_old

  !! Calculate U,V from div,rot

  use prm_phconst,only : ER_INV,ER
  use yderiv_dfs_old, only : yderiv_dfs_old__run1
  use com_dfs, only : MMAX,NNUM,NNUMHF1,NNUMHF2
  use com_dfs, only : jcn_dfs,            &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS,  N_TRUNC_M3_DFS,            &
   &                  N_L2TRUNC1_M0_DFS, N_L2TRUNC2_M0_DFS, N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS, &
   &                  N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS, N_L2TRUNC1_M3_DFS, N_L2TRUNC2_M3_DFS, &
   &                  NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS
  use pentadiagonal, only : pentadiagonal__allocate, pentadiagonal__ini, pentadiagonal__solve
  use tridiagonal, only : tridiagonal__allocate, tridiagonal__ini, tridiagonal__solve
  use gmean_to_zero_dfs, only : gmean_to_zero_dfs__run
  use e_time, only : e_time__start, e_time__end

  implicit none

  private
  public :: divrot2uv_dfs_old__run, divrot2uv_dfs_old__poisson
  
  type :: type_diagonal
    private
    real(8),allocatable :: d(:,:)
    integer,allocatable :: iww(:)
  end type
  type(type_diagonal),save,allocatable :: td1(:)
  type(type_diagonal),save,allocatable :: td2(:)

contains

  subroutine divrot2uv_dfs_old__run(qdiv,qrot,qu,qv) !INOUT,INOUT,OUT,OUT
    !
    complex(8),intent(INOUT) :: qdiv(NNUM,0:MMAX)  !! Divergence
    complex(8),intent(INOUT) :: qrot(NNUM,0:MMAX)  !! Rotation
    complex(8),intent(OUT) :: qu(NNUM,0:MMAX)   !! U
    complex(8),intent(OUT) :: qv(NNUM,0:MMAX)   !! V
    !
    complex(8) :: qpsi(NNUM,0:MMAX) !! (Stream function)/er**2
    complex(8) :: qchi(NNUM,0:MMAX) !! (Velocity potential)/er**2
    !
!    complex(8) :: qtmp(NUMMA_YZ_S,NUMJ_YZ_S)
!    real(8) :: gdata(NUMI_I,NUMFA_I,NUMJ_S)
    !
    integer :: m,j
    real(8) :: am,an
    integer,save :: nt=0
    !
    m=0
!    write(6,*) '3:qrot(:,m,1)=',qrot(:,m,1)

!#if defined(MODWV34)
!    call modwave_cos2lat__run(qrot) !IN,IN,OUT,OUT
!    call modwave_cos2lat__run(qdiv) !IN,IN,OUT,OUT
!#endif

    call divrot2uv_dfs_old__poisson( qdiv, qchi ) !INOUT,OUT
    call divrot2uv_dfs_old__poisson( qrot, qpsi ) !INOUT,OUT


!    m=0
!    write(6,*) '11:m,qchi(:,m)*ER*ER=',qchi(:,m)*ER*ER
!    write(6,*) '11:m,qpsi(:,m)*ER*ER=',qpsi(:,m)*ER*ER
!    write(0,*)

      !
      !! U = ( cos(lat)*@psi/@colat + @chi/@ramda )/er
      !! V = ( - cos(lat)*@chi/@colat + @psi/@ramda )/er
      !
      call yderiv_dfs_old__run1( 1, qpsi, qu ) !IN,IN,OUT
      call yderiv_dfs_old__run1( 1, qchi, qv ) !IN,IN,OUT
      !
!$OMP PARALLEL default(SHARED), private(j,m)
!$OMP DO schedule(STATIC)
      do m=0,MMAX
        do j=1,NNUM
          qu(j,m) = ( -qu(j,m) + qchi(j,m)*m*dcmplx(0.0d0,1.0d0) )*ER_INV
          qv(j,m) = (  qv(j,m) + qpsi(j,m)*m*dcmplx(0.0d0,1.0d0) )*ER_INV
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL


!    write(0,*) '1:mg,qu(mg,:,1)*ER=',mg_debug,qu(mg_debug,:,1)*ER
!    write(0,*) '1:mg,qv(mg,:,1)*ER=',mg_debug,qv(mg_debug,:,1)*ER
    !
  end subroutine divrot2uv_dfs_old__run


  subroutine divrot2uv_dfs_old__poisson( qvar1, qvar2 ) !IN,OUT
    !
    complex(8),intent(INOUT) :: qvar1(NNUM,0:MMAX)
    complex(8),intent(OUT) :: qvar2(NNUM,0:MMAX)
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
      allocate( a1(NNUMHF1+1,0:MMAX) )
      allocate( b1(NNUMHF1+1,0:MMAX) )
      allocate( c1(NNUMHF1+1,0:MMAX) )
      !
      allocate( a2(NNUMHF2+1,0:MMAX) )
      allocate( b2(NNUMHF2+1,0:MMAX) )
      allocate( c2(NNUMHF2+1,0:MMAX) )
      !
        call calc_abc_old      &
         &( a1, b1, c1,    &!OUT
         &  a2, b2, c2  )   !OUT
      !
    end if
    !
      call poisson2_core_old &
       &( a1, b1, c1,    &!IN
       &  a2, b2, c2,    &!IN
       &  qvar1,         &!IN
       &  qvar2      )    !OUT
    !
  end subroutine divrot2uv_dfs_old__poisson


  subroutine poisson2_core_old &
   &( a1, b1, c1,       &!IN
   &  a2, b2, c2,       &!IN
   &  qvar1,            &!IN
   &  qvar2      )       !OUT
    !
    real(8),intent(in) :: a1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: b1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: c1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: a2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: b2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: c2(NNUMHF2+1,0:MMAX)
    complex(8),intent(inout) :: qvar1(NNUM,0:MMAX)
    complex(8),intent(out) :: qvar2(NNUM,0:MMAX)
    !
    complex(8) :: g1(NNUMHF1+1)
    complex(8) :: g2(NNUMHF2+1)
    !
    complex(8) :: x1(NNUMHF1)
    complex(8) :: x2(NNUMHF2)
    !
    integer,parameter :: nnmax=8
    !
    integer :: m,l,l2,nn1,nn2
    complex(8) :: x1aa,x2aa,ww
    !
    !
    call gmean_to_zero_dfs__run(qvar1) !INOUT
    !
!$OMP PARALLEL default(SHARED), private(m,l,l2,g1,g2,x1,x2,nn1,nn2,x1aa,x2aa,ww)
!$OMP DO schedule(STATIC)
      do m=0,mmax
        if ( m == 0 ) then
          qvar1(N_TRUNC_M0_DFS+1:NNUM,m) = 0.0d0
          l2= 1
          l = 3     !! Attention
          g1(l2) = - qvar1(l-2,m)*2.0d0 + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          l = 2                    !! l=2
          g2(l2) = qvar1(l,m) - qvar1(l+2,m)

          l2= 2
          l = 5     !! Attention
          g1(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          l = 4
          g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          do l2=3,NNUMHF1-2
            l=l2*2+1                 !! l = 5,7,...,JMAX-3, Attention
            g1(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          end do
          do l2=3,NNUMHF2-2
            l=l2*2                    !! l = 4,6,...,JMAX-4
            g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          end do
          l2=NNUMHF1-1
          l=l2*2+1                 !! l = JMAX-1, Attention
          g1(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0
          l2=NNUMHF2-1
          l=l2*2                    !! l = JMAX-2
          g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          l2=NNUMHF1
          l=l2*2+1                     !! l = JMAX+1, Attention
          g1(l2) = - qvar1(l-2,m)
          l2=NNUMHF2
          l=l2*2                        !! l = JMAX
          g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0
          l2=NNUMHF2+1
          l=l2*2                        !! l = JMAX
          g2(l2) = - qvar1(l-2,m)
        else
          if ( m == 1 ) then
             qvar1(N_TRUNC_M1_DFS+1:NNUM,m) = 0.0d0     !!!!!!!!!!!!!!!
          else if ( mod(m,2) == 0 ) then
             qvar1(N_TRUNC_M2_DFS+1:NNUM,m) = 0.0d0     !!!!!!!!!!!!!!!
          else
             qvar1(N_TRUNC_M3_DFS+1:NNUM,m) = 0.0d0     !!!!!!!!!!!!!!!
          end if
          
          l2=1
          l =1
          g1(l2) = qvar1(l,m)*3.0d0 - qvar1(l+2,m)
          l =2
          g2(l2) = qvar1(l,m)*2.0d0 - qvar1(l+2,m)

          l2=2
          l =3
          g1(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          l =4
          g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          
          do l2=3,NNUMHF1-1
            l=l2*2-1                 !! l = 3,5,...,JMAX-3
            g1(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          end do
          do l2=3,NNUMHF2-1
            l=l2*2                   !! l = 4,6,...,JMAX-2
            g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0 - qvar1(l+2,m)
          end do
          l2=NNUMHF1
          l=l2*2-1                     !! l = JMAX-1
          g1(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0
          l2=NNUMHF2
          l=l2*2                        !! l = JMAX
          g2(l2) = - qvar1(l-2,m) + qvar1(l,m)*2.0d0
          !
          l2=NNUMHF1+1
          l=l2*2-1                     !! l = JMAX-1
          g1(l2) = - qvar1(l-2,m)
          l2=NNUMHF2+1
          l=l2*2                        !! l = JMAX
          g2(l2) = - qvar1(l-2,m)
        end if

          if ( m == 0 ) then
!xx            nn1 = NNUMHF -1
            nn1 = N_L2TRUNC1_M0_DFS-1  !! NNUMHF -1 !! Attention
            nn2 = N_L2TRUNC2_M0_DFS  !! NNUMHF - 1
          else if ( m == 1 ) then
            nn1 = N_L2TRUNC1_M1_DFS  !! NNUMHF -1
            nn2 = N_L2TRUNC2_M1_DFS  !! NNUMHF - 1
          else if ( mod(m,2) == 0 ) then
            nn1 = N_L2TRUNC1_M2_DFS  !! NNUMHF - 1
            nn2 = N_L2TRUNC2_M2_DFS  !! NNUMHF-2
          else
            nn1 = N_L2TRUNC1_M3_DFS
            nn2 = N_L2TRUNC2_M3_DFS
          end if

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

        if ( m == 0 ) then
          l=1
          qvar2(l,m) = 0.0d0        !! l=1
          do l2=1,NNUMHF1-1
            l = l2*2+1
            qvar2(l,m) = x1(l2)     !! l = 3,5,7,...,JMAX-1
          end do
          do l2=1,NNUMHF2-1
            l = l2*2
            qvar2(l,m) = x2(l2)     !! l = 2,4,6,...,JMAX-2
          end do
          l2=NNUMHF1
          l = l2*2
          qvar2(l,m) = x2(l2)       !! l = JMAX
          l2=NNUMHF2
          l = l2*2
          qvar2(l,m) = x2(l2)       !! l = JMAX
        else
          do l2=1,NNUMHF1
            l = l2*2-1
            qvar2(l,m) = x1(l2)     !! l   = 1,3,5,...,JMAX-1
          end do
          do l2=1,NNUMHF2
            l = l2*2-1
            qvar2(l+1,m) = x2(l2)     !! l+1 = 2,4,6,...,JMAX
          end do
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL

    call gmean_to_zero_dfs__run(qvar2) !INOUT
    !
  end subroutine poisson2_core_old


  subroutine calc_abc_old &!
   &( a1, b1, c1,     &!OUT
   &  a2, b2, c2  )    !OUT
    !
    real(8),intent(out) :: a1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: b1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: c1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: a2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: b2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: c2(NNUMHF2+1,0:MMAX)
    !
    integer :: m,l2,l,nn
    !
!$OMP PARALLEL default(SHARED), private(m,l2,l,nn)
!$OMP DO schedule(STATIC)
    do m=0,mmax
      if ( m == 0 ) then
        do l2=1,NNUMHF1+1           !! m = 0
!          l=l2*2-2                   !! l = 0,2,4,6,...,JMAX-2
!          a1(l2,m) = (l-1)*(l-2)
!          b1(l2,m) = -2*l*l
!          c1(l2,m) = (l+1)*(l+2)

          l=l2*2                   !! l = 2,4,6,...,JMAX-2,JMAX, Attention!!
          a1(l2,m) = (l-1)*(l-2)/ER**2
          b1(l2,m) = -2*l*l/ER**2
          c1(l2,m) = (l+1)*(l+2)/ER**2
        end do
        do l2=1,NNUMHF2+1           !! m = 0
          l=l2*2-1                      !! l = 1,3,5,...,JMAX-1
          a2(l2,m) = (l-1)*(l-2)/ER**2
          b2(l2,m) = -2*l*l/ER**2
          c2(l2,m) = (l+1)*(l+2)/ER**2
        end do
      else if ( mod(m,2) == 1 ) then   !! if wavenumber is odd
          do l2=1,NNUMHF1+1
            l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
            a1(l2,m) = (l-1)*(l-2)/ER**2
            b1(l2,m) = ( -2*l*l -4*m*m )/ER**2
            c1(l2,m) = (l+1)*(l+2)/ER**2
          end do
          do l2=1,NNUMHF2+1
            l=l2*2                    !! l = 2,4,6,...,JMAX
            a2(l2,m) = (l-1)*(l-2)/ER**2
            b2(l2,m) = ( -2*l*l -4*m*m )/ER**2
            c2(l2,m) = (l+1)*(l+2)/ER**2
          end do
      else                         !! if wavenumber is even
          do l2=1,NNUMHF1+1
            l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
            a1(l2,m) = l*(l-1)/ER**2
            b1(l2,m) = ( -2*l*l -4*m*m )/ER**2
            c1(l2,m) = l*(l+1)/ER**2
          end do
          do l2=1,NNUMHF2+1
            l=l2*2                   !! l = 2,4,6,...,JMAX
            a2(l2,m) = l*(l-1)/ER**2
            b2(l2,m) = ( -2*l*l -4*m*m )/ER**2
            c2(l2,m) = l*(l+1)/ER**2
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
    !! --------------------------------------
    !!     LU decomposition
    !! --------------------------------------
    !
!xx!$OMP PARALLEL default(SHARED), private(l2,m)
!xx!$OMP DO schedule(STATIC)
!    do m=ms,mmax
!        b1(1,m) = 1.0d0/b1(1,m)
!      do l2=2,NNUMHF
!          c1(l2-1,m) = c1(l2-1,m)*b1(l2-1,m)
!          b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
!      end do
!    end do
!xx!$OMP END DO
!xx!$OMP END PARALLEL
    !
!$OMP PARALLEL default(SHARED), private(m,l2,nn)
!$OMP DO schedule(STATIC)
  do m=0,mmax
      if ( m == 0 ) then
        nn = N_L2TRUNC1_M0_DFS-1  !! NNUMHF - 1, Attention !!!!!!!!!!!!!!!!
      else if ( m == 1 ) then
        nn = N_L2TRUNC1_M1_DFS
      else if ( mod(m,2) == 0 ) then
        nn = N_L2TRUNC1_M2_DFS
      else
        nn = N_L2TRUNC1_M3_DFS
      end if

!xx      if ( mg == 0 ) then
!        l2=NNUMHF-2
!        b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
!        c1(l2,m) = c1(l2,m)*b1(l2,m)
!        d1(l2,m) = ( d1(l2,m) - a1(l2,m)*d1(l2-1,m) )*b1(l2,m)
!        l2=NNUMHF-1
!        b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
!        c1(l2,m) = ( c1(l2,m)- a1(l2,m)*d1(l2-1,m) )*b1(l2,m)
!        l2=NNUMHF
!        b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
!xx      else
        b1(1,m) = 1.0d0/b1(1,m)
        c1(1,m) = c1(1,m)*b1(1,m)
        do l2=2,nn-2
          b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
          c1(l2,m) = c1(l2,m)*b1(l2,m)
        end do
        l2=nn-1
        if ( l2 >= 2 ) then
          b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
          c1(l2,m) = c1(l2,m)*b1(l2,m)
        end if
        l2=nn
        b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
        c1(l2,m) = 0.0d0
        do l2=nn+1,NNUMHF1+1
          a1(l2,m) = 0.0d0
          b1(l2,m) = 0.0d0
          c1(l2,m) = 0.0d0
        end do
!xx      end if
  end do
!$OMP END DO
!$OMP END PARALLEL
    !
!xx!$OMP PARALLEL default(SHARED), private(l2,m)
!xx!$OMP DO schedule(STATIC)
!    do m=0,mmax                   !! m = 0,1,...,mmax
!        b2(1,m) = 1.0d0/b2(1,m)
!      do l2=2,NNUMHF
!          c2(l2-1,m) = c2(l2-1,m)*b2(l2-1,m)
!          b2(l2,m) = 1.0d0/( b2(l2,m) - a2(l2,m)*c2(l2-1,m) )
!      end do
!    end do
!xx!$OMP END DO
!xx!$OMP END PARALLEL


!$OMP PARALLEL default(SHARED), private(m,l2,nn)
!$OMP DO schedule(STATIC)
  do m=0,mmax                   !! m = 0,1,...,mmax
      if ( m == 0 ) then
        nn = N_L2TRUNC2_M0_DFS  !! NNUMHF - 1
      else if ( m == 1 ) then
        nn = N_L2TRUNC2_M1_DFS
      else if ( mod(m,2) == 0 ) then
        nn = N_L2TRUNC2_M2_DFS
      else
        nn = N_L2TRUNC2_M3_DFS
      end if

      b2(1,m) = 1.0d0/b2(1,m)
      c2(1,m) = c2(1,m)*b2(1,m)
      do l2=2,nn-2
        b2(l2,m) = 1.0d0/( b2(l2,m) - a2(l2,m)*c2(l2-1,m) )
        c2(l2,m) = c2(l2,m)*b2(l2,m)
      end do
      l2=nn-1
      if ( l2 >= 2 ) then
        b2(l2,m) = 1.0d0/( b2(l2,m) - a2(l2,m)*c2(l2-1,m) )
        c2(l2,m) = c2(l2,m)*b2(l2,m)
      end if
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
  end subroutine calc_abc_old

end module divrot2uv_dfs_old
