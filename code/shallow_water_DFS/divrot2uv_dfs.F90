module divrot2uv_dfs

  !! Calculate U,V from div,rot

  use prm_phconst,only : ER_INV,ER
  use com_dfs, only : MMAX,NNUM,NNUMHF1,NNUMHF2
  use com_dfs, only : jcn_dfs,            &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS,  N_TRUNC_M3_DFS,            &
   &                  N_L2TRUNC1_M0_DFS, N_L2TRUNC2_M0_DFS, N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS, &
   &                  N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS, N_L2TRUNC1_M3_DFS, N_L2TRUNC2_M3_DFS, &
   &                  NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS
  use pentadiagonal, only : pentadiagonal__allocate, pentadiagonal__ini, pentadiagonal__solve
  use tridiagonal, only : tridiagonal__allocate, tridiagonal__ini, tridiagonal__solve
  use gmean_to_zero_dfs, only : gmean_to_zero_dfs__run
  use divrot2uv_dfs_old, only : divrot2uv_dfs_old__run, divrot2uv_dfs_old__poisson
  use e_time, only : e_time__start, e_time__end

  implicit none

  private
  public :: divrot2uv_dfs__run, divrot2uv_dfs__poisson
  
  type :: type_diagonal
    private
    real(8),allocatable :: d(:,:)
    integer,allocatable :: iww(:)
  end type
  type(type_diagonal),save,allocatable :: td1(:)
  type(type_diagonal),save,allocatable :: td2(:)

contains

  subroutine divrot2uv_dfs__run(qdiv,qrot,qu,qv) !INOUT,INOUT,OUT,OUT
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
    
    !! Complex: [cosine component] - dcmplx(0.0d0,1.0d0)*[sine component]
    
    call e_time__start(38,"divrot2uv_dfs__run")

    if ( jcn_dfs >= 1 ) then

      call divrot2uv_dfs__poisson( qdiv, qchi ) !INOUT,OUT
      call divrot2uv_dfs__poisson( qrot, qpsi ) !INOUT,OUT

!$OMP PARALLEL default(SHARED), private(m,j,nt,am,an)
!$OMP DO schedule(STATIC)
        do m=0,mmax
          am=m
          if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
            j=1
            an=j-1
            qu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*qchi(j,m) & !! n=0
             &        - ER_INV*0.25d0*an*qpsi(j+1,m)
            qv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)  & !! n=0
             &        + ER_INV*0.25d0*an*qchi(j+1,m)
            j=2
            an=j-1
            qu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)    & !! n=1
             &        + ER_INV*0.25d0*an*( 3.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
            qv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)    & !! n=1
             &        + ER_INV*0.25d0*an*( -3.0d0*qchi(j-1,m) + qchi(j+1,m) )
            j=3
            an=j-1
            qu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -qchi(j-2,m) + qchi(j,m) ) &
             &        + ER_INV*0.25d0*an*( 2.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
            qv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -qpsi(j-2,m) + qpsi(j,m) ) &
             &        + ER_INV*0.25d0*an*( -2.0d0*qchi(j-1,m) + qchi(j+1,m) )
            do j=4,NNUM-1
              an=j-1
              qu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -qchi(j-2,m) + qchi(j,m) ) &
               &        + ER_INV*0.25d0*an*( -qpsi(j-3,m) + 2.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
              qv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -qpsi(j-2,m) + qpsi(j,m) ) &
               &        + ER_INV*0.25d0*an*( qchi(j-3,m) -2.0d0*qchi(j-1,m) + qchi(j+1,m) )
            end do
            j=NNUM
            an=j-1
            qu(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -qchi(j-2,m) + qchi(j,m) ) &
             &        + ER_INV*0.25d0*an*( -qpsi(j-3,m) + 2.0d0*qpsi(j-1,m) )
            qv(j,m) = ER_INV*0.5d0*dcmplx(0.0d0,1.0d0)*am*( -qpsi(j-2,m) + qpsi(j,m) ) &
             &        + ER_INV*0.25d0*an*( qchi(j-3,m) -2.0d0*qchi(j-1,m) )
          else if ( m == 0 ) then
            do j=1,NNUM-1
              an=j
              qu(j,m) = -ER_INV*an*qpsi(j+1,m)
              qv(j,m) =  ER_INV*an*qchi(j+1,m)
            end do
            j=NNUM
            qu(j,m) = 0.0d0
            qv(j,m) = 0.0d0
          else if ( mod(m,2) == 1 ) then
            j=1
            an=j-1
            qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)  !! n=0
            qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)  !! n=0
            j=2
            an=j-1
            qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)           & !! n=1
             &        + ER_INV*0.5d0*an*( 2.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
            qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)           & !! n=1
             &        + ER_INV*0.5d0*an*( -2.0d0*qchi(j-1,m) + qchi(j+1,m) )
            do j=3,NNUM-1
              an=j-1
              qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)  &
               &        + ER_INV*0.5d0*an*( qpsi(j-1,m) - qpsi(j+1,m) )
              qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m) &
               &        + ER_INV*0.5d0*an*( -qchi(j-1,m) + qchi(j+1,m) )
            end do
            j=NNUM
            an=j-1
            qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)      &
             &        + ER_INV*0.5d0*an*qpsi(j-1,m) 
            qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)      &
             &        - ER_INV*0.5d0*an*qchi(j-1,m)
          else  !! mod(m,2) == 0
            j=1
            an=j
            qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)  &
             &        - ER_INV*0.5d0*an*qpsi(j+1,m)
            qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)  &
             &        + ER_INV*0.5d0*an*qchi(j+1,m)
            do j=2,NNUM-1
              an=j
              qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m) &
               &        + ER_INV*0.5d0*an*( qpsi(j-1,m) - qpsi(j+1,m) )
              qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m) &
               &        + ER_INV*0.5d0*an*( -qchi(j-1,m) + qchi(j+1,m) )
            end do
            j=NNUM
            an=j 
            qu(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qchi(j,m)              &
             &        + ER_INV*0.5d0*an*qpsi(j-1,m)
            qv(j,m) = ER_INV*dcmplx(0.0d0,1.0d0)*am*qpsi(j,m)             &
             &        - ER_INV*0.5d0*an*qchi(j-1,m)
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL
      !
    else
      call divrot2uv_dfs_old__run(qdiv,qrot,qu,qv) !INOUT,INOUT,OUT,OUT
      
    end if


!    write(0,*) '1:mg,qu(mg,:,1)*ER=',mg_debug,qu(mg_debug,:,1)*ER
!    write(0,*) '1:mg,qv(mg,:,1)*ER=',mg_debug,qv(mg_debug,:,1)*ER

    
    call e_time__end(38,"divrot2uv_dfs__run")
    !
  end subroutine divrot2uv_dfs__run


  subroutine divrot2uv_dfs__poisson( qvar1, qvar2 ) !IN,OUT
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
    if ( jcn_dfs >= 1 ) then
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
        if ( jcn_dfs >= 1 ) then
          call calc_abc
        end if
        !
      end if
      !
      call poisson_core           &
       &( qvar1,                   &!IN
       &  qvar2      )              !OUT
    else
      call divrot2uv_dfs_old__poisson( qvar1, qvar2 ) !IN,OUT
    end if
    !
  end subroutine divrot2uv_dfs__poisson


  subroutine poisson_core     &
   &( qvar1,                   &!IN
   &  qvar2      )              !OUT
    !
    complex(8),intent(inout) :: qvar1(NNUM,0:MMAX)
    complex(8),intent(out) :: qvar2(NNUM,0:MMAX)
    !
    complex(8) :: g1(NNUMHF1+1)
    complex(8) :: g2(NNUMHF2+1)
    complex(8) :: x1(NNUMHF1+1)
    complex(8) :: x2(NNUMHF2+1)
    complex(8) :: work(max(NNUMHF1,NNUMHF2)+1)
    !
    integer,parameter :: nnmax=8
    !
    integer :: m,l,ll,l2
    !
    !
    call gmean_to_zero_dfs__run(qvar1) !INOUT
    
!   write(0,*) 'poisson_core:mg,qvar1(mg,:,1)*ER*ER=',mg_debug,qvar1(mg_debug,:,1)*ER*ER
!   write(0,*)

    !
!$OMP PARALLEL default(SHARED), private(m,l,ll,l2,g1,g2,x1,x2,work)
!$OMP DO schedule(STATIC)
    do m=0,mmax
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        l2= 1
        l = 1
        g1(l2) = 10.0d0*qvar1(l,m) - 5.0d0*qvar1(l+2,m) + qvar1(l+4,m)  !! n=1
        l = 2
        g2(l2) = 5.0d0*qvar1(l,m) - 4.0d0*qvar1(l+2,m) + qvar1(l+4,m)   !! n=2

        l2= 2
        l = 3
        g1(l2) = - 5.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) &     !! n=3
         &       - 4.0d0*qvar1(l+2,m) + qvar1(l+4,m)
        if ( NNUMHF2 >= 4 ) then
          l = 4
          g2(l2) = - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) &
           &       - 4.0d0*qvar1(l+2,m) + qvar1(l+4,m)
        end if

        do l2=3,NNUMHF1-2
          l=l2*2-1
          g1(l2) = qvar1(l-4,m) - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) &
           &       - 4.0d0*qvar1(l+2,m) + qvar1(l+4,m)
        end do
        do l2=3,NNUMHF2-2
          l=l2*2
          g2(l2) = qvar1(l-4,m) - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) &
           &       - 4.0d0*qvar1(l+2,m) + qvar1(l+4,m)
        end do

        l2=NNUMHF1-1
        l=l2*2-1
        g1(l2) = qvar1(l-4,m) - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) &
         &       - 4.0d0*qvar1(l+2,m)
        l2=NNUMHF2-1
        l=l2*2
        if ( NNUMHF2 >= 4 ) then
          g2(l2) = qvar1(l-4,m) - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) &
           &       - 4.0d0*qvar1(l+2,m)
        else
          g2(l2) = - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m) - 4.0d0*qvar1(l+2,m)
        end if

        l2=NNUMHF1
        l=l2*2-1
        g1(l2) = qvar1(l-4,m) - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m)
        l2=NNUMHF2
        l=l2*2
        g2(l2) = qvar1(l-4,m) - 4.0d0*qvar1(l-2,m) + 6.0d0*qvar1(l,m)
        
      else if ( m == 0 ) then
        l2= 1
        l = 1
        g1(l2) = 0.0d0    !!! Attention, not used !!!
        l = 2
        g2(l2) = qvar1(l,m) - qvar1(l+2,m)       !! n=1

        l2= 2
        l = 3
        g1(l2) = -2.0d0*qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)  !! n=2
        l = 4
        g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)

        do l2=3,NNUMHF1-1
          l=l2*2-1
          g1(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        end do
        do l2=3,NNUMHF2-1
          l=l2*2
          g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        end do

        l2=NNUMHF1
        l=l2*2-1
        g1(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m)
        l2=NNUMHF2
        l=l2*2
        g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m)

      else if ( mod(m,2) == 1 ) then
        l2= 1
        l = 1
        g1(l2) = 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        l = 2
        g2(l2) = qvar1(l,m) - qvar1(l+2,m)                        !! n=1

        l2= 2
        l = 3
        g1(l2) = -2.0d0*qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)  !! n=2
        l = 4
        g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)

        do l2=3,NNUMHF1-1
          l=l2*2-1
          g1(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        end do
        do l2=3,NNUMHF2-1
          l=l2*2
          g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        end do

        l2=NNUMHF1
        l=l2*2-1
        g1(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m)
        l2=NNUMHF2
        l=l2*2
        g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m)

      else    !! mod(m,2) == 0
        l2= 1
        l = 1
        g1(l2) = 3.0d0*qvar1(l,m) - qvar1(l+2,m)
        l = 2
        g2(l2) = 2.0d0*qvar1(l,m) - qvar1(l+2,m)

        do l2=2,NNUMHF1-1
          l=l2*2-1
          g1(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        end do
        do l2=2,NNUMHF2-1
          l=l2*2
          g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m) - qvar1(l+2,m)
        end do

        l2=NNUMHF1
        l=l2*2-1
        g1(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m)
        l2=NNUMHF2
        l=l2*2
        g2(l2) = -qvar1(l-2,m) + 2.0d0*qvar1(l,m)
      end if

      !! Calculate x1,x2 from g1,g2        
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        call pentadiagonal__solve( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww, &
         &                         work, g1, x1 )
        call pentadiagonal__solve( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww, &
         &                         work, g2, x2 )
      else
        call tridiagonal__solve( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww, &
         &                       work, g1, x1 )
        call tridiagonal__solve( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww, &
         &                       work, g2, x2 )
      end if

      do l2=1,NNUMHF1
        l = l2*2-1
        qvar2(l,m) = x1(l2)     !! l   = 1,3,5,...,JMAX-1
      end do
      do l2=1,NNUMHF2
        l = l2*2
        qvar2(l,m) = x2(l2)     !! l+1 = 2,4,6,...,JMAX
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    call gmean_to_zero_dfs__run(qvar2) !INOUT
    !
  end subroutine poisson_core


  subroutine calc_abc
    !
    integer :: m,l2,nn
    real(8) :: am,al
        
    allocate( td1(0:MMAX) )
    allocate( td2(0:MMAX) )

!    write(0,*) '000000000000000000'

!$OMP PARALLEL default(SHARED), private(m,l2,nn,am,al)
!$OMP DO schedule(STATIC)
    do m=0,mmax
      am=m
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        call pentadiagonal__allocate( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call pentadiagonal__allocate( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      else
        call tridiagonal__allocate( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call tridiagonal__allocate( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        do l2=1,NNUMHF1
          al=l2*2-1
          td1(m)%d(l2,1) = -(al-2)*(al-1)/ER**2
          td1(m)%d(l2,2) = ( 4*al*al -6*al +4 +4*am*am )/ER**2
          td1(m)%d(l2,3) = ( -6*al*al -4 -8*am*am )/ER**2
          td1(m)%d(l2,4) = ( 4*al*al +6*al +4 +4*am*am )/ER**2
          td1(m)%d(l2,5) = -(al+2)*(al+1)/ER**2
        end do
        l2=1
        al=1
        td1(m)%d(l2,1) = -999.0d33
        td1(m)%d(l2,2) = -999.0d33
        td1(m)%d(l2,3) = ( -12 - 12*am*am )/ER**2  !! n=1
        l2=2
        al=3
        td1(m)%d(l2,1) = -999.0d33
        td1(m)%d(l2,2) = ( 24 +4*am*am )/ER**2     !! n=3
        
        do l2=1,NNUMHF2
          al=l2*2
          td2(m)%d(l2,1) = -(al-2)*(al-1)/ER**2
          td2(m)%d(l2,2) = ( 4*al*al -6*al +4 +4*am*am )/ER**2
          td2(m)%d(l2,3) = ( -6*al*al -4 -8*am*am )/ER**2
          td2(m)%d(l2,4) = ( 4*al*al +6*al +4 +4*am*am )/ER**2
          td2(m)%d(l2,5) = -(al+2)*(al+1)/ER**2
        end do
      else if ( m == 0 ) then
        do l2=1,NNUMHF1           !! m = 0
          al=l2*2-2                   !! l = 0,2,4,6,...,JMAX-2
          td1(m)%d(l2,1) = (al-1)*(al-2)/ER**2
          td1(m)%d(l2,2) = -2*al*al/ER**2
          td1(m)%d(l2,3) = (al+1)*(al+2)/ER**2
        end do
        l2=1
        td1(m)%d(l2,1) = 0.0d0
        td1(m)%d(l2,2) = 1.0d0      !!!!! Attention !!!!!
        td1(m)%d(l2,3) = 0.0d0      !!!!! Attention !!!!!
        
        do l2=1,NNUMHF2           !! m = 0
          al=l2*2-1                      !! l = 1,3,5,...,JMAX-1
          td2(m)%d(l2,1) = (al-1)*(al-2)/ER**2
          td2(m)%d(l2,2) = -2*al*al/ER**2
          td2(m)%d(l2,3) = (al+1)*(al+2)/ER**2
        end do

      else if ( mod(m,2) == 1 ) then
        do l2=1,NNUMHF1
          al=l2*2-2                 !! l = 1,3,5,...,JMAX-1
          td1(m)%d(l2,1) = (al-1)*al/ER**2
          td1(m)%d(l2,2) = ( -2*al*al -4*am*am )/ER**2
          td1(m)%d(l2,3) = (al+1)*al/ER**2
        end do
        l2=2
        al=2
        td1(m)%d(l2,1) = 2*(al-1)*al/ER**2   !! n=2
        
        do l2=1,NNUMHF2
          al=l2*2-1                    !! l = 2,4,6,...,JMAX
          td2(m)%d(l2,1) = (al-1)*al/ER**2
          td2(m)%d(l2,2) = ( -2*al*al -4*am*am )/ER**2
          td2(m)%d(l2,3) = (al+1)*al/ER**2
        end do

      else  !! mod(m,2) == 0
        do l2=1,NNUMHF1
          al=l2*2-1                 !! l = 1,3,5,...,JMAX-1
          td1(m)%d(l2,1) = (al-1)*al/ER**2
          td1(m)%d(l2,2) = ( -2*al*al -4*am*am )/ER**2
          td1(m)%d(l2,3) = (al+1)*al/ER**2
        end do
        do l2=1,NNUMHF2
          al=l2*2                      !! l = 2,4,6,...,JMAX
          td2(m)%d(l2,1) = (al-1)*al/ER**2
          td2(m)%d(l2,2) = ( -2*al*al -4*am*am )/ER**2
          td2(m)%d(l2,3) = (al+1)*al/ER**2
        end do
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        call pentadiagonal__ini( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call pentadiagonal__ini( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      else
        call tridiagonal__ini( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call tridiagonal__ini( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine calc_abc

end module divrot2uv_dfs
