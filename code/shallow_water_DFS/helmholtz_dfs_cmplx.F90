module helmholtz_dfs_cmplx

  !! Helmholtz equation

  use prm_phconst, only : ER
  use com_dfs, only : NUMNA_I=>NNUM, NNUMHF1, NNUMHF2, MMAX
  use com_dfs, only : jcn_dfs,                                                                    &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS, N_TRUNC_M3_DFS,             &
   &                  N_L2TRUNC1_M0_DFS, N_L2TRUNC2_M0_DFS, N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS, &
   &                  N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS, N_L2TRUNC1_M3_DFS, N_L2TRUNC2_M3_DFS, &
   &                  NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS
  use pentadiagonal_cmplx, only : pentadiagonal_cmplx__allocate, pentadiagonal_cmplx__ini, pentadiagonal_cmplx__solve
  use tridiagonal_cmplx, only : tridiagonal_cmplx__allocate, tridiagonal_cmplx__ini, tridiagonal_cmplx__solve
  use helmholtz_dfs_cmplx_old, only : helmholtz_dfs_cmplx_old__run1, helmholtz_dfs_cmplx_old__run2
  
  implicit none

  private
  public :: helmholtz_dfs_cmplx__run1, helmholtz_dfs_cmplx__run2
  
  type :: type_diagonal
    private
    complex(8),allocatable :: d(:,:)
    integer,allocatable :: iww(:)
  end type

contains
  
  subroutine helmholtz_dfs_cmplx__run1 &
   &( delt, gamma,        &!IN
   &  qvar         )       !INOUT
    !
    real(8),intent(in) :: delt
    complex(8),intent(in) :: gamma
    complex(8),intent(inout) :: qvar(NUMNA_I,0:MMAX)
    real(8),save :: delt_save = -999.0d0
    !
    type(type_diagonal),save,allocatable :: td1(:)
    type(type_diagonal),save,allocatable :: td2(:)
    !
    logical,save :: first_helmholtz_dfs_cmplx = .true.
    integer :: m
    !
    if ( jcn_dfs >= 1 ) then 
      if ( first_helmholtz_dfs_cmplx ) then
        !
        allocate( td1(0:MMAX) )
        allocate( td2(0:MMAX) )
       !$OMP PARALLEL default(SHARED), private(m)
       !$OMP DO schedule(STATIC)
        do m=0,MMAX
          if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
            call pentadiagonal_cmplx__allocate( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
            call pentadiagonal_cmplx__allocate( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
          else
            call tridiagonal_cmplx__allocate( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
            call tridiagonal_cmplx__allocate( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
          end if
        end do
       !$OMP END DO
       !$OMP END PARALLEL
      end if
      !
      first_helmholtz_dfs_cmplx = .false.
      !
      if ( delt_save /= delt ) then
        call calc_abcde_helmholtz &
         &( gamma, td1, td2 ) !IN,IO,IO
        delt_save = delt
      end if
    !
      call do_helmholtz &
       &( td1, td2, qvar            )         !INOUT
    else
      call helmholtz_dfs_cmplx_old__run1 &
       &( delt, gamma,        &!IN
       &  qvar         )       !INOUT
    end if
    !
  end subroutine helmholtz_dfs_cmplx__run1
  
  
  subroutine helmholtz_dfs_cmplx__run2 &
   &( delt, gamma,        &!IN
   &  qvar         )       !INOUT
    !
    real(8),intent(in) :: delt
    complex(8),intent(in) :: gamma
    complex(8),intent(inout) :: qvar(NUMNA_I,0:MMAX)
    real(8),save :: delt_save = -999.0d0
    !
    type(type_diagonal),save,allocatable :: td1(:)
    type(type_diagonal),save,allocatable :: td2(:)
    !
    logical,save :: first_helmholtz_dfs_cmplx = .true.
    integer :: m
    !
    if ( jcn_dfs >= 1 ) then
      if ( first_helmholtz_dfs_cmplx ) then
      !
        allocate( td1(0:MMAX) )
        allocate( td2(0:MMAX) )
       !$OMP PARALLEL default(SHARED), private(m)
       !$OMP DO schedule(STATIC)
        do m=0,MMAX
          if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
            call pentadiagonal_cmplx__allocate( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
            call pentadiagonal_cmplx__allocate( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
          else
            call tridiagonal_cmplx__allocate( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
            call tridiagonal_cmplx__allocate( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
          end if
        end do
       !$OMP END DO
       !$OMP END PARALLEL
      end if
      !
      first_helmholtz_dfs_cmplx = .false.
      !
      if ( delt_save /= delt ) then
        call calc_abcde_helmholtz &
         &( gamma, td1, td2 ) !IN,IO,IO
        delt_save = delt
      end if
      !
      call do_helmholtz              &
       &( td1, td2, qvar            )       !INOUT
    else
     call helmholtz_dfs_cmplx_old__run2 &
      &( delt, gamma,        &!IN
      &  qvar         )       !INOUT
    end if
    !
  end subroutine helmholtz_dfs_cmplx__run2
  
  
  subroutine do_helmholtz      &
   &( td1, td2, qvar            )     !INOUT
    !
    type(type_diagonal),intent(inout) :: td1(0:MMAX)
    type(type_diagonal),intent(inout) :: td2(0:MMAX)
    complex(8),intent(inout) :: qvar(NUMNA_I,0:MMAX)
    !
    complex(8) :: qtmp(NUMNA_I+6)
    complex(8) :: g1(NNUMHF1+1)
    complex(8) :: g2(NNUMHF1+1)   !! NNUMHF1 >= NNUMHF2
    complex(8) :: x1(NNUMHF1+1)
    complex(8) :: x2(NNUMHF1+1)   !! NNUMHF1 >= NNUMHF2
    complex(8) :: work(NNUMHF1+1) !! NNUMHF1 >= NNUMHF2
    !
    integer :: k,m,l,ll,l2,nn1,nn2
    complex(8) :: x1aa,x2aa,ww
    !
   !$OMP PARALLEL default(SHARED), private(k,m,l,ll,l2,nn1,nn2,x1aa,x2aa,ww,qtmp,g1,g2,x1,x2,work)
   !$OMP DO schedule(STATIC)
    do m=0,MMAX
      qtmp(1:NUMNA_I) = qvar(1:NUMNA_I,m)
      qtmp(NUMNA_I+1:NUMNA_I+6) = 0.0d0

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        l2= 1
        l = 1
        ll = 2
        g1(l2) = 10.0d0*qtmp(l) - 5.0d0*qtmp(l+2) + qtmp(l+4)
        g2(l2) = 5.0d0*qtmp(ll) - 4.0d0*qtmp(ll+2) + qtmp(ll+4)
  
        l2= 2
        l = 3
        ll = 4
        g1(l2) = -5.0d0*qtmp(l-2) + 6.0d0*qtmp(l) -4.0d0*qtmp(l+2) + qtmp(l+4)
        g2(l2) = -4.0d0*qtmp(ll-2) + 6.0d0*qtmp(ll) -4.0d0*qtmp(ll+2) + qtmp(ll+4)
      else if ( m == 0 .or. mod(m,2) == 1 ) then
        l2= 1
        l = 1
        ll = 2
        g1(l2) = 2.0d0*qtmp(l) - qtmp(l+2)
        g2(l2) = qtmp(ll) - qtmp(ll+2)
  
        l2= 2
        l = 3
        ll = 4
        g1(l2) = - 2.0d0*qtmp(l-2) + 2.0d0*qtmp(l) - qtmp(l+2)
        g2(l2) = - qtmp(ll-2) + 2.0d0*qtmp(ll) - qtmp(ll+2)
      else     !! mod(m,2) == 0
        l2= 1
        l = 1
        ll = 2
        g1(l2) = 3.0d0*qtmp(l) - qtmp(l+2)
        g2(l2) = 2.0d0*qtmp(ll) - qtmp(ll+2)
  
        l2= 2
        l = 3
        ll = 4
        g1(l2) = - qtmp(l-2) + 2.0d0*qtmp(l) - qtmp(l+2)
        g2(l2) = - qtmp(ll-2) + 2.0d0*qtmp(ll) - qtmp(ll+2)
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        do l2=3,NNUM1_M_DFS(m)
          l=l2*2-1
          ll=l2*2
          g1(l2) = qtmp(l-4) - 4.0d0*qtmp(l-2) + 6.0d0*qtmp(l) &
           &       - 4.0d0*qtmp(l+2) + qtmp(l+4)
          g2(l2) = qtmp(ll-4) - 4.0d0*qtmp(ll-2) + 6.0d0*qtmp(ll) &
           &       - 4.0d0*qtmp(ll+2) + qtmp(ll+4)
        end do
      else
        do l2=3,NNUM1_M_DFS(m)
          l=l2*2-1
          ll=l2*2
          g1(l2) = -qtmp(l-2) + 2.0d0*qtmp(l) - qtmp(l+2)
          g2(l2) = -qtmp(ll-2) + 2.0d0*qtmp(ll) - qtmp(ll+2)
        end do
      end if

      !! Calculate x1,x2 from g1,g2        
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        call pentadiagonal_cmplx__solve( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww, &
         &                         work, g1, x1 )
        call pentadiagonal_cmplx__solve( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww, &
         &                         work, g2, x2 )
      else
        call tridiagonal_cmplx__solve( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww, &
         &                       work, g1, x1 )
        call tridiagonal_cmplx__solve( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww, &
         &                       work, g2, x2 )
      end if

      do l2=1,NNUMHF1
        l = l2*2-1
        qvar(l,m) = x1(l2)     !! l   = 1,3,5,...,JMAX-1
      end do
      do l2=1,NNUMHF2
        l = l2*2
        qvar(l,m) = x2(l2)     !! l+1 = 2,4,6,...,JMAX
      end do

!      do l2=1,NNUMHF2
!        l = l2*2-1
!        qvar(l,m) = x1(l2)       !! l   = 1,3,5,...,JMAX-1
!        qvar(l+1,m) = x2(l2)     !! l+1 = 2,4,6,...,JMAX
!      end do
!      if ( NNUMHF1 > NNUMHF2 ) then   
!        l2=NNUMHF1
!        l = l2*2-1     
!        qvar(l,m) = x1(l2)     !! l   = 1,3,5,...,JMAX-1
!      end if
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    !
  end subroutine do_helmholtz
  
  
  subroutine calc_abcde_helmholtz  &
   &( gamma, td1, td2 ) !IN

    complex(8),intent(in) :: gamma
    type(type_diagonal),intent(inout) :: td1(0:MMAX)
    type(type_diagonal),intent(inout) :: td2(0:MMAX)
    
    real(8) :: aaa1(NNUMHF1+1)
    real(8) :: bbb1(NNUMHF1+1)
    real(8) :: ccc1(NNUMHF1+1)
    real(8) :: ddd1(NNUMHF1+1)
    real(8) :: eee1(NNUMHF1+1)
    real(8) :: aaa2(NNUMHF2+1)
    real(8) :: bbb2(NNUMHF2+1)
    real(8) :: ccc2(NNUMHF2+1)
    real(8) :: ddd2(NNUMHF2+1)
    real(8) :: eee2(NNUMHF2+1)
    complex(8) :: gamma2

    integer :: l2,m
    real(8) :: am,al

    gamma2 = gamma/ER**2
  
   !$OMP PARALLEL default(SHARED), private(l2,m,am,al,aaa1,bbb1,ccc1,ddd1,eee1,aaa2,bbb2,ccc2,ddd2,eee2)
   !$OMP DO schedule(STATIC)
    do m=0,MMAX
      am=m
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        do l2=1,NNUMHF1
          td1(m)%d(l2,1) = 1.0d0
          td1(m)%d(l2,2) = -4.0d0
          td1(m)%d(l2,3) = 6.0d0
          td1(m)%d(l2,4) = -4.0d0
          td1(m)%d(l2,5) = 1.0d0
        end do
        td1(m)%d(1,1) = -999.0d33
        td1(m)%d(1,2) = -999.0d33
        td1(m)%d(1,3) = 10.0d0
        td1(m)%d(1,4) = -5.0d0
        td1(m)%d(2,1)  = -999.0d33
        td1(m)%d(2,2) = -5.0d0
        
        do l2=1,NNUMHF2
          td2(m)%d(l2,1) = 1.0d0
          td2(m)%d(l2,2) = -4.0d0
          td2(m)%d(l2,3) = 6.0d0
          td2(m)%d(l2,4) = -4.0d0
          td2(m)%d(l2,5) = 1.0d0
        end do
        td2(m)%d(1,1) = -999.0d33
        td2(m)%d(1,2) = -999.0d33
        td2(m)%d(1,3) = 5.0d0
        
      else if ( m == 0 .or. mod(m,2) == 1 ) then
        do l2=1,NNUMHF1
          td1(m)%d(l2,1) = -1.0d0
          td1(m)%d(l2,2) = 2.0d0
          td1(m)%d(l2,3) = -1.0d0
        end do                    
        td1(m)%d(2,1) = -2.0d0
        
        do l2=1,NNUMHF2
          td2(m)%d(l2,1) = -1.0d0
          td2(m)%d(l2,2) = 2.0d0
          td2(m)%d(l2,3) = -1.0d0
        end do         
        td2(m)%d(1,1) = 999.0d33
        td2(m)%d(1,2) = 1.0d0
      else   !! mod(m,2) == 0
        do l2=1,NNUMHF1
          td1(m)%d(l2,1) = -1.0d0
          td1(m)%d(l2,2) = 2.0d0
          td1(m)%d(l2,3) = -1.0d0
        end do
        td1(m)%d(1,2) = 3.0d0
        
        do l2=1,NNUMHF2
          td2(m)%d(l2,1) = -1.0d0
          td2(m)%d(l2,2) = 2.0d0
          td2(m)%d(l2,3) = -1.0d0
        end do
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        do l2=1,NNUMHF1+1
          al=l2*2-1
          aaa1(l2) = -(al-2)*(al-1)
          bbb1(l2) = 4*al*al -6*al +4 +4*am*am
          ccc1(l2) = -6*al*al -4 -8*am*am
          ddd1(l2) = 4*al*al +6*al +4 +4*am*am
          eee1(l2) = -(al+2)*(al+1)
        end do
        l2=1
        al=1
        aaa1(l2) = -999.0d33
        bbb1(l2) = -999.0d33
        ccc1(l2) = -12 - 12*am*am
        l2=2
        al=3
        aaa1(l2) = -999.0d33
        bbb1(l2) = 24 +4*am*am
        
        do l2=1,NNUMHF2+1
          al=l2*2
          aaa2(l2) = -(al-2)*(al-1)
          bbb2(l2) = 4*al*al -6*al +4 +4*am*am
          ccc2(l2) = -6*al*al -4 -8*am*am
          ddd2(l2) = 4*al*al +6*al +4 +4*am*am
          eee2(l2) = -(al+2)*(al+1)
        end do
      else if ( m == 0 ) then
        do l2=1,NNUMHF1+1
          al=l2*2-2                    !! l = 0,2,4,6,...,JMAX-2
          aaa1(l2) = 0.0d0
          bbb1(l2) = (al-1)*(al-2)
          ccc1(l2) = -2*al*al
          ddd1(l2) = (al+1)*(al+2)
          eee1(l2) = 0.0d0
        end do
        
        do l2=1,NNUMHF2+1
          al=l2*2-1                   !! l = 1,3,5,...,JMAX-1
          aaa2(l2) = 0.0d0
          bbb2(l2) = (al-1)*(al-2)
          ccc2(l2) = -2*al*al
          ddd2(l2) = (al+1)*(al+2)
          eee2(l2) = 0.0d0
        end do
      else if ( mod(m,2) == 1 ) then
        do l2=1,NNUMHF1+1
          al=l2*2-2
          aaa1(l2) = 0.0d0
          bbb1(l2) = (al-1)*al
          ccc1(l2) = -2*al*al -4*am*am
          ddd1(l2) = (al+1)*al
          eee1(l2) = 0.0d0
        end do
        l2=2
        al=2
        bbb1(l2) = 4
        
        do l2=1,NNUMHF2+1
          al=l2*2-1
          aaa2(l2) = 0.0d0
          bbb2(l2) = (al-1)*al
          ccc2(l2) = -2*al*al -4*am*am
          ddd2(l2) = (al+1)*al
          eee2(l2) = 0.0d0
        end do
      else
        do l2=1,NNUMHF1+1
          al=l2*2-1
          aaa1(l2) = 0.0d0
          bbb1(l2) = (al-1)*al
          ccc1(l2) = -2*al*al -4*am*am
          ddd1(l2) = (al+1)*al
          eee1(l2) = 0.0d0
        end do
        
        do l2=1,NNUMHF2+1
          al=l2*2
          aaa2(l2) = 0.0d0
          bbb2(l2) = (al-1)*al
          ccc2(l2) = -2*al*al -4*am*am
          ddd2(l2) = (al+1)*al
          eee2(l2) = 0.0d0
        end do
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        do l2=1,NNUMHF1
          td1(m)%d(l2,1) = td1(m)%d(l2,1) - aaa1(l2)*gamma2
          td1(m)%d(l2,2) = td1(m)%d(l2,2) - bbb1(l2)*gamma2
          td1(m)%d(l2,3) = td1(m)%d(l2,3) - ccc1(l2)*gamma2
          td1(m)%d(l2,4) = td1(m)%d(l2,4) - ddd1(l2)*gamma2
          td1(m)%d(l2,5) = td1(m)%d(l2,5) - eee1(l2)*gamma2
        end do
        do l2=1,NNUMHF2
          td2(m)%d(l2,1) = td2(m)%d(l2,1) - aaa2(l2)*gamma2
          td2(m)%d(l2,2) = td2(m)%d(l2,2) - bbb2(l2)*gamma2
          td2(m)%d(l2,3) = td2(m)%d(l2,3) - ccc2(l2)*gamma2
          td2(m)%d(l2,4) = td2(m)%d(l2,4) - ddd2(l2)*gamma2
          td2(m)%d(l2,5) = td2(m)%d(l2,5) - eee2(l2)*gamma2
        end do
      else
        do l2=1,NNUMHF1
          td1(m)%d(l2,1) = td1(m)%d(l2,1) - bbb1(l2)*gamma2
          td1(m)%d(l2,2) = td1(m)%d(l2,2) - ccc1(l2)*gamma2
          td1(m)%d(l2,3) = td1(m)%d(l2,3) - ddd1(l2)*gamma2
        end do
        do l2=1,NNUMHF2
          td2(m)%d(l2,1) = td2(m)%d(l2,1) - bbb2(l2)*gamma2
          td2(m)%d(l2,2) = td2(m)%d(l2,2) - ccc2(l2)*gamma2
          td2(m)%d(l2,3) = td2(m)%d(l2,3) - ddd2(l2)*gamma2
        end do
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        call pentadiagonal_cmplx__ini( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call pentadiagonal_cmplx__ini( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      else
        call tridiagonal_cmplx__ini( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call tridiagonal_cmplx__ini( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      end if
    end do
   !$OMP END DO
   !$OMP END PARALLEL
  
  end subroutine calc_abcde_helmholtz

end module helmholtz_dfs_cmplx
