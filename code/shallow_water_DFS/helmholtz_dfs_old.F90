module helmholtz_dfs_old

  !! Helmholtz equation

  use prm_phconst, only : ER
  use com_dfs, only : NUMNA_I=>NNUM, NNUMHF1, NNUMHF2, MMAX
  use com_dfs, only : jcn_dfs,                                                                    &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS, N_TRUNC_M3_DFS,             &
   &                  N_L2TRUNC1_M0_DFS, N_L2TRUNC2_M0_DFS, N_L2TRUNC1_M1_DFS, N_L2TRUNC2_M1_DFS, &
   &                  N_L2TRUNC1_M2_DFS, N_L2TRUNC2_M2_DFS, N_L2TRUNC1_M3_DFS, N_L2TRUNC2_M3_DFS, &
   &                  NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS
  use pentadiagonal, only : pentadiagonal__allocate, pentadiagonal__ini, pentadiagonal__solve
  use tridiagonal, only : tridiagonal__allocate, tridiagonal__ini, tridiagonal__solve
  
  implicit none

  private
  public :: helmholtz_dfs_old__run, helmholtz_dfs_old__hdiff2
  
  type :: type_diagonal
    private
    real(8),allocatable :: d(:,:)
    integer,allocatable :: iww(:)
  end type

contains
  
  subroutine helmholtz_dfs_old__run &
   &( delt, gamma,        &!IN
   &  qvar         )       !INOUT
    !
    real(8),intent(in) :: delt
    real(8),intent(in) :: gamma
    complex(8),intent(inout) :: qvar(NUMNA_I,0:MMAX)
    real(8),save :: delt_save = -999.0d0
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
    logical,save :: first_helmholtz_dfs = .true.
    !
    if ( first_helmholtz_dfs ) then
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
      !
      first_helmholtz_dfs = .false.
      !
    end if
    !
    if ( delt_save /= delt ) then
      call calc_abc_helmholtz_old &
       &( gamma,          &!IN
       &  aa1, bb1, cc1,  &!OUT
       &  aa2, bb2, cc2,  &!OUT
       &  a1, b1, c1,     &!OUT
       &  a2, b2, c2   )   !OUT
      delt_save = delt
    end if
    !
    call do_helmholtz_old         &
     &( aa1, bb1, cc1,  &!IN
     &  aa2, bb2, cc2,  &!IN
     &  a1, b1, c1,     &!IN
     &  a2, b2, c2,     &!IN
     &  qvar      )      !INOUT
    !
  end subroutine helmholtz_dfs_old__run
  
  
  subroutine helmholtz_dfs_old__hdiff2 &
   &( delt, gamma,        &!IN
   &  qvar         )       !INOUT
    !
    real(8),intent(in) :: delt
    real(8),intent(in) :: gamma
    complex(8),intent(inout) :: qvar(NUMNA_I,0:MMAX)
    real(8),save :: delt_save = -999.0d0
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
    logical,save :: first_helmholtz_dfs = .true.
    !
    if ( first_helmholtz_dfs ) then
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
      !
      first_helmholtz_dfs = .false.
      !
    end if
    !
    if ( delt_save /= delt ) then
      call calc_abc_helmholtz_old &
       &( gamma,  &!IN
       &  aa1, bb1, cc1,  &!IN
       &  aa2, bb2, cc2,  &!IN
       &  a1, b1, c1,     &!OUT
       &  a2, b2, c2   )   !OUT
      delt_save = delt
    end if
    !
    call do_helmholtz_old         &
     &( aa1, bb1, cc1,  &!IN
     &  aa2, bb2, cc2,  &!IN
     &  a1, b1, c1,     &!IN
     &  a2, b2, c2,     &!IN
     &  qvar      )      !INOUT
    !
  end subroutine helmholtz_dfs_old__hdiff2
  
  
  subroutine do_helmholtz_old  &
   &( aa1, bb1, cc1, &!IN
   &  aa2, bb2, cc2, &!IN
   &  a1, b1, c1,    &!IN
   &  a2, b2, c2,    &!IN
   &  qvar      )     !INOUT
    !
    real(8),intent(in) :: aa1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: bb1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: cc1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: aa2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: bb2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: cc2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: a1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: b1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: c1(NNUMHF1+1,0:MMAX)
    real(8),intent(in) :: a2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: b2(NNUMHF2+1,0:MMAX)
    real(8),intent(in) :: c2(NNUMHF2+1,0:MMAX)
    complex(8),intent(inout) :: qvar(NUMNA_I,0:MMAX)
    !
    complex(8) :: g1(NNUMHF1+1)
    complex(8) :: g2(NNUMHF2+1)
    !
    complex(8) :: x1(NNUMHF1+1)
    complex(8) :: x2(NNUMHF2+1)
    !
    integer :: m,l,ll,l2
    !
   !$OMP PARALLEL default(SHARED), private(m,l,ll,l2,g1,g2,x1,x2)
   !$OMP DO schedule(STATIC)
    do m=0,MMAX
       l2= 1
       l = 1
      ll = 2
       g1(l2) = qvar(l,m)*bb1(l2,m) + qvar(l+2,m)*cc1(l2,m)
       g2(l2) = qvar(ll,m)*bb2(l2,m) + qvar(ll+2,m)*cc2(l2,m)

       do l2=2,NNUMHF1-1
         l=l2*2-1
         ll=l+1
         g1(l2) = qvar(l-2,m)*aa1(l2,m) + qvar(l,m)*bb1(l2,m) + qvar(l+2,m)*cc1(l2,m)
       end do

       do l2=2,NNUMHF2-1
       l=l2*2-1
       ll=l+1
         g2(l2) = qvar(ll-2,m)*aa2(l2,m) + qvar(ll,m)*bb2(l2,m) + qvar(ll+2,m)*cc2(l2,m)
       end do

 
       l2=NNUMHF1
       l=l2*2-1
       ll=l+1
       g1(l2) = qvar(l-2,m)*aa1(l2,m) + qvar(l,m)*bb1(l2,m)
 
       l2=NNUMHF2
       l=l2*2-1
       ll=l+1
       g2(l2) = qvar(ll-2,m)*aa2(l2,m) + qvar(ll,m)*bb2(l2,m)
 
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
 
 !      if ( mg == 0 .and. k == 1 ) then
 !        write(6,*) '1:x2(:,m)=',x2(:,m)
 !      end if
 
       do l2=1,NNUMHF1
       l = l2*2-1
         qvar(l  ,m) = x1(l2)     !! l   = 1,3,5,...,JMAX-1
       end do

       do l2=1,NNUMHF2
         l = l2*2-1
         qvar(l+1,m) = x2(l2)     !! l+1 = 2,4,6,...,JMAX
       end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    !
  end subroutine do_helmholtz_old
  
  
  subroutine calc_abc_helmholtz_old &
   &( gamma,          &!IN
   &  aa1, bb1, cc1,  &!OUT
   &  aa2, bb2, cc2,  &!OUT
   &  a1, b1, c1,     &!OUT
   &  a2, b2, c2   )   !OUT
    !
    real(8),intent(in) :: gamma
    real(8),intent(out) :: aa1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: bb1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: cc1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: aa2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: bb2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: cc2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: a1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: b1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: c1(NNUMHF1+1,0:MMAX)
    real(8),intent(out) :: a2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: b2(NNUMHF2+1,0:MMAX)
    real(8),intent(out) :: c2(NNUMHF2+1,0:MMAX)
    !
    real(8) :: aaa1(NNUMHF1+1)
    real(8) :: bbb1(NNUMHF1+1)
    real(8) :: ccc1(NNUMHF1+1)
    real(8) :: aaa2(NNUMHF2+1)
    real(8) :: bbb2(NNUMHF2+1)
    real(8) :: ccc2(NNUMHF2+1)
    real(8) :: gamma2
    !
    integer :: l,ll,l2,m
    !
    gamma2 = gamma/ER**2
    !
   !$OMP PARALLEL default(SHARED), private(l,ll,l2,m,aaa1,bbb1,ccc1,aaa2,bbb2,ccc2)
   !$OMP DO schedule(STATIC)
    do m=0,MMAX
        do l2=1,NNUMHF1+1
          aa1(l2,m) = - 1.0d0
          bb1(l2,m) = + 2.0d0
          cc1(l2,m) = - 1.0d0
        end do

        do l2=1,NNUMHF2+1
          aa2(l2,m) = - 1.0d0
          bb2(l2,m) = + 2.0d0
          cc2(l2,m) = - 1.0d0
        end do

        if ( m == 0 ) then
          bb2(1,m) = + 1.0d0
          aa1(2,m) = - 2.0d0
        else
          bb1(1,m) = + 3.0d0
        end if

        if ( m == 0 ) then          !! zonal wavenumber = 0
          do l2=1,NNUMHF1+1
            l=l2*2-2                   !! l = 0,2,4,6,...,JMAX-2,JMAX
            aaa1(l2) = (l-1)*(l-2)
            bbb1(l2) = -2*l*l
            ccc1(l2) = (l+1)*(l+2)
          end do
          do l2=1,NNUMHF2+1
            ll=l2*2-1                      !! l = 1,3,5,...,JMAX-1
            aaa2(l2) = (ll-1)*(ll-2)
            bbb2(l2) = -2*ll*ll
            ccc2(l2) = (ll+1)*(ll+2)
          end do
          
        else if ( mod(m,2) == 1 ) then   !! if zonal wavenumber is odd
            do l2=1,NNUMHF1+1
              l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
              aaa1(l2) = (l-1)*(l-2)
              bbb1(l2) = -2*l*l -4*m*m
              ccc1(l2) = (l+1)*(l+2)
            end do
            do l2=1,NNUMHF2+1
              ll=l2*2                   !! l = 2,4,6,...,JMAX
              aaa2(l2) = (ll-1)*(ll-2)
              bbb2(l2) = -2*ll*ll -4*m*m
              ccc2(l2) = (ll+1)*(ll+2)
            end do
        else                   !! if zonal wavenumber is even, not zero
            do l2=1,NNUMHF1+1
              l=l2*2-1                 !! l = 1,3,5,...,JMAX-1
              aaa1(l2) = l*(l-1)
              bbb1(l2) = -2*l*l -4*m*m
              ccc1(l2) = l*(l+1)
            end do
            do l2=1,NNUMHF2+1
              ll=l2*2                    !! l = 2,4,6,...,JMAX
              aaa2(l2) = ll*(ll-1)
              bbb2(l2) = -2*ll*ll -4*m*m
              ccc2(l2) = ll*(ll+1)
            end do
        end if

        do l2=1,NNUMHF1+1
          a1(l2,m) = aa1(l2,m) - aaa1(l2)*gamma2
          b1(l2,m) = bb1(l2,m) - bbb1(l2)*gamma2
          c1(l2,m) = cc1(l2,m) - ccc1(l2)*gamma2
        end do

        do l2=1,NNUMHF2+1
          a2(l2,m) = aa2(l2,m) - aaa2(l2)*gamma2
          b2(l2,m) = bb2(l2,m) - bbb2(l2)*gamma2
          c2(l2,m) = cc2(l2,m) - ccc2(l2)*gamma2
        end do

   
  !  if ( ML2GW(0,MYRANK_AGCM_KSET) == 0 ) then
  !    write(6,*) 'a1(:,0,1)=',a1(:,0,1)
  !    write(6,*) 'b1(:,0,1)=',b1(:,0,1)
  !    write(6,*) 'c1(:,0,1)=',c1(:,0,1)
  !    write(6,*) 'd1(:,0,1)=',d1(:,0,1)
  !  end if
  !  write(6,*) 'WAVE0A=',WAVE0A
    !
    !! --------------------------------------
    !!     LU decomposition
    !! --------------------------------------

        b1(1,m) = 1.0d0/b1(1,m)
        c1(1,m) = c1(1,m)*b1(1,m)

        do l2=2,max(N_L2TRUNC1_M0_DFS,N_L2TRUNC1_M1_DFS,N_L2TRUNC1_M2_DFS)
          b1(l2,m) = 1.0d0/( b1(l2,m) - a1(l2,m)*c1(l2-1,m) )
          c1(l2,m) = c1(l2,m)*b1(l2,m)
        end do

        do l2=min(N_L2TRUNC1_M0_DFS,N_L2TRUNC1_M1_DFS,N_L2TRUNC1_M2_DFS)+1,NNUMHF1+1
          if ( l2 >= NNUM1_M_DFS(m)+1 ) then
            a1(l2,m) = 0.0d0
            b1(l2,m) = 0.0d0
            c1(l2,m) = 0.0d0
          end if
        end do

        b2(1,m) = 1.0d0/b2(1,m)
        c2(1,m) = c2(1,m)*b2(1,m)

        do l2=2,max(N_L2TRUNC2_M0_DFS,N_L2TRUNC2_M1_DFS,N_L2TRUNC2_M2_DFS)
          b2(l2,m) = 1.0d0/( b2(l2,m) - a2(l2,m)*c2(l2-1,m) )
          c2(l2,m) = c2(l2,m)*b2(l2,m)
        end do

        do l2=min(N_L2TRUNC2_M0_DFS,N_L2TRUNC2_M1_DFS,N_L2TRUNC2_M2_DFS)+1,NNUMHF2+1
          if ( l2 >= NNUM2_M_DFS(m)+1 ) then
            a2(l2,m) = 0.0d0
            b2(l2,m) = 0.0d0
            c2(l2,m) = 0.0d0
          end if
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    !
  !  if ( ML2GW(0,MYRANK_AGCM_KSET) == 0 ) then
  !    write(6,*) '2:a1(:,0,1)=',a1(:,0,1)
  !    write(6,*) '2:b1(:,0,1)=',b1(:,0,1)
  !    write(6,*) '2:c1(:,0,1)=',c1(:,0,1)
  !    write(6,*) '2:d1(:,0,1)=',d1(:,0,1)
  !  end if
  
  end subroutine calc_abc_helmholtz_old
  
end module helmholtz_dfs_old
