module uv2divrot_dfs

  !! Calculate div,rot from U,V

  use prm_phconst,only : ER_INV,ER
  use com_dfs, only : jcn_dfs, MMAX, NNUM, NNUMHF1, NNUMHF2,                                      &
   &                  N_TRUNC_M0_DFS, N_TRUNC_M1_DFS, N_TRUNC_M2_DFS, N_TRUNC_M3_DFS,             &
   &                  NNUM1_M_DFS, NNUM2_M_DFS, NNUM_M_DFS
  use pentadiagonal_4, only : pentadiagonal_4__allocate, pentadiagonal_4__ini, pentadiagonal_4__solve
  use ninediagonal_4, only : ninediagonal_4__allocate, ninediagonal_4__ini, ninediagonal_4__solve
  use pentadiagonal, only : pentadiagonal__allocate, pentadiagonal__ini, pentadiagonal__solve
  use tridiagonal, only : tridiagonal__allocate, tridiagonal__ini, tridiagonal__solve
  use uv2divrot_dfs_old, only : uv2divrot_dfs_old__run, uv2divrot_dfs_old__laplacian
  use e_time, only : e_time__start, e_time__end

  implicit none

  private
  public :: uv2divrot_dfs__run, uv2divrot_dfs__uv2chipsi
  public :: uv2divrot_dfs__laplacian
  
  type :: type_diagonal
    private
    real(8),allocatable :: d(:,:)
    integer,allocatable :: iww(:)
  end type
  
  type(type_diagonal),save,allocatable :: td(:)

  type(type_diagonal),save,allocatable :: td1(:)
  type(type_diagonal),save,allocatable :: td2(:)

contains


  subroutine uv2divrot_dfs__run(qu,qv,qdiv,qrot) !IN,IN,OUT,OUT
    !
    complex(8),intent(IN) :: qu(NNUM,0:MMAX)   !! U
    complex(8),intent(IN) :: qv(NNUM,0:MMAX)   !! V
    complex(8),intent(OUT) :: qdiv(NNUM,0:MMAX) !! Divergence
    complex(8),intent(OUT) :: qrot(NNUM,0:MMAX) !! Rotation
    
    call e_time__start(37,"uv2divrot_dfs__run")

!    write(0,*)
!    write(0,*) 'cal_divrot: mg,qu(mg,:,1)*ER=',mg_debug,qu(mg_debug,:,1)*ER
!    write(0,*) 'cal_divrot: mg,qv(mg,:,1)*ER=',mg_debug,qv(mg_debug,:,1)*ER
    !
    if ( jcn_dfs >= 1 ) then
      call uv2divrot_dfs__uv2chipsi( qu, qv, qdiv, qrot )  !! u,v => chi,psi
      call uv2divrot_dfs__laplacian(qdiv)   !! chi => div
      call uv2divrot_dfs__laplacian(qrot)   !! psi => vor
      !
    else
      !
      call uv2divrot_dfs_old__run(qu,qv,qdiv,qrot) !IN,IN,OUT,OUT
      !
    endif
    
    call e_time__end(37,"uv2divrot_dfs__run")
    !
  end subroutine uv2divrot_dfs__run


  subroutine uv2divrot_dfs__uv2chipsi( qu, qv, qchi, qpsi )
    complex(8),intent(IN) :: qu(NNUM,0:MMAX)   !! U
    complex(8),intent(IN) :: qv(NNUM,0:MMAX)   !! V
    complex(8),intent(OUT) :: qchi(NNUM,0:MMAX) !! Velocity potential
    complex(8),intent(OUT) :: qpsi(NNUM,0:MMAX) !! Stream function
    !
    logical,save :: first_uv2chipsi = .true.
    !
    call e_time__start(36,"uv2_dfs__uv2chipsi")
    !
    if ( jcn_dfs >= 1 ) then
      if ( first_uv2chipsi ) then
        !
        first_uv2chipsi = .false.
        !
        call uv2chipsi_pre
        !
      end if

      call uv2chipsi_core             &
       &( qu, qv,                     &!IN
       &  qchi, qpsi )                 !OUT
    else
      write(6,*) "Error: uv2divrot_dfs__uv2chipsi: jcn_dfs == 0 is not supported in this subroutine."
      stop 999
    end if

    call e_time__end(36,"uv2_dfs_uv2chipsi")
    
    return
  end subroutine uv2divrot_dfs__uv2chipsi



  subroutine uv2chipsi_core       &
   &( qu, qv,                     &!IN
   &  qchi, qpsi )                 !OUT
    !
    complex(8),intent(IN) :: qu(NNUM,0:MMAX)   !! U
    complex(8),intent(IN) :: qv(NNUM,0:MMAX)   !! V
    complex(8),intent(OUT) :: qchi(NNUM,0:MMAX) !! Velocity potential
    complex(8),intent(OUT) :: qpsi(NNUM,0:MMAX) !! Stream function

    real(8) :: qu_c(0:NNUM)
    real(8) :: qu_s(0:NNUM)
    real(8) :: qv_c(0:NNUM)
    real(8) :: qv_s(0:NNUM)
    real(8) :: gg(NNUM+3,4)
    real(8) :: xx(NNUM+3,4)
    real(8) :: work(NNUM+3,4)
    real(8) :: w1,a1,b1,c1,d1,e1
    integer :: m,l,n
    !
    ! =============================================================
    !
!$OMP PARALLEL default(SHARED), private(m,l,n,w1,a1,b1,c1,d1,e1,gg,xx, work) &
!$OMP  &         private(qu_c,qu_s,qv_c,qv_s )
!$OMP DO schedule(STATIC)
      do m=0,mmax
        qu_c(0) = 0.0d0
        qu_s(0) = 0.0d0

        qv_c(0) = 0.0d0
        qv_s(0) = 0.0d0

        do l=1,NNUM
          if ( mod(m,2) == 0 ) then
            qu_c(l) = real(qu(l,m),kind=8)
            qu_s(l) = -dimag(qu(l,m))
            qv_c(l) = real(qv(l,m),kind=8)
            qv_s(l) = -dimag(qv(l,m))
          else
            qu_c(l) = real(qu(l,m),kind=8)
            qu_s(l) = -dimag(qu(l,m))
            qv_c(l) = real(qv(l,m),kind=8)
            qv_s(l) = -dimag(qv(l,m))
          end if
        end do
        
        if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
          l=1
          n=l
          b1 = -16.0d0*m
          c1 = -12.0d0*n
          d1 = 8.0d0*m
          e1 = 4.0d0*(n+2)
          gg(l,1) = + b1*qu_s(l) + c1*qv_c(l+1) &
           &         + d1*qu_s(l+2) + e1*qv_c(l+3)     !!(a) chic,psis,...
          gg(l,2) = - b1*qu_c(l) + c1*qv_s(l+1) &
           &         - d1*qu_c(l+2) + e1*qv_s(l+3)     !!(b) chis,psic,...
          gg(l,3) = - b1*qv_s(l) + c1*qu_c(l+1) &
           &         - d1*qv_s(l+2) + e1*qu_c(l+3)     !!(c) psic,chis,...
          gg(l,4) = - b1*qv_c(l) - c1*qu_s(l+1) &
           &         - d1*qv_c(l+2) - e1*qu_s(l+3)     !!(d) psis,chic,...
          do l=2,NNUM_M_DFS(m)
            n=l
            a1 = 4.0d0*(n-2)
            b1 = -8.0d0*m
            c1 = -8.0d0*n
            d1 = 8.0d0*m
            e1 = 4.0d0*(n+2)
            if ( mod(l,2) == 1 ) then
              gg(l,1) =  a1*qv_c(l-1) + b1*qu_s(l) + c1*qv_c(l+1) &
               &         + d1*qu_s(l+2) + e1*qv_c(l+3)     !!(a) chic,psis,...
              gg(l,2) =  a1*qv_s(l-1) - b1*qu_c(l) + c1*qv_s(l+1) &
               &         - d1*qu_c(l+2) + e1*qv_s(l+3)     !!(b) chis,psic,...
              gg(l,3) =  a1*qu_c(l-1) - b1*qv_s(l) + c1*qu_c(l+1) &
               &         - d1*qv_s(l+2) + e1*qu_c(l+3)     !!(c) psic,chis,...
              gg(l,4) = -a1*qu_s(l-1) - b1*qv_c(l) - c1*qu_s(l+1) &
               &         - d1*qv_c(l+2) - e1*qu_s(l+3)     !!(d) psis,chic,...
            else
              gg(l,1) = -a1*qu_s(l-1) - b1*qv_c(l) - c1*qu_s(l+1) &
               &         - d1*qv_c(l+2) - e1*qu_s(l+3)     !!(d) chic,psis,...
              gg(l,2) =  a1*qu_c(l-1) - b1*qv_s(l) + c1*qu_c(l+1) &
               &         - d1*qv_s(l+2) + e1*qu_c(l+3)     !!(c) chis,psic,...
              gg(l,3) =  a1*qv_s(l-1) - b1*qu_c(l) + c1*qv_s(l+1) &
               &         - d1*qu_c(l+2) + e1*qv_s(l+3)     !!(b) psic,chis,...
              gg(l,4) =  a1*qv_c(l-1) + b1*qu_s(l) + c1*qv_c(l+1) &
               &         + d1*qu_s(l+2) + e1*qv_c(l+3)     !!(a) psis,chic,...
            end if
          end do
        else if ( m == 0 ) then
          l=1
          n=l-1
          gg(l,1) = 0.0d0      !!(a) chic,psis,...
          gg(l,2) = 0.0d0      !!(b) chis,psic,...
          gg(l,3) = 0.0d0      !!(c) psic,chis,...
          gg(l,4) = 0.0d0      !!(d) psis,chic,...
          
          do l=2,NNUM_M_DFS(m)
            n=l-1
            w1=1.0d0/n
            if ( mod(l,2) == 1 ) then
              gg(l,1) = qv_c(l-1)*w1  !!(a) chic,psis,...
              gg(l,2) = 0.0d0      !!(b) chis,psic,...
              gg(l,3) = qu_c(l-1)*w1  !!(c) psic,chis,...
              gg(l,4) = 0.0d0      !!(d) psis,chic,...
            else
              gg(l,1) = 0.0d0      !!(d) chic,psis,...
              gg(l,2) = qu_c(l-1)*w1  !!(c) chis,psic,...
              gg(l,3) = 0.0d0      !!(b) psic,chis,...
              gg(l,4) = qv_c(l-1)*w1  !!(a) psis,chic,...
            end if
          end do

        else if ( mod(m,2) == 1 ) then
          l = 1
          n = l-1
          b1 = -8.0d0*m
          c1 = -4.0d0*(n+1)
          gg(l,1) = + b1*qu_s(l) + c1*qv_c(l+1)  !!(a) chic,psis,...
          gg(l,2) = - b1*qu_c(l) + c1*qv_s(l+1)  !!(b) chis,psic,...
          gg(l,3) = - b1*qv_s(l) + c1*qu_c(l+1)  !!(c) psic,chis,...
          gg(l,4) = - b1*qv_c(l) - c1*qu_s(l+1)  !!(d) psis,chic,...
          do l=2,NNUM_M_DFS(m)
            n = l-1
            a1 = 2.0d0*(n-1)
            b1 = -4.0d0*m
            c1 = -2.0d0*(n+1)
            if ( mod(l,2) == 1 ) then
              gg(l,1) =  a1*qv_c(l-1) + b1*qu_s(l) + c1*qv_c(l+1)  !!(a) chic,psis,...
              gg(l,2) =  a1*qv_s(l-1) - b1*qu_c(l) + c1*qv_s(l+1)  !!(b) chis,psic,...
              gg(l,3) =  a1*qu_c(l-1) - b1*qv_s(l) + c1*qu_c(l+1)  !!(c) psic,chis,...
              gg(l,4) = -a1*qu_s(l-1) - b1*qv_c(l) - c1*qu_s(l+1)  !!(d) psis,chic,...
            else
              gg(l,1) = -a1*qu_s(l-1) - b1*qv_c(l) - c1*qu_s(l+1)  !!(d) chic,psis,...
              gg(l,2) =  a1*qu_c(l-1) - b1*qv_s(l) + c1*qu_c(l+1)  !!(c) chis,psic,...
              gg(l,3) =  a1*qv_s(l-1) - b1*qu_c(l) + c1*qv_s(l+1)  !!(b) psic,chis,...
              gg(l,4) =  a1*qv_c(l-1) + b1*qu_s(l) + c1*qv_c(l+1)  !!(a) psis,chic,...
            end if
          end do

        else   !! mod(m,2) == 0
          do l=1,NNUM_M_DFS(m)
            n=l
            c1 = 2.0d0*(n-1)
            d1 = -4.0d0*m
            e1 = -2.0d0*(n+1)
            if ( mod(l,2) == 1 ) then
              gg(l,1) = + c1*qv_c(l-1) + d1*qu_s(l) + e1*qv_c(l+1)  !!(a) chic,psis,...
              gg(l,2) = + c1*qv_s(l-1) - d1*qu_c(l) + e1*qv_s(l+1)  !!(b) chis,psic,...
              gg(l,3) = + c1*qu_c(l-1) - d1*qv_s(l) + e1*qu_c(l+1)  !!(c) psic,chis,...
              gg(l,4) = - c1*qu_s(l-1) - d1*qv_c(l) - e1*qu_s(l+1)  !!(d) psis,chic,...
            else
              gg(l,1) = - c1*qu_s(l-1) - d1*qv_c(l) - e1*qu_s(l+1)  !!(d) chic,psis,...
              gg(l,2) = + c1*qu_c(l-1) - d1*qv_s(l) + e1*qu_c(l+1)  !!(c) chis,psic,...
              gg(l,3) = + c1*qv_s(l-1) - d1*qu_c(l) + e1*qv_s(l+1)  !!(b) psic,chis,...
              gg(l,4) = + c1*qv_c(l-1) + d1*qu_s(l) + e1*qv_c(l+1)  !!(a) psis,chic,...
            end if
          end do
        end if
              
        if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
          call ninediagonal_4__solve( NNUM_M_DFS(m), NNUM+3, td(m)%d, td(m)%iww, work, gg, xx ) !IN,IN,OOT,IO,OUT
        else if ( m == 0 ) then
          xx(:,:) = gg(:,:)
          xx(NNUM_M_DFS(m)+1:NNUM,:) = 0.0d0
        else
          call pentadiagonal_4__solve( NNUM_M_DFS(m), NNUM+3, td(m)%d, td(m)%iww, work, gg, xx ) !IN,IN,OUT,IO,OUT
        end if

        do l=1,NNUM
          if ( mod(l,2) == 1 ) then
            qchi(l,m) = ER*dcmplx( xx(l,1),-xx(l,2))
            qpsi(l,m) = ER*dcmplx(-xx(l,3),-xx(l,4))
          else
            qchi(l,m) = ER*dcmplx( xx(l,4),-xx(l,3))
            qpsi(l,m) = ER*dcmplx(-xx(l,2),-xx(l,1))
          end if
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    !
    return
  end subroutine uv2chipsi_core



  subroutine uv2chipsi_pre
    integer :: m,l,lmax
    real(8) :: am,an
        
    allocate( td(0:MMAX) )

!$OMP PARALLEL default(SHARED), private(m,l,lmax,am,an)
!$OMP DO schedule(STATIC)
    do m=0,mmax
      am=m
    
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then    
        call ninediagonal_4__allocate( NNUM_M_DFS(m), NNUM+3, td(m)%d, td(m)%iww ) !IN,IN,IO,OUT
      else if ( m == 0 ) then
      else
        call pentadiagonal_4__allocate( NNUM_M_DFS(m), NNUM+3, td(m)%d, td(m)%iww ) !IN,IN,IO,OUT
      end if     
      
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then    
        do l=1,NNUM
          an=l
          td(m)%d(l,1) = (an-2)**2
          td(m)%d(l,2) = 2.0d0*am
          td(m)%d(l,3) = -4.0d0*am**2 -4.0d0*an**2 + 8.0d0*an - 8.0d0
          td(m)%d(l,4) = -2.0d0*am
          td(m)%d(l,5) = 8.0d0*am**2 + 6.0d0*an**2 + 8.0d0
          td(m)%d(l,6) = -2.0d0*am
          td(m)%d(l,7) = -4.0d0*am**2 - 4.0d0*an**2 - 8.0d0*an - 8.0d0
          td(m)%d(l,8) = 2.0d0*am
          td(m)%d(l,9) = (an+2)**2
        end do
        
        td(m)%d(1,1) = -999.0d33
        td(m)%d(1,2) = -999.0d33
        td(m)%d(1,3) = -999.0d33
        td(m)%d(1,4) = -999.0d33
        td(m)%d(1,5) = 12.0d0*am**2 + 18.0d0  !! n=1
        td(m)%d(1,6) = -4.0d0*am              !! n=1
        td(m)%d(1,7) = -4.0d0*am**2 - 21.0d0  !! n=1
        
        td(m)%d(2,1) = -999.0d33
        td(m)%d(2,2) = -999.0d33
        td(m)%d(2,3) = -999.0d33
        td(m)%d(2,4) = -4.0d0*am              !! n=2
        
        td(m)%d(3,1) = -999.0d33
        td(m)%d(3,2) = -999.0d33
        td(m)%d(3,3) = -4.0d0*am**2 - 21.0d0  !! n=3
        
        td(m)%d(4,1) = -999.0d33
      else if ( m == 0 ) then
      else if ( mod(m,2) == 1 ) then        !! if wavenumber is 1
        do l=1,NNUM
          an=l-1
          td(m)%d(l,1) = -1.0d0*(an-1)**2
          td(m)%d(l,2) = -2.0d0*am
          td(m)%d(l,3) = 4.0d0*am**2 + 2.0d0*an**2 + 2.0d0
          td(m)%d(l,4) = -2.0d0*am
          td(m)%d(l,5) = -1.0d0*(an+1)**2
        end do
        td(m)%d(1,1) = -999.0d33
        td(m)%d(1,2) = -999.0d33
        td(m)%d(1,3) = 8.0d0*am**2 + 4.0d0   !! n=0
        td(m)%d(1,4) = -4.0d0*am             !! n=0
        td(m)%d(1,5) = -2.0d0                !! n=0
          
        td(m)%d(2,1) = -999.0d33
        td(m)%d(2,2) = -4.0d0*am             !! n=1
          
        td(m)%d(3,1) = -2.0d0                !! n=2
      else    !! mod(m,2) == 0
        do l=1,NNUM
          an=l
          td(m)%d(l,1) = -1.0d0*(an-1)**2
          td(m)%d(l,2) = -2.0d0*am
          td(m)%d(l,3) = 4.0d0*am**2 + 2.0d0*an**2 + 2.0d0
          td(m)%d(l,4) = -2.0d0*am
          td(m)%d(l,5) = -1.0d0*(an+1)**2
        end do
      end if

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then  
        call ninediagonal_4__ini( NNUM_M_DFS(m), NNUM+3, td(m)%d, td(m)%iww ) !IN,IN,IO,OUT
      else if ( m == 0 ) then
      else
        call pentadiagonal_4__ini( NNUM_M_DFS(m), NNUM+3, td(m)%d, td(m)%iww ) !IN,IN,IO,OUT
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL
    !
  end subroutine uv2chipsi_pre
  

  subroutine uv2divrot_dfs__laplacian( qvar ) !IN,OUT
    !
    complex(8),intent(INOUT) :: qvar(NNUM,0:MMAX)
    !
    real(8),save,allocatable :: a1(:,:)
    real(8),save,allocatable :: b1(:,:)
    real(8),save,allocatable :: c1(:,:)
    real(8),save,allocatable :: d1(:,:)
    real(8),save,allocatable :: e1(:,:)
    !
    real(8),save,allocatable :: a2(:,:)
    real(8),save,allocatable :: b2(:,:)
    real(8),save,allocatable :: c2(:,:)
    real(8),save,allocatable :: d2(:,:)
    real(8),save,allocatable :: e2(:,:)
    !
    logical,save :: first_laplacian = .true.
    !
    if ( jcn_dfs >= 1 ) then
      !
      if ( first_laplacian ) then
        !
        first_laplacian = .false.
        !
        allocate( a1(NNUMHF1+1,0:MMAX) )
        allocate( b1(NNUMHF1+1,0:MMAX) )
        allocate( c1(NNUMHF1+1,0:MMAX) )
        allocate( d1(NNUMHF1+1,0:MMAX) )
        allocate( e1(NNUMHF1+1,0:MMAX) )
        
        allocate( a2(NNUMHF2+1,0:MMAX) )
        allocate( b2(NNUMHF2+1,0:MMAX) )
        allocate( c2(NNUMHF2+1,0:MMAX) )
        allocate( d2(NNUMHF2+1,0:MMAX) )
        allocate( e2(NNUMHF2+1,0:MMAX) )
        
        call calc_abcde_lap         &
         &( a1, b1, c1, d1, e1,     &!OUT
         &  a2, b2, c2, d2, e2    )  !OUT
        !
      end if
    
      call laplacian_core          &
       &( a1, b1, c1, d1, e1,      &!IN
       &  a2, b2, c2, d2, e2,      &!IN
       &  qvar     )                !INOUT

    else
      call uv2divrot_dfs_old__laplacian( qvar ) !IN,OUT
    end if
    !
  end subroutine uv2divrot_dfs__laplacian



  subroutine laplacian_core    &!
   &( a1, b1, c1, d1, e1,      &!IN
   &  a2, b2, c2, d2, e2,      &!IN
   &  qvar     )                !INOUT
    !
    real(8),intent(in) :: a1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: b1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: c1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: d1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: e1(NNUMHF1+1,0:mmax)
    real(8),intent(in) :: a2(NNUMHF2+1,0:mmax)
    real(8),intent(in) :: b2(NNUMHF2+1,0:mmax)
    real(8),intent(in) :: c2(NNUMHF2+1,0:mmax)
    real(8),intent(in) :: d2(NNUMHF2+1,0:mmax)
    real(8),intent(in) :: e2(NNUMHF2+1,0:mmax)
    complex(8),intent(inout) :: qvar(NNUM,0:MMAX)
    !
    complex(8) :: g1(NNUMHF1+1)
    complex(8) :: g2(NNUMHF2+1)
    complex(8) :: x1(NNUMHF1+1)
    complex(8) :: x2(NNUMHF2+1)
    complex(8) :: work(max(NNUMHF1,NNUMHF2)+1)
    !
    integer :: m,mm,l,ll,l2
    !
    ! =============================================================

!$OMP PARALLEL default(SHARED), private(m,mm,l,ll,l2,g1,g2,x1,x2,work)
!$OMP DO schedule(STATIC)
    do m=0,mmax
        l2= 1
        l = 1
        g1(l2) = qvar(l,m)*c1(l2,m) + qvar(l+2,m)*d1(l2,m) + qvar(l+4,m)*e1(l2,m)
        ll = 2
        g2(l2) = qvar(ll,m)*c2(l2,m) + qvar(ll+2,m)*d2(l2,m) + qvar(ll+4,m)*e2(l2,m)

        l2= 2
        l = 3
        g1(l2) = qvar(l-2,m)*b1(l2,m) + qvar(l,m)*c1(l2,m) + qvar(l+2,m)*d1(l2,m) + qvar(l+4,m)*e1(l2,m)
        if ( NNUMHF2 >= 4 ) then
          ll = 4
          g2(l2) = qvar(ll-2,m)*b2(l2,m) + qvar(ll,m)*c2(l2,m) + qvar(ll+2,m)*d2(l2,m) + qvar(ll+4,m)*e2(l2,m)
        end if

        do l2=3,NNUMHF1-2
          l=l2*2-1
          g1(l2) = qvar(l-4,m)*a1(l2,m) + qvar(l-2,m)*b1(l2,m) + qvar(l,m)*c1(l2,m) &
           &       + qvar(l+2,m)*d1(l2,m) + qvar(l+4,m)*e1(l2,m)
        end do

        do l2=3,NNUMHF2-2
          ll=l2*2
          g2(l2) = qvar(ll-4,m)*a2(l2,m) + qvar(ll-2,m)*b2(l2,m) + qvar(ll,m)*c2(l2,m) &
           &       + qvar(ll+2,m)*d2(l2,m) + qvar(ll+4,m)*e2(l2,m)
        end do

        l2=NNUMHF1-1
        l=l2*2-1
        g1(l2) = qvar(l-4,m)*a1(l2,m) + qvar(l-2,m)*b1(l2,m) + qvar(l,m)*c1(l2,m) + qvar(l+2,m)*d1(l2,m)

        l2=NNUMHF1
        l=l2*2-1
        g1(l2) = qvar(l-4,m)*a1(l2,m) + qvar(l-2,m)*b1(l2,m) + qvar(l,m)*c1(l2,m)

        l2=NNUMHF2-1
        ll=l2*2
        if ( NNUMHF2 >= 4 ) then
          g2(l2) = qvar(ll-4,m)*a2(l2,m) + qvar(ll-2,m)*b2(l2,m) + qvar(ll,m)*c2(l2,m) + qvar(ll+2,m)*d2(l2,m)
        else
          g2(l2) = qvar(ll-2,m)*b2(l2,m) + qvar(ll,m)*c2(l2,m) + qvar(ll+2,m)*d2(l2,m)  
        end if
        l2=NNUMHF2
        ll=l2*2
        g2(l2) = qvar(ll-4,m)*a2(l2,m) + qvar(ll-2,m)*b2(l2,m) + qvar(ll,m)*c2(l2,m)

!        if ( m == 4 ) then
!          write(6,*) "m=",m
!          write(6,*) "qvar(:,m)=",qvar(:,m)
!          write(6,*) "g1(:)=",g1(:)
!          write(6,*) "g2(:)=",g2(:)
!        end if

      !! Calculate x1,x2 from g1,g2    
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        mm = 3
        call pentadiagonal__solve( NNUM1_M_DFS(mm), NNUMHF1+1, td1(mm)%d, td1(mm)%iww, &
         &                         work, g1, x1 )
        call pentadiagonal__solve( NNUM2_M_DFS(mm), NNUMHF2+1, td2(mm)%d, td2(mm)%iww, &
         &                         work, g2, x2 )
      else
        if ( m == 0 ) then
          mm = 0
        else if ( mod(m,2) == 1 ) then
          mm = 1
        else
          mm = 2
        end if
        call tridiagonal__solve( NNUM1_M_DFS(mm), NNUMHF1, td1(mm)%d, td1(mm)%iww, &
         &                       work, g1, x1 )
        call tridiagonal__solve( NNUM2_M_DFS(mm), NNUMHF2, td2(mm)%d, td2(mm)%iww, &
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
    end do     
!$OMP END DO
!$OMP END PARALLEL

    !
  end subroutine laplacian_core
  

  subroutine calc_abcde_lap   &!
   &( a1, b1, c1, d1, e1,     &!OUT
   &  a2, b2, c2, d2, e2    )  !OUT
    !
    real(8),intent(out) :: a1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: b1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: c1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: d1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: e1(NNUMHF1+1,0:mmax)
    real(8),intent(out) :: a2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: b2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: c2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: d2(NNUMHF2+1,0:mmax)
    real(8),intent(out) :: e2(NNUMHF2+1,0:mmax)
    !
    integer :: m,l2,nn
    real(8) :: am,al
        
    allocate( td1(0:3) )
    allocate( td2(0:3) )
    
  !$OMP PARALLEL default(SHARED), private(m,l2,nn)
   !$OMP DO schedule(STATIC)
    do m=0,3
      if ( m == 0 .or. m == 1 .or. m == 2 ) then
        call tridiagonal__allocate( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call tridiagonal__allocate( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      else if ( jcn_dfs >= 2 .and. m == 3 ) then
        call pentadiagonal__allocate( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call pentadiagonal__allocate( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      end if

      if ( m == 0 .or. m == 1 ) then
        do l2=1,NNUM1_M_DFS(m)
          td1(m)%d(l2,1) = -1.0d0
          td1(m)%d(l2,2) = 2.0d0
          td1(m)%d(l2,3) = -1.0d0
        end do
        td1(m)%d(2,1) = -2.0d0    !! n=2
        
        do l2=1,NNUM2_M_DFS(m)
          td2(m)%d(l2,1) = -1.0d0
          td2(m)%d(l2,2) = 2.0d0
          td2(m)%d(l2,3) = -1.0d0
        end do
        td2(m)%d(1,1) = 999.0d33
        td2(m)%d(1,2) = 1.0d0     !! n=1
      else if ( m == 2 ) then
        do l2=1,NNUM1_M_DFS(m)
          td1(m)%d(l2,1) = -1.0d0
          td1(m)%d(l2,2) = 2.0d0
          td1(m)%d(l2,3) = -1.0d0
        end do
        td1(m)%d(1,2) = 3.0d0
        
        do l2=1,NNUM2_M_DFS(m)
          td2(m)%d(l2,1) = -1.0d0
          td2(m)%d(l2,2) = 2.0d0
          td2(m)%d(l2,3) = -1.0d0
        end do
      else if ( jcn_dfs >= 2 .and. m == 3 ) then
        do l2=1,NNUM1_M_DFS(m)
          td1(m)%d(l2,1) = 1.0d0
          td1(m)%d(l2,2) = -4.0d0
          td1(m)%d(l2,3) = 6.0d0
          td1(m)%d(l2,4) = -4.0d0
          td1(m)%d(l2,5) = 1.0d0
        end do
        td1(m)%d(1,1) = -999.0d33
        td1(m)%d(1,2) = -999.0d33
        td1(m)%d(1,3) = 10.0d0     !! n=1
        td1(m)%d(1,4) = -5.0d0     !! n=1
        td1(m)%d(2,1) = -999.0d33
        td1(m)%d(2,2) = -5.0d0     !! n=3
        
        do l2=1,NNUM2_M_DFS(m)
          td2(m)%d(l2,1) = 1.0d0
          td2(m)%d(l2,2) = -4.0d0
          td2(m)%d(l2,3) = 6.0d0
          td2(m)%d(l2,4) = -4.0d0
          td2(m)%d(l2,5) = 1.0d0
        end do
        td2(m)%d(1,1) = -999.0d33
        td2(m)%d(1,2) = -999.0d33
        td2(m)%d(1,3) = 5.0d0      !! n=2
      end if
      
      if ( m == 0 .or. m == 1 .or. m == 2 ) then
        call tridiagonal__ini( NNUM1_M_DFS(m), NNUMHF1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call tridiagonal__ini( NNUM2_M_DFS(m), NNUMHF2, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      else if ( jcn_dfs >= 2 .and. m == 3 ) then
        call pentadiagonal__ini( NNUM1_M_DFS(m), NNUMHF1+1, td1(m)%d, td1(m)%iww ) !IN,IN,IO,OUT
        call pentadiagonal__ini( NNUM2_M_DFS(m), NNUMHF2+1, td2(m)%d, td2(m)%iww ) !IN,IN,IO,OUT
      end if   
    end do
   !$OMP END DO
  !$OMP END PARALLEL


  !$OMP PARALLEL default(SHARED), private(m,l2,am,al)
   !$OMP DO schedule(STATIC)
    do m=0,mmax
      am=m
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        do l2=1,NNUMHF1+1
          al=l2*2-1
          a1(l2,m) = -(al-2)*(al-1)/(ER**2)
          b1(l2,m) = ( 4*al*al -6*al +4 +4*am*am )/(ER**2)
          c1(l2,m) = ( -6*al*al -4 -8*am*am )/(ER**2)
          d1(l2,m) = ( 4*al*al +6*al +4 +4*am*am )/(ER**2)
          e1(l2,m) = -(al+2)*(al+1)/(ER**2)
        end do
        l2=1
        al=1
        a1(l2,m) = -999.0d33
        b1(l2,m) = -999.0d33
        c1(l2,m) = ( -12 - 12*am*am )/(ER**2)   !! n=1
        l2=2
        al=3
        a1(l2,m) = -999.0d33
        b1(l2,m) = ( 24 +4*am*am )/(ER**2)      !! n=3
        
        do l2=1,NNUMHF2+1
          al=l2*2
          a2(l2,m) = -(al-2)*(al-1)/(ER**2)
          b2(l2,m) = ( 4*al*al -6*al +4 +4*am*am )/(ER**2)
          c2(l2,m) = ( -6*al*al -4 -8*am*am )/(ER**2)
          d2(l2,m) = ( 4*al*al +6*al +4 +4*am*am )/(ER**2)
          e2(l2,m) = -(al+2)*(al+1)/(ER**2)
        end do
      else if ( m == 0 ) then
        do l2=1,NNUMHF1+1
          al=l2*2-2                   !! l = 0,2,4,6,...,JMAX-2
          a1(l2,m) = 0.0d0
          b1(l2,m) = (al-1)*(al-2)/(ER**2)
          c1(l2,m) = -2*al*al/(ER**2)
          d1(l2,m) = (al+1)*(al+2)/(ER**2)
          e1(l2,m) = 0.0d0
        end do
        
        do l2=1,NNUMHF2+1
          al=l2*2-1                      !! l = 1,3,5,...,JMAX-1
          a2(l2,m) = 0.0d0
          b2(l2,m) = (al-1)*(al-2)/(ER**2)
          c2(l2,m) = -2*al*al/(ER**2)
          d2(l2,m) = (al+1)*(al+2)/(ER**2)
          e2(l2,m) = 0.0d0
        end do
      else if ( mod(m,2) == 1 ) then         
        do l2=1,NNUMHF1+1
          al=l2*2-2                 !! l = 1,3,5,...,JMAX-1
          a1(l2,m) = 0.0d0
          b1(l2,m) = (al-1)*al/(ER**2)
          c1(l2,m) = ( -2*al*al -4*am*am )/(ER**2)
          d1(l2,m) = (al+1)*al/(ER**2)
          e1(l2,m) = 0.0d0
        end do
        l2=2
        al=2
        b1(l2,m) = 2*(al-1)*al/(ER**2)    !! n=2
        
        do l2=1,NNUMHF2+1
          al=l2*2-1                 !! l = 2,4,6,...,JMAX
          a2(l2,m) = 0.0d0
          b2(l2,m) = (al-1)*al/(ER**2)
          c2(l2,m) = ( -2*al*al -4*am*am )/(ER**2)
          d2(l2,m) = (al+1)*al/(ER**2)
          e2(l2,m) = 0.0d0
        end do
      else
        do l2=1,NNUMHF1+1
          al=l2*2-1                !! l = 1,3,5,...,JMAX-1
          a1(l2,m) = 0.0d0
          b1(l2,m) = (al-1)*al/(ER**2)
          c1(l2,m) = ( -2*al*al -4*am*am )/(ER**2)
          d1(l2,m) = (al+1)*al/(ER**2)
          e1(l2,m) = 0.0d0
        end do
        
        do l2=1,NNUMHF2+1
          al=l2*2                   !! l = 2,4,6,...,JMAX
          a2(l2,m) = 0.0d0
          b2(l2,m) = (al-1)*al/(ER**2)
          c2(l2,m) = ( -2*al*al -4*am*am )/(ER**2)
          d2(l2,m) = (al+1)*al/(ER**2)
          e2(l2,m) = 0.0d0
        end do
      end if
    end do
   !$OMP END DO
  !$OMP END PARALLEL

    !
  end subroutine calc_abcde_lap

end module uv2divrot_dfs
