module ninediagonal

  implicit none
  
  private
  public :: ninediagonal__allocate, ninediagonal__ini, ninediagonal__solve, ninediagonal__test
  
!  integer,parameter :: jcn_cyclic = 0  !! LU decomposition of nine-diagonal matrix
  integer,parameter :: jcn_cyclic = 1  !! LU decomposition of block tri-diagonal matrix
!  integer,parameter :: jcn_cyclic = 2  !! Cyclic reduction for block tri-diagonal matrix

!  call main_sub

contains


  subroutine ninediagonal__test

    integer,parameter :: imax=1000
    
    real(8),allocatable :: dat(:,:)
    integer,allocatable :: iww(:)
    complex(8) :: yy(imax+3)
    complex(8) :: yvar(imax+3)
    complex(8) :: xvar(imax+3)
    integer :: n   

    call ninediagonal__allocate &
     &( imax, imax+3, dat, iww) !IN,IN,IO,OUT
    
    dat(:,1) = 2.0d0
    dat(:,2) = 1.0d0
    dat(:,3) = 2.0d0
    dat(:,4) = 1.0d0
    dat(:,5) = 10.0d0
    dat(:,6) = 1.0d0
    dat(:,7) = 2.0d0
    dat(:,8) = 1.0d0
    dat(:,9) = 2.0d0

    call ninediagonal__ini     &
     &( imax, imax+3, dat, iww) !IN,IN,IO,OUT

    n=1
!    do n=1,100000

      yvar(1) = 16.0d0
      yvar(2) = 17.0d0
      yvar(3) = 19.0d0
      yvar(4) = 20.0d0
      yvar(5:imax-4) = 22.0d0
      yvar(imax-3)   = 20.0d0
      yvar(imax-2)   = 19.0d0
      yvar(imax-1)   = 17.0d0
      yvar(imax)     = 16.0d0
      
      call ninediagonal__solve    &
       &( imax, imax+3, dat, iww, &!IN
       &  yy, yvar, xvar )         !OUT,IO,OUT
      
!    end do
  
    write(6,*) "xvar(1:imax)=",xvar(1:imax)
  
  end subroutine ninediagonal__test
 
  
  subroutine ninediagonal__allocate( imax, ilen, dat, iww )
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
!xx    real(8),intent(out),pointer :: dat(:,:)
!xx    integer,intent(out),pointer :: iww(:)
    real(8),intent(out),allocatable :: dat(:,:)
    integer,intent(out),allocatable :: iww(:)
    
    integer :: nn_tmp, nn_cyclic

    if ( ilen < imax+3 ) then
      write(6,*) "Error: ninediagonal__allocate: ilen should be equal to or larger than imax+3."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
    
    if ( jcn_cyclic == 0 ) then
      allocate( dat(ilen,9) )
      allocate( iww(1) )
    
    else if ( jcn_cyclic == 1 ) then
      allocate( dat(ilen,12) )
      allocate( iww(2) )
      nn_cyclic = 0
      iww(1) = nn_cyclic
    
    else if ( jcn_cyclic >= 2 ) then
      nn_tmp = int(log(imax+0.001d0)/log(2.0d0))
      nn_cyclic = max( nn_tmp-2, 0 )
      allocate( dat(ilen,20) )
      allocate( iww(2+2*nn_cyclic) )
      
      iww(1) = nn_cyclic
      write(6,*) "nn_cyclic=",nn_cyclic
    end if
  
  end subroutine ninediagonal__allocate
 
 
  subroutine ninediagonal__ini                               &
   &( imax, ilen, dat, iww)

    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    real(8),intent(inout) :: dat(ilen,*)
    integer,intent(out)   :: iww(*)
  
    integer :: nn_tmp, nn_cyclic

    if ( ilen < imax+3 ) then
      write(6,*) "Error: ninediagonal__ini: ilen should be equal to or larger than imax+3."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
    
    if ( jcn_cyclic == 0 ) then
      call ninediagonal_simple_pre( imax,                                                         &!IN
       & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6), dat(1,7), dat(1,8), dat(1,9) )!IO
    
    else if ( jcn_cyclic == 1 ) then
      nn_cyclic = iww(1)
      call ninediagonal_cyclic_pre( nn_cyclic, imax,                                               &!IN
       & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6), dat(1,7), dat(1,8), dat(1,9), &!IO
       & dat(1,10), dat(1,11), dat(1,12), iww(2) )                                                  !OUT
       
    else !! jcn_cyclic >= 2
      nn_cyclic = iww(1)
      call ninediagonal_cyclic_pre( nn_cyclic, imax,                                               &!IN
       & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6), dat(1,7), dat(1,8), dat(1,9), &!IO
       & dat(1,10), dat(1,11), dat(1,12), iww(2),                                                  &!OUT
       & dat(1,13), dat(1,14), dat(1,15), dat(1,16), dat(1,17), dat(1,18), dat(1,19), dat(1,20),   &!OUT
       & iww(3+nn_cyclic) )                                                                         !OUT
    end if
    
  end subroutine ninediagonal__ini
   
  
  subroutine ninediagonal__solve             &
   &( imax, ilen, dat, iww,                  &
   &  yy, yvar, xvar )
   
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    real(8),intent(in) :: dat(ilen,20)
    integer,intent(in) :: iww(ilen)
    complex(8),intent(inout) :: yy(ilen)
    complex(8),intent(inout) :: yvar(ilen)
    complex(8),intent(inout) :: xvar(ilen)
      
    integer :: nn_cyclic, jmax2

    if ( ilen < imax+3 ) then
      write(6,*) "Error: ninediagonal__solve: ilen should be equal to or larger than imax+3."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
      
    yvar(imax+1:ilen) = 0.0d0
    if ( jcn_cyclic == 0 ) then
      call ninediagonal_simple_solve( imax,                                                        &
       & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6), dat(1,7), dat(1,8), dat(1,9), &!IN
       & yvar, xvar )                                                                               !IO,OUT
    else if ( jcn_cyclic == 1 ) then
      nn_cyclic = iww(1)
      jmax2 = iww(2)
      call ninediagonal_block_solve                                                                      &
       &( jmax2, dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6), dat(1,7), dat(1,8), dat(1,9), &!IN
       &  dat(1,10), dat(1,11), dat(1,12), yvar,                                                           &!IN
       &  xvar   )                                                                                          !OUT

    else      !! jcn_cyclic >= 2
       nn_cyclic = iww(1)
       call ninediagonal_cyclic_solve( nn_cyclic, imax, &
        & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6), dat(1,7), dat(1,8), dat(1,9), &!IN
        & dat(1,10), dat(1,11), dat(1,12), dat(1,13), dat(1,14), dat(1,15), dat(1,16), dat(1,17),   &!IN
        & dat(1,18), dat(1,19), dat(1,20), iww(2), iww(3+nn_cyclic),                                &!IN
        & yy, yvar, xvar )                                                                   !OUT,IO,OUT
    end if
    xvar(imax+1:ilen) = 0.0d0
    
  end subroutine ninediagonal__solve
  
  
  subroutine inverse_4x4( aa, bb, cc, dd )
     real(8),intent(inout) :: aa(4)
     real(8),intent(inout) :: bb(4)
     real(8),intent(inout) :: cc(4)
     real(8),intent(inout) :: dd(4)
     real(8) :: aa2(4)
     real(8) :: bb2(4)
     real(8) :: cc2(4)
     real(8) :: dd2(4)
     real(8) :: ainv,aaa
     integer :: jj,j
  
     aa2(:) = 0.0d0
     bb2(:) = 0.0d0
     cc2(:) = 0.0d0
     dd2(:) = 0.0d0
     aa2(1) = 1.0d0
     bb2(2) = 1.0d0
     cc2(3) = 1.0d0
     dd2(4) = 1.0d0
  
!     aa(:) = 1.0d0
!     bb(:) = 2.0d0
!     cc(:) = 3.0d0
!     dd(:) = 4.0d0
!     aa(1) = 10.0d0
!     bb(2) = 10.0d0
!     cc(3) = 10.0d0
!     dd(4) = 10.0d0
     
!     write(6,*) "0:aa(:)=",aa(:)
!     write(6,*) "0:bb(:)=",bb(:)
!     write(6,*) "0:cc(:)=",cc(:)
!     write(6,*) "0:dd(:)=",dd(:)    
     
     
!     write(6,*) "1:aa(:)=",aa(:)
!     write(6,*) "1:bb(:)=",bb(:)
!     write(6,*) "1:cc(:)=",cc(:)
!     write(6,*) "1:dd(:)=",dd(:)  
     
     do jj=1,4
        if ( jj == 1 ) then
           ainv = 1.0d0/aa(jj)
        else if ( jj == 2 ) then
           ainv = 1.0d0/bb(jj)
        else if ( jj == 3 ) then
           ainv = 1.0d0/cc(jj)
        else
           ainv = 1.0d0/dd(jj)
        end if
        aa(jj)  = aa(jj)*ainv
        bb(jj)  = bb(jj)*ainv
        cc(jj)  = cc(jj)*ainv
        dd(jj)  = dd(jj)*ainv
        aa2(jj) = aa2(jj)*ainv
        bb2(jj) = bb2(jj)*ainv
        cc2(jj) = cc2(jj)*ainv
        dd2(jj) = dd2(jj)*ainv
        do j=jj+1,4
           if ( jj == 1 ) then
              aaa = aa(j)
           else if ( jj == 2 ) then
              aaa = bb(j)
           else
              aaa = cc(j)
           end if
           aa(j) = aa(j) - aaa*aa(jj)
           bb(j) = bb(j) - aaa*bb(jj)
           cc(j) = cc(j) - aaa*cc(jj)
           dd(j) = dd(j) - aaa*dd(jj)
           aa2(j) = aa2(j) - aaa*aa2(jj)
           bb2(j) = bb2(j) - aaa*bb2(jj)
           cc2(j) = cc2(j) - aaa*cc2(jj)
           dd2(j) = dd2(j) - aaa*dd2(jj)
        end do
     end do  
     
     do jj=4,2,-1
        do j=jj-1,1,-1
           if ( jj == 4 ) then
              aaa = dd(j)
           else if ( jj == 3 ) then
              aaa = cc(j)
           else
              aaa = bb(j)
           end if
           aa(j) = aa(j) - aaa*aa(jj)
           bb(j) = bb(j) - aaa*bb(jj)
           cc(j) = cc(j) - aaa*cc(jj)
           dd(j) = dd(j) - aaa*dd(jj)
           aa2(j) = aa2(j) - aaa*aa2(jj)
           bb2(j) = bb2(j) - aaa*bb2(jj)
           cc2(j) = cc2(j) - aaa*cc2(jj)
           dd2(j) = dd2(j) - aaa*dd2(jj)
        end do
     end do
     
!     write(6,*) "aa(:)=",aa(:)
!     write(6,*) "bb(:)=",bb(:)
!     write(6,*) "cc(:)=",cc(:)
!     write(6,*) "dd(:)=",dd(:)
     
!     write(6,*) "aa2(:)=",aa2(:)
!     write(6,*) "bb2(:)=",bb2(:)
!     write(6,*) "cc2(:)=",cc2(:)
!     write(6,*) "dd2(:)=",dd2(:)
     
     aa(:) = aa2(:)
     bb(:) = bb2(:)
     cc(:) = cc2(:)
     dd(:) = dd2(:)
     
  end subroutine inverse_4x4
  

  subroutine multiply_4x4( aa, bb, cc, dd, ee, ff, gg, hh, &
   &                       pp, qq, rr, ss )
    real(8),intent(in) :: aa(4)
    real(8),intent(in) :: bb(4)
    real(8),intent(in) :: cc(4)
    real(8),intent(in) :: dd(4)
    real(8),intent(in) :: ee(4)
    real(8),intent(in) :: ff(4)
    real(8),intent(in) :: gg(4)
    real(8),intent(in) :: hh(4)
    real(8),intent(out) :: pp(4)
    real(8),intent(out) :: qq(4)
    real(8),intent(out) :: rr(4)
    real(8),intent(out) :: ss(4)
    integer :: j
  
    !! [p q r s;p q r s;p q r s;p q r s] 
    !! = [a b c d;a b c d;a b c d;a b c d]*[e f g h;e f g h;e f g h;e f g h]
    do j=1,4
       pp(j) = aa(j)*ee(1) + bb(j)*ee(2) + cc(j)*ee(3) + dd(j)*ee(4)
       qq(j) = aa(j)*ff(1) + bb(j)*ff(2) + cc(j)*ff(3) + dd(j)*ff(4)
       rr(j) = aa(j)*gg(1) + bb(j)*gg(2) + cc(j)*gg(3) + dd(j)*gg(4)
       ss(j) = aa(j)*hh(1) + bb(j)*hh(2) + cc(j)*hh(3) + dd(j)*hh(4)
    end do
     
  end subroutine multiply_4x4


  subroutine ninediagonal_simple_pre( jmax, aa, bb, cc, dd, ee, ff, gg, hh, ii )
    
    integer,intent(in) :: jmax
    real(8),intent(inout) :: aa(jmax)
    real(8),intent(inout) :: bb(jmax)
    real(8),intent(inout) :: cc(jmax)
    real(8),intent(inout) :: dd(jmax)
    real(8),intent(inout) :: ee(jmax)
    real(8),intent(inout) :: ff(jmax)
    real(8),intent(inout) :: gg(jmax)
    real(8),intent(inout) :: hh(jmax)
    real(8),intent(inout) :: ii(jmax)
    integer :: j
    
!    aa(:) = 0.0d0
!    bb(:) = 0.0d0
!    cc(:) = 0.0d0
!    dd(:) = 0.0d0
!    ee(:) = 0.0d0
!    ff(:) = 0.0d0
!    gg(:) = 0.0d0
!    hh(:) = 0.0d0
!    ii(:) = 0.0d0

    ee(1) = 1.0d0/ee(1)
    ff(1) = ff(1)*ee(1)
    gg(1) = gg(1)*ee(1)
    hh(1) = hh(1)*ee(1)
    ii(1) = ii(1)*ee(1)

    ee(2) = 1.0d0/(ee(2) - dd(2)*ff(1))
    ff(2) = ( ff(2) - dd(2)*gg(1) )*ee(2)
    gg(2) = ( gg(2) - dd(2)*hh(1) )*ee(2)
    hh(2) = ( hh(2) - dd(2)*ii(1) )*ee(2)
    ii(2) = ii(2)*ee(2)

    dd(3) = dd(3) - cc(3)*ff(1)
    ee(3) = 1.0d0/( ee(3) - cc(3)*gg(1) - dd(3)*ff(2) )
    ff(3) = ( ff(3) - cc(3)*hh(1) - dd(3)*gg(2) )*ee(3)
    gg(3) = ( gg(3) - cc(3)*ii(1) - dd(3)*hh(2) )*ee(3)
    hh(3) = ( hh(3)               - dd(3)*ii(2) )*ee(3)
    ii(3) = ( ii(3)                             )*ee(3)

    cc(4) = cc(4) - bb(4)*ff(1)
    dd(4) = dd(4) - bb(4)*gg(1) - cc(4)*ff(2)
    ee(4) = 1.0d0/( ee(4) - bb(4)*hh(1) - cc(4)*gg(2) - dd(4)*ff(3) )
    ff(4) = ( ff(4) - bb(4)*ii(1) - cc(4)*hh(2) - dd(4)*gg(3) )*ee(4)
    gg(4) = ( gg(4)               - cc(4)*ii(2) - dd(4)*hh(3) )*ee(4)
    hh(4) = ( hh(4)                             - dd(4)*ii(3) )*ee(4)
    ii(4) = ( ii(4)                                           )*ee(4)

    do j=5,jmax
       bb(j) = bb(j) - aa(j)*ff(j-4)
       cc(j) = cc(j) - aa(j)*gg(j-4) - bb(j)*ff(j-3)
       dd(j) = dd(j) - aa(j)*hh(j-4) - bb(j)*gg(j-3) - cc(j)*ff(j-2)
       ee(j) = 1.0d0/( ee(j) - aa(j)*ii(j-4) - bb(j)*hh(j-3) - cc(j)*gg(j-2) - dd(j)*ff(j-1) )
       ff(j) = ( ff(j) - bb(j)*ii(j-3) - cc(j)*hh(j-2) - dd(j)*gg(j-1) )*ee(j)
       gg(j) = ( gg(j)                 - cc(j)*ii(j-2) - dd(j)*hh(j-1) )*ee(j)
       hh(j) = ( hh(j)                                 - dd(j)*ii(j-1) )*ee(j)
       ii(j) = ( ii(j)                                                 )*ee(j)
       aa(j) = aa(j)*ee(j)
       bb(j) = bb(j)*ee(j)
       cc(j) = cc(j)*ee(j)
       dd(j) = dd(j)*ee(j)
    end do

  end subroutine ninediagonal_simple_pre
  
  
  
  subroutine ninediagonal_simple_solve( jmax, aa, bb, cc, dd, ee, ff, gg, hh, ii, yvar, xvar )
    integer,intent(in) :: jmax
    real(8),intent(in) :: aa(jmax)
    real(8),intent(in) :: bb(jmax)
    real(8),intent(in) :: cc(jmax)
    real(8),intent(in) :: dd(jmax)
    real(8),intent(in) :: ee(jmax)
    real(8),intent(in) :: ff(jmax)
    real(8),intent(in) :: gg(jmax)
    real(8),intent(in) :: hh(jmax)
    real(8),intent(in) :: ii(jmax)
    complex(8),intent(in) :: yvar(jmax)
    complex(8),intent(out) :: xvar(jmax)
    integer :: j
    
    !! ---  Forward substitution  ---
    xvar(1) = yvar(1)*ee(1)

    xvar(2) = ( yvar(2) - dd(2)*xvar(1) )*ee(2)
    
    xvar(3) = ( yvar(3) - dd(3)*xvar(2) - cc(3)*xvar(1) )*ee(3)

    xvar(4) = ( yvar(4) - dd(4)*xvar(3) - cc(4)*xvar(2) - bb(4)*xvar(1) )*ee(4)
    do j=5,jmax
       xvar(j) = yvar(j)*ee(j) - dd(j)*xvar(j-1) - cc(j)*xvar(j-2) - bb(j)*xvar(j-3) - aa(j)*xvar(j-4)
    end do

    !! ---  Backward substitution  ---
    j=jmax-1
    xvar(j) = xvar(j) - ff(j)*xvar(j+1)
    j=jmax-2
    xvar(j) = xvar(j) - ff(j)*xvar(j+1) - gg(j)*xvar(j+2)
    j=jmax-3
    xvar(j) = xvar(j) - ff(j)*xvar(j+1) - gg(j)*xvar(j+2) - hh(j)*xvar(j+3)
    do j=jmax-4,1,-1
       xvar(j) = xvar(j) - ff(j)*xvar(j+1) - gg(j)*xvar(j+2) - hh(j)*xvar(j+3) - ii(j)*xvar(j+4) 
    end do
    
  end subroutine ninediagonal_simple_solve

  
  subroutine ninediagonal_cyclic_pre                                               &
   &( nn_cyclic, imax, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, jmax2_iter, &
   &  pp, qq, rr, ss, tt, uu, vv, ww, jptr_iter )
  
    integer,intent(in) :: nn_cyclic
    integer,intent(in) :: imax
    real(8),intent(inout) :: aa(imax+3)
    real(8),intent(inout) :: bb(imax+3)
    real(8),intent(inout) :: cc(imax+3)
    real(8),intent(inout) :: dd(imax+3)
    real(8),intent(inout) :: ee(imax+3)
    real(8),intent(inout) :: ff(imax+3)
    real(8),intent(inout) :: gg(imax+3)
    real(8),intent(inout) :: hh(imax+3)
    real(8),intent(inout) :: ii(imax+3)
    real(8),intent(inout) :: jj(imax+3)
    real(8),intent(inout) :: kk(imax+3)
    real(8),intent(inout) :: ll(imax+3)
    integer,intent(out)  :: jmax2_iter(0:nn_cyclic)
    real(8),intent(out),optional :: pp(imax)
    real(8),intent(out),optional :: qq(imax)
    real(8),intent(out),optional :: rr(imax)
    real(8),intent(out),optional :: ss(imax)
    real(8),intent(out),optional :: tt(imax)
    real(8),intent(out),optional :: uu(imax)
    real(8),intent(out),optional :: vv(imax)
    real(8),intent(out),optional :: ww(imax)
    integer,intent(out),optional :: jptr_iter(nn_cyclic)
    real(8) :: aa1(imax)
    real(8) :: bb1(imax)
    real(8) :: cc1(imax)
    real(8) :: dd1(imax)
    real(8) :: ee1(imax)
    real(8) :: ff1(imax)
    real(8) :: gg1(imax)
    real(8) :: hh1(imax)
    real(8) :: ii1(imax)
    real(8) :: jj1(imax)
    real(8) :: kk1(imax)
    real(8) :: ll1(imax)
    real(8) :: aa2(imax)
    real(8) :: bb2(imax)
    real(8) :: cc2(imax)
    real(8) :: dd2(imax)
    real(8) :: ee2(imax)
    real(8) :: ff2(imax)
    real(8) :: gg2(imax)
    real(8) :: hh2(imax)
    real(8) :: ii2(imax)
    real(8) :: jj2(imax)
    real(8) :: kk2(imax)
    real(8) :: ll2(imax)
    integer :: imax0,jmax0,jmax1,jmax2,j,j2,j3,j4,n,jptr,jptr0
    
    imax0 = ( (imax-1)/4 + 1 )*4
    
    do j=imax+1,imax0
       aa(j) = 0.0d0
       bb(j) = 0.0d0
       cc(j) = 0.0d0
       dd(j) = 0.0d0
       ee(j) = 1.0d0
       ff(j) = 0.0d0
       gg(j) = 0.0d0
       hh(j) = 0.0d0
       ii(j) = 0.0d0
       jj(j) = 0.0d0
       kk(j) = 0.0d0
       ll(j) = 0.0d0
    end do
    
    do j=1,imax0,4
       ll(j) = 0.0d0
       kk(j) = 0.0d0
       jj(j) = 0.0d0
       j2=j+1
       ll(j2) = 0.0d0
       kk(j2) = 0.0d0
       jj(j2) = ii(j2)
       ii(j2) = hh(j2)
       hh(j2) = gg(j2)
       gg(j2) = ff(j2)
       ff(j2) = ee(j2)
       ee(j2) = dd(j2)
       dd(j2) = cc(j2)
       cc(j2) = bb(j2)
       bb(j2) = aa(j2)
       aa(j2) = 0.0d0
       j3=j+2
       ll(j3) = 0.0d0
       kk(j3) = ii(j3)
       jj(j3) = hh(j3)
       ii(j3) = gg(j3)
       hh(j3) = ff(j3)
       gg(j3) = ee(j3)
       ff(j3) = dd(j3)
       ee(j3) = cc(j3)
       dd(j3) = bb(j3)
       cc(j3) = aa(j3)
       bb(j3) = 0.0d0
       aa(j3) = 0.0d0
       j4=j+3
       ll(j4) = ii(j4)
       kk(j4) = hh(j4)
       jj(j4) = gg(j4)
       ii(j4) = ff(j4)
       hh(j4) = ee(j4)
       gg(j4) = dd(j4)
       ff(j4) = cc(j4)
       ee(j4) = bb(j4)
       dd(j4) = aa(j4)
       cc(j4) = 0.0d0
       bb(j4) = 0.0d0
       aa(j4) = 0.0d0
    end do
    
!   | ee ff gg hh | ii          |             |
!   | ee ff gg hh | ii jj       |             |
!   | ee ff gg hh | ii jj kk    |             |
!   | ee ff gg hh | ii jj kk ll |             |
!   |-----------------------------------------|
!   | aa bb cc dd | ee ff gg hh | ii          |
!   |    bb cc dd | ee ff gg hh | ii jj       |
!   |       cc dd | ee ff gg hh | ii jj kk    |
!   |          dd | ee ff gg hh | ii jj kk ll |
!   |-----------------------------------------|
!   |             | aa bb cc dd | ee ff gg hh |
!   |             |    bb cc dd | ee ff gg hh |
!   |             |       cc dd | ee ff gg hh |
!   |             |          dd | ee ff gg hh |
   
!   | aa bb cc dd |
!   | aa bb cc dd |
!   | aa bb cc dd |
!   | aa bb cc dd | = A

!   ! ee ff gg hh |
!   ! ee ff gg hh |
!   ! ee ff gg hh |
!   | ee ff gg hh | = B

!   ! ii jj kk ll |
!   ! ii jj kk ll |
!   ! ii jj kk ll |
!   | ii jj kk ll | = C
    
    jmax0 = imax0/4
    jmax2_iter(0) = jmax0
    jptr = 1
    
    do n=1,nn_cyclic
       jmax1 = (jmax0+1)/2
       jmax2 = jmax0 - jmax1
       jmax2_iter(n) = jmax2
       jptr_iter(n) = jptr
       if ( n == 1 ) then
          call set_array_cyclic                                                               &
           &( jmax0, jmax1, jmax2, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll,            &
           &  aa1(jptr), bb1(jptr), cc1(jptr), dd1(jptr), ee1(jptr), ff1(jptr),               &
           &  gg1(jptr), hh1(jptr), ii1(jptr), jj1(jptr), kk1(jptr), ll1(jptr),               &
           &  pp(jptr), qq(jptr), rr(jptr), ss(jptr), tt(jptr), uu(jptr), vv(jptr), ww(jptr), &
           &  aa2(jptr), bb2(jptr), cc2(jptr), dd2(jptr), ee2(jptr), ff2(jptr),               &
           &  gg2(jptr), hh2(jptr), ii2(jptr), jj2(jptr), kk2(jptr), ll2(jptr) )
       else
          call set_array_cyclic                                                                            &
           &( jmax0, jmax1, jmax2, aa2(jptr0), bb2(jptr0), cc2(jptr0), dd2(jptr0), ee2(jptr0), ff2(jptr0), &
           &  gg2(jptr0), hh2(jptr0), ii2(jptr0), jj2(jptr0), kk2(jptr0), ll2(jptr0),         &
           &  aa1(jptr), bb1(jptr), cc1(jptr), dd1(jptr), ee1(jptr), ff1(jptr),               &
           &  gg1(jptr), hh1(jptr), ii1(jptr), jj1(jptr), kk1(jptr), ll1(jptr),               &
           &  pp(jptr), qq(jptr), rr(jptr), ss(jptr), tt(jptr), uu(jptr), vv(jptr), ww(jptr), &
           &  aa2(jptr), bb2(jptr), cc2(jptr), dd2(jptr), ee2(jptr), ff2(jptr),               &
           &  gg2(jptr), hh2(jptr), ii2(jptr), jj2(jptr), kk2(jptr), ll2(jptr) )
       end if
       jmax0 = jmax2
       jptr0 = jptr
       jptr = jptr + jmax1*4
    end do
    
    if ( nn_cyclic > 0 ) then
       do j=1,jptr-1
          aa(j) = aa1(j)
          bb(j) = bb1(j)
          cc(j) = cc1(j)
          dd(j) = dd1(j)
          ee(j) = ee1(j)
          ff(j) = ff1(j)
          gg(j) = gg1(j)
          hh(j) = hh1(j)
          ii(j) = ii1(j)
          jj(j) = jj1(j)
          kk(j) = kk1(j)
          ll(j) = ll1(j)
       end do
    
       do j=jptr,jptr+jmax2*4-1
          aa(j)   = aa2(j-jmax1*4)
          bb(j)   = bb2(j-jmax1*4)
          cc(j)   = cc2(j-jmax1*4)
          dd(j)   = dd2(j-jmax1*4)
          ee(j)   = ee2(j-jmax1*4)
          ff(j)   = ff2(j-jmax1*4)
          gg(j)   = gg2(j-jmax1*4)
          hh(j)   = hh2(j-jmax1*4)
          ii(j)   = ii2(j-jmax1*4)
          jj(j)   = jj2(j-jmax1*4)
          kk(j)   = kk2(j-jmax1*4)
          ll(j)   = ll2(j-jmax1*4)
       end do
    end if

!    write(6,*) "1:aa(jptr:jptr+jmax2*4-1)=",aa(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:bb(jptr:jptr+jmax2*4-1)=",bb(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:cc(jptr:jptr+jmax2*4-1)=",cc(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:dd(jptr:jptr+jmax2*4-1)=",dd(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:ee(jptr:jptr+jmax2*4-1)=",ee(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:ff(jptr:jptr+jmax2*4-1)=",ff(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:gg(jptr:jptr+jmax2*4-1)=",gg(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:hh(jptr:jptr+jmax2*4-1)=",hh(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:ii(jptr:jptr+jmax2*4-1)=",ii(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:jj(jptr:jptr+jmax2*4-1)=",jj(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:kk(jptr:jptr+jmax2*4-1)=",kk(jptr:jptr+jmax0*4-1)
!    write(6,*) "1:ll(jptr:jptr+jmax2*4-1)=",ll(jptr:jptr+jmax0*4-1)

    call ninediagonal_block_solve_pre &
     &( jmax0, aa(jptr), bb(jptr), cc(jptr), dd(jptr), ee(jptr), ff(jptr), &
     &  gg(jptr), hh(jptr), ii(jptr), jj(jptr), kk(jptr), ll(jptr) )
  
  end subroutine ninediagonal_cyclic_pre
  
  
  subroutine set_array_cyclic                                              &
   &( jmax0, jmax1, jmax2, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, &
   &  aa1, bb1, cc1, dd1, ee1, ff1, gg1, hh1, ii1, jj1, kk1, ll1,          &
   &  pp2, qq2, rr2, ss2, tt2, uu2, vv2, ww2,                              &
   &  aa2, bb2, cc2, dd2, ee2, ff2, gg2, hh2, ii2, jj2, kk2, ll2 )
   
    integer,intent(in) :: jmax0, jmax1, jmax2
    real(8),intent(inout) :: aa(jmax0*4)
    real(8),intent(inout) :: bb(jmax0*4)
    real(8),intent(inout) :: cc(jmax0*4)
    real(8),intent(inout) :: dd(jmax0*4)
    real(8),intent(inout) :: ee(jmax0*4)
    real(8),intent(inout) :: ff(jmax0*4)
    real(8),intent(inout) :: gg(jmax0*4)
    real(8),intent(inout) :: hh(jmax0*4)
    real(8),intent(inout) :: ii(jmax0*4)
    real(8),intent(inout) :: jj(jmax0*4)
    real(8),intent(inout) :: kk(jmax0*4)
    real(8),intent(inout) :: ll(jmax0*4)
    
    real(8),intent(out) :: aa1(jmax1*4)
    real(8),intent(out) :: bb1(jmax1*4)
    real(8),intent(out) :: cc1(jmax1*4)
    real(8),intent(out) :: dd1(jmax1*4)
    real(8),intent(out) :: ee1(jmax1*4)
    real(8),intent(out) :: ff1(jmax1*4)
    real(8),intent(out) :: gg1(jmax1*4)
    real(8),intent(out) :: hh1(jmax1*4)
    real(8),intent(out) :: ii1(jmax1*4)
    real(8),intent(out) :: jj1(jmax1*4)
    real(8),intent(out) :: kk1(jmax1*4)
    real(8),intent(out) :: ll1(jmax1*4)
    
    real(8),intent(out) :: pp2(jmax2*4)
    real(8),intent(out) :: qq2(jmax2*4)
    real(8),intent(out) :: rr2(jmax2*4)
    real(8),intent(out) :: ss2(jmax2*4)
    real(8),intent(out) :: tt2(jmax2*4)
    real(8),intent(out) :: uu2(jmax2*4)
    real(8),intent(out) :: vv2(jmax2*4)
    real(8),intent(out) :: ww2(jmax2*4)
    
    real(8),intent(out) :: aa2(jmax2*4)
    real(8),intent(out) :: bb2(jmax2*4)
    real(8),intent(out) :: cc2(jmax2*4)
    real(8),intent(out) :: dd2(jmax2*4)
    real(8),intent(out) :: ee2(jmax2*4)
    real(8),intent(out) :: ff2(jmax2*4)
    real(8),intent(out) :: gg2(jmax2*4)
    real(8),intent(out) :: hh2(jmax2*4)
    real(8),intent(out) :: ii2(jmax2*4)
    real(8),intent(out) :: jj2(jmax2*4)
    real(8),intent(out) :: kk2(jmax2*4)
    real(8),intent(out) :: ll2(jmax2*4)
    
    real(8) :: w1(0:3)
    real(8) :: w2(0:3)
    real(8) :: w3(0:3)
    real(8) :: w4(0:3)
    integer :: j1,j,j0,k

!    write(6,*) "aa(1:jmax0*4)=",aa(1:jmax0*4)
!    write(6,*) "bb(1:jmax0*4)=",bb(1:jmax0*4)
!    write(6,*) "cc(1:jmax0*4)=",cc(1:jmax0*4)
!    write(6,*) "dd(1:jmax0*4)=",dd(1:jmax0*4)
!    write(6,*) "ee(1:jmax0*4)=",ee(1:jmax0*4)
!    write(6,*) "ff(1:jmax0*4)=",ff(1:jmax0*4)
!    write(6,*)

    do j0 = 1, jmax1
       j  = j0*8 - 3
       j1 = j0*4 - 3   
       
       !! [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1] 
       !! = [e f g h;e f g h;e f g h;e f g h]^(-1)(j0-1)
       do k=0,3
          ee1(j1+k) = ee(j-4+k)
          ff1(j1+k) = ff(j-4+k)
          gg1(j1+k) = gg(j-4+k)
          hh1(j1+k) = hh(j-4+k)
       end do
       call inverse_4x4( ee1(j1), ff1(j1), gg1(j1), hh1(j1) )
       
       !! [a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1] 
       !! = [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]*[a b c d;a b c d;a b c d;a b c d](j0-1)
       call multiply_4x4                                                            &
        & ( ee1(j1), ff1(j1), gg1(j1), hh1(j1), aa(j-4), bb(j-4), cc(j-4), dd(j-4), &
        &   aa1(j1), bb1(j1), cc1(j1), dd1(j1) )
       
       !! [i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1] 
       !! = [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]*[i j k l;i j k l;i j k l;i j k l](j0-1)
       call multiply_4x4                                                          &
        &( ee1(j1), ff1(j1), gg1(j1), hh1(j1), ii(j-4), jj(j-4), kk(j-4), ll(j-4), &
        &  ii1(j1), jj1(j1), kk1(j1), ll1(j1) )
    end do
    
!    write(6,*) "jmax1,jmax0=",jmax1,jmax0
!    write(6,*) "aa1(1:jmax1*2)=",aa1(1:jmax1*2)
!    write(6,*) "bb1(1:jmax1*2)=",bb1(1:jmax1*2)
!    write(6,*) "cc1(1:jmax1*2)=",cc1(1:jmax1*2)
!    write(6,*) "dd1(1:jmax1*2)=",dd1(1:jmax1*2)
!    write(6,*) "ee1(1:jmax1*2)=",ee1(1:jmax1*2)
!    write(6,*) "ff1(1:jmax1*2)=",ff1(1:jmax1*2)
!    stop 333
    
    do j0 = 1, jmax2
       j  = j0*8 - 3
       j1 = j0*4 - 3   

       if ( jmax1 == jmax2 .and. j0 == jmax2 ) then
       
          !! [p2 q2 r2 s2;p2 q2 r2 s2;p2 q2 r2 s2;p2 q2 r2 s2] 
          !! = [a b c d;a b c d;a b c d;a b c d]*[e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]
          call multiply_4x4                                                   &
           &( aa(j), bb(j), cc(j), dd(j), ee1(j1), ff1(j1), gg1(j1), hh1(j1), &
           &  pp2(j1), qq2(j1), rr2(j1), ss2(j1) )
           
          !! [t2 u2 v2 w2;t2 u2 v2 w2;t2 u2 v2 w2;t2 u2 v2 w2] = 0.0d0
          tt2(j1:j1+3) = 0.0d0
          uu2(j1:j1+3) = 0.0d0
          vv2(j1:j1+3) = 0.0d0
          ww2(j1:j1+3) = 0.0d0
       
          !! [a2 b2 c2 d2;a2 b2 c2 d2;a2 b2 c2 d2;a2 b2 c2 d2] 
          !! = - [a b c d;a b c d;a b c d;a b c d]*[a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1]
          call multiply_4x4                                                   &
           &( aa(j), bb(j), cc(j), dd(j), aa1(j1), bb1(j1), cc1(j1), dd1(j1), &
           &  aa2(j1), bb2(j1), cc2(j1), dd2(j1) )
          do k=0,3
             aa2(j1+k) = - aa2(j1+k)
             bb2(j1+k) = - bb2(j1+k)
             cc2(j1+k) = - cc2(j1+k)
             dd2(j1+k) = - dd2(j1+k)
          end do  
       
          !! [e2 f2 g2 h2;e2 f2 g2 h2;e2 f2 g2 h2;e2 f2 g2 h2]
          !! = [e f g h;e f g h;e f g h;e f g h]
          !!   - [a b c d;a b c d;a b c d;a b c d]*[i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1]
          !!   - [i j k l;i j k l;i j k l;i j k l]*[a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1](j0+1)
          call multiply_4x4                                                   &
           &( aa(j), bb(j), cc(j), dd(j), ii1(j1), jj1(j1), kk1(j1), ll1(j1), &
           &  ee2(j1), ff2(j1), gg2(j1), hh2(j1) )
          do k=0,3
             ee2(j1+k) = ee(j+k) - ee2(j1+k)
             ff2(j1+k) = ff(j+k) - ff2(j1+k)
             gg2(j1+k) = gg(j+k) - gg2(j1+k)
             hh2(j1+k) = hh(j+k) - hh2(j1+k)
          end do  
          
          !! [i2 j2 k2 l2;i2 j2 k2 l2;i2 j2 k2 l2;i2 j2 k2 l2] = 0.0
          ii2(j1:j1+3) = 0.0d0
          jj2(j1:j1+3) = 0.0d0
          kk2(j1:j1+3) = 0.0d0
          ll2(j1:j1+3) = 0.0d0
          
       else
       
          !! [p2 q2 r2 s2;p2 q2 r2 s2;p2 q2 r2 s2;p2 q2 r2 s2] 
          !! = [a b c d;a b c d;a b c d;a b c d]*[e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]
          call multiply_4x4                                                   &
           &( aa(j), bb(j), cc(j), dd(j), ee1(j1), ff1(j1), gg1(j1), hh1(j1), &
           &  pp2(j1), qq2(j1), rr2(j1), ss2(j1) )
       
          !! [t2 u2 v2 w2;t2 u2 v2 w2;t2 u2 v2 w2;t2 u2 v2 w2] 
          !! = [i j k l;i j k l;i j k l;i j k l]*[e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1](j0+1)
          call multiply_4x4                                                   &
           &( ii(j), jj(j), kk(j), ll(j), ee1(j1+4), ff1(j1+4), gg1(j1+4), hh1(j1+4), &
           &  tt2(j1), uu2(j1), vv2(j1), ww2(j1) )
       
          !! [a2 b2 c2 d2;a2 b2 c2 d2;a2 b2 c2 d2;a2 b2 c2 d2] 
          !! = - [a b c d;a b c d;a b c d;a b c d]*[a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1]
          call multiply_4x4                                                   &
           &( aa(j), bb(j), cc(j), dd(j), aa1(j1), bb1(j1), cc1(j1), dd1(j1), &
           &  aa2(j1), bb2(j1), cc2(j1), dd2(j1) )
          do k=0,3
             aa2(j1+k) = - aa2(j1+k)
             bb2(j1+k) = - bb2(j1+k)
             cc2(j1+k) = - cc2(j1+k)
             dd2(j1+k) = - dd2(j1+k)
          end do  
       
          !! [e2 f2 g2 h2;e2 f2 g2 h2;e2 f2 g2 h2;e2 f2 g2 h2]
          !! = [e f g h;e f g h;e f g h;e f g h]
          !!   - [a b c d;a b c d;a b c d;a b c d]*[i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1]
          !!   - [i j k l;i j k l;i j k l;i j k l]*[a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1](j0+1)
          call multiply_4x4                                                   &
           &( aa(j), bb(j), cc(j), dd(j), ii1(j1), jj1(j1), kk1(j1), ll1(j1), &
           &  ee2(j1), ff2(j1), gg2(j1), hh2(j1) )
          call multiply_4x4                                                   &
           &( ii(j), jj(j), kk(j), ll(j), aa1(j1+4), bb1(j1+4), cc1(j1+4), dd1(j1+4), &
           &  w1, w2, w3, w4 )
          do k=0,3
             ee2(j1+k) = ee(j+k) - ee2(j1+k) - w1(k)
             ff2(j1+k) = ff(j+k) - ff2(j1+k) - w2(k)
             gg2(j1+k) = gg(j+k) - gg2(j1+k) - w3(k)
             hh2(j1+k) = hh(j+k) - hh2(j1+k) - w4(k)
          end do  
          
          !! [i2 j2 k2 l2;i2 j2 k2 l2;i2 j2 k2 l2;i2 j2 k2 l2]
          !! = - [i j k l;i j k l;i j k l;i j k l]*[i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1](j1+1)
          call multiply_4x4                                                   &
           &( ii(j), jj(j), kk(j), ll(j), ii1(j1+4), jj1(j1+4), kk1(j1+4), ll1(j1+4), &
           &  ii2(j1), jj2(j1), kk2(j1), ll2(j1) )
          do k=0,3
             ii2(j1+k) = - ii2(j1+k)
             jj2(j1+k) = - jj2(j1+k)
             kk2(j1+k) = - kk2(j1+k)
             ll2(j1+k) = - ll2(j1+k)
          end do  
       end if
    end do


!    write(6,*) "aa2(1:jmax2*2)=",aa2(1:jmax2*2)
!    write(6,*) "bb2(1:jmax2*2)=",bb2(1:jmax2*2)
!    write(6,*) "cc2(1:jmax2*2)=",cc2(1:jmax2*2)
!    write(6,*) "dd2(1:jmax2*2)=",dd2(1:jmax2*2)
!    write(6,*) "ee2(1:jmax2*2)=",ee2(1:jmax2*2)
!    write(6,*) "ff2(1:jmax2*2)=",ff2(1:jmax2*2)


!     write(6,*) "ii2=",ii2
!     write(6,*) "jj2=",jj2
!     write(6,*) "kk2=",kk2
!     write(6,*) "ll2=",ll2
     
!    stop 333
    
!    write(6,*) "array(:)=",array(:)
!    write(6,*) "a1(:)=",a1(:)
!    write(6,*) "a2(:)=",a2(:)
!    write(6,*) "a3(:)=",a3(:)




!    stop 333

     
    
  end subroutine set_array_cyclic


  subroutine ninediagonal_block_solve_pre &
   &( jmax0, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll )
   
    integer,intent(in) :: jmax0
    real(8),intent(inout) :: aa(jmax0*4)
    real(8),intent(inout) :: bb(jmax0*4)
    real(8),intent(inout) :: cc(jmax0*4)
    real(8),intent(inout) :: dd(jmax0*4)
    real(8),intent(inout) :: ee(jmax0*4)
    real(8),intent(inout) :: ff(jmax0*4)
    real(8),intent(inout) :: gg(jmax0*4)
    real(8),intent(inout) :: hh(jmax0*4)
    real(8),intent(inout) :: ii(jmax0*4)
    real(8),intent(inout) :: jj(jmax0*4)
    real(8),intent(inout) :: kk(jmax0*4)
    real(8),intent(inout) :: ll(jmax0*4)
    real(8) :: w1(4)
    real(8) :: w2(4)
    real(8) :: w3(4)
    real(8) :: w4(4)
    integer :: j0,j

!    write(6,*) "jmax0=",jmax0

  !! LU factorization of block ninediagonal matrix

    j0 = 1
    j = j0*4 - 3

    !! [e f g h;e f g h;e f g h;e f g h] 
    !! = [e f g h;e f g h;e f g h;e f g h]^(-1)
    call inverse_4x4( ee(j), ff(j), gg(j), hh(j) )
    
!    write(6,*) "ee=",ee
!    write(6,*) "ff=",ff
!    write(6,*) "gg=",gg
!    write(6,*) "hh=",hh
       
    !! [i j k l;i j k l;i j k l;i j k l] 
    !! = [e f g h;e f g h;e f g h;e f g h]*[i j k l;i j k l;i j k l;i j k l]
    call multiply_4x4( ee(j), ff(j), gg(j), hh(j), ii(j), jj(j), kk(j), ll(j), &
     &             w1, w2, w3, w4 )
    ii(j:j+3) = w1(:)
    jj(j:j+3) = w2(:)
    kk(j:j+3) = w3(:)
    ll(j:j+3) = w4(:)

    do j0 = 2, jmax0
       j = j0*4 - 3
       
       !! [e f g h;e f g h;e f g h;e f g h] = [e f g h;e f g h;e f g h;e f g h]
       !! - [a b c d;a b c d;a b c d;a b c d]*[i j k l;i j k l;i j k l;i j k l](j0-1)
       call multiply_4x4( aa(j), bb(j), cc(j), dd(j), ii(j-4), jj(j-4), kk(j-4), ll(j-4), &
        &             w1, w2, w3, w4 )
       ee(j:j+3) = ee(j:j+3) - w1(:)
       ff(j:j+3) = ff(j:j+3) - w2(:)
       gg(j:j+3) = gg(j:j+3) - w3(:)
       hh(j:j+3) = hh(j:j+3) - w4(:)

       !! [e f g h;e f g h;e f g h;e f g h] 
       !! = [e f g h;e f g h;e f g h;e f g h]^(-1)
       call inverse_4x4( ee(j), ff(j), gg(j), hh(j) )
       
       !! [a b c d;a b c d;a b c d;a b c d] 
       !! = [e f g h;e f g h;e f g h;e f g h]*[a b c d;a b c d;a b c d;a b c d]
       call multiply_4x4( ee(j), ff(j), gg(j), hh(j), aa(j), bb(j), cc(j), dd(j),  &
        &             w1, w2, w3, w4 )
       aa(j:j+3) = w1(:)
       bb(j:j+3) = w2(:)
       cc(j:j+3) = w3(:)
       dd(j:j+3) = w4(:)
       
       !! [i j k l;i j k l;i j k l;i j k l]
       !! = [e f g h;e f g h;e f g h;e f g h]*[i j k l;i j k l;i j k l;i j k l]
       call multiply_4x4( ee(j), ff(j), gg(j), hh(j), ii(j), jj(j), kk(j), ll(j),  &
        &             w1, w2, w3, w4 )
       ii(j:j+3) = w1(:)
       jj(j:j+3) = w2(:)
       kk(j:j+3) = w3(:)
       ll(j:j+3) = w4(:)
    end do
    
  end subroutine ninediagonal_block_solve_pre


  subroutine ninediagonal_block_solve                       &
   &( jmax0, aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, &
   &  yvar, xvar )
   
    integer,intent(in) :: jmax0
    real(8),intent(in) :: aa(jmax0*4)
    real(8),intent(in) :: bb(jmax0*4)
    real(8),intent(in) :: cc(jmax0*4)
    real(8),intent(in) :: dd(jmax0*4)
    real(8),intent(in) :: ee(jmax0*4)
    real(8),intent(in) :: ff(jmax0*4)
    real(8),intent(in) :: gg(jmax0*4)
    real(8),intent(in) :: hh(jmax0*4)
    real(8),intent(in) :: ii(jmax0*4)
    real(8),intent(in) :: jj(jmax0*4)
    real(8),intent(in) :: kk(jmax0*4)
    real(8),intent(in) :: ll(jmax0*4)
    complex(8),intent(in) :: yvar(jmax0*4)
    complex(8),intent(out) :: xvar(jmax0*4)
    integer :: j0,j1,j


!    write(6,*) "jmax0=",jmax0

  !! Forward substitution

    j0 = 1
    j1 = j0*4 - 3
    !! [xvar;xvar;var;xvar] = [e f g h;e f g h;e f g h;e f g h]*[yvar;yvar;yvar;yvar]
    do j= j1, j1+3
       xvar(j) = ee(j)*yvar(j1) + ff(j)*yvar(j1+1) + gg(j)*yvar(j1+2) + hh(j)*yvar(j1+3)
    end do

    do j0 = 2, jmax0
       j1 = j0*4 - 3
       !! [xvar;xvar;xvar;xvar] = [e f g h;e f g h;e f g h;e f g h]*[yvar;yvar;yvar;yvar]
       !!                         - [a b c d;a b c d;a b c d;a b c d]*[xvar;xvar;xvar;xvar](j0-1)
       do j= j1, j1+3
          xvar(j) = ee(j)*yvar(j1) + ff(j)*yvar(j1+1) + gg(j)*yvar(j1+2) + hh(j)*yvar(j1+3)   &
           &        - aa(j)*xvar(j1-4) - bb(j)*xvar(j1-3) - cc(j)*xvar(j1-2) - dd(j)*xvar(j1-1) 
       end do
    end do

  !! Backward substitution

    do j0 = jmax0-1,1,-1
       j1 = j0*4 - 3
       !! [xvar;xvar;xvar;xvar] = [xvar;xvar;xvar;xvar] - [i j k l;i j k l;i j k l;i j k l]*[xvar;xvar;xvar;xvar](j0+1)
       do j= j0*4-3, j0*4
          xvar(j) = xvar(j) - ii(j)*xvar(j1+4) - jj(j)*xvar(j1+5) - kk(j)*xvar(j1+6) - ll(j)*xvar(j1+7) 
       end do
    end do
    
  end subroutine ninediagonal_block_solve
  

  subroutine ninediagonal_cyclic_solve                               &
   &( nn_cyclic, imax, aa1, bb1, cc1, dd1, ee1, ff1, gg1, hh1, ii1, jj1, kk1, ll1,   &
   &  pp2, qq2, rr2, ss2, tt2, uu2, vv2, ww2, jmax2_iter, jptr_iter, &
   &  y1, yvar, xvar )
   
    integer,intent(in) :: nn_cyclic
    integer,intent(in) :: imax
    real(8),intent(in) :: aa1(imax+3)
    real(8),intent(in) :: bb1(imax+3)
    real(8),intent(in) :: cc1(imax+3)
    real(8),intent(in) :: dd1(imax+3)
    real(8),intent(in) :: ee1(imax+3)
    real(8),intent(in) :: ff1(imax+3)
    real(8),intent(in) :: gg1(imax+3)
    real(8),intent(in) :: hh1(imax+3)
    real(8),intent(in) :: ii1(imax+3)
    real(8),intent(in) :: jj1(imax+3)
    real(8),intent(in) :: kk1(imax+3)
    real(8),intent(in) :: ll1(imax+3)
    real(8),intent(in) :: pp2(imax)
    real(8),intent(in) :: qq2(imax)
    real(8),intent(in) :: rr2(imax)
    real(8),intent(in) :: ss2(imax)
    real(8),intent(in) :: tt2(imax)
    real(8),intent(in) :: uu2(imax)
    real(8),intent(in) :: vv2(imax)
    real(8),intent(in) :: ww2(imax)
    integer,intent(in) :: jmax2_iter(0:nn_cyclic)
    integer,intent(in) :: jptr_iter(nn_cyclic)
    complex(8),intent(inout) :: y1(imax+3)
    complex(8),intent(inout) :: yvar(imax+3)
    complex(8),intent(inout) :: xvar(imax+3)
    integer :: jmax0,jmax1,jmax2,n,jptr
    
    !! For 'call ninediagonal_block_solve' when nn_cyclic = 0
    n=0
    jmax2 = jmax2_iter(0)
    jmax1 = 0
    jptr  = 1
    
    do n=1,nn_cyclic
      jmax0 = jmax2_iter(n-1)
      jmax2 = jmax2_iter(n)
      jmax1 = jmax0 - jmax2
      jptr = jptr_iter(n)
       
      call ninediagonal_cyclic_forward                                                                  &
       &( jmax0, jmax1, jmax2, ee1(jptr), ff1(jptr), gg1(jptr), hh1(jptr),                              &
       &  pp2(jptr), qq2(jptr), rr2(jptr), ss2(jptr), tt2(jptr), uu2(jptr), vv2(jptr), ww2(jptr), yvar, &
       &  y1(jptr), xvar  )
      yvar = xvar
    end do

    call ninediagonal_block_solve                                                          &
     &( jmax2, aa1(jptr+jmax1*4), bb1(jptr+jmax1*4), cc1(jptr+jmax1*4), dd1(jptr+jmax1*4),  &
     &  ee1(jptr+jmax1*4), ff1(jptr+jmax1*4), gg1(jptr+jmax1*4), hh1(jptr+jmax1*4),         &
     &  ii1(jptr+jmax1*4), jj1(jptr+jmax1*4), kk1(jptr+jmax1*4), ll1(jptr+jmax1*4), yvar,   &
     &  xvar   )

    do n=nn_cyclic,1,-1
      jmax0 = jmax2_iter(n-1)
      jmax2 = jmax2_iter(n)
      jmax1 = jmax0 - jmax2
      jptr  = jptr_iter(n)
      call ninediagonal_cyclic_backward                                        &
       &( jmax0, jmax1, jmax2, aa1(jptr), bb1(jptr), cc1(jptr), dd1(jptr),    &
       &  ii1(jptr), jj1(jptr), kk1(jptr), ll1(jptr), y1(jptr), xvar,         &
       &  yvar )
      xvar = yvar
    end do
   
  end subroutine ninediagonal_cyclic_solve


  subroutine ninediagonal_cyclic_forward  &
   &( jmax0, jmax1, jmax2, ee1, ff1, gg1, hh1, pp2, qq2, rr2, ss2, tt2, uu2, vv2, ww2, yvar, &
   &  y1, yvar2  )
    integer,intent(in) :: jmax0
    integer,intent(in) :: jmax1
    integer,intent(in) :: jmax2
    real(8),intent(in) :: ee1(jmax1*4)
    real(8),intent(in) :: ff1(jmax1*4)
    real(8),intent(in) :: gg1(jmax1*4)
    real(8),intent(in) :: hh1(jmax1*4)
    real(8),intent(in) :: pp2(jmax2*4)
    real(8),intent(in) :: qq2(jmax2*4)
    real(8),intent(in) :: rr2(jmax2*4)
    real(8),intent(in) :: ss2(jmax2*4)
    real(8),intent(in) :: tt2(jmax2*4)
    real(8),intent(in) :: uu2(jmax2*4)
    real(8),intent(in) :: vv2(jmax2*4)
    real(8),intent(in) :: ww2(jmax2*4)
    complex(8),intent(in) :: yvar(jmax0*4)
    complex(8),intent(out) :: y1(jmax1*4)
    complex(8),intent(out) :: yvar2(jmax2*4)
    integer :: j1, j, j0, k
    
    
!    write(6,*) "jmax0,jmax1,jmax2=",jmax0,jmax1,jmax2
    
    if ( jmax1 > jmax2 ) then
       do j0 = 1, jmax2
          j  = j0*8 - 3
          j1 = j0*4 - 3
          
          !! [y1;y1;y1;y1] = [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]*[yvar;yvar;yvar;yvar](j0-1)
          do k=0,3
             y1(j1+k) = ee1(j1+k)*yvar(j-4) + ff1(j1+k)*yvar(j-3) + gg1(j1+k)*yvar(j-2) + hh1(j1+k)*yvar(j-1)
          end do
          
          !! [yvar2;yvar2;yvar2;yvar2] = [yvar;yvar;yvar;yvar]
          !!  - [pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2]*[yvar;yvar;yvar;yvar](jj-1)
          !!  - [tt2 uu2 vv2 ww2;tt2 uu2 vv2 ww2;tt2 uu2 vv2 ww2;tt2 uu2 vv2 ww2]*[yvar;yvar;yvar;yvar](jj+1)
          do k=0,3
             yvar2(j1+k) = yvar(j+k) - pp2(j1+k)*yvar(j-4) - qq2(j1+k)*yvar(j-3) - rr2(j1+k)*yvar(j-2) - ss2(j1+k)*yvar(j-1) &
              &                      - tt2(j1+k)*yvar(j+4) - uu2(j1+k)*yvar(j+5) - vv2(j1+k)*yvar(j+6) - ww2(j1+k)*yvar(j+7)
          end do
       end do
       
       j0 = jmax1
       j  = j0*8 - 3
       j1 = j0*4 - 3
          
       !! [y1;y1;y1;y1] = [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]*[yvar;yvar;yvar;yvar](j0-1)
       do k=0,3
          y1(j1+k) = ee1(j1+k)*yvar(j-4) + ff1(j1+k)*yvar(j-3) + gg1(j1+k)*yvar(j-2) + hh1(j1+k)*yvar(j-1)
       end do
    else
       do j0 = 1, jmax2-1
          j  = j0*8 - 3
          j1 = j0*4 - 3
          
          !! [y1;y1;y1;y1] = [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]*[yvar;yvar;yvar;yvar](j0-1)
          do k=0,3
             y1(j1+k) = ee1(j1+k)*yvar(j-4) + ff1(j1+k)*yvar(j-3) + gg1(j1+k)*yvar(j-2) + hh1(j1+k)*yvar(j-1)
          end do
          
          !! [yvar2;yvar2;yvar2;yvar2] = [yvar;yvar;yvar;yvar]
          !!  - [pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2]*[yvar;yvar;yvar;yvar](jj-1)
          !!  - [tt2 uu2 vv2 ww2;tt2 uu2 vv2 ww2;tt2 uu2 vv2 ww2;tt2 uu2 vv2 ww2]*[yvar;yvar;yvar;yvar](jj+1)
          do k=0,3
             yvar2(j1+k) = yvar(j+k) - pp2(j1+k)*yvar(j-4) - qq2(j1+k)*yvar(j-3) - rr2(j1+k)*yvar(j-2) - ss2(j1+k)*yvar(j-1) &
              &                      - tt2(j1+k)*yvar(j+4) - uu2(j1+k)*yvar(j+5) - vv2(j1+k)*yvar(j+6) - ww2(j1+k)*yvar(j+7)
          end do
       end do
       j0 = jmax2
       j  = j0*8 - 3
       j1 = j0*4 - 3
          
       !! [y1;y1;y1;y1] = [e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1;e1 f1 g1 h1]*[yvar;yvar;yvar;yvar](j0-1)
       do k=0,3
          y1(j1+k) = ee1(j1+k)*yvar(j-4) + ff1(j1+k)*yvar(j-3) + gg1(j1+k)*yvar(j-2) + hh1(j1+k)*yvar(j-1)
       end do   
          
       !! [yvar2;yvar2;yvar2;yvar2] = [yvar;yvar;yvar;yvar]
       !!  - [pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2;pp2 qq2 rr2 ss2]*[yvar;yvar;yvar;yvar](jj-1)
       do k=0,3
          yvar2(j1+k) = yvar(j+k) - pp2(j1+k)*yvar(j-4) - qq2(j1+k)*yvar(j-3) - rr2(j1+k)*yvar(j-2) - ss2(j1+k)*yvar(j-1)
       end do
    end if
    
!    write(6,*)    
!    write(6,*) "jmax1=",jmax1
!    write(6,*) "yvar(:)=",yvar(:)
!    write(6,*) "y1(:)=",y1(:)
!    write(6,*) "yvar2(:)=",yvar2(:)
!    write(6,*) "pp2=",pp2
!    write(6,*) "qq2=",qq2
!    write(6,*) "rr2=",rr2
!    write(6,*) "ss2=",ss2
    
  end subroutine ninediagonal_cyclic_forward


  subroutine ninediagonal_cyclic_backward                                       &
   &( jmax0, jmax1, jmax2, aa1, bb1, cc1, dd1, ii1, jj1, kk1, ll1, y1, xvar2,  &
   &  xvar )
    integer,intent(in) :: jmax0
    integer,intent(in) :: jmax1
    integer,intent(in) :: jmax2
    real(8),intent(in) :: aa1(jmax1*4)
    real(8),intent(in) :: bb1(jmax1*4)
    real(8),intent(in) :: cc1(jmax1*4)
    real(8),intent(in) :: dd1(jmax1*4)
    real(8),intent(in) :: ii1(jmax1*4)
    real(8),intent(in) :: jj1(jmax1*4)
    real(8),intent(in) :: kk1(jmax1*4)
    real(8),intent(in) :: ll1(jmax1*4)
    complex(8),intent(in) :: y1(jmax1*4)
    complex(8),intent(in) :: xvar2(jmax2*4)
    complex(8),intent(out) :: xvar(jmax0*4)
    integer :: j1, j, j0, k
    
    j0 = 1
    j  = j0*8 - 3
    j1 = j0*4 - 3
    
       !! [xvar;xvar;xvar;xvar](j0-1) = [y1;y1;y1;y1] 
       !!        - [i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1]*[xvar2;xvar2;xvar2;xvar2]
       do k=0,3
           xvar(j-4+k) = y1(j1+k) - ii1(j1+k)*xvar2(j1)   - jj1(j1+k)*xvar2(j1+1) &
            &                     - kk1(j1+k)*xvar2(j1+2) - ll1(j1+k)*xvar2(j1+3)
       end do
       
       !! [xvar;xvar;xvar;xvar] = [xvar2;xvar2;xvar2;xvar2]
       xvar(j:j+3)   = xvar2(j1:j1+3)

    do j0 = 2, jmax2
       j  = j0*8 - 3
       j1 = j0*4 - 3
             
       !! [xvar;xvar;xvar;xvar](j0-1) = [y1;y1;y1;y1] 
       !!  - [a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1]*[xvar2,xvar2;xvar2;xvar2](j0-1)
       !!  - [i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1;i1 j1 k1 l1]*[xvar2;xvar2;xvar2;xvar2]
       do k=0,3
           xvar(j-4+k) = y1(j1+k) - aa1(j1+k)*xvar2(j1-4) - bb1(j1+k)*xvar2(j1-3) &
            &                     - cc1(j1+k)*xvar2(j1-2) - dd1(j1+k)*xvar2(j1-1) &
            &                     - ii1(j1+k)*xvar2(j1)   - jj1(j1+k)*xvar2(j1+1) &
            &                     - kk1(j1+k)*xvar2(j1+2) - ll1(j1+k)*xvar2(j1+3)
       end do
       
       !! [xvar;xvar;xvar;xvar] = [xvar2;xvar2;xvar2;xvar2]
       xvar(j:j+3)   = xvar2(j1:j1+3)
    end do
    
    if ( jmax1 > jmax2 ) then
       j0 = jmax1
       j  = j0*8 - 3
       j1 = j0*4 - 3
       
       !! [xvar;xvar;xvar;xvar](jj-1) = [y1;y1;y1;y1] 
       !!  - [a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1;a1 b1 c1 d1]*[xvar2,xvar2;xvar2;xvar2](j0-1)
       do k=0,3
           xvar(j-4+k) = y1(j1+k) - aa1(j1+k)*xvar2(j1-4) - bb1(j1+k)*xvar2(j1-3) &
            &                     - cc1(j1+k)*xvar2(j1-2) - dd1(j1+k)*xvar2(j1-1)
       end do
    end if
    
  end subroutine ninediagonal_cyclic_backward

end module ninediagonal

