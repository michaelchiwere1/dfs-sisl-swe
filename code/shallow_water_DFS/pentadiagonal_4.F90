module pentadiagonal_4

  implicit none
  
  private
  public :: pentadiagonal_4__allocate, pentadiagonal_4__ini, pentadiagonal_4__solve, pentadiagonal_4__test
  
!  integer,parameter :: jcn_cyclic = 0  !! LU decomposition of nine-diagonal matrix
  integer,parameter :: jcn_cyclic = 1  !! LU decomposition of block tri-diagonal matrix
!  integer,parameter :: jcn_cyclic = 2  !! Cyclic reduction for block tri-diagonal matrix

!xx  call main_sub
  
contains


  subroutine pentadiagonal_4__test

    integer,parameter :: imax=1000
    
!xx    real(8),pointer :: dat(:,:)
!xx    integer,pointer :: iww(:)
    real(8),allocatable :: dat(:,:)
    integer,allocatable :: iww(:)
    real(8) :: yy(imax+1,4)
    real(8) :: yvar(imax+1,4)
    real(8) :: xvar(imax+1,4)
    integer :: n

    call pentadiagonal_4__allocate( imax, imax+1, dat, iww )

    dat(:,1) = 2.0d0
    dat(:,2) = 1.0d0
    dat(:,3) = 3.0d0
    dat(:,4) = 1.0d0
    dat(:,5) = 2.0d0

    call pentadiagonal_4__ini( imax, imax+1, dat, iww )

    n=1
!    do n=1,1000
      yvar(1,:) = 6.0d0
      yvar(2,:) = 7.0d0
      yvar(3:imax-2,:) = 9.0d0
      yvar(imax-1,:) = 7.0d0
      yvar(imax,:) = 6.0d0
      
      call pentadiagonal_4__solve   &
       &( imax, imax+1, dat, iww,   &
       &  yy, yvar, xvar )
!    end do
  
    write(6,*) "xvar(1:imax,3)=",xvar(1:imax,3)
  
  end subroutine pentadiagonal_4__test
  
  
  subroutine pentadiagonal_4__allocate( imax, ilen, dat, iww )
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
!xx    real(8),intent(out),pointer :: dat(:,:)
!xx    integer,intent(out),pointer :: iww(:)
    real(8),intent(out),allocatable :: dat(:,:)
    integer,intent(out),allocatable :: iww(:)
    
    integer :: nn_tmp, nn_cyclic

    if ( ilen < imax ) then
      write(6,*) "Error: tridiagonal__ini: ilen should be equal to or larger than imax."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
    
    if ( jcn_cyclic == 0 ) then
      allocate( dat(ilen,5) )
      allocate( iww(1) )
    
    else if ( jcn_cyclic == 1 ) then
      allocate( dat(ilen,6) )
      allocate( iww(2) )
      iww(1) = 0
    
    else if ( jcn_cyclic >= 2 ) then
      if ( jcn_cyclic == 1 ) then
        nn_cyclic = 0
      else
        nn_tmp = int(log(imax+0.001d0)/log(2.0d0))
        nn_cyclic = max( nn_tmp-2, 0 )
      end if   
      allocate( dat(ilen,10) )
      allocate( iww(2+2*nn_cyclic) )
      
      iww(1) = nn_cyclic
      write(6,*) "nn_cyclic=",nn_cyclic
    end if
  
  end subroutine pentadiagonal_4__allocate


  subroutine pentadiagonal_4__ini( imax, ilen, dat, iww )
  
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    real(8),intent(inout) :: dat(ilen,*)
    integer,intent(out) :: iww(*)
  
    integer :: nn_tmp, nn_cyclic

    if ( ilen < imax+1 ) then
      write(6,*) "Error: pentadiagonal_4__ini: ilen should be equal to or larger than imax+1."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
    
    if ( jcn_cyclic == 0 ) then
      call pentadiagonal_4_simple_pre( imax,  dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5) )
    else
      if ( jcn_cyclic == 1 ) then
        nn_cyclic = 0
        call pentadiagonal_4_cyclic_pre( nn_cyclic, imax,       &!IN
         & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5),    &!INOUT
         & dat(1,6), iww(2)    )                                 !OUT
      else
        nn_cyclic = iww(1)
        call pentadiagonal_4_cyclic_pre( nn_cyclic, imax,                               &!IN
         & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5),                            &!INOUT
         & dat(1,6), iww(2), dat(1,7), dat(1,8), dat(1,9), dat(1,10), iww(3+nn_cyclic) ) !OUT
      end if
    end if

  end subroutine pentadiagonal_4__ini


  subroutine pentadiagonal_4__solve  &
   &( imax, ilen, dat, iww,        &
   &  yy, yvar, xvar )
    
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    real(8),intent(in) :: dat(ilen,10)
    integer,intent(in) :: iww(ilen)
    real(8),intent(inout) :: yy(ilen)
    real(8),intent(inout) :: yvar(ilen,4)
    real(8),intent(inout) :: xvar(ilen,4)
      
    integer :: nn_cyclic,jmax2

    if ( ilen < imax+1 ) then
      write(6,*) "Error: pentadiagonal_4__solve: ilen should be equal to or larger than imax+1."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if

    yvar(imax+1:ilen,:) = 0.0d0
    if ( jcn_cyclic == 0 ) then
      call pentadiagonal_4_simple_solve( imax,                 &
       &  dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5),    &
       &  yvar(1,1), yvar(1,2), yvar(1,3), yvar(1,4),          &!IN
       &  xvar(1,1), xvar(1,2), xvar(1,3), xvar(1,4)   )        !OUT
    else if ( jcn_cyclic == 1 ) then
      jmax2 = iww(2)
      call pentadiagonal_4_block_solve                                        &!IN
       &( jmax2, dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6),  &!IN
       &  yvar(1,1), yvar(1,2), yvar(1,3), yvar(1,4),                         &!IN
       &  xvar(1,1), xvar(1,2), xvar(1,3), xvar(1,4)   )                       !OUT
    else !! jcn_cyclic == 2
      nn_cyclic = iww(1)
      call pentadiagonal_4_cyclic_solve( nn_cyclic, imax,                   &!IN
       & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), dat(1,6),        &!IN
       & dat(1,7), dat(1,8), dat(1,9), dat(1,10), iww(2), iww(3+nn_cyclic), &!IN
       & yy, yvar, xvar )                                            !OUT,IO,OUT
    end if
    xvar(imax+1:ilen,:) = 0.0d0

  end subroutine pentadiagonal_4__solve


  subroutine pentadiagonal_4_simple_pre( jmax, aa, bb, cc, dd, ee )
    
    integer,intent(in) :: jmax
    real(8),intent(inout) :: aa(jmax)
    real(8),intent(inout) :: bb(jmax)
    real(8),intent(inout) :: cc(jmax)
    real(8),intent(inout) :: dd(jmax)
    real(8),intent(inout) :: ee(jmax)
    integer :: j

    aa(1) = 0.0d0
    bb(1) = 0.0d0
    cc(1) = 1.0d0/cc(1)
    dd(1) = dd(1)*cc(1)
    ee(1) = ee(1)*cc(1)
    aa(2) = 0.0d0
    cc(2) = 1.0d0/( cc(2) - bb(2)*dd(1) )
    dd(2) = ( dd(2) - bb(2)*ee(1) )*cc(2)
    ee(2) = ee(2)*cc(2)
    do j=3,jmax
       bb(j) = bb(j) - aa(j)*dd(j-2)
       cc(j) = 1.0d0/( cc(j) - aa(j)*ee(j-2) - bb(j)*dd(j-1) )
       dd(j) = ( dd(j) - bb(j)*ee(j-1) )*cc(j)
       ee(j) = ee(j)*cc(j)
    end do
    ee(jmax-1) = 0.0d0
    dd(jmax) = 0.0d0
    ee(jmax) = 0.0d0

  end subroutine pentadiagonal_4_simple_pre
  
  
  
  subroutine pentadiagonal_4_simple_solve                    &
   &( jmax, aa, bb, cc, dd, ee, yvar1, yvar2, yvar3, yvar4,  &
   & xvar1, xvar2, xvar3, xvar4 )
    integer,intent(in) :: jmax
    real(8),intent(in) :: aa(jmax)
    real(8),intent(in) :: bb(jmax)
    real(8),intent(in) :: cc(jmax)
    real(8),intent(in) :: dd(jmax)
    real(8),intent(in) :: ee(jmax)
    real(8),intent(in) :: yvar1(jmax)
    real(8),intent(in) :: yvar2(jmax)
    real(8),intent(in) :: yvar3(jmax)
    real(8),intent(in) :: yvar4(jmax)
    real(8),intent(out) :: xvar1(jmax)
    real(8),intent(out) :: xvar2(jmax)
    real(8),intent(out) :: xvar3(jmax)
    real(8),intent(out) :: xvar4(jmax)
    integer :: j

    !! Forward substitution
    xvar1(1) = yvar1(1)*cc(1)
    xvar2(1) = yvar2(1)*cc(1)
    xvar3(1) = yvar3(1)*cc(1)
    xvar4(1) = yvar4(1)*cc(1)
    xvar1(2) = ( yvar1(2) - bb(2)*xvar1(1) )*cc(2)
    xvar2(2) = ( yvar2(2) - bb(2)*xvar2(1) )*cc(2)
    xvar3(2) = ( yvar3(2) - bb(2)*xvar3(1) )*cc(2)
    xvar4(2) = ( yvar4(2) - bb(2)*xvar4(1) )*cc(2)
    do j=3,jmax
      xvar1(j) = ( yvar1(j) - bb(j)*xvar1(j-1) - aa(j)*xvar1(j-2) )*cc(j)
      xvar2(j) = ( yvar2(j) - bb(j)*xvar2(j-1) - aa(j)*xvar2(j-2) )*cc(j)
      xvar3(j) = ( yvar3(j) - bb(j)*xvar3(j-1) - aa(j)*xvar3(j-2) )*cc(j)
      xvar4(j) = ( yvar4(j) - bb(j)*xvar4(j-1) - aa(j)*xvar4(j-2) )*cc(j)
    end do
  
    !! Backward substitution
    j = jmax - 1
    xvar1(j) = xvar1(j) - dd(j)*xvar1(j+1)
    xvar2(j) = xvar2(j) - dd(j)*xvar2(j+1)
    xvar3(j) = xvar3(j) - dd(j)*xvar3(j+1)
    xvar4(j) = xvar4(j) - dd(j)*xvar4(j+1)
    do j=jmax-2,1,-1
      xvar1(j) = xvar1(j) - dd(j)*xvar1(j+1) - ee(j)*xvar1(j+2)
      xvar2(j) = xvar2(j) - dd(j)*xvar2(j+1) - ee(j)*xvar2(j+2)
      xvar3(j) = xvar3(j) - dd(j)*xvar3(j+1) - ee(j)*xvar3(j+2)
      xvar4(j) = xvar4(j) - dd(j)*xvar4(j+1) - ee(j)*xvar4(j+2)
    end do
    
  end subroutine pentadiagonal_4_simple_solve

  
  subroutine pentadiagonal_4_cyclic_pre( nn_cyclic, imax, aa, bb, cc, dd, ee, ff, jmax2_iter, &
   &                                   pp, qq, rr, ss, jptr_iter )
  
    integer,intent(in) :: nn_cyclic
    integer,intent(in) :: imax
    real(8),intent(inout) :: aa(imax+1)
    real(8),intent(inout) :: bb(imax+1)
    real(8),intent(inout) :: cc(imax+1)
    real(8),intent(inout) :: dd(imax+1)
    real(8),intent(inout) :: ee(imax+1)
    real(8),intent(out)   :: ff(imax+1)
    integer,intent(out)   :: jmax2_iter(0:nn_cyclic)
    real(8),intent(out),optional :: pp(imax)
    real(8),intent(out),optional :: qq(imax)
    real(8),intent(out),optional :: rr(imax)
    real(8),intent(out),optional :: ss(imax)
    integer,intent(out),optional :: jptr_iter(nn_cyclic)
    real(8) :: aa1(imax)
    real(8) :: bb1(imax)
    real(8) :: cc1(imax)
    real(8) :: dd1(imax)
    real(8) :: ee1(imax)
    real(8) :: ff1(imax)
    real(8) :: aa2(imax)
    real(8) :: bb2(imax)
    real(8) :: cc2(imax)
    real(8) :: dd2(imax)
    real(8) :: ee2(imax)
    real(8) :: ff2(imax)
    integer :: imax0,jmax0,jmax1,jmax2,j,n,jptr,jptr0
    
    imax0 = ( (imax-1)/2 + 1 )*2

    do j=imax+1,imax0
       aa(j) = 0.0d0
       bb(j) = 0.0d0
       cc(j) = 1.0d0
       dd(j) = 0.0d0
       ee(j) = 0.0d0
    end do
    ff(:) = 0.0d0
    
    do j=2,imax0,2
       ff(j) = ee(j)
       ee(j) = dd(j)
       dd(j) = cc(j)
       cc(j) = bb(j)
       bb(j) = aa(j)
       aa(j) = 0.0d0
    end do
    aa(1) = 0.0d0
    bb(1) = 0.0d0
    aa(2) = 0.0d0
    bb(2) = 0.0d0
    ee(imax0-1) = 0.0d0
    ff(imax0-1) = 0.0d0
    ee(imax0)   = 0.0d0
    ff(imax0)   = 0.0d0
    
!   | cc dd | ee    |       |       |
!   | cc dd | ee ff |       |       |
!   |-------------------------------|
!   | aa bb | cc dd | ee    |       |
!   |    bb | cc dd | ee ff |       |
!   |-------------------------------|
!   |       | aa bb | cc dd | ee    |
!   |       |    bb | cc dd | ee ff |
!   |-------------------------------|
!   |       |       | aa bb | cc dd |
!   |       |       |    bb | cc dd |
   
!   | aa bb |
!   | aa bb | = A

!   ! cc dd |
!   | cc dd | = B

!   ! ee ff |
!   | ee ff | = C
    
    jmax0 = imax0/2
    jmax2_iter(0) = jmax0
    jptr = 1
    
    if ( nn_cyclic > 0 ) then
       do n=1,nn_cyclic
          jmax1 = (jmax0+1)/2
          jmax2 = jmax0 - jmax1
          jmax2_iter(n) = jmax2
          jptr_iter(n) = jptr
          if ( n == 1 ) then
             call set_array_cyclic                                                 &
              &( jmax0, jmax1, jmax2, aa, bb, cc, dd, ee, ff,                      &
              &  aa1(jptr), bb1(jptr), cc1(jptr), dd1(jptr), ee1(jptr), ff1(jptr), &
              &  pp(jptr), qq(jptr), rr(jptr), ss(jptr),                           &
              &  aa2(jptr), bb2(jptr), cc2(jptr), dd2(jptr), ee2(jptr), ff2(jptr) )
          else
             call set_array_cyclic                                                                            &
              &( jmax0, jmax1, jmax2, aa2(jptr0), bb2(jptr0), cc2(jptr0), dd2(jptr0), ee2(jptr0), ff2(jptr0), &
              &  aa1(jptr), bb1(jptr), cc1(jptr), dd1(jptr), ee1(jptr), ff1(jptr),                            &
              &  pp(jptr), qq(jptr), rr(jptr), ss(jptr),                                                      &
              &  aa2(jptr), bb2(jptr), cc2(jptr), dd2(jptr), ee2(jptr), ff2(jptr) )
          end if
          jmax0 = jmax2
          jptr0 = jptr
          jptr = jptr + jmax1*2
       end do

       do j=1,jptr-1
          aa(j) = aa1(j)
          bb(j) = bb1(j)
          cc(j) = cc1(j)
          dd(j) = dd1(j)
          ee(j) = ee1(j)
          ff(j) = ff1(j)
       end do
    
       do j=jptr,jptr+jmax2*2-1
          aa(j)   = aa2(j-jmax1*2)
          bb(j)   = bb2(j-jmax1*2)
          cc(j)   = cc2(j-jmax1*2)
          dd(j)   = dd2(j-jmax1*2)
          ee(j)   = ee2(j-jmax1*2)
          ff(j)   = ff2(j-jmax1*2)
       end do
    end if

!    write(6,*) "1:aa(jptr:jptr+jmax2*2-1)=",aa(jptr:jptr+jmax0*2-1)
!    write(6,*) "1:bb(jptr:jptr+jmax2*2-1)=",bb(jptr:jptr+jmax0*2-1)
!    write(6,*) "1:cc(jptr:jptr+jmax2*2-1)=",cc(jptr:jptr+jmax0*2-1)
!    write(6,*) "1:dd(jptr:jptr+jmax2*2-1)=",dd(jptr:jptr+jmax0*2-1)
!    write(6,*) "1:ee(jptr:jptr+jmax2*2-1)=",ee(jptr:jptr+jmax0*2-1)
!    write(6,*) "1:ff(jptr:jptr+jmax2*2-1)=",ff(jptr:jptr+jmax0*2-1)

    call pentadiagonal_4_block_solve_pre( jmax0, aa(jptr), bb(jptr), cc(jptr), dd(jptr), ee(jptr), ff(jptr) )

!    write(6,*) "2:aa(jptr:jptr+jmax2*2-1)=",aa(jptr:jptr+jmax0*2-1)
!    write(6,*) "2:bb(jptr:jptr+jmax2*2-1)=",bb(jptr:jptr+jmax0*2-1)
!    write(6,*) "2:cc(jptr:jptr+jmax2*2-1)=",cc(jptr:jptr+jmax0*2-1)
!    write(6,*) "2:dd(jptr:jptr+jmax2*2-1)=",dd(jptr:jptr+jmax0*2-1)
!    write(6,*) "2:ee(jptr:jptr+jmax2*2-1)=",ee(jptr:jptr+jmax0*2-1)
!    write(6,*) "2:ff(jptr:jptr+jmax2*2-1)=",ff(jptr:jptr+jmax0*2-1)
  
  end subroutine pentadiagonal_4_cyclic_pre
  
  
  subroutine set_array_cyclic( jmax0, jmax1, jmax2, aa, bb, cc, dd, ee, ff,      &
   &                           aa1, bb1, cc1, dd1, ee1, ff1, pp2, qq2, rr2, ss2, &
   &                           aa2, bb2, cc2, dd2, ee2, ff2 )
   
    integer,intent(in) :: jmax0, jmax1, jmax2
    real(8),intent(inout) :: aa(jmax0*2)
    real(8),intent(inout) :: bb(jmax0*2)
    real(8),intent(inout) :: cc(jmax0*2)
    real(8),intent(inout) :: dd(jmax0*2)
    real(8),intent(inout) :: ee(jmax0*2)
    real(8),intent(inout) :: ff(jmax0*2)
    real(8),intent(out) :: aa1(jmax1*2)
    real(8),intent(out) :: bb1(jmax1*2)
    real(8),intent(out) :: cc1(jmax1*2)
    real(8),intent(out) :: dd1(jmax1*2)
    real(8),intent(out) :: ee1(jmax1*2)
    real(8),intent(out) :: ff1(jmax1*2)
    real(8),intent(out) :: pp2(jmax2*2)
    real(8),intent(out) :: qq2(jmax2*2)
    real(8),intent(out) :: rr2(jmax2*2)
    real(8),intent(out) :: ss2(jmax2*2)
    real(8),intent(out) :: aa2(jmax2*2)
    real(8),intent(out) :: bb2(jmax2*2)
    real(8),intent(out) :: cc2(jmax2*2)
    real(8),intent(out) :: dd2(jmax2*2)
    real(8),intent(out) :: ee2(jmax2*2)
    real(8),intent(out) :: ff2(jmax2*2)
    real(8) :: div
    integer :: j1,j,jj

!    write(6,*) "aa(1:jmax0*2)=",aa(1:jmax0*2)
!    write(6,*) "bb(1:jmax0*2)=",bb(1:jmax0*2)
!    write(6,*) "cc(1:jmax0*2)=",cc(1:jmax0*2)
!    write(6,*) "dd(1:jmax0*2)=",dd(1:jmax0*2)
!    write(6,*) "ee(1:jmax0*2)=",ee(1:jmax0*2)
!    write(6,*) "ff(1:jmax0*2)=",ff(1:jmax0*2)
!    write(6,*)

    do jj = 1, jmax1
       j  = jj*4 - 1
       j1 = jj*2 - 1 
       
       !! [cc1,dd1;cc1,dd1] = [cc,dd;cc,dd]^(-1)(jj-1)
       div = 1.0d0/( cc(j-2)*dd(j-1) - dd(j-2)*cc(j-1) )
       cc1(j1)   =  dd(j-1)*div
       dd1(j1)   = -dd(j-2)*div
       cc1(j1+1) = -cc(j-1)*div
       dd1(j1+1) =  cc(j-2)*div

       !! [aa1,bb1;aa1,bb1] = [cc1,dd1;cc1,dd1]*[aa,bb;aa,bb]
       aa1(j1)   = cc1(j1)*aa(j-2)   + dd1(j1)*aa(j-1)
       bb1(j1)   = cc1(j1)*bb(j-2)   + dd1(j1)*bb(j-1)
       aa1(j1+1) = cc1(j1+1)*aa(j-2) + dd1(j1+1)*aa(j-1)
       bb1(j1+1) = cc1(j1+1)*bb(j-2) + dd1(j1+1)*bb(j-1)
       
       !! [ee1,ff1;ee1,ff1] = [cc1,dd1;cc1,dd1]*[ee,ff;ee,ff]
       ee1(j1)   = cc1(j1)*ee(j-2)   + dd1(j1)*ee(j-1)
       ff1(j1)   = cc1(j1)*ff(j-2)   + dd1(j1)*ff(j-1)
       ee1(j1+1) = cc1(j1+1)*ee(j-2) + dd1(j1+1)*ee(j-1)
       ff1(j1+1) = cc1(j1+1)*ff(j-2) + dd1(j1+1)*ff(j-1)
    end do
    
!    write(6,*) "jmax1,jmax0=",jmax1,jmax0
!    write(6,*) "aa1(1:jmax1*2)=",aa1(1:jmax1*2)
!    write(6,*) "bb1(1:jmax1*2)=",bb1(1:jmax1*2)
!    write(6,*) "cc1(1:jmax1*2)=",cc1(1:jmax1*2)
!    write(6,*) "dd1(1:jmax1*2)=",dd1(1:jmax1*2)
!    write(6,*) "ee1(1:jmax1*2)=",ee1(1:jmax1*2)
!    write(6,*) "ff1(1:jmax1*2)=",ff1(1:jmax1*2)
!    stop 333
    
    do jj = 1, jmax2
       j  = jj*4 - 1
       j1 = jj*2 - 1

       if ( jmax1 == jmax2 .and. jj == jmax2 ) then
          !! [pp2,qq2;pp2,qq2] = [aa,bb;aa,bb]*[cc1,dd1;cc1,dd1]
          pp2(j1)   = aa(j)*cc1(j1)   + bb(j)*cc1(j1+1)
          qq2(j1)   = aa(j)*dd1(j1)   + bb(j)*dd1(j1+1)
          pp2(j1+1) = aa(j+1)*cc1(j1) + bb(j+1)*cc1(j1+1)
          qq2(j1+1) = aa(j+1)*dd1(j1) + bb(j+1)*dd1(j1+1)
          
          !! [rr2,ss2;rr2,ss2] = 0.0d0
          rr2(j1)   = 0.0d0
          ss2(j1)   = 0.0d0
          rr2(j1+1) = 0.0d0
          ss2(j1+1) = 0.0d0
          
          !! [aa2,bb2;aa2,bb2] = -[aa,bb;aa,bb]*[aa1,bb1;aa1,bb1]
          aa2(j1)   = - aa(j)*aa1(j1)   - bb(j)*aa1(j1+1)
          bb2(j1)   = - aa(j)*bb1(j1)   - bb(j)*bb1(j1+1)
          aa2(j1+1) = - aa(j+1)*aa1(j1) - bb(j+1)*aa1(j1+1)
          bb2(j1+1) = - aa(j+1)*bb1(j1) - bb(j+1)*bb1(j1+1)
          
          !! [cc2,dd2;cc2,dd2] = [cc,dd;cc,dd] - [aa,bb;aa,bb]*[ee1,ff1;ee1,ff1]
          cc2(j1)   = cc(j)   - aa(j)*ee1(j1)     - bb(j)*ee1(j1+1)
          dd2(j1)   = dd(j)   - aa(j)*ff1(j1)     - bb(j)*ff1(j1+1)
          cc2(j1+1) = cc(j+1) - aa(j+1)*ee1(j1)   - bb(j+1)*ee1(j1+1)
          dd2(j1+1) = dd(j+1) - aa(j+1)*ff1(j1)   - bb(j+1)*ff1(j1+1)
          
          !! [ee2,ff2;ee2,ff2] = 0.0d0
          ee2(j1)   = 0.0d0
          ff2(j1)   = 0.0d0
          ee2(j1+1) = 0.0d0
          ff2(j1+1) = 0.0d0
       else
          !! [pp2,qq2;pp2,qq2] = [aa,bb;aa,bb]*[cc1,dd1;cc1,dd1]
          pp2(j1)   = aa(j)*cc1(j1)   + bb(j)*cc1(j1+1)
          qq2(j1)   = aa(j)*dd1(j1)   + bb(j)*dd1(j1+1)
          pp2(j1+1) = aa(j+1)*cc1(j1) + bb(j+1)*cc1(j1+1)
          qq2(j1+1) = aa(j+1)*dd1(j1) + bb(j+1)*dd1(j1+1)
          
          !! [rr2,ss2;rr2,ss2] = [ee,ff,ee,ff]*[cc1,dd1;cc1,dd1](jj+1)
          rr2(j1)   = ee(j)*cc1(j1+2)   + ff(j)*cc1(j1+3)
          ss2(j1)   = ee(j)*dd1(j1+2)   + ff(j)*dd1(j1+3)
          rr2(j1+1) = ee(j+1)*cc1(j1+2) + ff(j+1)*cc1(j1+3)
          ss2(j1+1) = ee(j+1)*dd1(j1+2) + ff(j+1)*dd1(j1+3)
          
          !! [aa2,bb2;aa2,bb2] = -[aa,bb;aa,bb]*[aa1,bb1;aa1,bb1]
          aa2(j1)   = - aa(j)*aa1(j1)   - bb(j)*aa1(j1+1)
          bb2(j1)   = - aa(j)*bb1(j1)   - bb(j)*bb1(j1+1)
          aa2(j1+1) = - aa(j+1)*aa1(j1) - bb(j+1)*aa1(j1+1)
          bb2(j1+1) = - aa(j+1)*bb1(j1) - bb(j+1)*bb1(j1+1)
          
          !! [cc2,dd2;cc2,dd2] = [cc,dd;cc,dd] - [aa,bb;aa,bb]*[ee1,ff1;ee1,ff1]
          !!                             - [ee,ff;ee,ff]*[aa1,bb1;aa1,bb1](jj+1)
          cc2(j1)   = cc(j)   - aa(j)*ee1(j1)     - bb(j)*ee1(j1+1) &
           &                  - ee(j)*aa1(j1+2)   - ff(j)*aa1(j1+3)
          dd2(j1)   = dd(j)   - aa(j)*ff1(j1)     - bb(j)*ff1(j1+1) &
           &                  - ee(j)*bb1(j1+2)   - ff(j)*bb1(j1+3)
          cc2(j1+1) = cc(j+1) - aa(j+1)*ee1(j1)   - bb(j+1)*ee1(j1+1) &
           &                  - ee(j+1)*aa1(j1+2) - ff(j+1)*aa1(j1+3)
          dd2(j1+1) = dd(j+1) - aa(j+1)*ff1(j1)   - bb(j+1)*ff1(j1+1) &
           &                  - ee(j+1)*bb1(j1+2) - ff(j+1)*bb1(j1+3)
          
          !! [ee2,ff2;ee2,ff2] = -[ee,ff;ee,ff]*[ee1,ff1;ee1,ff1](jj+1)
          ee2(j1)   = - ee(j)*ee1(j1+2)   - ff(j)*ee1(j1+3)
          ff2(j1)   = - ee(j)*ff1(j1+2)   - ff(j)*ff1(j1+3)
          ee2(j1+1) = - ee(j+1)*ee1(j1+2) - ff(j+1)*ee1(j1+3)
          ff2(j1+1) = - ee(j+1)*ff1(j1+2) - ff(j+1)*ff1(j1+3)
       end if
    end do


!    write(6,*) "aa2(1:jmax2*2)=",aa2(1:jmax2*2)
!    write(6,*) "bb2(1:jmax2*2)=",bb2(1:jmax2*2)
!    write(6,*) "cc2(1:jmax2*2)=",cc2(1:jmax2*2)
!    write(6,*) "dd2(1:jmax2*2)=",dd2(1:jmax2*2)
!    write(6,*) "ee2(1:jmax2*2)=",ee2(1:jmax2*2)
!    write(6,*) "ff2(1:jmax2*2)=",ff2(1:jmax2*2)
!    stop 333
    
!    write(6,*) "array(:)=",array(:)
!    write(6,*) "a1(:)=",a1(:)
!    write(6,*) "a2(:)=",a2(:)
!    write(6,*) "a3(:)=",a3(:)

!    stop 333
    
  end subroutine set_array_cyclic


  subroutine pentadiagonal_4_block_solve_pre( jmax0, aa, bb, cc, dd, ee, ff )
   
    integer,intent(in) :: jmax0
    real(8),intent(inout) :: aa(jmax0*2)
    real(8),intent(inout) :: bb(jmax0*2)
    real(8),intent(inout) :: cc(jmax0*2)
    real(8),intent(inout) :: dd(jmax0*2)
    real(8),intent(inout) :: ee(jmax0*2)
    real(8),intent(inout) :: ff(jmax0*2)
    real(8) :: div,c1,c2,d1,d2,e1,e2,f1,f2
    integer :: j0,j

  !! LU factorization of block pentadiagonal_4 matrix

    j0 = 1
    j = j0*2 - 1

    aa(j) = 0.0d0
    bb(j) = 0.0d0

    !! [cc,dd;cc,dd] = [cc,dd;cc,dd]^(-1)
    c1 = cc(j)
    d1 = dd(j)
    c2 = cc(j+1)
    d2 = dd(j+1)
    div = 1.0d0/( c1*d2 - d1*c2 )
    cc(j)   =  d2*div
    dd(j)   = -d1*div
    cc(j+1) = -c2*div
    dd(j+1) =  c1*div
       
    !! [ee,ff;ee,ff] = [cc,dd;cc,dd]*[ee,ff;ee,ff]
    e1 = cc(j)*ee(j)   + dd(j)*ee(j+1)
    f1 = cc(j)*ff(j)   + dd(j)*ff(j+1)
    e2 = cc(j+1)*ee(j) + dd(j+1)*ee(j+1)
    f2 = cc(j+1)*ff(j) + dd(j+1)*ff(j+1)
    ee(j)   = e1
    ff(j)   = f1
    ee(j+1) = e2
    ff(j+1) = f2

    do j0 = 2, jmax0
       j = j0*2 - 1

       !! [cc,dd;cc,dd] = [cc,dd;cc,dd] - [aa,bb;aa,bb]*[ee,ff;ee,ff](j0-1)
       c1 = cc(j)   - aa(j)*ee(j-2)   - bb(j)*ee(j-1)
       d1 = dd(j)   - aa(j)*ff(j-2)   - bb(j)*ff(j-1)
       c2 = cc(j+1) - aa(j+1)*ee(j-2) - bb(j+1)*ee(j-1)
       d2 = dd(j+1) - aa(j+1)*ff(j-2) - bb(j+1)*ff(j-1)

       !! [cc,dd;cc,dd] = [cc,dd;cc,dd]^(-1)
       div = 1.0d0/( c1*d2 - d1*c2 )
       cc(j)   =  d2*div
       dd(j)   = -d1*div
       cc(j+1) = -c2*div
       dd(j+1) =  c1*div

       !! [aa,bb;aa,bb] = [cc,dd;cc,dd]*[aa,bb;aa,bb]
       e1 = cc(j)*aa(j)   + dd(j)*aa(j+1)
       f1 = cc(j)*bb(j)   + dd(j)*bb(j+1)
       e2 = cc(j+1)*aa(j) + dd(j+1)*aa(j+1)
       f2 = cc(j+1)*bb(j) + dd(j+1)*bb(j+1)
       aa(j)   = e1
       bb(j)   = f1
       aa(j+1) = e2
       bb(j+1) = f2
       
       !! [ee,ff;ee,ff] = [cc,dd;cc,dd]*[ee,ff;ee,ff]
       e1 = cc(j)*ee(j)   + dd(j)*ee(j+1)
       f1 = cc(j)*ff(j)   + dd(j)*ff(j+1)
       e2 = cc(j+1)*ee(j) + dd(j+1)*ee(j+1)
       f2 = cc(j+1)*ff(j) + dd(j+1)*ff(j+1)
       ee(j)   = e1
       ff(j)   = f1
       ee(j+1) = e2
       ff(j+1) = f2
    end do
    
  end subroutine pentadiagonal_4_block_solve_pre


  subroutine pentadiagonal_4_block_solve                         &
   &( jmax0, aa, bb, cc, dd, ee, ff, yvar1, yvar2, yvar3, yvar4, &
   &  xvar1, xvar2, xvar3, xvar4 )
   
    integer,intent(in) :: jmax0
    real(8),intent(in) :: aa(jmax0*2)
    real(8),intent(in) :: bb(jmax0*2)
    real(8),intent(in) :: cc(jmax0*2)
    real(8),intent(in) :: dd(jmax0*2)
    real(8),intent(in) :: ee(jmax0*2)
    real(8),intent(in) :: ff(jmax0*2)
    real(8),intent(in) :: yvar1(jmax0*2)
    real(8),intent(in) :: yvar2(jmax0*2)
    real(8),intent(in) :: yvar3(jmax0*2)
    real(8),intent(in) :: yvar4(jmax0*2)
    real(8),intent(out) :: xvar1(jmax0*2)
    real(8),intent(out) :: xvar2(jmax0*2)
    real(8),intent(out) :: xvar3(jmax0*2)
    real(8),intent(out) :: xvar4(jmax0*2)
    integer :: j0,j

  !! Forward substitution

    j0 = 1
    j = j0*2 - 1
    !! [xvar;xvar] = [cc,dd;cc,dd]*[yvar;yvar]
    xvar1(j)   = cc(j)*yvar1(j)   + dd(j)*yvar1(j+1)
    xvar2(j)   = cc(j)*yvar2(j)   + dd(j)*yvar2(j+1)
    xvar3(j)   = cc(j)*yvar3(j)   + dd(j)*yvar3(j+1)
    xvar4(j)   = cc(j)*yvar4(j)   + dd(j)*yvar4(j+1)
    xvar1(j+1) = cc(j+1)*yvar1(j) + dd(j+1)*yvar1(j+1)
    xvar2(j+1) = cc(j+1)*yvar2(j) + dd(j+1)*yvar2(j+1)
    xvar3(j+1) = cc(j+1)*yvar3(j) + dd(j+1)*yvar3(j+1)
    xvar4(j+1) = cc(j+1)*yvar4(j) + dd(j+1)*yvar4(j+1)

    do j0 = 2, jmax0
       j = j0*2 - 1
       !! [xvar;xvar] = [cc,dd;cc,dd]*[yvar;yvar] - [aa,bb;aa,bb]*[xvar;xvar](j0-1)
       xvar1(j)   = cc(j)*yvar1(j)   + dd(j)*yvar1(j+1)   - aa(j)*xvar1(j-2)   - bb(j)*xvar1(j-1)
       xvar2(j)   = cc(j)*yvar2(j)   + dd(j)*yvar2(j+1)   - aa(j)*xvar2(j-2)   - bb(j)*xvar2(j-1)
       xvar3(j)   = cc(j)*yvar3(j)   + dd(j)*yvar3(j+1)   - aa(j)*xvar3(j-2)   - bb(j)*xvar3(j-1)
       xvar4(j)   = cc(j)*yvar4(j)   + dd(j)*yvar4(j+1)   - aa(j)*xvar4(j-2)   - bb(j)*xvar4(j-1)
       xvar1(j+1) = cc(j+1)*yvar1(j) + dd(j+1)*yvar1(j+1) - aa(j+1)*xvar1(j-2) - bb(j+1)*xvar1(j-1)
       xvar2(j+1) = cc(j+1)*yvar2(j) + dd(j+1)*yvar2(j+1) - aa(j+1)*xvar2(j-2) - bb(j+1)*xvar2(j-1)
       xvar3(j+1) = cc(j+1)*yvar3(j) + dd(j+1)*yvar3(j+1) - aa(j+1)*xvar3(j-2) - bb(j+1)*xvar3(j-1)
       xvar4(j+1) = cc(j+1)*yvar4(j) + dd(j+1)*yvar4(j+1) - aa(j+1)*xvar4(j-2) - bb(j+1)*xvar4(j-1)
    end do

  !! Backward substitution

    do j0 = jmax0-1,1,-1
       j = j0*2 - 1
       !! [xvar;xvar] = [xvar;xvar] - [cc,dd;cc,dd]*[xvar;xvar](j0+1)
       xvar1(j)   = xvar1(j)   - ee(j)*xvar1(j+2)   - ff(j)*xvar1(j+3)
       xvar2(j)   = xvar2(j)   - ee(j)*xvar2(j+2)   - ff(j)*xvar2(j+3)
       xvar3(j)   = xvar3(j)   - ee(j)*xvar3(j+2)   - ff(j)*xvar3(j+3)
       xvar4(j)   = xvar4(j)   - ee(j)*xvar4(j+2)   - ff(j)*xvar4(j+3)
       xvar1(j+1) = xvar1(j+1) - ee(j+1)*xvar1(j+2) - ff(j+1)*xvar1(j+3)
       xvar2(j+1) = xvar2(j+1) - ee(j+1)*xvar2(j+2) - ff(j+1)*xvar2(j+3)
       xvar3(j+1) = xvar3(j+1) - ee(j+1)*xvar3(j+2) - ff(j+1)*xvar3(j+3)
       xvar4(j+1) = xvar4(j+1) - ee(j+1)*xvar4(j+2) - ff(j+1)*xvar4(j+3)
    end do
    
  end subroutine pentadiagonal_4_block_solve
  

  subroutine pentadiagonal_4_cyclic_solve( nn_cyclic, imax,                   &
   & aa1, bb1, cc1, dd1, ee1, ff1, pp2, qq2, rr2, ss2, jmax2_iter, jptr_iter, &
   & y1, yvar, xvar )
   integer,intent(in) :: nn_cyclic
   integer,intent(in) :: imax
    real(8),intent(in) :: aa1(imax)
    real(8),intent(in) :: bb1(imax)
    real(8),intent(in) :: cc1(imax)
    real(8),intent(in) :: dd1(imax)
    real(8),intent(in) :: ee1(imax)
    real(8),intent(in) :: ff1(imax)
    real(8),intent(in) :: pp2(imax)
    real(8),intent(in) :: qq2(imax)
    real(8),intent(in) :: rr2(imax)
    real(8),intent(in) :: ss2(imax)
    integer,intent(in) :: jmax2_iter(0:nn_cyclic)
    integer,intent(in) :: jptr_iter(nn_cyclic)
    real(8),intent(inout) :: y1(imax,4)
    real(8),intent(inout) :: yvar(imax+1,4)
    real(8),intent(inout) :: xvar(imax+1,4)
    integer :: jmax0,jmax1,jmax2,n,jptr,k
 
    if ( nn_cyclic == 0 ) then
      jptr  = 1
      jmax1 = 0
      jmax2 = jmax2_iter(0)
    end if

    do n=1,nn_cyclic
      jmax0 = jmax2_iter(n-1)
      jmax2 = jmax2_iter(n)
      jmax1 = jmax0 - jmax2
      jptr = jptr_iter(n)
      do k=1,4
        call pentadiagonal_4_cyclic_forward                        &
         &( jmax0, jmax1, jmax2, cc1(jptr), dd1(jptr),             &
         &  pp2(jptr), qq2(jptr), rr2(jptr), ss2(jptr), yvar(1,k), &
         &  y1(jptr,k), xvar(1,k)  ) 
        yvar(:,k) = xvar(:,k)
      end do
    end do

    call pentadiagonal_4_block_solve                                     &
     &( jmax2, aa1(jptr+jmax1*2), bb1(jptr+jmax1*2), cc1(jptr+jmax1*2),  &
     &  dd1(jptr+jmax1*2), ee1(jptr+jmax1*2), ff1(jptr+jmax1*2),         &
     &  yvar(1,1), yvar(1,2), yvar(1,3), yvar(1,4),                      &
     &  xvar(1,1), xvar(1,2), xvar(1,3), xvar(1,4)  )

    do n=nn_cyclic,1,-1
      jmax0 = jmax2_iter(n-1)
      jmax2 = jmax2_iter(n)
      jmax1 = jmax0 - jmax2
      jptr = jptr_iter(n)
      do k=1,4
        call pentadiagonal_4_cyclic_backward              &
         &( jmax0, jmax1, jmax2, aa1(jptr), bb1(jptr),    &
         &  ee1(jptr), ff1(jptr), y1(jptr,k), xvar(1,k),  &
         &  yvar(1,k) )
        xvar(:,k) = yvar(:,k)
      end do
    end do
   
  end subroutine pentadiagonal_4_cyclic_solve


  subroutine pentadiagonal_4_cyclic_forward  &
   &( jmax0, jmax1, jmax2, cc1, dd1, pp2, qq2, rr2, ss2, yvar, &
   &  y1, yvar2  )
    integer,intent(in) :: jmax0
    integer,intent(in) :: jmax1
    integer,intent(in) :: jmax2
    real(8),intent(in) :: cc1(jmax1*2)
    real(8),intent(in) :: dd1(jmax1*2)
    real(8),intent(in) :: pp2(jmax2*2)
    real(8),intent(in) :: qq2(jmax2*2)
    real(8),intent(in) :: rr2(jmax2*2)
    real(8),intent(in) :: ss2(jmax2*2)
    real(8),intent(in) :: yvar(jmax0*2)
    real(8),intent(out) :: y1(jmax1*2)
    real(8),intent(out) :: yvar2(jmax2*2)
    integer :: j1, j, jj
    
    if ( jmax1 > jmax2 ) then
       do jj = 1, jmax2
          j  = jj*4 - 1
          j1 = jj*2 - 1
          
          !! [y1,y1]T = [cc1,dd1;cc1,dd1]*[yvar,yvar]T(jj-1)
          y1(j1)   = cc1(j1)*yvar(j-2)   + dd1(j1)*yvar(j-1)
          y1(j1+1) = cc1(j1+1)*yvar(j-2) + dd1(j1+1)*yvar(j-1)
          
          !! [yvar2,yvar2]T = [yvar,yvar]T - [pp2,qq2;pp2,qq2]*[yvar,yvar]T(jj-1)
          !!                               - [rr2,ss2;rr2,ss2]*[yvar,yvar]T(jj+1)
          yvar2(j1)   = yvar(j) - pp2(j1)*yvar(j-2) - qq2(j1)*yvar(j-1) &
           &                    - rr2(j1)*yvar(j+2) - ss2(j1)*yvar(j+3)
          yvar2(j1+1) = yvar(j+1) - pp2(j1+1)*yvar(j-2) - qq2(j1+1)*yvar(j-1) &
           &                      - rr2(j1+1)*yvar(j+2) - ss2(j1+1)*yvar(j+3)
       end do
       
       jj = jmax1
       j  = jj*4 - 1
       j1 = jj*2 - 1
          
       !! [y1,y1]T = [cc1,dd1;cc1,dd1]*[yvar,yvar]T(jj-1)
       y1(j1)   = cc1(j1)*yvar(j-2)   + dd1(j1)*yvar(j-1)
       y1(j1+1) = cc1(j1+1)*yvar(j-2) + dd1(j1+1)*yvar(j-1)
    else
       do jj = 1, jmax2-1
          j  = jj*4 - 1
          j1 = jj*2 - 1
          
          !! [y1,y1]T = [cc1,dd1;cc1,dd1]*[yvar,yvar]T(jj-1)
          y1(j1)   = cc1(j1)*yvar(j-2)   + dd1(j1)*yvar(j-1)
          y1(j1+1) = cc1(j1+1)*yvar(j-2) + dd1(j1+1)*yvar(j-1)
          
          !! [yvar2,yvar2]T = [yvar,yvar]T - [pp2,qq2;pp2,qq2]*[yvar,yvar]T(jj-1)
          !!                               - [rr2,ss2;rr2,ss2]*[yvar,yvar]T(jj+1)
          yvar2(j1)   = yvar(j) - pp2(j1)*yvar(j-2) - qq2(j1)*yvar(j-1) &
           &                    - rr2(j1)*yvar(j+2) - ss2(j1)*yvar(j+3)
          yvar2(j1+1) = yvar(j+1) - pp2(j1+1)*yvar(j-2) - qq2(j1+1)*yvar(j-1) &
           &                      - rr2(j1+1)*yvar(j+2) - ss2(j1+1)*yvar(j+3)
       end do
       jj = jmax2
       j  = jj*4 - 1
       j1 = jj*2 - 1
          
       !! [y1,y1]T = [cc1,dd1;cc1,dd1]*[yvar,yvar]T(jj-1)
       y1(j1)   = cc1(j1)*yvar(j-2)   + dd1(j1)*yvar(j-1)
       y1(j1+1) = cc1(j1+1)*yvar(j-2) + dd1(j1+1)*yvar(j-1)
          
       !! [yvar2,yvar2]T = [yvar,yvar]T - [pp2,qq2;pp2,qq2]*[yvar,yvar]T(jj-1)
       yvar2(j1)   = yvar(j)   - pp2(j1)*yvar(j-2)   - qq2(j1)*yvar(j-1)
       yvar2(j1+1) = yvar(j+1) - pp2(j1+1)*yvar(j-2) - qq2(j1+1)*yvar(j-1)
    end if
    
!    write(6,*) "jmax1=",jmax1
!    write(6,*) "yvar(:)=",yvar(:)
!    write(6,*) "y1(:)=",y1(:)
!    write(6,*) "yvar2(:)=",yvar2(:)
    
  end subroutine pentadiagonal_4_cyclic_forward


  subroutine pentadiagonal_4_cyclic_backward                   &
   &( jmax0, jmax1, jmax2, aa1, bb1, ee1, ff1, y1, xvar2,  &
   &  xvar )
    integer,intent(in) :: jmax0
    integer,intent(in) :: jmax1
    integer,intent(in) :: jmax2
    real(8),intent(in) :: aa1(jmax1*2)
    real(8),intent(in) :: bb1(jmax1*2)
    real(8),intent(in) :: ee1(jmax1*2)
    real(8),intent(in) :: ff1(jmax1*2)
    real(8),intent(in) :: y1(jmax1*2)
    real(8),intent(in) :: xvar2(jmax2*2)
    real(8),intent(out) :: xvar(jmax0*2)
    integer :: j1, j, jj
    
    jj = 1
    j  = jj*4 - 1
    j1 = jj*2 - 1
    
       !! [xvar,xvar]T(jj-1) = [y1,y1]T - [ee1,ff1;ee2,ff2]*[xvar2,xvar2]T
       xvar(j-2) = y1(j1) - ee1(j1)*xvar2(j1)   - ff1(j1)*xvar2(j1+1)
       xvar(j-1) = y1(j1+1) - ee1(j1+1)*xvar2(j1)   - ff1(j1+1)*xvar2(j1+1)
       
       !! [xvar,xvar]T = [xvar2,xvar2]T 
       xvar(j)   = xvar2(j1)
       xvar(j+1) = xvar2(j1+1)

    do jj = 2, jmax2
       j  = jj*4 - 1
       j1 = jj*2 - 1
             
       !! [xvar,xvar]T(jj-1) = [y1,y1]T - [aa1,bb1;aa1,bb2]*[xvar2,xvar2]T(jj-1)
       !!                         - [ee1,ff1;ee2,ff2]*[xvar2,xvar2]T
       xvar(j-2) = y1(j1) - aa1(j1)*xvar2(j1-2) - bb1(j1)*xvar2(j1-1) &
        &                 - ee1(j1)*xvar2(j1)   - ff1(j1)*xvar2(j1+1)
       xvar(j-1) = y1(j1+1) - aa1(j1+1)*xvar2(j1-2) - bb1(j1+1)*xvar2(j1-1) &
        &                   - ee1(j1+1)*xvar2(j1)   - ff1(j1+1)*xvar2(j1+1)
       
       !! [xvar,xvar]T = [xvar2,xvar2]T 
       xvar(j)   = xvar2(j1)
       xvar(j+1) = xvar2(j1+1)
    end do
    
    if ( jmax1 > jmax2 ) then
       jj = jmax1
       j  = jj*4 - 1
       j1 = jj*2 - 1
       
       !! [xvar,xvar]T(jj-1) = [y1,y1]T - [aa1,bb1;aa1,bb2]*[xvar2,xvar2]T(jj-1)
       xvar(j-2) = y1(j1) - aa1(j1)*xvar2(j1-2) - bb1(j1)*xvar2(j1-1)
       xvar(j-1) = y1(j1+1) - aa1(j1+1)*xvar2(j1-2) - bb1(j1+1)*xvar2(j1-1)
    end if
    
  end subroutine pentadiagonal_4_cyclic_backward

end module pentadiagonal_4

