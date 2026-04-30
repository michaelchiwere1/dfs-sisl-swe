module tridiagonal_cmplx

  implicit none
  
  private
  public :: tridiagonal_cmplx__allocate, tridiagonal_cmplx__ini, tridiagonal_cmplx__solve, tridiagonal_cmplx__test

!  integer,parameter :: jcn_cyclic = 0  !! LU decomposition of tri-diagonal matrix
  integer,parameter :: jcn_cyclic = 1  !! Modified LU decomposition of tri-diagonal matrix
!  integer,parameter :: jcn_cyclic = 2  !! Cyclic reduction for block tri-diagonal matrix

!xx  call main_sub

contains

  subroutine tridiagonal_cmplx__test
  
    integer,parameter :: imax=1000
    integer,parameter :: ilen=imax+1

!    complex(8) :: dat(ilen,5)
!    integer :: iww(ilen)
    complex(8),allocatable :: dat(:,:)
    integer,allocatable :: iww(:)
    complex(8) :: yy(ilen)
    complex(8) :: yvar(ilen)
    complex(8) :: xvar(ilen)
    integer :: n
    
    call tridiagonal_cmplx__allocate( imax, ilen, dat, iww )!IN,IN,IO,OUT
    
    dat(:,1) = 1.0d0
    dat(:,2) = 2.0d0
    dat(:,3) = 1.0d0
    
    call tridiagonal_cmplx__ini( imax, ilen, dat, iww )!IN,IN,IO,OUT
    
    n=1
!    do n=1,100000
      yvar(1) = 3.0d0
      yvar(2:imax-1) = 4.0d0
      yvar(imax) = 3.0d0
  
      call tridiagonal_cmplx__solve( imax, ilen, dat, iww, &!IN
       &                       yy, yvar, xvar )       !OUT,IO,OUT
      
!    end do
  
    write(6,*) "xvar(1:imax)=",xvar(1:imax)
  
  end subroutine tridiagonal_cmplx__test
  
  
  subroutine tridiagonal_cmplx__allocate( imax, ilen, dat, iww )
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    complex(8),intent(out),allocatable :: dat(:,:)
    integer,intent(out),allocatable :: iww(:)
    
    integer :: nn_tmp, nn_cyclic

    if ( ilen < imax ) then
      write(6,*) "Error: tridiagonal_cmplx__ini: ilen should be equal to or larger than imax."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
    
    if ( jcn_cyclic == 0 ) then
      allocate( dat(ilen,3) )
      allocate( iww(0) )
    
    else if ( jcn_cyclic == 1 ) then
      allocate( dat(ilen,4) )
      allocate( iww(0) )
    
    else if ( jcn_cyclic >= 2 ) then
      nn_tmp = int(log(imax+0.001d0)/log(2.0d0))
      nn_cyclic = max( nn_tmp-2, 0 )
      allocate( dat(ilen,5) )
      allocate( iww(2+2*nn_cyclic) )
      
      iww(1) = nn_cyclic
      write(6,*) "nn_cyclic=",nn_cyclic
    end if
  
  end subroutine tridiagonal_cmplx__allocate
  
  
  subroutine tridiagonal_cmplx__ini( imax, ilen, dat, iww )
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    complex(8),intent(inout) :: dat(ilen,5)
    integer,intent(out)   :: iww(ilen)
    
    integer :: nn_tmp, nn_cyclic

    if ( ilen < imax ) then
      write(6,*) "Error: tridiagonal_cmplx__ini: ilen should be equal to or larger than imax."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
    
    if ( jcn_cyclic == 0 ) then
      call simple_pre( imax, dat(1,1), dat(1,2), dat(1,3))
    else if ( jcn_cyclic == 1 ) then
      call simple_mod_pre( imax, dat(1,1), dat(1,2), dat(1,3), dat(1,4) )
    else if ( jcn_cyclic >= 2 ) then
      nn_cyclic = iww(1)
      call cyclic_pre( nn_cyclic, imax,      &!IN
       & dat(1,1), dat(1,2), dat(1,3),                   &!INOUT
       & dat(1,4), dat(1,5), iww(2), iww(3+nn_cyclic) )   !OUT
      write(6,*) "nn_cyclic=",nn_cyclic
    end if

  end subroutine tridiagonal_cmplx__ini
  
  
  subroutine tridiagonal_cmplx__solve( imax, ilen, dat, iww, &
   &                             yy, yvar, xvar )
    integer,intent(in) :: imax
    integer,intent(in) :: ilen
    complex(8),intent(in) :: dat(ilen,5)
    integer,intent(in) :: iww(ilen)
    complex(8),intent(inout) :: yy(ilen)
    complex(8),intent(inout) :: yvar(ilen)
    complex(8),intent(out)   :: xvar(ilen)
    
    integer :: nn_cyclic

    yvar(imax+1:ilen) = 0.0d0
    if ( ilen < imax ) then
      write(6,*) "Error: tridiagonal_cmplx__solve: ilen should be equal to or larger than imax."
      write(6,*) "       imax, ilen = ", imax, ilen
      stop 333
    end if
  
    if ( jcn_cyclic >= 2 ) then
      nn_cyclic = iww(1)
      call cyclic_solve( nn_cyclic, imax,                                &!IN
       & dat(1,1), dat(1,2), dat(1,3), dat(1,4), dat(1,5), iww(2), iww(3+nn_cyclic), &!IN
       & yy, yvar, xvar )                                                      !OUT,IO,OUT
    else if ( jcn_cyclic == 1 ) then
      call simple_mod_solve( imax, dat(1,1), dat(1,2), dat(1,3), dat(1,4), yy, yvar, xvar )
    else
      call simple_solve( imax, dat(1,1), dat(1,2), dat(1,3), yvar, xvar )
    end if
    xvar(imax+1:ilen) = 0.0d0
  
  end subroutine tridiagonal_cmplx__solve


  subroutine simple_pre( jmax, aa, bb, cc )
    
    integer,intent(in) :: jmax
    complex(8),intent(inout) :: aa(jmax)
    complex(8),intent(inout) :: bb(jmax)
    complex(8),intent(inout) :: cc(jmax)
    integer :: j

    aa(1) = 0.0d0
    bb(1) = 1.0d0/bb(1)
    cc(1) = cc(1)*bb(1)
    do j=2,jmax
       bb(j) = 1.0d0/( bb(j) - aa(j)*cc(j-1) )
       aa(j) = aa(j)*bb(j)
       cc(j) = cc(j)*bb(j)
    end do
    cc(jmax) = 0.0d0
    
!    write(6,*) "aa(:)=",aa(:)
!    write(6,*) "1/bb(:)=",1.0d0/bb(:)
!    write(6,*) "cc(:)=",cc(:)

  end subroutine simple_pre
  
  
  
  subroutine simple_solve( jmax, aa, bb, cc, yvar, xvar )
    integer,intent(in) :: jmax
    complex(8),intent(in) :: aa(jmax)
    complex(8),intent(in) :: bb(jmax)
    complex(8),intent(in) :: cc(jmax)
    complex(8),intent(in) :: yvar(jmax)
    complex(8),intent(out) :: xvar(jmax)
    integer :: j

    !! Forward substitution
    xvar(1) = yvar(1)*bb(1)
    do j=2,jmax
       xvar(j) = bb(j)*yvar(j) - aa(j)*xvar(j-1)
    end do
    
    !! Backward substitution
    do j=jmax-1,1,-1
       xvar(j) = xvar(j) - cc(j)*xvar(j+1)
    end do
    
  end subroutine simple_solve
  
  

  subroutine simple_mod_pre( jmax, aa, bb, cc, dd )
    
    integer,intent(in) :: jmax
    complex(8),intent(inout) :: aa(jmax)
    complex(8),intent(inout) :: bb(jmax)
    complex(8),intent(inout) :: cc(jmax)
    complex(8),intent(inout) :: dd(jmax)
    integer :: j,j2,j2max

    call simple_pre( jmax, aa, bb, cc )
    
    do j2=1,(jmax-1)/2
      j=j2*2
      dd(j2)  = aa(j+1)*bb(j)
      aa(j+1) = -aa(j+1)*aa(j)
    end do
    
    do j2=1,(jmax-1)/2
      j=j2*2-1
      dd(jmax/2+j2) = cc(j)
      cc(j) = -cc(j)*cc(j+1)
    end do
    
!    write(6,*) "aa(:)=",aa(:)
!    write(6,*) "1/bb(:)=",1.0d0/bb(:)
!    write(6,*) "cc(:)=",cc(:)

  end subroutine simple_mod_pre
  

  
  subroutine simple_mod_solve( jmax, aa, bb, cc, dd, work, yvar, xvar )
    integer,intent(in) :: jmax
    complex(8),intent(in) :: aa(jmax)
    complex(8),intent(in) :: bb(jmax)
    complex(8),intent(in) :: cc(jmax)
    complex(8),intent(in) :: dd(jmax)
    complex(8),intent(out) :: work(jmax)
    complex(8),intent(in) :: yvar(jmax)
    complex(8),intent(out) :: xvar(jmax)
    integer :: j,j2

!    !! Forward substitution
!    work(1) = yvar(1)*bb(1)
!    do j=2,jmax
!       work(j) = bb(j)*yvar(j) - aa(j)*work(j-1)
!    end do

    !! Forward substitution
    work(1) = yvar(1)*bb(1)
    do j2=1,(jmax-1)/2
       j=j2*2
       work(j)   = bb(j)*yvar(j) - aa(j)*work(j-1)
       work(j+1) = bb(j+1)*yvar(j+1) - dd(j2)*yvar(j) - aa(j+1)*work(j-1)
    end do
    if ( mod(jmax,2) == 0 ) then
      j=jmax
       work(j)   = bb(j)*yvar(j) - aa(j)*work(j-1)
    end if
    
!    !! Backward substitution
!    xvar(jmax) = work(jmax)
!    do j=jmax-1,1,-1
!       xvar(j) = work(j) - cc(j)*xvar(j+1)
!    end do
    
    !! Backward substitution
    xvar(jmax) = work(jmax)
    if ( mod(jmax,2) == 0 ) then
      j=jmax-1
      xvar(j) = work(j) - cc(j)*xvar(j+1)
    end if
    do j2=(jmax-1)/2,1,-1
      j=j2*2-1
      xvar(j+1) = work(j+1) - cc(j+1)*xvar(j+2)
      xvar(j)   = work(j) - dd(jmax/2+j2)*work(j+1) - cc(j)*xvar(j+2)
    end do
   
  end subroutine simple_mod_solve
  
  
  
  subroutine cyclic_pre( nn_cyclic, imax, aa, bb, cc, qq, rr, jmax2_iter, jptr_iter )
  
    integer,intent(inout) :: nn_cyclic
    integer,intent(in) :: imax
    complex(8),intent(inout) :: aa(imax)
    complex(8),intent(inout) :: bb(imax)
    complex(8),intent(inout) :: cc(imax)
    complex(8),intent(out)   :: qq(imax)
    complex(8),intent(out)   :: rr(imax)
    integer,intent(out)  :: jmax2_iter(0:nn_cyclic)
    integer,intent(out)  :: jptr_iter(nn_cyclic)
    complex(8) :: aa1(imax)
    complex(8) :: bb1(imax)
    complex(8) :: cc1(imax)
    complex(8) :: aa2(imax)
    complex(8) :: bb2(imax)
    complex(8) :: cc2(imax)
    integer :: jmax0,jmax1,jmax2,j,n,jptr,jptr0
    
!    aa(1) = 0.0d0
!    cc(imax) = 0.0d0
    qq(:)=0.0d0
    rr(:)=0.0d0
    
    jmax0 = imax
    jmax2_iter(0) = jmax0
    jptr = 1

    if ( nn_cyclic > 0 ) then

       do n=1,nn_cyclic
       
          jmax1 = (jmax0+1)/2
          jmax2 = jmax0 - jmax1
          jmax2_iter(n) = jmax2
          jptr_iter(n) = jptr
       
          if ( n == 1 ) then
             call set_array_cyclic                                      &
              &( jmax0, jmax1, jmax2, aa, bb, cc,                       &
              &  aa1(jptr), bb1(jptr), cc1(jptr),                       &
              &  qq(jptr), rr(jptr), aa2(jptr), bb2(jptr), cc2(jptr) )
          else
             call set_array_cyclic                                        &
              &( jmax0, jmax1, jmax2, aa2(jptr0), bb2(jptr0), cc2(jptr0), &
              &  aa1(jptr), bb1(jptr), cc1(jptr),                         &
              &  qq(jptr), rr(jptr), aa2(jptr), bb2(jptr), cc2(jptr) )
          end if
          jmax0 = jmax2
          jptr0 = jptr
          jptr = jptr + jmax1
       end do
    
    
       do j=1,jptr-1
          aa(j) = aa1(j)
          bb(j) = bb1(j)
          cc(j) = cc1(j)
       end do
       
       do j=jptr,jptr+jmax2-1
          aa(j) = aa2(j-jmax1)
          bb(j) = bb2(j-jmax1)
          cc(j) = cc2(j-jmax1)
       end do
    end if

    call simple_pre( jmax0, aa(jptr), bb(jptr), cc(jptr) )
    
  
  end subroutine cyclic_pre
  
  
  subroutine set_array_cyclic( jmax0, jmax1, jmax2, aa, bb, cc,        &
   &                          aa1, bb1, cc1, qq2, rr2, aa2, bb2, cc2 )
   
    integer,intent(in) :: jmax0, jmax1, jmax2
    complex(8),intent(inout) :: aa(jmax0)
    complex(8),intent(inout) :: bb(jmax0)
    complex(8),intent(inout) :: cc(jmax0)
    complex(8),intent(out) :: aa1(jmax1)
    complex(8),intent(out) :: bb1(jmax1)
    complex(8),intent(out) :: cc1(jmax1)
    complex(8),intent(out) :: qq2(jmax2)
    complex(8),intent(out) :: rr2(jmax2)
    complex(8),intent(out) :: aa2(jmax2)
    complex(8),intent(out) :: bb2(jmax2)
    complex(8),intent(out) :: cc2(jmax2)
    integer :: j1,j

    do j1 = 1, jmax1
       j = j1*2
       bb1(j1) = 1.0d0/bb(j-1)
       aa1(j1) = aa(j-1)*bb1(j1)
       cc1(j1) = cc(j-1)*bb1(j1)
    end do
    
!    write(6,*) "jmax1,jmax0=",jmax1,jmax0
!    write(6,*) "aa(:)=",aa(:)
!    write(6,*) "aa2(1:jmax1)=",aa2(1:jmax1)
!    stop 333
    
    do j1 = 1, jmax2
       j = j1*2
       if ( j == jmax0 ) then
          qq2(j1) = aa(j)*bb1(j1)
          rr2(j1) = 0.0d0
          aa2(j1) = - aa(j)*aa1(j1)
          bb2(j1) = bb(j) - aa(j)*cc1(j1)
          cc2(j1) = 0.0d0
       else
          qq2(j1) = aa(j)*bb1(j1)
          rr2(j1) = cc(j)*bb1(j1+1)
          aa2(j1) = - aa(j)*aa1(j1)
          bb2(j1) = bb(j) - aa(j)*cc1(j1) - cc(j)*aa1(j1+1)
          cc2(j1) = - cc(j)*cc1(j1+1)
       end if
    end do
    
!    write(6,*) "array(:)=",array(:)
!    write(6,*) "a1(:)=",a1(:)
!    write(6,*) "a2(:)=",a2(:)
!    write(6,*) "a3(:)=",a3(:)
    
  end subroutine set_array_cyclic
  

  subroutine cyclic_solve( nn_cyclic, imax, aa1, bb1, cc1, qq2, rr2, jmax2_iter, jptr_iter, &
   &                                   y1, yvar, xvar )
    integer,intent(in) :: nn_cyclic
    integer,intent(in) :: imax
    complex(8),intent(in) :: aa1(imax)
    complex(8),intent(in) :: bb1(imax)
    complex(8),intent(in) :: cc1(imax)
    complex(8),intent(in) :: qq2(imax)
    complex(8),intent(in) :: rr2(imax)
    integer,intent(in) :: jmax2_iter(0:nn_cyclic)
    integer,intent(in) :: jptr_iter(nn_cyclic)
    complex(8),intent(out) :: y1(imax)
    complex(8),intent(inout) :: yvar(imax)
    complex(8),intent(out) :: xvar(imax)
    integer :: jmax0,jmax1,jmax2,n,jptr
    
    do n=1,nn_cyclic
       jmax0 = jmax2_iter(n-1)
       jmax2 = jmax2_iter(n)
       jmax1 = jmax0 - jmax2
       jptr = jptr_iter(n)
       
       call cyclic_forward                                 &
        &( jmax0, jmax1, jmax2, bb1(jptr), qq2(jptr), rr2(jptr), yvar, &
        &  y1(jptr), xvar  )
        
       yvar = xvar

!       n=n2*2
!       jmax0 = jmax_iter(n-1)
!       jmax1 = jmax_iter(n)
!       jmax2 = jmax0 - jmax1
!       jptr = jptr_iter(n)
!       call cyclic_forward                          &
!        &( jmax0, jmax1, jmax2, bb1(jptr), qq2(jptr), rr2(jptr), xvar, &
!        &  y1(jptr), yvar  )
    end do
 
    if ( nn_cyclic == 0 ) then
       jptr  = 1
       jmax1 = 0
       jmax2 = jmax2_iter(0)
    end if
    
    call simple_solve                                       &
     &( jmax2, aa1(jptr+jmax1), bb1(jptr+jmax1), cc1(jptr+jmax1), yvar, &
     &  xvar )
    
    do n=nn_cyclic,1,-1
       jmax0 = jmax2_iter(n-1)
       jmax2 = jmax2_iter(n)
       jmax1 = jmax0 - jmax2
       jptr = jptr_iter(n)
       call cyclic_backward                  &
        &( jmax0, jmax1, jmax2, aa1(jptr), cc1(jptr), y1(jptr), xvar, &
        &  yvar )
        
       xvar = yvar
        
!       n=n2*2-1
!       jmax0 = jmax_iter(n-1)
!       jmax1 = jmax_iter(n)
!       jmax2 = jmax0 - jmax1
!       jptr = jptr_iter(n)
!       call cyclic_backward                  &
!        &( jmax0, jmax1, jmax2, aa1(jptr), cc1(jptr), y1(jptr), yvar, &
!        &  xvar )
    end do

  end subroutine cyclic_solve


  subroutine cyclic_forward  &
   &( jmax0, jmax1, jmax2, bb1, qq2, rr2, yvar, &
   &  y1, yvar2  )
    integer,intent(in) :: jmax0
    integer,intent(in) :: jmax1
    integer,intent(in) :: jmax2
    complex(8),intent(in) :: bb1(jmax1)
    complex(8),intent(in) :: qq2(jmax2)
    complex(8),intent(in) :: rr2(jmax2)
    complex(8),intent(in) :: yvar(jmax0)
    complex(8),intent(out) :: y1(jmax1)
    complex(8),intent(out) :: yvar2(jmax2)
    integer :: j1, j
    
    if ( jmax1 > jmax2 ) then
       do j1 = 1, jmax2
          j = j1*2
          y1(j1) = yvar(j-1)*bb1(j1)
          yvar2(j1) = yvar(j) - yvar(j-1)*qq2(j1) - yvar(j+1)*rr2(j1)
       end do
       j1 = jmax1
       j = j1*2
       y1(j1) = yvar(j-1)*bb1(j1)
    else
       do j1 = 1, jmax2-1
          j = j1*2
          y1(j1) = yvar(j-1)*bb1(j1)
          yvar2(j1) = yvar(j) - yvar(j-1)*qq2(j1) - yvar(j+1)*rr2(j1)
       end do
       j1 = jmax2
       j = j1*2
       y1(j1) = yvar(j-1)*bb1(j1)
       yvar2(j1) = yvar(j) - yvar(j-1)*qq2(j1)
    end if
    
!    write(6,*) "imax1=",imax1
!    write(6,*) "a1(:)=",a1(:)
!    write(6,*) "yvar(:)=",yvar(:)
!    write(6,*) "y1(:)=",y1(:)
!    write(6,*) "yvar2(:)=",yvar2(:)
    
  end subroutine cyclic_forward


  subroutine cyclic_backward  &
   &( jmax0, jmax1, jmax2, aa1, cc1, y1, xvar2,  &
   &  xvar )
    integer,intent(in) :: jmax0
    integer,intent(in) :: jmax1
    integer,intent(in) :: jmax2
    complex(8),intent(in) :: aa1(jmax1)
    complex(8),intent(in) :: cc1(jmax1)
    complex(8),intent(in) :: y1(jmax1)
    complex(8),intent(in) :: xvar2(jmax2)
    complex(8),intent(out) :: xvar(jmax0)
    integer :: j1, j
    
    j1 = 1
    j = j1*2
    xvar(j-1) = y1(j1) - xvar2(j1)*cc1(j1)
    xvar(j)   = xvar2(j1)
    do j1 = 2, jmax2
       j = j1*2
       xvar(j-1) = y1(j1) - xvar2(j1-1)*aa1(j1) - xvar2(j1)*cc1(j1)
       xvar(j)   = xvar2(j1)
    end do
    if ( jmax1 > jmax2 ) then
       j1 = jmax1
       j = j1*2
       xvar(j-1) = y1(j1) - xvar2(j1-1)*aa1(j1)
    end if
    
  end subroutine cyclic_backward

end module tridiagonal_cmplx

