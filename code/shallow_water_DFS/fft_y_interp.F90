module fft_y_interp

  use com_dfs, only : jcn_dfs, jcn_grid, NNUM, N_TRUNC_M0_DFS, N_TRUNC_M1_DFS,     &
   &   NNUMHF1, NNUMHF2, NNUM1_M_DFS, NNUM2_M_DFS
  use prm_phconst, only : PI
  use tridiagonal, only : tridiagonal__allocate, tridiagonal__ini, tridiagonal__solve
  use pentadiagonal, only : pentadiagonal__allocate, pentadiagonal__ini, pentadiagonal__solve
  use e_time, only : e_time__start, e_time__end
   
!--------------------------------------------------------------
  implicit none
!--------------------------------------------------------------
  private
  public :: fft_y__ini_interp, fft_y__g2w_interp, fft_y__g2w_uv_interp, &
   &        fft_y__w2g_interp, fft_y__w2g_dy_interp, fft_y__w2g_uv_interp
  public :: fft_y__g2w_cheong, fft_y__w2g_cheong, fft_y__w2g_dy_cheong
  public :: fft_y__g2w_cheonglike, fft_y__w2g_cheonglike
  public :: fft_y__g2w_orszag, fft_y__w2g_orszag
  public :: grid2sin, sin2grid, grid2cos, cos2grid
  public :: sinsin2none
!
!  integer,parameter :: nw = 12500
!  real(8),save :: w(0:nw)
!
!  integer,save :: ip(0:100)
!
  integer,save :: IMAX
  integer,save :: JMAX
!  integer,save :: NNUMHF1
!  integer,save :: NNUMHF2
  integer,save :: MMAX

  real(8),allocatable :: DLAT
  real(8),allocatable :: ALAT(:)
  real(8),allocatable :: COSLAT(:)
  real(8),allocatable :: COSLAT_INV(:)

  real(8),allocatable :: WSAVE_SIN(:)
  real(8),allocatable :: WSAVE_COS(:)
!$OMP threadprivate(WSAVE_SIN,WSAVE_COS)

  type :: type_diagonal
    private
    real(8),allocatable :: d(:,:)
    integer,allocatable :: iww(:)
  end type
  type(type_diagonal),allocatable :: md1(:)
  type(type_diagonal),allocatable :: md2(:)

  real(8),save,allocatable :: aa1(:,:)
  real(8),save,allocatable :: bb1(:,:)
  real(8),save,allocatable :: cc1(:,:)
  real(8),save,allocatable :: dd1(:,:)
  real(8),save,allocatable :: ee1(:,:)

  real(8),save,allocatable :: aa2(:,:)
  real(8),save,allocatable :: bb2(:,:)
  real(8),save,allocatable :: cc2(:,:)
  real(8),save,allocatable :: dd2(:,:)
  real(8),save,allocatable :: ee2(:,:)
!
!********************************************************************
contains
!********************************************************************


subroutine fft_y__ini_interp(imax_in,jmax_in,mmax_in,alat_out,weight_out)
!
  integer,intent(in) :: imax_in
  integer,intent(in) :: jmax_in
  integer,intent(in) :: mmax_in
  real(8),intent(out) :: alat_out(jmax_in)
  real(8),intent(out) :: weight_out(jmax_in)
  real(8) :: work(jmax_in+2)
  real(8) :: www(jmax_in+2)
  real(8) :: w, an
  integer :: m, j, j2
!
  write(6,*) 'fft_y: initialization start.'
!
  IMAX = imax_in
  JMAX = jmax_in
!  NNUMHF1 = (JMAX+1)/2
!  NNUMHF2 = JMAX/2  !! NNUMHF1+NNUMHF2=JMAX
  MMAX = mmax_in
  
  DLAT = PI/(jmax-jcn_grid)
  
  allocate( ALAT(JMAX) )
  allocate( COSLAT(JMAX) )
  allocate( COSLAT_INV(JMAX) )

!$OMP PARALLEL default(SHARED), private(j)
 !$OMP DO schedule(STATIC)
  do j=1,(JMAX+1)/2
     ALAT(j) = DLAT*( (JMAX+1)*0.5d0 - j ) !## pi/2 >= alat >= -pi/2
  end do
 !$OMP END DO
 !$OMP DO schedule(STATIC)
  do j=1,JMAX/2
     ALAT(JMAX+1-j) = -ALAT(j)
  end do
 !$OMP END DO
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
     alat_out(j) = ALAT(j)
     COSLAT(j) = cos(ALAT(j))
     COSLAT_INV(j) = 1.0d0/COSLAT(j)
  end do
 !$OMP END DO
!$OMP END PARALLEL
  
!$OMP PARALLEL default(SHARED)
  allocate( WSAVE_SIN(3*JMAX+15) )
  allocate( WSAVE_COS(3*JMAX+15) )
  
  if ( jcn_grid == 0 ) then
    call dsinqi( JMAX, WSAVE_SIN )
    call dcosqi( JMAX, WSAVE_COS )
  else if ( jcn_grid == 1 ) then
    call dsinti( JMAX-2, WSAVE_SIN )
    call dcosti( JMAX, WSAVE_COS )
  else if ( jcn_grid == -1 ) then
    call dsinti( JMAX, WSAVE_SIN )
    call dcosti( JMAX+2, WSAVE_COS )
  else
    write(0,*) "Error: jcn_grid should be 0 or 1 or -1."
    stop 999
  end if
!$OMP END PARALLEL
!
!$OMP PARALLEL default(SHARED), private(j,j2,an,w,work,www)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
    work(:) = 0.0d0
    work(j) = 1.0d0
    call grid2cos( work, www )
!    write(6,*) "j,data(:)=",j,data(:)
    w = 0.0d0
!xx    do j2=1,JMAX,2
    do j2=1,NNUM,2
      an=j2-1
      w = w + work(j2)/(1.0d0-an*an)
    end do
    weight_out(j) = w
  end do
 !$OMP END DO
 !$OMP DO schedule(STATIC)
  do j=1,JMAX/2
    w = ( weight_out(j) + weight_out(JMAX-j+1) )*0.5d0
    weight_out(j) = w
    weight_out(JMAX-j+1) = w
  end do
 !$OMP END DO
!$OMP END PARALLEL

!  write(6,*) "JMAX,sum(weight_out(1:JMAX))=",JMAX,sum(weight_out(1:JMAX))

  weight_out(:) = weight_out(:)/sum(weight_out(:))  
!  weight_out(:) = 2.0d0*weight_out(:)/sum(weight_out(:))  

  write(6,*) "JMAX,sum(weight_out(1:JMAX))=",JMAX,sum(weight_out(1:JMAX))
    
  call none2sinsin_pre
    
!
end subroutine fft_y__ini_interp


!********************************************************************


subroutine fft_y__g2w_interp &
 &( data,             &!IN
 &  qdata      )       !INOUT
!
  real(8),intent(inout) :: data(2,0:imax/2-1,jmax)
  complex(8),intent(out) :: qdata(NNUM,0:MMAX)
!
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(21,"fft_y__g2w")
!
  if ( jcn_dfs >= 1 ) then
  
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,JMAX
           work1(j) = data(1,m,j)
           work2(j) = data(2,m,j)
        end do
        if ( m == 0 ) then             !! m = 0
           call grid2cos( work1, www )
           work2(:) = 0.0d0
        else if ( mod(m,2) == 1 ) then !! Odd m
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        else                           !! Even m
!           call grid2cos( work1 )
!           call grid2cos( work2 )
           call grid2cos_polar0( work1, www )
           call grid2cos_polar0( work2, www )
        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
     call none2sinsin( qdata )
     
  else
  
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        if ( m == 0 ) then             !! m = 0
           do j=1,JMAX
              work1(j) = data(1,m,j)
           end do
           call grid2cos( work1, www )
           work2(:) = 0.0d0
        else if ( mod(m,2) == 1 ) then !! Odd m
           do j=1,JMAX
              work1(j) = data(1,m,j)
              work2(j) = data(2,m,j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        else                           !! Even m   
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat_inv(j)
              work2(j) = data(2,m,j)*coslat_inv(j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL

  end if
  
  call e_time__end(21,"fft_y__g2w")
!
end subroutine fft_y__g2w_interp


!********************************************************************

subroutine fft_y__g2w_cheong &
 &( data,             &!IN
 &  qdata      )       !INOUT
!
  real(8),intent(inout) :: data(2,0:imax/2-1,jmax)
  complex(8),intent(out) :: qdata(NNUM,0:MMAX)
!
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(21,"fft_y__g2w")
  
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        if ( m == 0 ) then             !! m = 0
           do j=1,JMAX
              work1(j) = data(1,m,j)
           end do
           call grid2cos( work1, www )
           work2(:) = 0.0d0
        else if ( mod(m,2) == 1 ) then !! Odd m
           do j=1,JMAX
              work1(j) = data(1,m,j)
              work2(j) = data(2,m,j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        else                           !! Even m   
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat_inv(j)
              work2(j) = data(2,m,j)*coslat_inv(j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  
  call e_time__end(21,"fft_y__g2w")
!
end subroutine fft_y__g2w_cheong


!********************************************************************


subroutine fft_y__g2w_cheonglike &
 &( data,             &!IN
 &  qdata      )       !INOUT
!
  real(8),intent(inout) :: data(2,0:imax/2-1,jmax)
  complex(8),intent(out) :: qdata(NNUM,0:MMAX)
!
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(21,"fft_y__g2w")

!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        if ( m == 0 ) then              !! m = 0
           do j=1,JMAX
              work1(j) = data(1,m,j)
           end do
           call grid2cos( work1, www )
           work2(:) = 0.0d0
        else if ( m == 1 ) then         !! m = 1
           do j=1,JMAX
              work1(j) = data(1,m,j)
              work2(j) = data(2,m,j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        else if ( mod(m,2) == 0 ) then  !! Even m   
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat_inv(j)
              work2(j) = data(2,m,j)*coslat_inv(j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        else
           do j=1,JMAX
!              work1(j) = data(1,m,j)*coslat_inv(j)*coslat_inv(j)
!              work2(j) = data(2,m,j)*coslat_inv(j)*coslat_inv(j)
              work1(j) = data(1,m,j)
              work2(j) = data(2,m,j)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  
  call e_time__end(21,"fft_y__g2w")
!
end subroutine fft_y__g2w_cheonglike



!********************************************************************


subroutine fft_y__g2w_orszag &
 &( data,             &!IN
 &  qdata      )       !INOUT
!
  real(8),intent(inout) :: data(2,0:imax/2-1,jmax)
  complex(8),intent(out) :: qdata(NNUM,0:MMAX)
!
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(21,"fft_y__g2w")

!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        if ( m == 0 ) then             !! m = 0
           do j=1,JMAX
              work1(j) = data(1,m,j)
           end do
           call grid2cos( work1, www )
           work2(:) = 0.0d0
        else if ( m == 1 ) then !! m=1
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat_inv(j)
              work2(j) = data(2,m,j)*coslat_inv(j)
           end do
           call grid2cos( work1, www )
           call grid2cos( work2, www )
        else if ( mod(m,2) == 0 ) then !! Even m
           do j=1,JMAX
              work1(j) = data(1,m,j)
              work2(j) = data(2,m,j)
           end do
           call grid2cos( work1, www )
           call grid2cos( work2, www )
           work1(NNUM-1) = -sum(work1(NNUM-3:1:-2))
           work1(NNUM)   = -sum(work1(NNUM-2:1:-2))
           work2(NNUM-1) = -sum(work2(NNUM-3:1:-2))
           work2(NNUM)   = -sum(work2(NNUM-2:1:-2))
        else                           !! Odd m   
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat_inv(j)
              work2(j) = data(2,m,j)*coslat_inv(j)
           end do
           call grid2cos( work1, www )
           call grid2cos( work2, www )
           work1(NNUM-1) = -sum(work1(NNUM-3:1:-2))
           work1(NNUM)   = -sum(work1(NNUM-2:1:-2))
           work2(NNUM-1) = -sum(work2(NNUM-3:1:-2))
           work2(NNUM)   = -sum(work2(NNUM-2:1:-2))

!           write(6,*) 'work1(1:NNUM)=',work1(1:NNUM)
!           write(6,*) 'work1(NNUM-1:1:-2)=',work1(NNUM-1:1:-2)
!           write(6,*) 'work1(NNUM  :1:-2)=',work1(NNUM  :1:-2)
!           write(6,*) 'sum(work1(NNUM-1:1:-2))=',sum(work1(NNUM-1:1:-2))
!           write(6,*) 'sum(work1(NNUM  :1:-2))=',sum(work1(NNUM  :1:-2))
!           write(6,*) 

        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  
  call e_time__end(21,"fft_y__g2w")
!
end subroutine fft_y__g2w_orszag


!********************************************************************


subroutine fft_y__g2w_uv_interp &
 &( data,             &!IN
 &  qdata      )       !INOUT
!
  real(8),intent(inout) :: data(2,0:IMAX/2-1,JMAX)
  complex(8),intent(out) :: qdata(NNUM,0:MMAX)
!
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(22,"fft_y__g2w_uv")
!
  if ( jcn_dfs >= 1 ) then
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,JMAX
           work1(j) = data(1,m,j)
           work2(j) = data(2,m,j)
        end do
        if ( m == 0 ) then             !! m = 0
           call grid2sin( work1, www )
           work2(:) = 0.0d0
        else if ( m == 1 ) then
           call grid2cos( work1, www )
           call grid2cos( work2, www )
        else if ( mod(m,2) == 1 ) then !! Odd m
!           call grid2cos( work1 )
!           call grid2cos( work2 )
           call grid2cos_polar0( work1, www )
           call grid2cos_polar0( work2, www )
        else                           !! Even m
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  else
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        if ( m == 0 ) then             !! m = 0
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat_inv(j) !! U/cos2(lat)
           end do
           call grid2cos( work1, www )
           work2(:) = 0.0d0
        else if ( mod(m,2) == 1 ) then !! Odd m
           do j=1,JMAX
              work1(j) = data(1,m,j)*coslat(j)  !! U (=u*cos(lat))
              work2(j) = data(2,m,j)*coslat(j)  !! U (=u*cos(lat))
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        else                           !! Even m
           do j=1,JMAX
              work1(j) = data(1,m,j)  !! U/cos(lat) (=u)
              work2(j) = data(2,m,j)  !! U/cos(lat) (=u)
           end do
           call grid2sin( work1, www )
           call grid2sin( work2, www )
        end if
        do j=1,NNUM
           qdata(j,m) = cmplx( work1(j), work2(j), kind=8 )
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  end if
  
  call e_time__end(22,"fft_y__g2w_uv")
!
end subroutine fft_y__g2w_uv_interp


!********************************************************************


subroutine fft_y__w2g_interp &
 &( qdata,            &!IN
 &  data      )        !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  complex(8) :: qwork(NNUM)
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(23,"fft_y__w2g")
!
  if ( jcn_dfs >= 1 ) then
     
!     write(6,*) "qwork(:,4)=",qwork(:,4)
     
!$OMP PARALLEL default(SHARED), private(m,j,qwork,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        call sinsin2none( m, qdata(:,m), qwork )
        do j=1,NNUM
           work1(j) = real( qwork(j), kind=8 )
           work2(j) = dimag( qwork(j) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
        else if ( mod(m,2) == 1 ) then !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
        else                           !! Even m
           call cos2grid( work1, www )
           call cos2grid( work2, www )
        end if
        do j=1,JMAX
           data(1,m,j) = work1(j)
           data(2,m,j) = work2(j)
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  else
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,NNUM
           work1(j) = real( qdata(j,m), kind=8 )
           work2(j) = dimag( qdata(j,m) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
!           data(1,m,1:JMAX) = work1(1:JMAX)
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = 0.0d0
           end do
        else if ( mod(m,2) == 1 ) then !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat(j)
              data(2,m,j) = work2(j)*coslat(j)
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
  end if
  
  call e_time__end(23,"fft_y__w2g")
!
end subroutine fft_y__w2g_interp


!********************************************************************


subroutine fft_y__w2g_cheong &
 &( qdata,            &!IN
 &  data      )        !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  complex(8) :: qwork(NNUM)
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(23,"fft_y__w2g")
!
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,NNUM
           work1(j) = real( qdata(j,m), kind=8 )
           work2(j) = dimag( qdata(j,m) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
!           data(1,m,1:JMAX) = work1(1:JMAX)
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = 0.0d0
           end do
        else if ( mod(m,2) == 1 ) then !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat(j)
              data(2,m,j) = work2(j)*coslat(j)
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
  
  call e_time__end(23,"fft_y__w2g")
!
end subroutine fft_y__w2g_cheong


!********************************************************************


subroutine fft_y__w2g_cheonglike &
 &( qdata,            &!IN
 &  data      )        !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  complex(8) :: qwork(NNUM,0:MMAX)
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(23,"fft_y__w2g")
!
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,NNUM
           work1(j) = real( qdata(j,m), kind=8 )
           work2(j) = dimag( qdata(j,m) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
!           data(1,m,1:JMAX) = work1(1:JMAX)
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = 0.0d0
           end do
        else if ( m == 1 ) then        !! m = 1
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        else if ( mod(m,2) == 0 ) then !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat(j)
              data(2,m,j) = work2(j)*coslat(j)
           end do
        else                           !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
!              data(1,m,j) = work1(j)*coslat(j)*coslat(j)
!              data(2,m,j) = work2(j)*coslat(j)*coslat(j)
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
  
  call e_time__end(23,"fft_y__w2g")
!
end subroutine fft_y__w2g_cheonglike


!********************************************************************


subroutine fft_y__w2g_orszag &
 &( qdata,            &!IN
 &  data      )        !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  complex(8) :: qwork(NNUM,0:MMAX)
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
  
  call e_time__start(23,"fft_y__w2g")
!
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,NNUM
           work1(j) = real( qdata(j,m), kind=8 )
           work2(j) = dimag( qdata(j,m) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
!           data(1,m,1:JMAX) = work1(1:JMAX)
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = 0.0d0
           end do
        else if ( m == 1 ) then        !! m = 1
           call cos2grid( work1, www )
           call cos2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat(j)
              data(2,m,j) = work2(j)*coslat(j)
           end do
        else if ( mod(m,2) == 0 ) then !! Even m
           call cos2grid( work1, www )
           call cos2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        else                           !! Odd m
           call cos2grid( work1, www )
           call cos2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat(j)
              data(2,m,j) = work2(j)*coslat(j)
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
  
  call e_time__end(23,"fft_y__w2g")
!
end subroutine fft_y__w2g_orszag


!********************************************************************


subroutine fft_y__w2g_dy_interp &
 &( qdata,            &!IN
 &  data      )        !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  complex(8) :: qwork(NNUM)
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j,n,nn

  logical,save :: first_togrid_dy_1 = .true.
  real(8),save,allocatable :: a(:,:)
  real(8),save,allocatable :: b(:,:)
  complex(8) :: aay
    
  call e_time__start(24,"fft_y__w2g_dy")

    if ( first_togrid_dy_1 ) then

      if ( jcn_dfs == 0 ) then
        allocate( a(0:NNUM-1,0:mmax) )
        allocate( b(0:NNUM-1,0:mmax) )

!$OMP PARALLEL default(SHARED), private(m,n,nn)
 !$OMP DO schedule(STATIC,1)
        do m=0,mmax
          n=0
          if ( m == 0 ) then
             nn = n  !! Wave number
             b(n,m) = (nn+1)*0.5d0
          else if ( mod(m,2) == 0 ) then
             nn = n+1   !! Wave number
             b(n,m) = nn*0.5d0
          else
             nn = n+1   !! Wave number
             b(n,m) = (nn+1)*0.5d0
          end if
!          if ( m == 0 ) then
!          else
!             b(n,m) = b(n,m)*D2
!          end if

          do n=1,NNUM-2
             if ( m == 0 ) then
               nn = n  !! Wave number
               a(n,m) = -(nn-1)*0.5d0
               b(n,m) = (nn+1)*0.5d0
             else if ( mod(m,2) == 0 ) then
               nn = n+1   !! Wave number
               a(n,m) = -nn*0.5d0
               b(n,m) = nn*0.5d0
             else
               nn = n+1   !! Wave number
               a(n,m) = -(nn-1)*0.5d0
               b(n,m) = (nn+1)*0.5d0
             end if
!             a(n,m) = a(n,m)*D2
!             b(n,m) = b(n,m)*D2
          end do

          n=NNUM-1
          if ( m == 0 ) then
            nn = n  !! Wave number
            a(n,m) = -(nn-1)*0.5d0
          else if ( mod(m,2) == 0 ) then
            nn = n+1   !! Wave number
            a(n,m) = -nn*0.5d0
          else
            nn = n+1   !! Wave number
            a(n,m) = -(nn-1)*0.5d0
          end if
!          if ( m == 0 ) then
!             a(n,m) = a(n,m)*D2
!          else
!          end if
        end do
 !$OMP END DO
!$OMP END PARALLEL
      end if

      first_togrid_dy_1 = .false.
    end if
  
  if ( jcn_dfs >= 1 ) then
     
!$OMP PARALLEL default(SHARED), private(m,j,qwork,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        call yderiv( m, qdata(:,m), qwork )
        do j=1,NNUM
           work1(j) = real( qwork(j), kind=8 )
           work2(j) = dimag( qwork(j) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call sin2grid( work1, www )
           work2(:) = 0.0d0
        else if ( mod(m,2) == 1 ) then !! Odd m
           call cos2grid( work1, www )
           call cos2grid( work2, www )
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
        end if
        do j=1,JMAX
           data(1,m,j) = work1(j)
           data(2,m,j) = work2(j)
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL

  else
  
!$OMP PARALLEL default(SHARED), private(m,n,j,qwork,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
       qwork(:) = qdata(:,m)
          j=1
          n=j-1
          aay = b(n,m)*qwork(j+1)
          work1(j) = real (aay,kind=8)
          work2(j) = dimag(aay)
          do j=2,NNUM-1
             n=j-1
             aay = a(n,m)*qwork(j-1) + b(n,m)*qwork(j+1)
             work1(j) = real (aay,kind=8)
             work2(j) = dimag(aay)
          end do
          j=NNUM
          n=j-1
          aay = a(n,m)*qwork(j-1)
          work1(j) = real (aay,kind=8)
          work2(j) = dimag(aay)
          do j=NNUM+1,JMAX
             work1(j) = 0.0d0
             work2(j) = 0.0d0
          end do

        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
!           data(1,m,:) = work1(:)
           do j=1,JMAX
!              data(1,m,j) = work1(j)
              data(1,m,j) = work1(j)*coslat_inv(j)
              data(2,m,j) = 0.0d0
           end do
        else if ( mod(m,2) == 1 ) then !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
!              data(1,m,j) = work1(j)
!              data(2,m,j) = work2(j)
              data(1,m,j) = work1(j)*coslat_inv(j)
              data(2,m,j) = work2(j)*coslat_inv(j)
           end do
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
!              data(1,m,j) = work1(j)*coslat(j)
!              data(2,m,j) = work2(j)*coslat(j)
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
  end if
!
  call e_time__end(24,"fft_y__w2g_dy")
!
end subroutine fft_y__w2g_dy_interp


!********************************************************************

subroutine fft_y__w2g_dy_cheong &
 &( qdata,            &!IN
 &  data      )        !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  complex(8) :: qwork(NNUM)
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j,n,nn

  logical,save :: first_togrid_dy_1 = .true.
  real(8),save,allocatable :: a(:,:)
  real(8),save,allocatable :: b(:,:)
  complex(8) :: aay
    
  call e_time__start(24,"fft_y__w2g_dy")

    if ( first_togrid_dy_1 ) then

!xx      if ( jcn_dfs == 0 ) then
        allocate( a(0:NNUM-1,0:mmax) )
        allocate( b(0:NNUM-1,0:mmax) )

!$OMP PARALLEL default(SHARED), private(m,n,nn)
 !$OMP DO schedule(STATIC,1)
        do m=0,mmax
          n=0
          if ( m == 0 ) then
             nn = n  !! Wave number
             b(n,m) = (nn+1)*0.5d0
          else if ( mod(m,2) == 0 ) then
             nn = n+1   !! Wave number
             b(n,m) = nn*0.5d0
          else
             nn = n+1   !! Wave number
             b(n,m) = (nn+1)*0.5d0
          end if
!          if ( m == 0 ) then
!          else
!             b(n,m) = b(n,m)*D2
!          end if

          do n=1,NNUM-2
             if ( m == 0 ) then
               nn = n  !! Wave number
               a(n,m) = -(nn-1)*0.5d0
               b(n,m) = (nn+1)*0.5d0
             else if ( mod(m,2) == 0 ) then
               nn = n+1   !! Wave number
               a(n,m) = -nn*0.5d0
               b(n,m) = nn*0.5d0
             else
               nn = n+1   !! Wave number
               a(n,m) = -(nn-1)*0.5d0
               b(n,m) = (nn+1)*0.5d0
             end if
!             a(n,m) = a(n,m)*D2
!             b(n,m) = b(n,m)*D2
          end do

          n=NNUM-1
          if ( m == 0 ) then
            nn = n  !! Wave number
            a(n,m) = -(nn-1)*0.5d0
          else if ( mod(m,2) == 0 ) then
            nn = n+1   !! Wave number
            a(n,m) = -nn*0.5d0
          else
            nn = n+1   !! Wave number
            a(n,m) = -(nn-1)*0.5d0
          end if
!          if ( m == 0 ) then
!             a(n,m) = a(n,m)*D2
!          else
!          end if
        end do
 !$OMP END DO
!$OMP END PARALLEL
!xx      end if

      first_togrid_dy_1 = .false.
    end if
  
!$OMP PARALLEL default(SHARED), private(m,n,j,qwork,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
       qwork(:) = qdata(:,m)
          j=1
          n=j-1
          aay = b(n,m)*qwork(j+1)
          work1(j) = real (aay,kind=8)
          work2(j) = dimag(aay)
          do j=2,NNUM-1
             n=j-1
             aay = a(n,m)*qwork(j-1) + b(n,m)*qwork(j+1)
             work1(j) = real (aay,kind=8)
             work2(j) = dimag(aay)
          end do
          j=NNUM
          n=j-1
          aay = a(n,m)*qwork(j-1)
          work1(j) = real (aay,kind=8)
          work2(j) = dimag(aay)
          do j=NNUM+1,JMAX
             work1(j) = 0.0d0
             work2(j) = 0.0d0
          end do

        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
!           data(1,m,:) = work1(:)
           do j=1,JMAX
!              data(1,m,j) = work1(j)
              data(1,m,j) = work1(j)*coslat_inv(j)
              data(2,m,j) = 0.0d0
           end do
        else if ( mod(m,2) == 1 ) then !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
!              data(1,m,j) = work1(j)
!              data(2,m,j) = work2(j)
              data(1,m,j) = work1(j)*coslat_inv(j)
              data(2,m,j) = work2(j)*coslat_inv(j)
           end do
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
!              data(1,m,j) = work1(j)*coslat(j)
!              data(2,m,j) = work2(j)*coslat(j)
              data(1,m,j) = work1(j)
              data(2,m,j) = work2(j)
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
!
  call e_time__end(24,"fft_y__w2g_dy")
!
end subroutine fft_y__w2g_dy_cheong


!********************************************************************


subroutine fft_y__w2g_uv_interp &
 &( qdata,               &!IN
 &  data      )           !INOUT
!
  complex(8),intent(in) :: qdata(NNUM,0:MMAX)
  real(8),intent(out) :: data(2,0:IMAX/2-1,JMAX)
!
  real(8) :: work1(JMAX+2)
  real(8) :: work2(JMAX+2)
  real(8) :: www(JMAX+2)
  integer :: m,j
!
  call e_time__start(25,"fft_y__w2g_uv")
  
  if ( jcn_dfs >= 1 ) then
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,NNUM
           work1(j) = real( qdata(j,m), kind=8 )
           work2(j) = dimag( qdata(j,m) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call sin2grid( work1, www )
           work2(:) = 0.0d0
        else if ( mod(m,2) == 1 ) then !! Odd m
           call cos2grid( work1, www )
           call cos2grid( work2, www )
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
        end if
        do j=1,JMAX
           data(1,m,j) = work1(j)
           data(2,m,j) = work2(j)
        end do
     end do
 !$OMP END DO
!$OMP END PARALLEL
  else
!$OMP PARALLEL default(SHARED), private(m,j,work1,work2,www)
 !$OMP DO schedule(STATIC,1)
     do m=0,MMAX
        do j=1,NNUM
           work1(j) = real( qdata(j,m), kind=8 )
           work2(j) = dimag( qdata(j,m) )
        end do
        do j=NNUM+1,JMAX
           work1(j) = 0.0d0
           work2(j) = 0.0d0
        end do
        if ( m == 0 ) then             !! m = 0
           call cos2grid( work1, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat(j) !! u (=U/cos(lat))
              data(2,m,j) = 0.0d0
           end do
        else if ( mod(m,2) == 1 ) then !! Odd m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)*coslat_inv(j) !! u (=U/cos(lat))
              data(2,m,j) = work2(j)*coslat_inv(j) !! u (=U/cos(lat))
           end do
        else                           !! Even m
           call sin2grid( work1, www )
           call sin2grid( work2, www )
           do j=1,JMAX
              data(1,m,j) = work1(j)       !! u (=U/cos(lat))
              data(2,m,j) = work2(j)       !! u (=U/cos(lat))
           end do
        end if
     end do
 !$OMP END DO
!$OMP END PARALLEL
  end if
  
  call e_time__end(25,"fft_y__w2g_uv")
!
end subroutine fft_y__w2g_uv_interp


!********************************************************************


subroutine grid2cos &
 &( data, work2      )            !INOUT
!
  real(8),intent(inout) :: data(JMAX+2)
  real(8),intent(out) :: work2(JMAX+2)
  integer :: j
! 
  if ( jcn_grid == -1 ) then
     work2(1:JMAX) = data(1:JMAX)*coslat(1:JMAX)/(JMAX+1)
     call dsint( JMAX, work2, WSAVE_SIN ) !! Sin transform
     data(JMAX+2) = 0.0d0
     data(JMAX+1) = 0.0d0
     data(JMAX)   = work2(JMAX)*2.0d0
     data(JMAX-1) = work2(JMAX-1)*2.0d0
     do j=JMAX-2,2,-1
        data(j) = work2(j)*2.0d0 + data(j+2)
     end do
     j=1
     data(j) = work2(j) + data(j+2)*0.5d0
  else if ( jcn_grid == 1 ) then
     call dcost( JMAX, data(:), WSAVE_COS )
     data(1)=data(1)/(2.0d0*(JMAX-1))
     data(2:JMAX-1)=data(2:JMAX-1)/(JMAX-1)
     data(JMAX)=data(JMAX)/(2.0d0*(JMAX-1))
  else
     call dcosqb( JMAX, data(:), WSAVE_COS )
     data(1)=data(1)/(4.0d0*JMAX)
     data(2:JMAX)=data(2:JMAX)/(2.0d0*JMAX)
  end if
!
end subroutine grid2cos


!********************************************************************


subroutine grid2cos_polar0 &
 &( data, work2     )            !INOUT
!
  real(8),intent(inout) :: data(JMAX+2)
  real(8),intent(out) :: work2(JMAX+2)
  integer :: j
! 
  if ( jcn_grid == -1 ) then
     work2(2:JMAX+1) = data(1:JMAX)
     work2(1)      = 0.0d0              !! Zero at the poles
     work2(JMAX+2) = 0.0d0              !! Zero at the poles
     call dcost( JMAX+2, work2, WSAVE_COS )
     data(1)=work2(1)/(2.0d0*(JMAX+1))
     data(2:JMAX+1)=work2(2:JMAX+1)/(JMAX+1)
     data(JMAX+2)=work2(JMAX+2)/(2.0d0*(JMAX+1))
  else if ( jcn_grid == 1 ) then
     call dcost( JMAX, data, WSAVE_COS )
     data(1)=data(1)/(2.0d0*(JMAX-1))
     data(2:JMAX-1)=data(2:JMAX-1)/(JMAX-1)
     data(JMAX)=data(JMAX)/(2.0d0*(JMAX-1))
  else
     call dcosqb( JMAX, data, WSAVE_COS )
     data(1)=data(1)/(4.0d0*JMAX)
     data(2:JMAX)=data(2:JMAX)/(2.0d0*JMAX)
  end if
!
end subroutine grid2cos_polar0


!********************************************************************


subroutine grid2sin &
 &( data, www     )            !INOUT
!
  real(8),intent(inout) :: data(JMAX+2)
  real(8),intent(out) :: www(JMAX+2)

  if ( jcn_grid == -1 ) then
     www(1:JMAX) = data(1:JMAX)
     call dsint( JMAX, www, WSAVE_SIN ) !! Sin transform
     data(1:JMAX)=www(1:JMAX)/(JMAX+1)
     data(JMAX+1) = 0.0d0
  else if ( jcn_grid == 1 ) then
     www(1:JMAX-2) = data(2:JMAX-1)
     call dsint( JMAX-2, www, WSAVE_SIN ) !! Sin transform
     data(1:JMAX-2)=www(1:JMAX-2)/(JMAX-1)
     data(JMAX-1) = 0.0d0
     data(JMAX) = 0.0d0
  else
     call dsinqb( JMAX, data(:), WSAVE_SIN ) !! Sin transform
     data(1:JMAX-1)=data(1:JMAX-1)/(2.0d0*JMAX)
     data(JMAX)=data(JMAX)/(4.0d0*jmax)
  end if

end subroutine grid2sin


!********************************************************************


subroutine cos2grid &
 &( data, www     )            !INOUT
!
  real(8),intent(inout) :: data(JMAX+2)
  real(8),intent(out) :: www(JMAX+2)
  integer :: j
! 
  if ( jcn_grid == -1 ) then
     www(1) = data(1)
     www(2:JMAX)  = data(2:JMAX)*0.5d0
     www(JMAX+1) = 0.0d0
     www(JMAX+2) = 0.0d0
     call dcost( JMAX+2, www(:), WSAVE_COS ) !! Cos transform
     data(1:JMAX) = www(2:JMAX+1)
  else if ( jcn_grid == 1 ) then
     data(2:JMAX-1)=data(2:JMAX-1)*0.5d0
     call dcost( JMAX, data(:), WSAVE_COS ) !! Cos transform
  else
     data(2:JMAX) = data(2:JMAX)*0.5d0
     call dcosqf( JMAX, data(:), WSAVE_COS ) !! Cos transform
  end if
!
end subroutine cos2grid


!********************************************************************


subroutine sin2grid &
 &( data, www     )            !INOUT
!
  real(8),intent(inout) :: data(JMAX)
  real(8),intent(out) :: www(JMAX+2)
!
  if ( jcn_grid == -1 ) then  
    www(1:JMAX) = data(1:JMAX)*0.5d0
    call dsint( JMAX, www, WSAVE_SIN ) !! Sin transform
    data(1:JMAX) = www(1:JMAX)
  else if ( jcn_grid == 1 ) then  
    www(1:JMAX-2) = data(1:JMAX-2)*0.5d0
    call dsint( JMAX-2, www, WSAVE_SIN ) !! Sin transform
    data(2:JMAX-1) = www(1:JMAX-2)
    data(1)    = 0.0d0
    data(JMAX) = 0.0d0
  else
    data(1:JMAX-1) = data(1:JMAX-1)*0.5d0
    call dsinqf( JMAX, data, WSAVE_SIN ) !! Cos transform
  end if
!
end subroutine sin2grid


!********************************************************************


  subroutine sinsin2none( m, wdata, wdata2 )
    !
    integer,intent(in) :: m
    complex(kind=8),dimension(NNUM),intent(IN) :: wdata
    complex(kind=8),dimension(NNUM),intent(OUT) :: wdata2
    !
    integer :: j

      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        j = 1
        wdata2(j) = 0.25d0*( 3.0d0*wdata(j) - wdata(j+2) )
        j = 2
        wdata2(j) = 0.25d0*( 2.0d0*wdata(j) - wdata(j+2) )
        do j=3,NNUM-2
          wdata2(j) = 0.25d0*( -wdata(j-2) + 2.0d0*wdata(j) - wdata(j+2) )
        end do
        do j=NNUM-1,NNUM
          wdata2(j) = 0.25d0*( -wdata(j-2) + 2.0d0*wdata(j) )
        end do
      else if ( m == 0 ) then  
        wdata2(1:NNUM) = wdata(1:NNUM)
      else if ( mod(m,2) == 1 ) then
        j = 1
        wdata2(j) = wdata(j) - 0.5d0*wdata(j+2)
        do j=2,NNUM-2
          wdata2(j) = 0.5d0*( wdata(j) - wdata(j+2) )
        end do
        do j=NNUM-1,NNUM
          wdata2(j) = 0.5d0*wdata(j)
        end do
      else
        do j=1,2
          wdata2(j) = 0.5d0*wdata(j)
        end do
        do j=3,NNUM
          wdata2(j) = 0.5d0*( -wdata(j-2) + wdata(j) )
        end do
      end if

    return
  end subroutine sinsin2none


!********************************************************************


  subroutine none2sinsin_pre

!    type(type_multidiagonal) :: t
    integer :: l2,m,nn
    
    allocate( md1(0:MMAX) )
    allocate( md2(0:MMAX) )
    
!$OMP PARALLEL default(SHARED), private(l2,m,nn)

   !$OMP DO schedule(STATIC)
    do m=1,3
!      t => md(m)
      if ( m == 1 .or. m == 2 ) then
!        allocate( md1(m)%d(NNUMHF1,5) )
!        allocate( md1(m)%iww(NNUMHF1) )
!        allocate( md2(m)%d(NNUMHF2,5) )
!        allocate( md2(m)%iww(NNUMHF2) )
        call tridiagonal__allocate( NNUM1_M_DFS(m), NNUMHF1, md1(m)%d, md1(m)%iww ) !IN,IN,OUT,OUT
        call tridiagonal__allocate( NNUM2_M_DFS(m), NNUMHF2, md2(m)%d, md2(m)%iww ) !IN,IN,OUT,OUT
      else if ( m == 3 .and. jcn_dfs >= 2 ) then
!        allocate( md1(m)%d(NNUMHF1+1,10) )
!        allocate( md1(m)%iww(NNUMHF1+1) )
!        allocate( md2(m)%d(NNUMHF2+1,10) )
!        allocate( md2(m)%iww(NNUMHF2+1) )
        call pentadiagonal__allocate( NNUM1_M_DFS(m), NNUMHF1+1, md1(m)%d, md1(m)%iww ) !IN,IN,OUT,OUT
        call pentadiagonal__allocate( NNUM2_M_DFS(m), NNUMHF2+1, md2(m)%d, md2(m)%iww ) !IN,IN,OUT,OUT
      end if
      
      if ( m == 1 ) then
        do l2=1,NNUM1_M_DFS(m)
          md1(m)%d(l2,1) = -1.0d0
          md1(m)%d(l2,2) = 2.0d0
          md1(m)%d(l2,3) = -1.0d0
        end do
        md1(m)%d(2,1) = -2.0d0  !! n=2
        
        do l2=1,NNUM2_M_DFS(m)
          md2(m)%d(l2,1) = -1.0d0
          md2(m)%d(l2,2) = 2.0d0
          md2(m)%d(l2,3) = -1.0d0
        end do
        md2(m)%d(1,2) = 1.0d0  !! n=1
      else if ( m == 2 ) then
        do l2=1,NNUM1_M_DFS(m)
          md1(m)%d(l2,1) = -1.0d0
          md1(m)%d(l2,2) = 2.0d0
          md1(m)%d(l2,3) = -1.0d0
        end do
        md1(m)%d(1,1) = -999.0d33
        md1(m)%d(1,2) = 3.0d0  !! n=1
        
        do l2=1,NNUM2_M_DFS(m)
          md2(m)%d(l2,1) = -1.0d0
          md2(m)%d(l2,2) = 2.0d0
          md2(m)%d(l2,3) = -1.0d0
        end do
      else if ( m == 3 .and. jcn_dfs >= 2 ) then
        do l2=1,NNUM1_M_DFS(m)
          md1(m)%d(l2,1) = 1.0d0
          md1(m)%d(l2,2) = -4.0d0
          md1(m)%d(l2,3) = 6.0d0
          md1(m)%d(l2,4) = -4.0d0
          md1(m)%d(l2,5) = 1.0d0
        end do
        
        md1(m)%d(1,1) = -999.0d33
        md1(m)%d(1,2) = -999.0d33
        md1(m)%d(1,3) = 10.0d0   !! n=1
        md1(m)%d(1,4) = -5.0d0   !! n=1

        md1(m)%d(2,1) = -999.0d33
        md1(m)%d(2,2) = -5.0d0   !! n=3
        
        do l2=1,NNUM2_M_DFS(m)
          md2(m)%d(l2,1) = 1.0d0
          md2(m)%d(l2,2) = -4.0d0
          md2(m)%d(l2,3) = 6.0d0
          md2(m)%d(l2,4) = -4.0d0
          md2(m)%d(l2,5) = 1.0d0
        end do
        
        md2(m)%d(1,1) = -999.0d33
        md2(m)%d(1,2) = -999.0d33
        md2(m)%d(1,3) = 5.0d0    !! n=2

        md2(m)%d(2,1) = -999.0d33
      end if
      
      if ( m == 1 .or. m == 2 ) then
        call tridiagonal__ini( NNUM1_M_DFS(m), NNUMHF1, md1(m)%d, md1(m)%iww ) !IN,IN,OUT,OUT
        call tridiagonal__ini( NNUM2_M_DFS(m), NNUMHF2, md2(m)%d, md2(m)%iww ) !IN,IN,OUT,OUT
      else if ( m == 3 .and. jcn_dfs >= 2 ) then
        call pentadiagonal__ini( NNUM1_M_DFS(m), NNUMHF1+1, md1(m)%d, md1(m)%iww ) !IN,IN,OUT,OUT
        call pentadiagonal__ini( NNUM2_M_DFS(m), NNUMHF2+1, md2(m)%d, md2(m)%iww ) !IN,IN,OUT,OUT
      end if

    end do
   !$OMP END DO

!$OMP END PARALLEL

    return
  end subroutine none2sinsin_pre


!********************************************************************


  subroutine none2sinsin( wdata )
    !
    complex(kind=8),dimension(NNUM,0:MMAX),intent(INOUT) :: wdata
    
    logical,save :: first_none2sinsin = .true.
    real(8),save,allocatable :: a1(:,:)
    real(8),save,allocatable :: b1(:,:)
    real(8),save,allocatable :: c1(:,:)
    real(8),save,allocatable :: a2(:,:)
    real(8),save,allocatable :: b2(:,:)
    real(8),save,allocatable :: c2(:,:)
    !
    complex(8) :: g1(NNUMHF1+1)
    complex(8) :: g2(NNUMHF2+1)
    complex(8) :: x1(NNUMHF1+1)
    complex(8) :: x2(NNUMHF2+1)
    complex(8) :: work(max(NNUMHF1,NNUMHF2)+1)
    integer :: m,mm,l,ll,l2

    call e_time__start(28,"none2sinsin")

!$OMP PARALLEL default(SHARED), private(m,mm,l,ll,l2,g1,g2,x1,x2,work)
 !$OMP DO schedule(STATIC)
    do m=1,mmax  !! m=0 is not necessary
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        mm = 3

        l2=1
        g1(1) = 12.0d0*wdata(1,m) - 4.0d0*wdata(1+2,m)
        g2(1) = 8.0d0*wdata(2,m) - 4.0d0*wdata(2+2,m)
        do l2=2,NNUM1_M_DFS(m)
          l=l2*2-1
          g1(l2) = -4.0d0*wdata(l-2,m) + 8.0d0*wdata(l,m) - 4.0d0*wdata(l+2,m)
        end do
        do l2=2,NNUM2_M_DFS(m)
          ll=l2*2
          g2(l2) = -4.0d0*wdata(ll-2,m) + 8.0d0*wdata(ll,m) - 4.0d0*wdata(ll+2,m)
        end do
      
      else if ( mod(m,2) == 1 ) then
        mm = 1

        l2=1
        g1(1) = 2.0d0*wdata(1,m)
        g2(1) = 2.0d0*wdata(2,m)
        do l2=2,NNUM1_M_DFS(m)
          l=l2*2-1
          g1(l2) = -2.0d0*wdata(l-2,m) + 2.0d0*wdata(l,m)
        end do
        do l2=2,NNUM2_M_DFS(m)
          ll=l2*2
          g2(l2) = -2.0d0*wdata(ll-2,m) + 2.0d0*wdata(ll,m)
        end do
        
      else
        mm = 2

        l2=1
        g1(1) = 4.0d0*wdata(1,m) - 2.0d0*wdata(1+2,m)
        g2(1) = 2.0d0*wdata(2,m) - 2.0d0*wdata(2+2,m)
        do l2=2,NNUM1_M_DFS(m)
          l=l2*2-1
          g1(l2) = 2.0d0*wdata(l,m) - 2.0d0*wdata(l+2,m)
        end do
        do l2=2,NNUM2_M_DFS(m)
          ll=l2*2
          g2(l2) = 2.0d0*wdata(ll,m) - 2.0d0*wdata(ll+2,m)
        end do
       
      end if
        
      !! Calculate x1,x2 from g1,g2        
      if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
        call pentadiagonal__solve( NNUM1_M_DFS(mm), NNUMHF1+1, md1(mm)%d, md1(mm)%iww, &
         &                         work, g1, x1 )
        call pentadiagonal__solve( NNUM2_M_DFS(mm), NNUMHF2+1, md2(mm)%d, md2(mm)%iww, &
         &                         work, g2, x2 )
      else
        call tridiagonal__solve( NNUM1_M_DFS(mm), NNUMHF1, md1(mm)%d, md1(mm)%iww, &
         &                       work, g1, x1 )
        call tridiagonal__solve( NNUM2_M_DFS(mm), NNUMHF2, md2(mm)%d, md2(mm)%iww, &
         &                       work, g2, x2 )
      end if

      do l2=1,NNUMHF1
        l = l2*2-1
        wdata(l,m) = x1(l2)
      end do
      do l2=1,NNUMHF2
        l = l2*2
        wdata(l,m) = x2(l2)
      end do
    end do
 !$OMP END DO
!$OMP END PARALLEL

    call e_time__end(28,"none2sinsin")
    !
    return
  end subroutine none2sinsin


!********************************************************************


  subroutine yderiv(m,qdata,qdata_dy) !IN,IN,OUT
    !
    integer,intent(IN) :: m
    complex(8),intent(IN) :: qdata(NNUM) 
    complex(8),intent(OUT) :: qdata_dy(NNUM)
    integer :: j,n

    if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
      j=1
      n=j-1
      qdata_dy(j) = -(-0.25d0*n*qdata(j+1))
      j=2
      n=j-1
      qdata_dy(j) = -0.25d0*n*( 3.0d0*qdata(j-1) - qdata(j+1) )
      j=3
      n=j-1
      qdata_dy(j) = -0.25d0*n*( 2.0d0*qdata(j-1) - qdata(j+1) )
      do j=4,NNUM-1
        n=j-1
        qdata_dy(j) = -0.25d0*n*( -qdata(j-3) + 2.0d0*qdata(j-1) - qdata(j+1) )
      end do
      n=NNUM
      n=j-1
      qdata_dy(j) = -0.25d0*n*( -qdata(j-3) + 2.0d0*qdata(j-1) )
          
    else if ( m == 0 ) then
      do j=1,NNUM-1
        n=j
        qdata_dy(j) = -(-n*qdata(j+1))
      end do
      j=NNUM
      qdata_dy(j) = 0.0d0
    else if ( mod(m,2) == 1 ) then
      j=1
      n=j-1
      qdata_dy(j) = 0.0d0
      j=2
      n=j-1
      qdata_dy(j) = -0.5d0*n*( 2.0d0*qdata(j-1) - qdata(j+1) )
      do j=3,NNUM-1
        n=j-1
        qdata_dy(j) = -0.5d0*n*( qdata(j-1) - qdata(j+1) )
      end do
      j=NNUM
      n=j-1
      qdata_dy(j) = -0.5d0*n*qdata(j-1)
    else
      j=1
      n=j
      qdata_dy(j) = 0.5d0*n*qdata(j+1)
      do j=2,NNUM-1
        n=j
        qdata_dy(j) = -0.5d0*n*( qdata(j-1) - qdata(j+1) )
      end do
      j=NNUM
      n=j 
      qdata_dy(j) = -0.5d0*n*qdata(j-1)
    end if
            
  end subroutine yderiv


!********************************************************************

end module fft_y_interp

