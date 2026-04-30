module fft

  use e_time, only : e_time__start, e_time__end
!--------------------------------------------------------------
  implicit none
!--------------------------------------------------------------
  private
  public :: fft_run_togrid, fft_run_towave
!  public :: fft_run_grid2sin, fft_run_grid2cos
!  public :: fft_run_sin2grid, fft_run_cos2grid
!--------------------------------------------------------------
!  integer,parameter :: ntrigs = 10000
!  real(8),save :: trigs(ntrigs)
!!
!  integer,save :: ifax(10)
!!
  integer,save :: imax_save = -1
!  integer,save :: jump
!!
!  integer,parameter :: inc = 1
!
  real(8),allocatable :: wsave(:)
!$OMP threadprivate(WSAVE)
!
!********************************************************************
contains
!********************************************************************


subroutine initialize(imax) !IN
!
  integer,intent(in) :: imax
!
  write(6,*) 'fft: initialization start.'

  imax_save = imax
   
!$OMP PARALLEL default(SHARED)
  allocate( wsave(2*imax+15) )
  call drffti(imax,wsave)
!$OMP END PARALLEL
!
end subroutine initialize


!********************************************************************


subroutine fft_run_togrid &
 &( imax, lot,            &!IN
 &  a,                    &!INOUT
 &  no_truncation    )     !IN,optional
!
  integer,intent(in) :: imax
  integer,intent(in) :: lot
  real(8),intent(inout) :: a(imax,lot)
  logical,intent(in),optional :: no_truncation
  real   (kind=8) :: aa(imax)
  integer :: i,j
  
  
  call e_time__start(31,"fft_run_togrid")
!
  if ( imax /= imax_save ) then
    call initialize(imax) !IN
  end if

!$OMP PARALLEL default(SHARED), private(i,j,aa)
 !$OMP DO schedule(STATIC)
  do j=1,lot
     aa(1)=a(1,j)
     do i=2,imax-1
        aa(i)=a(i+1,j)
     end do
     aa(imax)=0.0d0
     call drfftb(imax, aa, wsave )
     a(:,j)=aa(:)
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call e_time__end(31,"fft_run_togrid")
!
end subroutine fft_run_togrid


!********************************************************************


subroutine fft_run_towave &
 &( imax, lot,            &!IN
 &  a,                    &!INOUT
 &  no_truncation    )     !IN,optional
!
  integer,intent(in) :: imax
  integer,intent(in) :: lot
  real(8),intent(inout) :: a(imax,lot)
  logical,intent(in),optional :: no_truncation
  real   (kind=8) :: aa(imax)
  integer :: i,j
!
  call e_time__start(32,"fft_run_towave")
  
  if ( imax /= imax_save ) then
    call initialize(imax) !IN
  end if
!
!$OMP PARALLEL default(SHARED), private(i,j,aa)
 !$OMP DO schedule(STATIC)
  do j=1,lot
     aa(:)=a(:,j)/imax
     call drfftf(imax, aa, wsave )
     a(1,j)=aa(1)
     a(2,j)=0.0d0
     do i=3,imax
        a(i,j)=aa(i-1)
     end do
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call e_time__end(32,"fft_run_towave")
!
end subroutine fft_run_towave


end module fft
