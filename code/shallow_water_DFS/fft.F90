module fft

  use e_time, only : e_time__start, e_time__end
!--------------------------------------------------------------
  implicit none
!--------------------------------------------------------------
  private
  public :: fft__ini, fft__togrid, fft__towave

  integer,save :: IMAX
  integer,save :: JMAX

  real(8),allocatable :: WSAVE(:)
!$OMP threadprivate(WSAVE)

  integer,allocatable :: MC(:)

!********************************************************************
contains
!********************************************************************


subroutine fft__ini(jcn_zonalfilter,m0,imax_in,jmax_in,mmax,coslat) !IN
!
  integer,intent(in) :: jcn_zonalfilter
  integer,intent(in) :: m0
  integer,intent(in) :: imax_in
  integer,intent(in) :: jmax_in
  integer,intent(in) :: mmax
  real(8),intent(in) :: coslat(jmax_in)
!
  integer :: j,mm
!
!xx  integer,parameter :: m0 = 20
!
  write(6,*) 'fft: initialization start.'

  IMAX = imax_in   
  JMAX = jmax_in
   
 !$OMP PARALLEL default(SHARED)
  allocate( wsave(2*IMAX+15) )
  call drffti(IMAX,wsave)
 !$OMP END PARALLEL

  allocate( mc(JMAX) )
  if ( jcn_zonalfilter == 1 ) then
   !$OMP PARALLEL default(SHARED), private(j,mm)
   !$OMP DO schedule(STATIC)
    do j=1,JMAX
      mm = m0 + int( mmax*coslat(j) )
      MC(j) = min(mm,mmax)
    end do
   !$OMP END DO
   !$OMP END PARALLEL
  else
    MC(:) = mmax
  end if
!
end subroutine fft__ini


!********************************************************************


subroutine fft__togrid( a )   !INOUT

  real(8),intent(inout) :: a(IMAX,JMAX)
  real   (kind=8) :: aa(imax)
  integer :: i,j,ie
  
  call e_time__start(31,"fft__togrid")

!$OMP PARALLEL default(SHARED), private(i,j,ie,aa)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
     aa(1)=a(1,j)
!     do i=2,imax-1
!        aa(i)=a(i+1,j)
!     end do
!     aa(imax)=0.0d0
     ie=MC(j)*2+1
     do i=2,ie
        aa(i)=a(i+1,j)
     end do
     aa(ie+1:IMAX)=0.0d0
     call drfftb( IMAX, aa, WSAVE )
     a(:,j)=aa(:)
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call e_time__end(31,"fft__togrid")
!
end subroutine fft__togrid


!********************************************************************


subroutine fft__towave( a )!INOUT
!
  real(8),intent(inout) :: a(IMAX,JMAX)
  real   (kind=8) :: aa(IMAX)
  integer :: i,j,ie
!
  call e_time__start(32,"fft__towave")
!
!$OMP PARALLEL default(SHARED), private(i,j,ie,aa)
 !$OMP DO schedule(STATIC)
  do j=1,JMAX
     aa(:)=a(:,j)/IMAX
     call drfftf(IMAX, aa, WSAVE )
     a(1,j)=aa(1)
     a(2,j)=0.0d0
!     do i=3,imax
!        a(i,j)=aa(i-1)
!     end do
     ie=MC(j)*2+2
     do i=3,ie
        a(i,j)=aa(i-1)
     end do
     a(ie+1:IMAX,j) = 0.0d0
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call e_time__end(32,"fft__towave")
!
end subroutine fft__towave


end module fft
