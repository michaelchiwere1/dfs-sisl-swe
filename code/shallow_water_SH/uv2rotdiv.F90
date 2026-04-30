module uv2rotdiv
!
  use e_time, only : e_time__start, e_time__end
!
!  use prm, only : mnwav, mend1, nend1, jend1, kmx2
!
!  use prm_ptr_qdata, only :  &
!   & kqdmax,  &!## array size of qu
!   & kqrot,   &!## relative rotation
!   & kqdiv,   &!## divergence
!   & kqu,     &!## 1/(a*cos2(lat))*U ( U=u*cos(lat) )
!   & kqv       !## 1/(a*cos2(lat))*V ( V=v*cos(lat) )
!
!  use com_dynamics, only : mwave
!
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  private
  public :: uv2rotdiv_run
!-----------------------------------------------------------------------
  logical,save :: initialized = .false.
  
  integer(8),save :: MTRNC
  integer(8),save :: MNWAV
  integer(8),save :: MNWAV_UV
  
  real(8),save,allocatable :: emn0(:) !## (n-1)*epsmn
  real(8),save,allocatable :: emn1(:) !## (n+1)*epsmn
!
!**********************************************************************
contains
!**********************************************************************


subroutine initialize &
 &( mnwav_in, mnum, mnstart )  !IN
!
  integer,intent(in) :: mnwav_in
  integer,intent(in) :: mnum
  integer,intent(in) :: mnstart(0:mnum)
!
  real(8) :: am,an
  integer :: m,mn
!
! =============================================
!
  initialized = .true.
!
  MNWAV = mnwav_in
  MNWAV_UV = mnwav+mnum
  MTRNC = mnum - 1
  
  allocate( emn0(mnwav) )
  allocate( emn1(mnwav) )
  
!$OMP PARALLEL default(SHARED), private(m,mn,am,an)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
    do mn=MNSTART(m),MNSTART(m+1)-1
      am=m
      an=m+mn-MNSTART(m)
      emn0(mn) = an*sqrt(((an+1)**2-am**2)/(4*(an+1)**2-1))
      emn1(mn) = (an+1)*sqrt((an**2-am**2)/(4*an**2-1))
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
end subroutine initialize


!**********************************************************************


subroutine uv2rotdiv_run  &
 &( mnwav, mnum,         &!IN
 &  mnstart, mnstart_uv,  &!IN
 &  qu, qv,               &!IN
 &  qrot, qdiv   )         !INOUT 1/(a*cos2(lat))*(U,Vv) -> rot, div
!
  integer,intent(in) :: mnwav, mnum
  integer,intent(in) :: mnstart(0:mnum)
  integer,intent(in) :: mnstart_uv(0:mnum)
  real(8),intent(inout) :: qu (2,mnwav+mnum) !## 1/(a*cos2(lat))*U
  real(8),intent(inout) :: qv (2,mnwav+mnum) !## 1/(a*cos2(lat))*V
  real(8),intent(inout) :: qrot (2,mnwav)
  real(8),intent(inout) :: qdiv (2,mnwav)
!
  integer :: m
  integer :: mn, mn_uv
!
! =====================================================================

  call e_time__start(37,"uv2rotdiv_run")
!
  if ( .not. initialized ) then
    call initialize    &
     &( mnwav, mnum, mnstart )  !IN
  end if
!

!$OMP PARALLEL default(SHARED), private(m,mn,mn_uv)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
    mn=MNSTART(m)
    mn_uv = mn - MNSTART(m) + MNSTART_UV(m)
    qrot(1,mn  ) = -m*qv(2,mn_uv) - emn0(mn)*qu(1,mn_uv+1)
    qrot(2,mn  ) =  m*qv(1,mn_uv) - emn0(mn)*qu(2,mn_uv+1)
    qdiv(1,mn  ) = -m*qu(2,mn_uv) + emn0(mn)*qv(1,mn_uv+1)
    qdiv(2,mn  ) =  m*qu(1,mn_uv) + emn0(mn)*qv(2,mn_uv+1)
    do mn=MNSTART(m)+1,MNSTART(m+1)-1
      mn_uv = mn - MNSTART(m) + MNSTART_UV(m)
      qrot(1,mn  ) = -m*qv(2,mn_uv) - emn0(mn)*qu(1,mn_uv+1) + emn1(mn)*qu(1,mn_uv-1)
      qrot(2,mn  ) =  m*qv(1,mn_uv) - emn0(mn)*qu(2,mn_uv+1) + emn1(mn)*qu(2,mn_uv-1)
      qdiv(1,mn  ) = -m*qu(2,mn_uv) + emn0(mn)*qv(1,mn_uv+1) - emn1(mn)*qv(1,mn_uv-1)
      qdiv(2,mn  ) =  m*qu(1,mn_uv) + emn0(mn)*qv(2,mn_uv+1) - emn1(mn)*qv(2,mn_uv-1)
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL

!  write(6,*) "2: qu(1,1:100)=",qu(1,mnwav-100:mnwav)
!  write(6,*) "2: qrot(1,1:100)=",qdiv(1,mnwav-100:mnwav)


!  stop 333
! end if
  
!#endif

!  stop 333

  call e_time__end(37,"uv2rotdiv_run")
  
end subroutine uv2rotdiv_run


end module uv2rotdiv

