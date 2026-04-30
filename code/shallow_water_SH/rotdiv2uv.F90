module rotdiv2uv
!
  use e_time, only : e_time__start, e_time__end
  
!  use prm, only : kmx2, mnwav, mend1, nend1, jend1
!  use prm_ptr_qdata, only :  &
!   & kqdmax,  &!## array size of qu
!   & kqrot,   &!## relative rotation
!   & kqdiv,   &!## divergence
!   & kqu,     &!## U = u*cos(lat)
!   & kqv       !## V = v*cos(lat)
!  use prm_phconst, only : er
!
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  private
  public :: rotdiv2uv_run
!-----------------------------------------------------------------------
!
  logical,save :: initialized = .false.
!
  integer,save :: MTRNC

  real(8),save,allocatable :: eps0 (:) !## a/n*sqrt((n2-m2)/(4*n2-1))
  real(8),save,allocatable :: eps1 (:)       !## a*m/(n*(n+1))
!
!***********************************************************************
  contains
!***********************************************************************


subroutine initialize &
 &( mnwav, mnum, mnstart_uv,    &!IN
 &  er   )             !IN
!
  integer,intent(in) :: mnwav, mnum
  integer,intent(in) :: mnstart_uv(0:mnum)
  real(8),intent(in) :: er
!
  integer :: m,mn_uv
  real(8) :: am,an,an2
!
! ===================================================
!
  initialized = .true.
!
  allocate( eps0 (mnwav+mnum) )
  allocate( eps1 (mnwav+mnum) )
!
  MTRNC=mnum-1

!$OMP PARALLEL default(SHARED), private(m,mn_uv,am,an,an2)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
    do mn_uv=MNSTART_UV(m),MNSTART_UV(m+1)-1
      am=m
      an=m+mn_uv-MNSTART_UV(m)
      an2=an*an
      if(an.eq.0.0) then
        eps0(mn_uv)=0.0
        eps1(mn_uv)=0.0
      else
        eps0(mn_uv)=er*sqrt((an2-am*am)/(4.0*an2*an2-an2))
                          !## eps0 = a/n*sqrt((n2-m2)/(4*n2-1))
        eps1(mn_uv)=er*am/(an*an+an)  !## a*m/(n*(n+1))
      end if
    end do
  end do
 !$OMP END DO
!$OMP END PARALLEL
!
end subroutine initialize


!**********************************************************************


subroutine rotdiv2uv_run  &
 &( mnwav, mnum,         &!IN
 &  mnstart, mnstart_uv,  &!IN
 &  er,                   &!IN
 &  qrot, qdiv,           &!IN
 &  qu, qv             )   !INOUT
!
  integer,intent(in) :: mnwav, mnum
  integer,intent(in) :: mnstart(0:mnum)
  integer,intent(in) :: mnstart_uv(0:mnum)
  real(8),intent(in) :: er
  real(8),intent(in) :: qrot (2,mnwav)
  real(8),intent(in) :: qdiv (2,mnwav)
  real(8),intent(out) :: qu (2,mnwav+mnum) !## U
  real(8),intent(out) :: qv (2,mnwav+mnum) !## V(kmx2,mnwav)
!
!  integer :: l,m,nmax,n,k,ll,lx
!
  integer(8) :: mn, mn_uv
  integer(8) :: m
!
!=======================================================================

  call e_time__start(38,"rotdiv2uv_run")
!
  if ( .not. initialized ) then
    call initialize              &!## Set eps0, eps1
     &( mnwav, mnum, mnstart_uv, &!IN
     &  er   )          !IN 
  end if

!$OMP PARALLEL default(SHARED), private(m,mn,mn_uv)
 !$OMP DO schedule(DYNAMIC)
  do m=0,MTRNC
     mn_uv = MNSTART_UV(m)
     mn    = MNSTART(m)
!     if ( mn+1 <= MNSTART(m+1)-1 ) then
     if ( m /= MTRNC) then
        qu(1,mn_uv) =  eps1(mn_uv)*qdiv(2,mn) + eps0(mn_uv+1)*qrot(1,mn+1)
        qu(2,mn_uv) = -eps1(mn_uv)*qdiv(1,mn) + eps0(mn_uv+1)*qrot(2,mn+1)
        qv(1,mn_uv) =  eps1(mn_uv)*qrot(2,mn) - eps0(mn_uv+1)*qdiv(1,mn+1)
        qv(2,mn_uv) = -eps1(mn_uv)*qrot(1,mn) - eps0(mn_uv+1)*qdiv(2,mn+1)
     else
        qu(1,mn_uv) =  eps1(mn_uv)*qdiv(2,mn)
        qu(2,mn_uv) = -eps1(mn_uv)*qdiv(1,mn)
        qv(1,mn_uv) =  eps1(mn_uv)*qrot(2,mn)
        qv(2,mn_uv) = -eps1(mn_uv)*qrot(1,mn)
     end if
     do mn_uv=MNSTART_UV(m)+1,MNSTART_UV(m+1)-3
        mn = mn_uv - MNSTART_UV(m) + MNSTART(m)
        qu(1,mn_uv) =  eps1(mn_uv)*qdiv(2,mn) + eps0(mn_uv+1)*qrot(1,mn+1) - eps0(mn_uv)*qrot(1,mn-1)
        qu(2,mn_uv) = -eps1(mn_uv)*qdiv(1,mn) + eps0(mn_uv+1)*qrot(2,mn+1) - eps0(mn_uv)*qrot(2,mn-1)
        qv(1,mn_uv) =  eps1(mn_uv)*qrot(2,mn) - eps0(mn_uv+1)*qdiv(1,mn+1) + eps0(mn_uv)*qdiv(1,mn-1)
        qv(2,mn_uv) = -eps1(mn_uv)*qrot(1,mn) - eps0(mn_uv+1)*qdiv(2,mn+1) + eps0(mn_uv)*qdiv(2,mn-1)
     end do
     mn_uv = MNSTART_UV(m+1) - 2
!     if ( mn_uv >= MNSTART_UV(m)+1 ) then
     if ( m /= MTRNC) then
        mn = mn_uv - MNSTART_UV(m) + MNSTART(m)
        qu(1,mn_uv) =  eps1(mn_uv)*qdiv(2,mn) - eps0(mn_uv)*qrot(1,mn-1)
        qu(2,mn_uv) = -eps1(mn_uv)*qdiv(1,mn) - eps0(mn_uv)*qrot(2,mn-1)
        qv(1,mn_uv) =  eps1(mn_uv)*qrot(2,mn) + eps0(mn_uv)*qdiv(1,mn-1)
        qv(2,mn_uv) = -eps1(mn_uv)*qrot(1,mn) + eps0(mn_uv)*qdiv(2,mn-1)
     end if
     mn_uv = MNSTART_UV(m+1) - 1
!     if ( mn_uv >= MNSTART_UV(m)+1 ) then
        mn = mn_uv - MNSTART_UV(m) + MNSTART(m)
        qu(1,mn_uv) = -eps0(mn_uv)*qrot(1,mn-1)
        qu(2,mn_uv) = -eps0(mn_uv)*qrot(2,mn-1)
        qv(1,mn_uv) =  eps0(mn_uv)*qdiv(1,mn-1)
        qv(2,mn_uv) =  eps0(mn_uv)*qdiv(2,mn-1)
!     end if
  end do
 !$OMP END DO
!$OMP END PARALLEL


!  do m=0,MTRNC
!     write(6,*) "m=",m
!     write(6,*) "qu(1,MNSTART_UV(m):MNSTART_UV(m+1)-1)        =",qu(1,MNSTART_UV(m):MNSTART_UV(m+1)-1)
!     write(6,*) "qu0(1,MNSTART(m):MNSTART(m+1)-1),qucor(1,m+1)=",qu0(1,MNSTART(m):MNSTART(m+1)-1),qucor(1,m+1)
!  end do

  call e_time__end(38,"rotdiv2uv_run")

  return
end subroutine rotdiv2uv_run


!**************************************************************************

end module rotdiv2uv
