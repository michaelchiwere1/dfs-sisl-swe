module xderiv
!
  use e_time, only : e_time__start, e_time__end
  
  implicit none
!
!**************************************************************
contains
!**************************************************************

subroutine xderiv_run         &
 &( imax, jblk2, mnum, &!IN
 &  gd  ,                     &!IN
 &  gdx           )            !OUT
!----------------------------------------------------------
  integer,intent(in)  :: imax ,jblk2 ,mnum
  real(8),intent(out) :: gdx(imax,jblk2)
  real(8),intent(in)  :: gd(imax,jblk2)
!
  integer :: imin,j,m
!----------------------------------------------------------

  call e_time__start(35,"xderiv_run")

  imin=2*mnum+1
!$OMP PARALLEL default(SHARED), private(j,m)
 !$OMP DO schedule(STATIC)
  do j=1,jblk2
     do m=1,mnum
        gdx(2*m-1,j)=-float(m-1)*gd(2*m  ,j)
        gdx(2*m  ,j)= float(m-1)*gd(2*m-1,j)
     end do
     do m=imin,imax
        gdx(m,j)=0.d0
     end do
  end do
 !$OMP END DO
!$OMP END PARALLEL

  call e_time__end(35,"xderiv_run")
!
end subroutine xderiv_run

!**************************************************************

end module xderiv
