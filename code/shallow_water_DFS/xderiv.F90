module xderiv
!
  use com_dfs, only : jcn_grid, imax, jmax, mmax
  use e_time, only : e_time__start, e_time__end
!
  implicit none
!
!**************************************************************
contains
!**************************************************************

subroutine xderiv_run     &
 &( gd, gdy, coslat_inv,  &!IN
 &  gdx           )        !OUT
!----------------------------------------------------------
  real(8),intent(in)  :: gd(imax,jmax)
  real(8),intent(in)  :: gdy(imax,jmax)
  real(8),intent(in)  :: coslat_inv(jmax)
  real(8),intent(out) :: gdx(imax,jmax)
!
  integer :: imin,j,m,i
!----------------------------------------------------------

  call e_time__start(35,"xderiv_run")

  imin=2*mmax+3
  if ( jcn_grid == 1 ) then
  
!$OMP PARALLEL default(SHARED), private(j,m,i)
   !$OMP DO schedule(STATIC)
    do j=2,jmax-1
      do m=0,mmax
        gdx(2*m+1,j)=-float(m)*gd(2*m+2,j)*coslat_inv(j)
        gdx(2*m+2,j)= float(m)*gd(2*m+1,j)*coslat_inv(j)
      end do
      do i=imin,imax
        gdx(i,j)=0.d0
      end do
    end do
   !$OMP END DO NOWAIT

   !$OMP SINGLE
    j=1
    gdx(:,j)=0.0d0
    m=1
    gdx(2*m+1,j)=  gdy(2*m+2,j)
    gdx(2*m+2,j)= -gdy(2*m+1,j)
   !$OMP END SINGLE NOWAIT
    
   !$OMP SINGLE
    j=jmax
    gdx(:,j)=0.0d0
    m=1
    gdx(2*m+1,j)= -gdy(2*m+2,j)
    gdx(2*m+2,j)=  gdy(2*m+1,j)
   !$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

  else
  
!$OMP PARALLEL default(SHARED), private(j,m,i)
   !$OMP DO schedule(STATIC)
    do j=1,jmax
      do m=0,mmax
        gdx(2*m+1,j)=-float(m)*gd(2*m+2,j)*coslat_inv(j)
        gdx(2*m+2,j)= float(m)*gd(2*m+1,j)*coslat_inv(j)
      end do
      do i=imin,imax
        gdx(i,j)=0.d0
      end do
    end do
   !$OMP END DO NOWAIT
!$OMP END PARALLEL
    
  end if
  
  call e_time__end(35,"xderiv_run")
  !
end subroutine xderiv_run

!**************************************************************

end module xderiv
