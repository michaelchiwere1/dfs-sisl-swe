module hdiff4_dfs

  use com_dfs, only : MMAX, NNUM
  use helmholtz_dfs_cmplx, only : helmholtz_dfs_cmplx__run1, helmholtz_dfs_cmplx__run2

  implicit none

  private
  public :: hdiff4_dfs__run

contains

  subroutine hdiff4_dfs__run  &!biharmonic spectral filter
   &( delt, gamma,           &!IN
   &  qvar         )          !INOUT
    !
    real(8),intent(in) :: delt
    real(8),intent(in) :: gamma
    complex(8),intent(inout) :: qvar(NNUM,0:MMAX)
    !
    complex(8) :: gamma1
    complex(8) :: gamma2
    !
    gamma1 = cmplx( 0.0d0, -sqrt(gamma), kind=8 )
    gamma2 = cmplx( 0.0d0, sqrt(gamma), kind=8 )
    !
    call helmholtz_dfs_cmplx__run1     &!
     &( delt, gamma1,          &!IN
     &  qvar        )           !INOUT
    !
    call helmholtz_dfs_cmplx__run2     &!
     &( delt, gamma2,          &!IN
     &  qvar        )           !INOUT
    !
  end subroutine hdiff4_dfs__run

end module hdiff4_dfs
