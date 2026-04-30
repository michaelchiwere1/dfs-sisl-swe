module gmean_to_zero_dfs

  use com_dfs,only : MMAX,NNUM

  implicit none

  public :: gmean_to_zero_dfs__run

contains

  subroutine gmean_to_zero_dfs__run(qdata) !INOUT
    !
    complex(8),intent(INOUT) :: qdata(0:NNUM-1,0:MMAX)
    complex(8) :: w
    integer :: n
    !
    w = 0.0d0
    do n=2,NNUM-1,2
       w = w - qdata(n,0)/(1.0d0-n*n)
    end do
    qdata(0,0) = w
    !
  end subroutine gmean_to_zero_dfs__run

end module gmean_to_zero_dfs

