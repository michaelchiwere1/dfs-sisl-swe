        module interp
        use, intrinsic:: iso_fortran_env, only: error_unit
        use ISO_C_BINDING
        use omp_lib
       ! use finufft_mod
        implicit none
        include 'finufft.fh'




        contains

        ! the routine for computing DFS interpolant using NUFFT
        subroutine dfsInterp(coeffs, S, N1, N2, th, lb)
                integer(8), intent(in) :: N1, N2
                complex*16, intent(in) :: coeffs(:)
                real(8), intent(out) :: S(:)
                complex*16, allocatable :: fk(:)
                real(8), intent(in)     :: lb(:), th(:)

                integer :: ier, iflag, ntrans, ttype, dim
                integer(8) :: M, plan
                integer(8), allocatable :: n_modes(:)
               ! real(8), allocatable :: xx(:), yy(:)
                real(8) :: tol
                real(8), pointer :: dummy => null()
                type(nufft_opts), pointer :: defopts => null()

                ! initialize values
                M = size(lb)
                allocate(n_modes(3))
                n_modes(1) = N1
                n_modes(2) = N2
                iflag = 1
                ntrans = 2
                ttype = 2
                dim = 2
                tol = 1.0D-14
                allocate(fk(M))
c evaluating nufft using the guru interface

c               call finufft_makeplan(ttype,dim,n_modes,iflag,
c     +          ntrans,tol,plan,defopts,ier)
c
c                call finufft_setpts(plan,M, th, lb,dummy,dummy, 
c     +           dummy,dummy,dummy,ier)

c               call finufft_execute(plan, fk, coeffs, ier)

c              call finufft_destroy(plan, ier)
              call finufft2d2(M,th,lb,fk,iflag,tol
     +                    ,N1,N2,coeffs,defopts,ier)

                S = dble(fk)
                return
                end subroutine dfsInterp

        end module interp
                                                                                                                                                                                         1,1           Top

