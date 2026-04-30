        module interp
        use, intrinsic:: iso_fortran_env, only: error_unit
        use ISO_C_BINDING
        use omp_lib
       ! use finufft_mod
        implicit none
        include 'finufft.fh'
        
        



        contains

        ! the routine for computing DFS interpolant using NUFFT
        subroutine dfsInterp(coeffs, fk, N1, N2, th, lb)
            complex*16, intent(in) :: coeffs(:)
            real*8, intent(in) :: lb(:), th(:)
            integer(8), intent(in) :: N1, N2
            complex*16, intent(out), allocatable :: fk(:)
            integer :: ier, iflag, ttype, dim, ntrans
            integer*8 :: M, plan
            integer*8, allocatable :: n_modes(:)
            real(8) :: tol

c           to pass null pointers for unused arguments
            real*8, pointer :: dummy => null()

c           create the options struct
            type(nufft_opts), pointer :: defopts => null()

c           allocate the output array
            M = size(lb)
            allocate(fk(M))

c           set mandatory parameters to FINUFFT guru interface
            !ttype = 2
            !dim = 2
            iflag = 1
c            allocate(n_modes(3))
c            n_modes(1) = N1
c            n_modes(2) = N2
            ntrans = 2
            tol = 1D-14
c            print*, tol
c           use default opts
c            call finufft_makeplan(ttype, dim, n_modes, iflag, 
c     $         ntrans, tol, plan, defopts, ier)

c          arguments 6-9 ignored for type 2
c           call finufft_setpts(plan, M, th, lb, dummy, dummy,
c     $         dummy, dummy, dummy, ier)

c          evaluate, reads coeffs, writes fk and ier (states)
c           call finufft_execute(plan, fk, coeffs, ier)

c           call finufft_destroy(plan, ier)
          ! deallocate(n_modes)
        call finufft2d2(M,th,lb,fk,
     +     iflag,tol,N1,N2,coeffs,defopts,ier)
           return 
        end subroutine dfsInterp
        


        end module interp
