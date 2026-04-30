module vals
    use, intrinsic :: iso_fortran_env, only: error_unit
    use iso_c_binding
    use omp_lib
    implicit none

    include 'fftw3.f03'
    contains
!
    ! Computes the fftshift
    subroutine fftshift(data)
        implicit none

        complex(8), intent(inout) :: data(:)
        integer :: c, k, M 
        complex(8) :: tmp

        M = size(data)
        c = M / 2
        
        if(mod(M, 2) == 0) then
                !$OMP PARALLEL DO private(k, tmp)
                do k = 1, c 
                        tmp = data(k)
                        data(k) = data(k+c)
                        data(k+c) = tmp
                       
                end do
                !$OMP END PARALLEL DO
        else
                tmp = data(1)
                !$OMP PARALLEL DO private(k)
                do k = 1, c 
                        data(k) = data(c+k+1)
                        data(c+k+1) = data(k+1)
                end do
                !$OMP END PARALLEL DO
                data(c+1) = tmp
        end if
end subroutine fftshift

!
        subroutine shiftfftw(gridType, n, m, cn)
            integer, intent(in) :: m, n, gridType
            complex*16, intent(inout) :: cn(:)
            complex*16, allocatable :: ck(:)

            real(8), parameter :: PI = 4.0*DATAN(1.D0)
            complex*16, parameter :: i1 = dcmplx(0.0D0,1.0D0)
            real(8) :: tmp, tp
            integer :: k,j
!
            if(gridType .eq. 0) then
                tmp = (2d0*n - 2d0)/(4d0*n)
                tp = PI*tmp + PI
                allocate(ck(2* n))
                !$OMP PARALLEL DO private(j,k)
                do k=-n,n-1
                    j = k + n + 1
                    cn(j) = cdexp(i1*k*tp) /(2.0D0*n*m)
                   ! cn(j) = dcmplx(dcos(k*tp), dsin(k*tp)) / (2D0*n*m)
                end do
                !$OMP END PARALLEL DO         
!
            else if(gridType .eq. 1) then
                allocate(ck(2*(n-1)))
                tp =  3*PI/2
                !$OMP PARALLEL DO private(j,k)
                do k = -n+1, n-2
                    j = k + n
                    cn(j) = cdexp(i1*tp*k)  / (2*(n - 1)*m) 
                   ! cn(j) = dcmplx(dcos(k*tp), dsin(k*tp)) !/ dble(2*(n-1)*m) 
                end do
                !$OMP END PARALLEL DO
            end if
           ! cn = spread(ck, 2, m)
            deallocate(ck)
        end subroutine
!
!   ! computes the 2D Fourier coefficients on the sphere using DFS
        subroutine vals2coeffs(F, cn, plan, tr)
            complex*16, intent(inout) :: F(:,:)
            complex*16, intent(in) :: cn(:)
            integer(8) :: plan
            integer :: tr
!           
            integer ::  k, j
            integer :: m, n
            complex*16, allocatable:: tmp(:,:)
            
!
            n = size(F,1)
            m = size(F,2)

            call dfftw_execute_dft(plan,F,F)

!
            allocate(tmp(n,m))
            tmp(1:n/2,:) = F(n/2+1:n,:) 
            tmp(n/2+1:n,:) = F(1:n/2,:) 
            !$OMP PARALLEL DO private(j,k)
            do k = 1, m
                do j=1,n
                    tmp(j,k) =  tmp(j,k)  * cn(j)
                end do
            end do
            !$OMP END PARALLEL DO
            
            F(:,1:m/2) = tmp(:,m/2+1:m) 
            F(:,m/2+1:m) = tmp(:,1:m/2) 
            F(1:tr, :) = 0.0d0; F(:, 1:tr) = 0.0d0;
            F(n-tr+2:n, :) = 0.0d0; F(:, m-tr+2:m) = 0.0d0;
            deallocate(tmp)
        end subroutine vals2coeffs
!
        subroutine doubleUp(F, dF, gridType)
            real(8), intent(in) :: F(:,:)
            complex*16, intent(inout) :: dF(:,:)
            integer, intent(in) :: gridType
!
            integer :: m, n
            n = size(F, 1)
            m = size(F, 2)
            if(gridType .eq. 0) then
                dF(1:n,1:m/2) = F(:,m/2+1:m)
                dF(1:n,m/2+1:m) = F(1:n,1:m/2)
                dF(n+1:2*n,:) = F(n:1:-1,:)
            else if (gridType .eq. 1) then
                dF(1:n-1,1:m/2) = F(1:n-1,m/2+1:m)
                dF(1:n-1,m/2+1:m) = F(1:n-1,1:m/2)
                dF(n:2*(n-1),:) = F(n:2:-1,:)
               
            end if
        end subroutine doubleUp
!
        subroutine vals2coeffs2(F, cn)
            complex*16, intent(inout) :: F(:,:)
            complex*16, intent(in) :: cn(:)
            integer(8) :: plan
            complex(8), allocatable :: tmpin(:)
            integer :: m, n, k, j, iret
            complex*16, allocatable:: tmp(:,:)
!
            n = size(F,1)
            m = size(F,2)
            allocate(tmp(n,m))
            allocate(tmpin(n))
            call dfftw_init_threads(iret)
            call dfftw_plan_with_nthreads(8)
            call dfftw_plan_dft_1d(plan, n, tmpin, tmpin, FFTW_FORWARD, FFTW_ESTIMATE)
            do k = 1, m 
                tmpin = F(:, k)
                call dfftw_execute_dft(plan, tmpin, tmpin)
                call fftshift(tmpin)
                !$OMP PARALLEL DO private(j)
                do j = 1, n
                    F(j, k) = (tmpin(j)/m) * cn(j)
                        
                end do
                !$OMP END PARALLEL DO
            end do
            call  dfftw_destroy_plan(plan)
            deallocate(tmpin)
        !
            allocate(tmpin(m))
            call dfftw_init_threads(iret)
            call dfftw_plan_with_nthreads(8)
            call dfftw_plan_dft_1d(plan, n, tmpin, tmpin, FFTW_FORWARD, FFTW_ESTIMATE)
            do k = 1, n
                tmpin = F(k, :)
                call dfftw_execute_dft(plan, tmpin, tmpin)
                call fftshift(tmpin)
               !$OMP PARALLEL DO private(j)
                do j = 1, m
                         F(k, j) = tmpin(j) / n
                end do
                !$OMP END PARALLEL DO
        end do
        call  dfftw_destroy_plan(plan)
        !
        !tmp(1:n/2,:) = F(n/2+1:n,:) 
        !tmp(n/2+1:n,:) = F(1:n/2,:) 
        !F(:,1:m/2) = tmp(:,m/2+1:m)  
        !F(:,m/2+1:m) = tmp(:,1:m/2) 
        deallocate(tmpin)

        end subroutine vals2coeffs2


end module vals