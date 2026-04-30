module vals
        
        use, intrinsic:: iso_fortran_env, only: error_unit
        use ISO_C_BINDING
        use omp_lib
        implicit none

        include 'fftw3.f03'
contains
!         computes the flupuid in fortran
        subroutine flupud(A, F)
                implicit none

                ! dummy arguments
                real(8), intent(in) :: A(:,:)
                real(8), intent(out):: F(:,:)
                integer :: M, N, k, j
                !real(8) :: start, finish, time1, time2

                !real(8), allocatable::B(:,:)
                
                M = size(A,1)
                N = size(A,2)
               ! allocate(B(M,N))
               ! allocate(F(M,N))
                !start = omp_get_wtime()
               ! B = A(M:1:-1,:);
                !finish = omp_get_wtime()

               ! print*, "Time matrix:", finish - start

                ! using openMP
               ! time1 = omp_get_wtime()
                !$OMP PARALLEL DO private(j,k)
                do k = 1, N
                        
                        do j = 1, M
                                F(j,k) = A(M-j+1, k)
                        
                        end do
                        
                end do
                !$OMP END PARALLEL DO
               ! time2 = omp_get_wtime()

                !print*, "Time omp", time2 - time1

               ! print*, "Error", norm2(B - F)
            
        end subroutine flupud

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

        ! computes the ifft
        subroutine ifftshift(data)
                implicit none
                complex(8), intent(inout) :: data(:)
                integer :: k, c, M 
                complex(8) :: tmp

                M = size(data)
                c = M / 2
               
                if (mod(M,2) == 0) then
                        !$OMP PARALLEL  DO private(k, tmp)
                        do k = 1, c
                                tmp = data(k)
                                data(k) = data(k+c)
                                data(k+c) = tmp
                                
                        end do
                        !$OMP END PARALLEL DO
                else
                        tmp = data(M)
                        !$OMspread(cn, 2, m)P PARALLEL DO private(k)
                        do k = c, 1, -1
                                data(c+k+1) = data(k)
                                data(k) = data(c+k)
                        end do
                       !$OMP END PARALLEL DO
                        data(c+1) = tmp
                end if
        end subroutine ifftshift

!       computes the 2D Fourier coefficients on the sphere using DFS
        subroutine vals2coeffs(fjk, coeffs, m, n, gridtype, nthreads)
                implicit none
                integer, intent(in) :: m, n, gridtype, nthreads
                real(8), intent(in) :: fjk(:,:)
                complex(8), allocatable, intent(out) :: coeffs(:,:)

                complex(8), parameter :: i1 = dcmplx(0.0d0, 1.0d0)   ! sqrt(-1)  
                Real(8), parameter:: PI=4.D0*DATAN(1.D0)
                real(8), allocatable :: F(:,:), ft(:,:)
                real(8), allocatable :: S1(:,:), S2(:,:), tpin(:)
                complex(8), allocatable :: tmpin(:), cn(:), cm(:), tmpout(:), ck(:,:)
                integer :: j, k, nt,  iret
                integer(8) :: plan, planb, tn
                 real(8) :: tmp

                if(gridtype .eq. 0) then  
                        allocate(F(2*n, m))
                        allocate(ft(size(fjk,1), size(fjk,2)))
                        allocate(S1(n,m/2))
                        allocate(S2(n,m/2))

                        !call flupud(fjk,ft)
                        !ft = fjk(n:1:-1,:)
                        !ft = fjk
                       ! call flupud(ft(:,m/2+1:m), S1)
                        F(1:n, 1:m/2) = fjk(:,m/2+1:m)
                        !call flupud(ft(1:n,1:m/2), S2)
                        F(1:n, m/2+1:m) = fjk(1:n,1:m/2)
                        F(n+1:2*n, :) = fjk(n:1:-1,:)
                        
                        deallocate(ft)
                        deallocate(S1)
                        deallocate(S2)
                        nt = 2*n
                        !Write(*,*) ((F(k,j),j= 1,m),k= 1,nt)


                        ! computing coefficients for interpolating in [-pi/2, 3pi/2]
                        allocate(cn(nt))
                        tmp = (nt-2d0) / (2*nt)
                        if(mod(nt, 2)  .ne. 0) then
                                !$OMP PARALLEL Do private(j, k)
                                do k=-(nt-1)/2, (nt-1)/2
                                        j = k + (nt-1)/2 + 1
                                        cn(j) = (cdexp(i1*pi*k*tmp)  * (-1.0d0)**(k))
                                       
                                end do
                                !$OMP END PARALLEL DO
                        else
                               !$OMP PARALLEL DO private(j,k)
                                do k=-(nt/2), (nt/2) - 1
                                        j = k + (nt/2)+1
                                        cn(j)=cdexp(i1*pi*k*tmp) * (-1.0d0)**(k)
                                end do
                                !$OMP END PARALLEL DO
                        end if

                else if (gridtype .eq. 1) then
                        allocate(F(2*(n-1), m))
                        allocate(ft(size(fjk,1), size(fjk,2)))
                        allocate(S1(n-1,m/2))
                        allocate(S2(n-1,m/2))
                        ft = fjk(n:1:-1,:)
                        !ft = fjk
                        call flupud(ft(2:n,m/2+1:m), S1)
                        F(1:n-1, 1:m/2) = S1
                        call flupud(ft(2:n,1:m/2), S2)
                        F(1:n-1, m/2+1:m) = S2
                        F(n:2*(n-1), :) = ft(1:n-1,:)

                        nt = 2*(n-1)
                
                        deallocate(ft)
                        deallocate(S1)
                        deallocate(S2)
                        !computing coefficients for interpolating in [-pi/2, 3pi/2]
                       allocate(cn(1:nt))
                        if(mod(nt, 2)  .ne. 0) then
                                !$OMP PARALLEL Do private(j, k)
                                do k=-(nt-1)/2, (nt-1)/2
                                        j = k + (nt-1)/2 + 1
                                        cn(j) = (cdexp(i1*pi*k/2.0d0) * (-1.0d0)**k)
                                end do
                                !$OMP END PARALLEL DO
                        else
                                !$OMP PARALLEL DO private(j,k)
                                do k=-(nt/2), (nt/2)-1
                                        j = k + (nt/2)+1
                                       cn(j) =  cdexp(i1*pi*k/2.0d0) * ((-1.0d0)**(k))
                                       
                                       
                                end do
                                !$OMP END PARALLEL DO
                        end if

                        

               end if

                 ! compute the coefficients in theta direction
                allocate(coeffs(nt, m))
                allocate(ck(nt,m))
                allocate(tmpin(nt))
                allocate(tmpout(nt))
           
                call dfftw_init_threads(iret)
                call dfftw_plan_with_nthreads(nthreads)
                call dfftw_plan_dft_1d(plan, nt, tmpin, tmpout, FFTW_FORWARD, FFTW_ESTIMATE)
                do k = 1, m 
                        tmpin = F(:, k)
                        call dfftw_execute_dft(plan, tmpin, tmpout)
                  !      call fftshift(tmpout)
                        !$OMP PARALLEL DO private(j)
                        do j = 1, nt
                                coeffs(j, k) = tmpout(j)
                        
                        end do
                       !$OMP END PARALLEL DO
                end do
                call  dfftw_destroy_plan(plan)
                ck(1:nt/2,:) = coeffs(nt/2+1:nt,:)
                ck(nt/2+1:nt,:  ) = coeffs(1:nt/2,:)
                ck = ck * spread(cn, 2, m)/nt
                
                deallocate(tmpin)
                deallocate(tmpout)

               
  !              compute coefficients in phi direction
               ! allocate(cm(m))
               ! if(mod(m, 2)  .ne. 0) then
               !         !$OMP PARALLEL Do private(j, k)
               !         do k=-(m-1)/2, (m-1)/2
               !                 j = k + (m-1)/2 + 1
               !                 cm(j) =  (-1.0d0)**k
                !        end do
               !         !$OMP END PARALLEL DO
               ! else
                !        !$OMP PARALLEL DO private(j,k)
               !         do k=-(m/2), (m/2)-1
               !                 j = k + (m/2) + 1
               !                 cm(j) =  (-1.0d0)**k
               !         end do
               !         !$OMP END PARALLEL DO
               ! end if
        

                allocate(tmpin(m))
                allocate(tmpout(m))
                        

                call dfftw_init_threads(iret)
                call dfftw_plan_with_nthreads(nthreads)
                call dfftw_plan_dft_1d(planb, m, tmpin, tmpout, FFTW_FORWARD, FFTW_ESTIMATE)
                do k = 1, nt
                        tmpin = ck(k, :)
                        call dfftw_execute_dft(planb, tmpin, tmpout)
                       ! call fftshift(tmpout)
                       !$OMP PARALLEL DO private(j)
                        do j = 1, m
                                 ck(k, j) = tmpout(j)
                        end do
                        !$OMP END PARALLEL DO
                end do
                call  dfftw_destroy_plan(planb)
                deallocate(tmpin)
                deallocate(tmpout)
                coeffs(:,1:m/2) = ck(:,m/2+1:m)
                coeffs(:,m/2+1:m) = ck(:,1:m/2)
                coeffs = coeffs / m
                
               ! deallocate(cm)  
                deallocate(cn)
        end subroutine vals2coeffs


        subroutine maxerror(F, S, maxerr)
                real(8), intent(in) :: F(:), S(:)
                real(8), intent(out) :: maxerr
                real(8) :: maxdiff, maxf
                integer :: n, i
                maxf = maxval(dabs(F))
                
                maxdiff = 0.0d0
                n = size(S)
                do i=1, n 
                        
                        maxdiff = dmax1(dabs(S(i)-F(i)), maxdiff)
                end do
                
                maxerr = maxdiff / maxf
        end subroutine maxerror

end module vals
