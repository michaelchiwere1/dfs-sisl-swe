module new_dfs
    use fft
    use fft_y
    use uv2divrot_dfs
    use divrot2uv_dfs
!---------------------------------------------
    implicit none 

!--------------------------------------------------------------------  

 contains
 subroutine newdfs_complex_h(mydata, coeffs, NNMAX, cn)
    implicit none
    integer,            intent(in)    :: NNMAX
    complex(8),         intent(in)    :: cn(:)
    real(8),            intent(inout) :: mydata(:,:)
    complex(8),         intent(inout) :: coeffs(:)

    integer :: MNMAX, mm, nn
    complex(8), allocatable :: Tc(:,:), Ts(:,:), Tcoeff(:,:), myqdata(:,:)
    real(8),    allocatable :: tmp1(:), tmp2(:)
    complex(8), allocatable :: qwork(:)
    complex(8) :: I1
    integer :: need_size

    I1 = (0.0d0, 1.0d0)

    ! ---- 1) FFT in longitude
    call fft__towave(mydata)

    ! ---- 2) Temporary coefficients in latitude
    MNMAX = NNMAX
    allocate(myqdata(NNMAX+1, 0:MNMAX))
    call fft_y__g2w(mydata, myqdata)

    ! ---- 3) Convert to complex basis (Tc, Ts)
    allocate(Tc(2*NNMAX+1, MNMAX+2))
    allocate(Ts(2*NNMAX+1, MNMAX+1))
    Tc = (0.0d0, 0.0d0)
    Ts = (0.0d0, 0.0d0)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mm, nn, tmp1, tmp2, qwork)
    allocate(tmp1(NNMAX+1), tmp2(NNMAX+1), qwork(NNMAX+1))

!$OMP DO SCHEDULE(STATIC)
    do mm = 0, MNMAX
        call sinsin2none(mm, myqdata(:, mm), qwork)

        do nn = 1, NNMAX+1
            tmp1(nn) = real(qwork(nn), kind=8)
            tmp2(nn) = aimag(qwork(nn))
        end do

        if (mm == 0) then
            Tc(NNMAX:1:-1,       mm+1) = 0.5d0 * tmp1(2:NNMAX+1)
            Tc(NNMAX+2:2*NNMAX+1,mm+1) = 0.5d0 * tmp1(2:NNMAX+1)
            Tc(NNMAX+1,          mm+1) = tmp1(1)

        else if (mod(mm, 2) == 1) then
            Tc(NNMAX:1:-1,       mm+1) = -0.5d0*I1 * tmp1(1:NNMAX)
            Tc(NNMAX+2:2*NNMAX+1,mm+1) =  0.5d0*I1 * tmp1(1:NNMAX)

            Ts(NNMAX:1:-1,       mm+1) = -0.5d0*I1 * tmp2(1:NNMAX)
            Ts(NNMAX+2:2*NNMAX+1,mm+1) =  0.5d0*I1 * tmp2(1:NNMAX)

        else
            Tc(NNMAX:1:-1,       mm+1) = 0.5d0 * tmp1(2:NNMAX+1)
            Tc(NNMAX+2:2*NNMAX+1,mm+1) = 0.5d0 * tmp1(2:NNMAX+1)
            Tc(NNMAX+1,          mm+1) = tmp1(1)

            Ts(NNMAX:1:-1,       mm+1) = 0.5d0 * tmp2(2:NNMAX+1)
            Ts(NNMAX+2:2*NNMAX+1,mm+1) = 0.5d0 * tmp2(2:NNMAX+1)
            Ts(NNMAX+1,          mm+1) = tmp2(1)
        end if
    end do
!$OMP END DO

    deallocate(tmp1, tmp2, qwork)
!$OMP END PARALLEL

    ! ---- 4) Combine longitude direction into Tcoeff
    allocate(Tcoeff(2*NNMAX+1, 2*MNMAX+1))

!$OMP PARALLEL DO PRIVATE(nn, mm) SCHEDULE(STATIC)
    do nn = 1, 2*NNMAX+1
        Tcoeff(nn, MNMAX+1) = Tc(nn, 1) * cn(nn)
        do mm = 2, MNMAX+1
            Tcoeff(nn, mm+MNMAX   ) = (Tc(nn, mm) + I1*Ts(nn, mm)) * cn(nn)  ! plus
            Tcoeff(nn, MNMAX+2-mm ) = (Tc(nn, mm) - I1*Ts(nn, mm)) * cn(nn)  ! minus
        end do
    end do
!$OMP END PARALLEL DO

    ! ---- 5) Pack to 1D
    need_size = (2*NNMAX+1)*(2*MNMAX+1)
    if (size(coeffs) /= need_size) then
        stop 'newdfs_complex_h: coeffs has incorrect size.'
    end if
    coeffs = reshape(Tcoeff, [need_size])

    ! ---- 6) Cleanup
    deallocate(Tc, Ts, myqdata, Tcoeff)
end subroutine newdfs_complex_h
!
 !! computing shifts
subroutine fshisfts(cn, NNMAX, gridtype)
    implicit none
    !---- arguments
    integer,            intent(in)    :: NNMAX, gridtype
    complex(kind(0d0)), intent(inout) :: cn(:)

    !---- locals
    integer, parameter :: dp = kind(0d0)
    integer :: N, j, k
    real(dp),    parameter :: myPI = 4.0_dp*atan(1.0_dp)
    complex(dp), parameter :: I1   = cmplx(0.0_dp, 1.0_dp, kind=dp)
    complex(dp) :: tmp

    !---- basic size check (cn must be indexed j = k+NNMAX+1 for k=-NNMAX..NNMAX)
    if (size(cn) < 2*NNMAX + 1) then
        stop 'fshisfts: cn is too small; needs size >= 2*NNMAX+1'
    end if

    N = NNMAX + 2

    if (gridtype == 11) then
        ! tmp = ((2*N - 2)*pi*i)/(4*N)   (cast N to real to avoid integer division)
        tmp = ((2.0_dp*real(N,dp) - 2.0_dp) * myPI * I1) / (4.0_dp*real(N,dp))

!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(STATIC)
        do k = -NNMAX, NNMAX
            j     = k + NNMAX + 1
            cn(j) = exp(real(k,dp) * tmp) * ((-1.0_dp)**k)
        end do
!$OMP END PARALLEL DO

    else

!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(STATIC)
        do k = -NNMAX, NNMAX
            j     = k + NNMAX + 1
            cn(j) = exp(I1 * myPI * real(k,dp) / 2.0_dp) * ((-1.0_dp)**k)
        end do
!$OMP END PARALLEL DO

    end if
end subroutine fshisfts
!   
subroutine newdfs_complex_uv(u, v, coeffs_u, coeffs_v, NNMAX, cn)
    implicit none
    !---- arguments
    integer,            intent(in)    :: NNMAX
    real(8),            intent(inout) :: u(:,:), v(:,:)
    complex(8),         intent(inout) :: coeffs_u(:), coeffs_v(:)
    complex(8),         intent(in)    :: cn(:)     ! used as multiplier -> must be IN, not OUT

    !---- constants
    real(8), parameter :: ER     = 6.37122d6      ! Earth radius [m]
    real(8), parameter :: ER_INV = 1.0d0/ER
    complex(8), parameter :: I1  = cmplx(0.0d0, 1.0d0, kind=8)

    !---- locals
    integer :: mm, nn, MNMAX, MNUM
    integer :: m, j, mmax, jcn_dfs
    real(8) :: am, an

    ! work arrays
    real(8),    allocatable :: u_c(:), u_s(:), v_c(:), v_s(:)
    complex(8), allocatable :: qu(:,:), qv(:,:), Tc(:,:), Ts(:,:), Tcoeff(:,:)
    complex(8), allocatable :: Vcoeff(:,:), Vc(:,:), Vs(:,:)
    complex(8), allocatable :: qpsi(:,:), qchi(:,:)   ! stream function, velocity potential

    integer :: need_size

    !------------------------------
    ! 1) FFT in longitude (Eq. 56)
    !------------------------------
    call fft__towave(u)
    call fft__towave(v)

    MNMAX   = NNMAX
    MNUM    = NNMAX + 1
    mmax    = NNMAX
    jcn_dfs = 2

    !----------------------------------------------
    ! 2) Temporary coefficients in latitude (Eq.57)
    !----------------------------------------------
    allocate(qu(NNMAX+1, 0:MNMAX), qv(NNMAX+1, 0:MNMAX))
    allocate(qchi(NNMAX+1, 0:MNMAX), qpsi(NNMAX+1, 0:MNMAX))

    call fft_y__g2w_uv(u, qu)
    call fft_y__g2w_uv(v, qv)

    !--------------------------------------------------------
    ! 3) Velocity potential (chi) and streamfunction (psi)
    !    from u,v (Eqs. 63–64)
    !--------------------------------------------------------
    call uv2divrot_dfs__uv2chipsi(qu, qv, qchi, qpsi)

    !--------------------------------------------------------
    ! 4) Fourier coeffs of velocity from chi/psi (Eq. 53 etc.)
    !--------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m, j, am, an)
!$OMP DO SCHEDULE(STATIC)
    do m = 0, mmax
        am = real(m, kind=8)

        if ( jcn_dfs >= 2 .and. m >= 3 .and. mod(m,2) == 1 ) then
            j = 1;  an = j - 1.0d0
            qu(j,m) = ER_INV*0.5d0*I1*am*qchi(j,m) - ER_INV*0.25d0*an*qpsi(j+1,m)     ! n=0
            qv(j,m) = ER_INV*0.5d0*I1*am*qpsi(j,m) + ER_INV*0.25d0*an*qchi(j+1,m)     ! n=0

            j = 2;  an = j - 1.0d0
            qu(j,m) = ER_INV*0.5d0*I1*am*qchi(j,m) + ER_INV*0.25d0*an*( 3.0d0*qpsi(j-1,m) - qpsi(j+1,m) ) ! n=1
            qv(j,m) = ER_INV*0.5d0*I1*am*qpsi(j,m) + ER_INV*0.25d0*an*( -3.0d0*qchi(j-1,m) + qchi(j+1,m) ) ! n=1

            j = 3;  an = j - 1.0d0
            qu(j,m) = ER_INV*0.5d0*I1*am*( -qchi(j-2,m) + qchi(j,m) ) + ER_INV*0.25d0*an*( 2.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
            qv(j,m) = ER_INV*0.5d0*I1*am*( -qpsi(j-2,m) + qpsi(j,m) ) + ER_INV*0.25d0*an*( -2.0d0*qchi(j-1,m) + qchi(j+1,m) )

            do j = 4, MNUM-1
                an = j - 1.0d0
                qu(j,m) = ER_INV*0.5d0*I1*am*( -qchi(j-2,m) + qchi(j,m) ) &
                        + ER_INV*0.25d0*an*( -qpsi(j-3,m) + 2.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
                qv(j,m) = ER_INV*0.5d0*I1*am*( -qpsi(j-2,m) + qpsi(j,m) ) &
                        + ER_INV*0.25d0*an*(  qchi(j-3,m) - 2.0d0*qchi(j-1,m) + qchi(j+1,m) )
            end do

            j = MNUM;  an = j - 1.0d0
            qu(j,m) = ER_INV*0.5d0*I1*am*( -qchi(j-2,m) + qchi(j,m) ) &
                    + ER_INV*0.25d0*an*( -qpsi(j-3,m) + 2.0d0*qpsi(j-1,m) )
            qv(j,m) = ER_INV*0.5d0*I1*am*( -qpsi(j-2,m) + qpsi(j,m) ) &
                    + ER_INV*0.25d0*an*(  qchi(j-3,m) - 2.0d0*qchi(j-1,m) )

        else if ( m == 0 ) then
            do j = 1, MNUM-1
                an = real(j, kind=8)
                qu(j,m) = -ER_INV*an*qpsi(j+1,m)
                qv(j,m) =  ER_INV*an*qchi(j+1,m)
            end do
            j = MNUM
            qu(j,m) = 0.0d0
            qv(j,m) = 0.0d0

        else if ( mod(m,2) == 1 ) then
            j = 1; an = j - 1.0d0
            qu(j,m) = ER_INV*I1*am*qchi(j,m)   ! n=0
            qv(j,m) = ER_INV*I1*am*qpsi(j,m)   ! n=0

            j = 2; an = j - 1.0d0
            qu(j,m) = ER_INV*I1*am*qchi(j,m) + ER_INV*0.5d0*an*(  2.0d0*qpsi(j-1,m) - qpsi(j+1,m) )
            qv(j,m) = ER_INV*I1*am*qpsi(j,m) + ER_INV*0.5d0*an*( -2.0d0*qchi(j-1,m) + qchi(j+1,m) )

            do j = 3, MNUM-1
                an = j - 1.0d0
                qu(j,m) = ER_INV*I1*am*qchi(j,m) + ER_INV*0.5d0*an*(  qpsi(j-1,m) - qpsi(j+1,m) )
                qv(j,m) = ER_INV*I1*am*qpsi(j,m) + ER_INV*0.5d0*an*( -qchi(j-1,m) + qchi(j+1,m) )
            end do

            j = MNUM; an = j - 1.0d0
            qu(j,m) = ER_INV*I1*am*qchi(j,m) + ER_INV*0.5d0*an*qpsi(j-1,m)
            qv(j,m) = ER_INV*I1*am*qpsi(j,m) - ER_INV*0.5d0*an*qchi(j-1,m)

        else   ! mod(m,2) == 0
            j = 1; an = real(j,kind=8)
            qu(j,m) = ER_INV*I1*am*qchi(j,m) - ER_INV*0.5d0*an*qpsi(j+1,m)
            qv(j,m) = ER_INV*I1*am*qpsi(j,m) + ER_INV*0.5d0*an*qchi(j+1,m)

            do j = 2, MNUM-1
                an = real(j,kind=8)
                qu(j,m) = ER_INV*I1*am*qchi(j,m) + ER_INV*0.5d0*an*(  qpsi(j-1,m) - qpsi(j+1,m) )
                qv(j,m) = ER_INV*I1*am*qpsi(j,m) + ER_INV*0.5d0*an*( -qchi(j-1,m) + qchi(j+1,m) )
            end do

            j = MNUM; an = real(j,kind=8)
            qu(j,m) = ER_INV*I1*am*qchi(j,m) + ER_INV*0.5d0*an*qpsi(j-1,m)
            qv(j,m) = ER_INV*I1*am*qpsi(j,m) - ER_INV*0.5d0*an*qchi(j-1,m)
        end if
    end do
!$OMP END DO
!$OMP END PARALLEL

    !-------------------------------------------
    ! 5) Convert to complex basis (Tc/Ts, Vc/Vs)
    !-------------------------------------------
    allocate(Tc(2*NNMAX+1, MNMAX+1), Ts(2*NNMAX+1, MNMAX+1))
    allocate(Vc(2*NNMAX+1, MNMAX+1), Vs(2*NNMAX+1, MNMAX+1))
    Tc = (0.0d0, 0.0d0); Ts = (0.0d0, 0.0d0)
    Vc = (0.0d0, 0.0d0); Vs = (0.0d0, 0.0d0)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mm, nn, u_c, v_c, u_s, v_s)
    allocate(u_c(NNMAX+1), v_c(NNMAX+1), u_s(NNMAX+1), v_s(NNMAX+1))
!$OMP DO SCHEDULE(STATIC)
    do mm = 0, MNMAX
        do nn = 1, NNMAX+1
            u_c(nn) = real( qu(nn, mm), kind=8 )
            u_s(nn) = aimag(qu(nn, mm))
            v_c(nn) = real( qv(nn, mm), kind=8 )
            v_s(nn) = aimag(qv(nn, mm))
        end do

        if (mm == 0) then
            Tc(NNMAX:1:-1,        mm+1) = -0.5d0*I1*u_c(1:NNMAX)
            Tc(NNMAX+2:2*NNMAX+1, mm+1) =  0.5d0*I1*u_c(1:NNMAX)

            Vc(NNMAX:1:-1,        mm+1) = -0.5d0*I1*v_c(1:NNMAX)
            Vc(NNMAX+2:2*NNMAX+1, mm+1) =  0.5d0*I1*v_c(1:NNMAX)

        else if (mod(mm, 2) == 1) then
            Tc(NNMAX:1:-1,        mm+1) = 0.5d0 * u_c(2:NNMAX+1)
            Tc(NNMAX+2:2*NNMAX+1, mm+1) = 0.5d0 * u_c(2:NNMAX+1)
            Tc(NNMAX+1,           mm+1) = u_c(1)

            Vc(NNMAX:1:-1,        mm+1) = 0.5d0 * v_c(2:NNMAX+1)
            Vc(NNMAX+2:2*NNMAX+1, mm+1) = 0.5d0 * v_c(2:NNMAX+1)
            Vc(NNMAX+1,           mm+1) = v_c(1)

            Ts(NNMAX:1:-1,        mm+1) = 0.5d0 * u_s(2:NNMAX+1)
            Ts(NNMAX+2:2*NNMAX+1, mm+1) = 0.5d0 * u_s(2:NNMAX+1)
            Ts(NNMAX+1,           mm+1) = u_s(1)

            Vs(NNMAX:1:-1,        mm+1) = 0.5d0 * v_s(2:NNMAX+1)
            Vs(NNMAX+2:2*NNMAX+1, mm+1) = 0.5d0 * v_s(2:NNMAX+1)
            Vs(NNMAX+1,           mm+1) = v_s(1)

        else
            Tc(NNMAX:1:-1,        mm+1) = -0.5d0*I1*u_c(1:NNMAX)
            Tc(NNMAX+2:2*NNMAX+1, mm+1) =  0.5d0*I1*u_c(1:NNMAX)

            Vc(NNMAX:1:-1,        mm+1) = -0.5d0*I1*v_c(1:NNMAX)
            Vc(NNMAX+2:2*NNMAX+1, mm+1) =  0.5d0*I1*v_c(1:NNMAX)

            Ts(NNMAX:1:-1,        mm+1) = -0.5d0*I1*u_s(1:NNMAX)
            Ts(NNMAX+2:2*NNMAX+1, mm+1) =  0.5d0*I1*u_s(1:NNMAX)

            Vs(NNMAX:1:-1,        mm+1) = -0.5d0*I1*v_s(1:NNMAX)
            Vs(NNMAX+2:2*NNMAX+1, mm+1) =  0.5d0*I1*v_s(1:NNMAX)
        end if
    end do
!$OMP END DO
    deallocate(u_c, v_c, u_s, v_s)
!$OMP END PARALLEL

    !--------------------------------------------
    ! 6) Combine longitude direction into coeffs
    !--------------------------------------------
    allocate(Tcoeff(2*NNMAX+1, 2*MNMAX+1))
    allocate(Vcoeff(2*NNMAX+1, 2*MNMAX+1))

!$OMP PARALLEL DO PRIVATE(nn, mm) SCHEDULE(STATIC)
    do nn = 1, 2*NNMAX+1
        Tcoeff(nn, MNMAX+1) = Tc(nn, 1) * cn(nn)
        Vcoeff(nn, MNMAX+1) = Vc(nn, 1) * cn(nn)
        do mm = 2, NNMAX+1
            Tcoeff(nn, mm+MNMAX   ) = (Tc(nn, mm) + I1*Ts(nn, mm)) * cn(nn)  ! plus
            Tcoeff(nn, MNMAX+2-mm ) = (Tc(nn, mm) - I1*Ts(nn, mm)) * cn(nn)  ! minus

            Vcoeff(nn, mm+MNMAX   ) = (Vc(nn, mm) + I1*Vs(nn, mm)) * cn(nn)  ! plus
            Vcoeff(nn, MNMAX+2-mm ) = (Vc(nn, mm) - I1*Vs(nn, mm)) * cn(nn)  ! minus
        end do
    end do
!$OMP END PARALLEL DO

    !--------------------------------------------
    ! 7) Pack to 1D (with basic size checks)
    !--------------------------------------------
    need_size = (2*NNMAX+1)*(2*MNMAX+1)
    if (size(coeffs_u) /= need_size) then
        stop 'newdfs_complex_uv: coeffs_u has incorrect size.'
    end if
    if (size(coeffs_v) /= need_size) then
        stop 'newdfs_complex_uv: coeffs_v has incorrect size.'
    end if
    if (size(cn) < 2*NNMAX+1) then
        stop 'newdfs_complex_uv: cn is too small.'
    end if

    coeffs_u = reshape(Tcoeff, [need_size])
    coeffs_v = reshape(Vcoeff, [need_size])

    !--------------------------------------------
    ! 8) Cleanup
    !--------------------------------------------
    deallocate(Tc, Ts, Vc, Vs, Tcoeff, Vcoeff)
    deallocate(qu, qv, qchi, qpsi)
end subroutine newdfs_complex_uv
!  
end module new_dfs