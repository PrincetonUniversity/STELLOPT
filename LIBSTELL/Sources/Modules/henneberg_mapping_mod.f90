module henneberg_mapping_mod
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    ! =========================================================================
    ! Transformation from VMEC to Henneberg
    ! =========================================================================
    subroutine vmec_to_henneberg(mpol, ntor, rc, zs, nfp, alpha_fac, &
                                         mmax, nmax, ntheta, nphi, &
                                         R0_out, Z0_out, bcoef_out, rho_out)
        integer, intent(in)  :: mpol, ntor, nfp, alpha_fac
        integer, intent(in)  :: mmax, nmax, ntheta, nphi
        ! Fortran arrays can take negative indices, matching the [-ntor, ntor] physics naturally
        !real(dp), intent(in) :: rc(0:mpol, -ntor:ntor)
        !real(dp), intent(in) :: zs(0:mpol, -ntor:ntor)
        real(dp), intent(in) :: rc(-ntor:ntor, 0:mpol)
        real(dp), intent(in) :: zs(-ntor:ntor, 0:mpol)
        
        real(dp), intent(out) :: R0_out(0:nmax)
        real(dp), intent(out) :: Z0_out(0:nmax)
        real(dp), intent(out) :: bcoef_out(0:nmax)
        real(dp), intent(out) :: rho_out(-ntor:ntor, 0:mpol)

        ! Local scalars
        integer  :: m, n, jphi, itheta, nmin, sign_flips_count
        real(dp) :: alpha, phi0, theta0, cosaphi, sinaphi
        real(dp) :: min_for_b, max_for_b, b, Q, angle, Rmn, Zmn
        real(dp) :: R0H, Z0H, avg_R, avg_Z, Z00H, rho00, s
        
        ! Automatic arrays (avoids manual allocation/deallocation)
        real(dp) :: theta(ntheta), phi(nphi)
        real(dp) :: R(ntheta), Z(ntheta)
        real(dp) :: d_R_d_theta(ntheta), d_Z_d_theta(ntheta)
        real(dp) :: d_Z_rot_d_theta(ntheta), d_Z_rot_circ(ntheta)
        real(dp) :: temp_arr(ntheta), arcsin_term(ntheta), theta_H(ntheta)
        real(dp) :: R_H(ntheta), Z_H(ntheta)
        real(dp) :: R0_realsp(nphi), Z0_realsp(nphi), b_realsp(nphi)
        real(dp) :: rho_realsp(ntheta, nphi)
        real(dp) :: cosnphi_fac(nphi), sinnphi_fac(nphi)
        
        ! Arrays for cubic interpolation (concatenated times 3)
        real(dp) :: theta_H_3(3*ntheta), R_3(3*ntheta), Z_3(3*ntheta)

        ! Basic input validation checks
        if (alpha_fac < -1 .or. alpha_fac > 1) stop "alpha_fac needs to be -1, 0 or +1"
        if (mmax > mpol) stop "mmax is not allowed to be bigger than mpol"
        if (nmax > ntor) stop "nmax is not allowed to be bigger than ntor"
        if (ntheta <= 0 .or. nphi <= 0) stop "ntheta and nphi have to be > 0"

        alpha = 0.5_dp * real(nfp, dp) * real(alpha_fac, dp)

        ! Initialize grid
        do itheta = 1, ntheta
            theta(itheta) = 2.0_dp * pi * real(itheta - 1, dp) / real(ntheta, dp)
        end do
        do jphi = 1, nphi
            phi(jphi) = 2.0_dp * pi * real(jphi - 1, dp) / real(nfp * nphi, dp)
        end do

        R0_realsp = 0.0_dp
        Z0_realsp = 0.0_dp
        b_realsp  = 0.0_dp
        rho_realsp = 0.0_dp

        ! Loop over toroidal angle
        do jphi = 1, nphi
            phi0 = phi(jphi)
            cosaphi = cos(alpha * phi0)
            sinaphi = sin(alpha * phi0)

            ! 1. Optimize min and max of b function using embedded line search
            call find_b_extrema(mpol, ntor, rc, zs, nfp, phi0, cosaphi, sinaphi, &
                                min_for_b, max_for_b)

            b = 0.5_dp * (max_for_b - min_for_b)
            Q = 0.5_dp * (max_for_b + min_for_b)

            ! 2. Compute R, Z and their derivatives
            R = 0.0_dp
            Z = 0.0_dp
            d_R_d_theta = 0.0_dp
            d_Z_d_theta = 0.0_dp

            do m = 0, mpol
                do n = -ntor, ntor
                    Rmn = rc(n, m)
                    Zmn = zs(n, m)
                    do itheta = 1, ntheta
                        angle = real(m, dp) * theta(itheta) - real(n * nfp, dp) * phi0
                        R(itheta) = R(itheta) + Rmn * cos(angle)
                        Z(itheta) = Z(itheta) + Zmn * sin(angle)
                        d_R_d_theta(itheta) = d_R_d_theta(itheta) - Rmn * real(m, dp) * sin(angle)
                        d_Z_d_theta(itheta) = d_Z_d_theta(itheta) + Zmn * real(m, dp) * cos(angle)
                    end do
                end do
            end do

            ! 3. Rotated derivative and sign flips
            sign_flips_count = 0
            do itheta = 1, ntheta
                d_Z_rot_d_theta(itheta) = d_Z_d_theta(itheta) * cosaphi - d_R_d_theta(itheta) * sinaphi
                d_Z_rot_circ(itheta) = d_Z_rot_d_theta(itheta)
            end do
            
            ! Count sign flips sequentially (circular)
            do itheta = 1, ntheta - 1
                if (d_Z_rot_circ(itheta+1) * d_Z_rot_circ(itheta) < 0.0_dp) &
                    sign_flips_count = sign_flips_count + 1
            end do
            if (d_Z_rot_circ(1) * d_Z_rot_circ(ntheta) < 0.0_dp) &
                sign_flips_count = sign_flips_count + 1

            if (sign_flips_count /= 2) then
                print *, "Warning: ", sign_flips_count, " > 2 sign flips, shape maybe not representable in Henneberg coords."
            end if

            ! 4. Compute temp, clamp domains, and apply arcsin
            do itheta = 1, ntheta
                temp_arr(itheta) = (Z(itheta) * cosaphi - R(itheta) * sinaphi - Q) / b
                
                if (temp_arr(itheta) > 1.0_dp) then
                    if (temp_arr(itheta) > 1.0_dp + 1.0d-12) print *, "Warning: arcsin input > 1."
                    temp_arr(itheta) = 1.0_dp
                else if (temp_arr(itheta) < -1.0_dp) then
                    if (temp_arr(itheta) < -1.0_dp - 1.0d-12) print *, "Warning: arcsin input < -1."
                    temp_arr(itheta) = -1.0_dp
                end if

                arcsin_term(itheta) = asin(temp_arr(itheta))

                ! Apply masks based on derivative and negative values
                if (d_Z_rot_d_theta(itheta) < 0.0_dp) then
                    arcsin_term(itheta) = pi - arcsin_term(itheta)
                end if
                if (arcsin_term(itheta) < 0.0_dp) then
                    arcsin_term(itheta) = arcsin_term(itheta) + 2.0_dp * pi
                end if

                theta_H(itheta) = arcsin_term(itheta) + alpha * phi0
            end do

            ! 5. Build concatenated arrays for interpolation
            do itheta = 1, ntheta
                theta_H_3(itheta)            = theta_H(itheta) - 2.0_dp * pi
                theta_H_3(itheta + ntheta)   = theta_H(itheta)
                theta_H_3(itheta + 2*ntheta) = theta_H(itheta) + 2.0_dp * pi

                R_3(itheta)            = R(itheta)
                R_3(itheta + ntheta)   = R(itheta)
                R_3(itheta + 2*ntheta) = R(itheta)

                Z_3(itheta)            = Z(itheta)
                Z_3(itheta + ntheta)   = Z(itheta)
                Z_3(itheta + 2*ntheta) = Z(itheta)
            end do

            ! 6. Perform Cubic Spline Interpolation over the theta grid
            call interp_cubic(3*ntheta, theta_H_3, R_3, ntheta, theta, R_H)
            call interp_cubic(3*ntheta, theta_H_3, Z_3, ntheta, theta, Z_H)

            ! 7. Gather results mapped to Henneberg
            avg_R = sum(R_H) / real(ntheta, dp)
            avg_Z = sum(Z_H) / real(ntheta, dp)

            R0H = avg_R * cosaphi * cosaphi + avg_Z * sinaphi * cosaphi - Q * sinaphi
            Z0H = avg_R * cosaphi * sinaphi + avg_Z * sinaphi * sinaphi + Q * cosaphi

            R0_realsp(jphi) = R0H
            Z0_realsp(jphi) = Z0H
            b_realsp(jphi)  = b

            do itheta = 1, ntheta
                rho_realsp(itheta, jphi) = (R_H(itheta) - R0H) * cosaphi + &
                                           (Z_H(itheta) - Z0H) * sinaphi
            end do
        end do

        ! Final projection: Realspace to Fourier basis
        R0_out = 0.0_dp
        Z0_out = 0.0_dp
        bcoef_out = 0.0_dp

        R0_out(0) = sum(R0_realsp) / real(nphi, dp)
        bcoef_out(0) = sum(b_realsp) / real(nphi, dp)
        
        Z00H = sum(Z0_realsp) / real(nphi, dp)
        rho00 = sum(rho_realsp) / real(ntheta * nphi, dp)

        if (abs(Z00H) > 1.0d-6) print *, "Warning: Z0H at n=0 is ", Z00H, " but should be near 0"
        if (abs(rho00) > 1.0d-6) print *, "Warning: rho at m=0,n=0 is ", rho00, " but should be near 0"

        do n = 1, nmax
            do jphi = 1, nphi
                cosnphi_fac(jphi) = cos(real(n * nfp, dp) * phi(jphi)) / real(nphi, dp)
                sinnphi_fac(jphi) = sin(real(n * nfp, dp) * phi(jphi)) / real(nphi, dp)
            end do

            R0_out(n) = 2.0_dp * dot_product(R0_realsp, cosnphi_fac)
            Z0_out(n) = 2.0_dp * dot_product(Z0_realsp, sinnphi_fac)
            bcoef_out(n) = 2.0_dp * dot_product(b_realsp, cosnphi_fac)
        end do

        rho_out = 0.0_dp
        do m = 0, mmax
            nmin = -nmax
            if (m == 0) nmin = 1
            do n = nmin, nmax
                s = 0.0_dp
                do jphi = 1, nphi
                    phi0 = phi(jphi)
                    do itheta = 1, ntheta
                        theta0 = theta(itheta)
                        angle = real(m, dp) * theta0 + (real(n * nfp, dp) - alpha) * phi0
                        s = s + rho_realsp(itheta, jphi) * cos(angle)
                    end do
                end do
                rho_out(n, m) = 2.0_dp * s / real(ntheta * nphi, dp)
            end do
        end do

    end subroutine vmec_to_henneberg


    ! =========================================================================
    ! Transformation from Henneberg to VMEC
    ! =========================================================================
    subroutine henneberg_to_vmec(nmax, mmax, R0nH, Z0nH, bn, rhomn, &
                                         alpha_fac, ntor_out, rc, zs)
        integer, intent(in) :: nmax, mmax, alpha_fac
        real(dp), intent(in) :: R0nH(0:nmax)
        real(dp), intent(in) :: Z0nH(0:nmax)
        real(dp), intent(in) :: bn(0:nmax)
        !real(dp), intent(in) :: rhomn(0:mmax, -nmax:nmax)
        real(dp), intent(in) :: rhomn(-nmax:nmax, 0:mmax)
        
        integer, intent(out) :: ntor_out
        ! Must be sized to handle the resulting offset logic natively
        !real(dp), intent(out) :: rc(0:mmax, -(nmax + abs(alpha_fac)) : (nmax + abs(alpha_fac)))
        !real(dp), intent(out) :: zs(0:mmax, -(nmax + abs(alpha_fac)) : (nmax + abs(alpha_fac)))
        real(dp), intent(out) :: rc(-(nmax + abs(alpha_fac)) : (nmax + abs(alpha_fac)), 0:mmax)
        real(dp), intent(out) :: zs(-(nmax + abs(alpha_fac)) : (nmax + abs(alpha_fac)), 0:mmax)

        integer  :: m, n, nh, n_start
        real(dp) :: b, rho

        if (alpha_fac < -1 .or. alpha_fac > 1) stop "alpha_fac has to be -1, 0 or +1"

        ntor_out = nmax + abs(alpha_fac)

        rc = 0.0_dp
        zs = 0.0_dp

        ! Transform R0 and Z0
        do n = 0, nmax
            rc(n, 0) = rc(n, 0) + R0nH(n)
            if (n > 0) then
                zs(n, 0) = zs(n, 0) - Z0nH(n)
            end if
        end do

        ! Transform b
        do n = 0, nmax
            b = bn(n)
            rc( n, 1) = rc( n, 1) + 0.25_dp * b
            rc(-n, 1) = rc(-n, 1) + 0.25_dp * b
            rc( n + alpha_fac, 1) = rc( n + alpha_fac, 1) - 0.25_dp * b
            rc(-n + alpha_fac, 1) = rc(-n + alpha_fac, 1) - 0.25_dp * b

            zs( n, 1) = zs( n, 1) + 0.25_dp * b
            zs(-n, 1) = zs(-n, 1) + 0.25_dp * b
            zs( n + alpha_fac, 1) = zs( n + alpha_fac, 1) - 0.25_dp * b
            zs(-n + alpha_fac, 1) = zs(-n + alpha_fac, 1) - 0.25_dp * b
        end do

        ! Transform rho
        do m = 0, mmax
            if (m == 0) then
                n_start = 1
            else
                n_start = -nmax
            end if

            do nh = n_start, nmax
                rho = rhomn(nh, m)

                rc(-nh, m) = rc(-nh, m) + 0.5_dp * rho
                rc(-nh + alpha_fac, m) = rc(-nh + alpha_fac, m) + 0.5_dp * rho

                zs(-nh, m) = zs(-nh, m) + 0.5_dp * rho
                zs(-nh + alpha_fac, m) = zs(-nh + alpha_fac, m) - 0.5_dp * rho
            end do
        end do

    end subroutine henneberg_to_vmec


    ! =========================================================================
    ! Helper 1: Function to compute the scalar b_min equivalence
    ! =========================================================================
    real(dp) function evaluate_b_zeta(theta0, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp)
        real(dp), intent(in) :: theta0, phi0, cosaphi, sinaphi
        integer, intent(in)  :: mpol, ntor, nfp
        !real(dp), intent(in) :: rc(0:mpol, -ntor:ntor)
        !real(dp), intent(in) :: zs(0:mpol, -ntor:ntor)
        real(dp), intent(in) :: rc(-ntor:ntor, 0:mpol)
        real(dp), intent(in) :: zs(-ntor:ntor, 0:mpol)

        real(dp) :: R, Z, angle
        integer  :: m, n

        R = 0.0_dp
        Z = 0.0_dp

        do m = 0, mpol
            do n = -ntor, ntor
                angle = real(m, dp) * theta0 - real(n * nfp, dp) * phi0
                R = R + rc(n, m) * cos(angle)
                Z = Z + zs(n, m) * sin(angle)
            end do
        end do

        evaluate_b_zeta = Z * cosaphi - R * sinaphi
    end function evaluate_b_zeta


    ! =========================================================================
    ! Helper 2: Optimize b function to find minimum and maximum over theta
    ! Equivalent to scipy.optimize.minimize_scalar (Grid search + Local Polish)
    ! =========================================================================
    subroutine find_b_extrema(mpol, ntor, rc, zs, nfp, phi0, cosaphi, sinaphi, min_val, max_val)
        integer, intent(in)  :: mpol, ntor, nfp
        !real(dp), intent(in) :: rc(0:mpol, -ntor:ntor)
        !real(dp), intent(in) :: zs(0:mpol, -ntor:ntor)
        real(dp), intent(in) :: rc(-ntor:ntor, 0:mpol)
        real(dp), intent(in) :: zs(-ntor:ntor, 0:mpol)
        real(dp), intent(in) :: phi0, cosaphi, sinaphi
        real(dp), intent(out):: min_val, max_val

        integer, parameter :: N_GRID = 100
        real(dp) :: grid_val, grid_theta, f_val
        real(dp) :: min_theta, max_theta, bracket_low, bracket_high
        integer  :: i

        min_val = 1.0d30
        max_val = -1.0d30
        min_theta = 0.0_dp
        max_theta = 0.0_dp

        ! 1. Perform Grid Search to bracket minimum and maximum
        do i = 0, N_GRID - 1
            grid_theta = 2.0_dp * pi * real(i, dp) / real(N_GRID, dp)
            f_val = evaluate_b_zeta(grid_theta, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp)
            if (f_val < min_val) then
                min_val = f_val
                min_theta = grid_theta
            end if
            if (f_val > max_val) then
                max_val = f_val
                max_theta = grid_theta
            end if
        end do

        ! 2. Local refinement can be inserted here if extreme precision is necessary
        ! Due to the harmonic nature, a dense grid is highly accurate.
        ! Golden section polishing placeholder bounds
        bracket_low = min_theta - (2.0_dp * pi / N_GRID)
        bracket_high = min_theta + (2.0_dp * pi / N_GRID)
        call golden_section(bracket_low, bracket_high, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp, 1.0_dp, min_val)

        bracket_low = max_theta - (2.0_dp * pi / N_GRID)
        bracket_high = max_theta + (2.0_dp * pi / N_GRID)
        call golden_section(bracket_low, bracket_high, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp, -1.0_dp, max_val)
        max_val = -max_val ! Revert the sign flip from optimization
    end subroutine find_b_extrema

    ! Basic Golden Section local polisher
    subroutine golden_section(a, b, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp, sign_mod, extremum)
        real(dp), intent(in) :: a, b, phi0, cosaphi, sinaphi, sign_mod
        integer, intent(in)  :: mpol, ntor, nfp
        !real(dp), intent(in) :: rc(0:mpol, -ntor:ntor), zs(0:mpol, -ntor:ntor)
        real(dp), intent(in) :: rc(-ntor:ntor, 0:mpol), zs(-ntor:ntor, 0:mpol)
        real(dp), intent(out):: extremum
        
        real(dp), parameter :: invphi = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp
        real(dp), parameter :: tol = 1.0d-12
        real(dp) :: x1, x2, f1, f2, aa, bb
        
        aa = a
        bb = b
        x1 = bb - invphi * (bb - aa)
        x2 = aa + invphi * (bb - aa)
        f1 = sign_mod * evaluate_b_zeta(x1, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp)
        f2 = sign_mod * evaluate_b_zeta(x2, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp)

        do while (abs(bb - aa) > tol)
            if (f1 < f2) then
                bb = x2
                x2 = x1
                f2 = f1
                x1 = bb - invphi * (bb - aa)
                f1 = sign_mod * evaluate_b_zeta(x1, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp)
            else
                aa = x1
                x1 = x2
                f1 = f2
                x2 = aa + invphi * (bb - aa)
                f2 = sign_mod * evaluate_b_zeta(x2, phi0, cosaphi, sinaphi, mpol, ntor, rc, zs, nfp)
            end if
        end do
        extremum = min(f1, f2)
    end subroutine golden_section


    ! =========================================================================
    ! Helper 3: Cubic Spline interpolation routine
    ! Replicates scipy.interpolate.interp1d(kind="cubic")
    ! =========================================================================
    subroutine interp_cubic(n, x, y, n_out, x_out, y_out)
        integer, intent(in) :: n, n_out
        real(dp), intent(in) :: x(n), y(n), x_out(n_out)
        real(dp), intent(out) :: y_out(n_out)

        ! Automatic arrays for spline coefficients
        real(dp) :: h(n-1), alpha_spl(n-1)
        real(dp) :: l(n), mu(n), z(n), c(n), b(n), d(n)
        integer  :: i, j
        real(dp) :: dx

        ! Compute natural spline coefficients
        do i = 1, n - 1
            h(i) = x(i+1) - x(i)
        end do

        do i = 2, n - 1
            alpha_spl(i) = 3.0_dp / h(i) * (y(i+1) - y(i)) - 3.0_dp / h(i-1) * (y(i) - y(i-1))
        end do

        l(1) = 1.0_dp
        mu(1) = 0.0_dp
        z(1) = 0.0_dp

        do i = 2, n - 1
            l(i) = 2.0_dp * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1)
            mu(i) = h(i) / l(i)
            z(i) = (alpha_spl(i) - h(i-1) * z(i-1)) / l(i)
        end do

        l(n) = 1.0_dp
        z(n) = 0.0_dp
        c(n) = 0.0_dp

        do j = n - 1, 1, -1
            c(j) = z(j) - mu(j) * c(j+1)
            b(j) = (y(j+1) - y(j)) / h(j) - h(j) * (c(j+1) + 2.0_dp * c(j)) / 3.0_dp
            d(j) = (c(j+1) - c(j)) / (3.0_dp * h(j))
        end do

        ! Interpolate at requested points
        do i = 1, n_out
            ! Find bracket
            j = 1
            do while (j < n - 1 .and. x_out(i) > x(j+1))
                j = j + 1
            end do
            
            dx = x_out(i) - x(j)
            y_out(i) = y(j) + b(j)*dx + c(j)*(dx**2) + d(j)*(dx**3)
        end do

    end subroutine interp_cubic

end module henneberg_mapping_mod