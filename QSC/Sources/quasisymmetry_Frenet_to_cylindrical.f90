subroutine quasisymmetry_Frenet_to_cylindrical_linear

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: nu_1s, nu_1c, r2_source_1, r2_source_2
  real(dp), dimension(:), allocatable :: d_X1c_d_zeta, d_Y1c_d_zeta, d_Y1s_d_zeta

  R1c = (-binormal_cylindrical(:,3) * X1c_untwisted + normal_cylindrical(:,3) * Y1c_untwisted) * d_l_d_phi / R0
  R1s = (-binormal_cylindrical(:,3) * X1s_untwisted + normal_cylindrical(:,3) * Y1s_untwisted) * d_l_d_phi / R0
  Z1c = ( binormal_cylindrical(:,1) * X1c_untwisted - normal_cylindrical(:,1) * Y1c_untwisted) * d_l_d_phi / R0
  Z1s = ( binormal_cylindrical(:,1) * X1s_untwisted - normal_cylindrical(:,1) * Y1s_untwisted) * d_l_d_phi / R0

  if (.not. order_r_squared) return

  ! See the note "20190215-02 Frenet to cylindrical transformation to next order.docx" for derivation of the equations that follow:
  allocate(nu_1s(N_phi))
  allocate(nu_1c(N_phi))
  allocate(r2_source_1(N_phi))
  allocate(r2_source_2(N_phi))
  allocate(d_X1c_d_zeta(N_phi))
  allocate(d_Y1c_d_zeta(N_phi))
  allocate(d_Y1s_d_zeta(N_phi))
  nu_1s = B0_over_abs_G0 / d_l_d_phi * (R1s * R0p + z1s * z0p)
  nu_1c = B0_over_abs_G0 / d_l_d_phi * (R1c * R0p + z1c * z0p)
  d_X1c_d_zeta = matmul(d_d_zeta,X1c)
  d_Y1c_d_zeta = matmul(d_d_zeta,Y1c)
  d_Y1s_d_zeta = matmul(d_d_zeta,Y1s)

  ! r^2 terms that are independent of theta
  r2_source_1 = (0.25d+0) * (nu_1s * nu_1s + nu_1c * nu_1c) * abs_G0_over_B0 * abs_G0_over_B0 * curvature + (0.5d+0) * d_X1c_d_zeta * nu_1c &
       - (0.5d+0) * (Y1s * nu_1s + Y1c * nu_1c) * abs_G0_over_B0 * torsion + X20
  r2_source_2 = (0.5d+0) * abs_G0_over_B0 * torsion * X1c * nu_1c + (0.5d+0) * (d_Y1s_d_zeta * nu_1s + d_Y1c_d_zeta * nu_1c) + Y20
  R20             = (-binormal_cylindrical(:,3) * r2_source_1 + normal_cylindrical(:,3) * r2_source_2) * d_l_d_phi / R0
  z20_cylindrical = ( binormal_cylindrical(:,1) * r2_source_1 - normal_cylindrical(:,1) * r2_source_2) * d_l_d_phi / R0

  ! r^2 terms proportional to sin(2 theta)
  r2_source_1 = (0.25d+0) * (nu_1s * nu_1c + nu_1c * nu_1s) * abs_G0_over_B0 * abs_G0_over_B0 * curvature + (0.5d+0) * d_X1c_d_zeta * nu_1s &
       - (0.5d+0) * (Y1s * nu_1c + Y1c * nu_1s) * abs_G0_over_B0 * torsion + X2s
  r2_source_2 = (0.5d+0) * abs_G0_over_B0 * torsion * X1c * nu_1s + (0.5d+0) * (d_Y1s_d_zeta * nu_1c + d_Y1c_d_zeta * nu_1s) + Y2s
  R2s             = (-binormal_cylindrical(:,3) * r2_source_1 + normal_cylindrical(:,3) * r2_source_2) * d_l_d_phi / R0
  z2s_cylindrical = ( binormal_cylindrical(:,1) * r2_source_1 - normal_cylindrical(:,1) * r2_source_2) * d_l_d_phi / R0

  ! r^2 terms proportional to cos(2 theta)
  r2_source_1 = (0.25d+0) * (-nu_1s * nu_1s + nu_1c * nu_1c) * abs_G0_over_B0 * abs_G0_over_B0 * curvature + (0.5d+0) * d_X1c_d_zeta * nu_1c &
       - (0.5d+0) * (-Y1s * nu_1s + Y1c * nu_1c) * abs_G0_over_B0 * torsion + X2c
  r2_source_2 = (0.5d+0) * abs_G0_over_B0 * torsion * X1c * nu_1c + (0.5d+0) * (-d_Y1s_d_zeta * nu_1s + d_Y1c_d_zeta * nu_1c) + Y2c
  R2c             = (-binormal_cylindrical(:,3) * r2_source_1 + normal_cylindrical(:,3) * r2_source_2) * d_l_d_phi / R0
  z2c_cylindrical = ( binormal_cylindrical(:,1) * r2_source_1 - normal_cylindrical(:,1) * r2_source_2) * d_l_d_phi / R0


  print *,"nu_1s:",nu_1s
  print *,"nu_1c:",nu_1c
  print *,"r2_source_1:",r2_source_1
  print *,"r2_source_2:",r2_source_2
  print *,"d_X1c_d_zeta",d_X1c_d_zeta
  print *,"d_Y1c_d_zeta",d_Y1c_d_zeta
  print *,"d_Y1s_d_zeta",d_Y1s_d_zeta
  print *,"R20:",R20
  print *,"R2s:",R2s
  print *,"R2c:",R2c
  print *,"z20_cylindrical:",z20_cylindrical

  deallocate(nu_1s, nu_1c, r2_source_1, r2_source_2, d_X1c_d_zeta, d_Y1c_d_zeta, d_Y1s_d_zeta)

end subroutine quasisymmetry_Frenet_to_cylindrical_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module quasisymmetry_Frenet_to_cylindrical_mod

  ! This module is an implementation of the algorithm described in section 4.2 of
  ! the Landreman-Sengupta-Plunk Journal of Plasma Physics (2019) paper, to
  ! generate a finite-r surface from the near-axis expansion.

  ! The MATLAB version of this algorithm can be found in
  ! m20180501_01_quasisymmetricEquilibria_Frenet_altTransformation.m

  implicit none

  private

  public :: quasisymmetry_Frenet_to_cylindrical_nonlinear

contains

  subroutine quasisymmetry_Frenet_to_cylindrical_nonlinear

    use quasisymmetry_variables
    use vmec_input, only: vmec_nfp => nfp, lasym, mpol, ntor, RBC, RBS, ZBC, ZBS
    use vparams, only: ntord, mpol1d
    use quasisymmetry_splines

    implicit none

    integer :: N_theta, j_theta, N_phi_conversion, j_phi, i, m, n, nmin
    real(dp) :: costheta, sintheta, cos2theta, sin2theta, final_R, final_z
    real(dp), dimension(:,:), allocatable :: R_2D, z_2D, phi0_2D
    real(dp), dimension(:), allocatable :: theta, phi_conversion
    real(dp), dimension(:), allocatable :: X_at_this_theta, Y_at_this_theta, Z_at_this_theta
    type(periodic_spline) :: normal_R_spline, normal_phi_spline, normal_z_spline, binormal_R_spline, binormal_phi_spline, binormal_z_spline
    type(periodic_spline) :: tangent_R_spline, tangent_phi_spline, tangent_z_spline
    type(periodic_spline) :: X_spline, Y_spline, Z_spline, R0_spline, z0_spline
    real(dp) :: rootSolve_abserr, rootSolve_relerr, phi0_rootSolve_min, phi0_rootSolve_max
    real(dp) :: phi0_solution, phi_target, factor, factor2, angle, cosangle, sinangle
    integer :: fzeroFlag

    !----------------------------------------------

    N_theta = 20 ! Probably other values would work
    N_phi_conversion = N_phi
    allocate(theta(N_theta))
    allocate(phi_conversion(N_phi_conversion))
    allocate(phi0_2D(N_theta,N_phi_conversion))
    allocate(R_2D(N_theta,N_phi_conversion))
    allocate(z_2D(N_theta,N_phi_conversion))
    theta = [( 2*pi*i/N_theta, i=0,N_theta-1 )]
    phi_conversion = [( 2*pi*i/(nfp*N_phi_conversion), i=0,N_phi_conversion-1 )]

    allocate(X_at_this_theta(N_phi))
    allocate(Y_at_this_theta(N_phi))
    allocate(Z_at_this_theta(N_phi))

    call new_periodic_spline(N_phi, phi, R0, 2*pi/nfp, R0_spline)
    call new_periodic_spline(N_phi, phi, z0, 2*pi/nfp, z0_spline)
    call new_periodic_spline(N_phi, phi, normal_cylindrical(:,1), 2*pi/nfp, normal_R_spline)
    call new_periodic_spline(N_phi, phi, normal_cylindrical(:,2), 2*pi/nfp, normal_phi_spline)
    call new_periodic_spline(N_phi, phi, normal_cylindrical(:,3), 2*pi/nfp, normal_z_spline)
    call new_periodic_spline(N_phi, phi, binormal_cylindrical(:,1), 2*pi/nfp, binormal_R_spline)
    call new_periodic_spline(N_phi, phi, binormal_cylindrical(:,2), 2*pi/nfp, binormal_phi_spline)
    call new_periodic_spline(N_phi, phi, binormal_cylindrical(:,3), 2*pi/nfp, binormal_z_spline)
    if (order_r_squared) then
       call new_periodic_spline(N_phi, phi, tangent_cylindrical(:,1), 2*pi/nfp, tangent_R_spline)
       call new_periodic_spline(N_phi, phi, tangent_cylindrical(:,2), 2*pi/nfp, tangent_phi_spline)
       call new_periodic_spline(N_phi, phi, tangent_cylindrical(:,3), 2*pi/nfp, tangent_z_spline)
    end if

    rootSolve_abserr = 1.0e-10_dp
    rootSolve_relerr = 1.0e-10_dp

    do j_theta = 1, N_theta
       costheta = cos(theta(j_theta))
       sintheta = sin(theta(j_theta))
       X_at_this_theta = r * (X1c_untwisted * costheta + X1s_untwisted * sintheta)
       Y_at_this_theta = r * (Y1c_untwisted * costheta + Y1s_untwisted * sintheta)
       if (order_r_squared) then
          cos2theta = cos(2*theta(j_theta))
          sin2theta = sin(2*theta(j_theta))
          X_at_this_theta = X_at_this_theta + r*r*(X20_untwisted + X2c_untwisted * cos2theta + X2s_untwisted * sin2theta)
          Y_at_this_theta = Y_at_this_theta + r*r*(Y20_untwisted + Y2c_untwisted * cos2theta + Y2s_untwisted * sin2theta)
          Z_at_this_theta =                   r*r*(Z20_untwisted + Z2c_untwisted * cos2theta + Z2s_untwisted * sin2theta)
          call new_periodic_spline(N_phi, phi, Z_at_this_theta, 2*pi/nfp, Z_spline)
       end if
       call new_periodic_spline(N_phi, phi, X_at_this_theta, 2*pi/nfp, X_spline)
       call new_periodic_spline(N_phi, phi, Y_at_this_theta, 2*pi/nfp, Y_spline)
       do j_phi = 1, N_phi_conversion
          ! Solve for the phi0 such that r0 + X1 n + Y1 b has the desired phi

          phi_target = phi_conversion(j_phi)
          phi0_rootSolve_min = phi_target - 0.3
          phi0_rootSolve_max = phi_target + 0.3
          call quasisymmetry_fzero(Frenet_to_cylindrical_residual, phi0_rootSolve_min, phi0_rootSolve_max, phi_target, &
               rootSolve_relerr, rootSolve_abserr, fzeroFlag)
          ! Note: fzero returns its answer in phi0_rootSolve_min
          phi0_solution = phi0_rootSolve_min
          if (fzeroFlag == 4) then
             print *,"j_theta=",j_theta," j_phi=",j_phi
             stop "ERROR: fzero returned error 4: no sign change in residual"
          else if (fzeroFlag > 2) then
             print *,"WARNING in cosm: fzero returned an error code:",fzeroFlag
          end if

          call Frenet_to_cylindrical_1_point(phi0_solution, final_R, final_z)
          R_2D(j_theta,j_phi) = final_R
          z_2D(j_theta,j_phi) = final_z
          phi0_2D(j_theta,j_phi) = phi0_solution
       end do
       call delete_periodic_spline(X_spline)
       call delete_periodic_spline(Y_spline)
       if (order_r_squared) then
          call delete_periodic_spline(Z_spline)
       end if
    end do
    
!!$    print *,"N_theta:"
!!$    print *,N_theta
!!$    print *,"N_zeta:"
!!$    print *,N_phi_conversion
!!$    print *,"phi0_2D:"
!!$    do j_theta = 1,N_theta
!!$       !print "(*(f7.3))",phi0_2D(j_theta,:)
!!$       print "(*(f20.15))",phi0_2D(j_theta,:)
!!$    end do
!!$    print *,"R_2D:"
!!$    do j_theta = 1,N_theta
!!$       !print "(*(f7.3))",R_2D(j_theta,:)
!!$       print "(*(f20.15))",R_2D(j_theta,:)
!!$    end do
!!$    print *,"z_2D:"
!!$    do j_theta = 1,N_theta
!!$       !print "(*(f7.3))",Z_2D(j_theta,:)
!!$       print "(*(f20.15))",Z_2D(j_theta,:)
!!$    end do

    call delete_periodic_spline(R0_spline)
    call delete_periodic_spline(z0_spline)
    call delete_periodic_spline(normal_R_spline)
    call delete_periodic_spline(normal_phi_spline)
    call delete_periodic_spline(normal_z_spline)
    call delete_periodic_spline(binormal_R_spline)
    call delete_periodic_spline(binormal_phi_spline)
    call delete_periodic_spline(binormal_z_spline)
    if (order_r_squared) then
       call delete_periodic_spline(tangent_R_spline)
       call delete_periodic_spline(tangent_phi_spline)
       call delete_periodic_spline(tangent_z_spline)
    end if
    deallocate(X_at_this_theta, Y_at_this_theta, Z_at_this_theta)

    ! Fourier transform the result.
    ! This is not a rate-limiting step, so for clarity of code, we don't bother with an FFT.
    mpol = min(N_theta          / 2, mpol1d)
    ntor = min(N_phi_conversion / 2, ntord)
    mpol_nonzero = mpol
    RBC = 0
    RBS = 0
    ZBC = 0
    ZBS = 0
    factor = (2.0d+0) / (N_theta * N_phi_conversion)
    do j_phi = 1, N_phi_conversion
       do j_theta = 1, N_theta
          do m = 0, mpol
             nmin = -ntor
             if (m==0) nmin = 1
             do n = nmin, ntor
                angle = m * theta(j_theta) - n * nfp * phi_conversion(j_phi)
                sinangle = sin(angle)
                cosangle = cos(angle)
                factor2 = factor
                ! The next 2 lines ensure inverse Fourier transform(Fourier transform) = identity
                if (mod(N_theta,         2) == 0 .and.     m  == (N_theta/2))          factor2 = factor2 / 2
                if (mod(N_phi_conversion,2) == 0 .and. abs(n) == (N_phi_conversion/2)) factor2 = factor2 / 2
                RBC(n,m) = RBC(n,m) + R_2D(j_theta, j_phi) * cosangle * factor2
                RBS(n,m) = RBS(n,m) + R_2D(j_theta, j_phi) * sinangle * factor2
                ZBC(n,m) = ZBC(n,m) + Z_2D(j_theta, j_phi) * cosangle * factor2
                ZBS(n,m) = ZBS(n,m) + Z_2D(j_theta, j_phi) * sinangle * factor2
             end do
          end do
       end do
    end do
    RBC(0,0) = sum(R_2D) / (N_theta * N_phi_conversion)
    ZBC(0,0) = sum(Z_2D) / (N_theta * N_phi_conversion)
    
    if (.not. lasym) then
       RBS = 0
       ZBC = 0
    end if

    deallocate(theta,phi_conversion,R_2D,z_2D)

  contains

    function Frenet_to_cylindrical_residual(phi0)
      ! Given a point on the axis with toroidal angle phi0, compute phi for the associated point at r>0,
      ! and find the difference between this phi and the target value of phi.

      implicit none

      real(dp) :: Frenet_to_cylindrical_residual, phi0
      real(dp) :: total_x, total_y, R0_at_phi0, X_at_phi0, Y_at_phi0, Z_at_phi0
      real(dp) :: cosphi0, sinphi0
      real(dp) :: normal_R, normal_phi, normal_x, normal_y
      real(dp) :: binormal_R, binormal_phi, binormal_x, binormal_y
      real(dp) :: tangent_R, tangent_phi, tangent_x, tangent_y

      sinphi0 = sin(phi0)
      cosphi0 = cos(phi0)
      R0_at_phi0   = periodic_splint(phi0,R0_spline)
      X_at_phi0   = periodic_splint(phi0,X_spline)
      Y_at_phi0   = periodic_splint(phi0,Y_spline)
      normal_R     = periodic_splint(phi0,normal_R_spline)
      normal_phi   = periodic_splint(phi0,normal_phi_spline)
      binormal_R   = periodic_splint(phi0,binormal_R_spline)
      binormal_phi = periodic_splint(phi0,binormal_phi_spline)

      normal_x   =   normal_R * cosphi0 -   normal_phi * sinphi0
      normal_y   =   normal_R * sinphi0 +   normal_phi * cosphi0
      binormal_x = binormal_R * cosphi0 - binormal_phi * sinphi0
      binormal_y = binormal_R * sinphi0 + binormal_phi * cosphi0

      total_x = R0_at_phi0 * cosphi0 + X_at_phi0 * normal_x + Y_at_phi0 * binormal_x
      total_y = R0_at_phi0 * sinphi0 + X_at_phi0 * normal_y + Y_at_phi0 * binormal_y

      if (order_r_squared) then
         Z_at_phi0    = periodic_splint(phi0,Z_spline)
         tangent_R    = periodic_splint(phi0,tangent_R_spline)
         tangent_phi  = periodic_splint(phi0,tangent_phi_spline)

         tangent_x = tangent_R * cosphi0 - tangent_phi * sinphi0
         tangent_y = tangent_R * sinphi0 + tangent_phi * cosphi0

         total_x = total_x + Z_at_phi0 * tangent_x
         total_y = total_y + Z_at_phi0 * tangent_y
      end if

      Frenet_to_cylindrical_residual = atan2(total_y, total_x) - phi_target
      ! We expect the residual to be less than pi in absolute value, so if it is not, the reason must be the branch cut:
      if (Frenet_to_cylindrical_residual >  pi) Frenet_to_cylindrical_residual = Frenet_to_cylindrical_residual - 2*pi
      if (Frenet_to_cylindrical_residual < -pi) Frenet_to_cylindrical_residual = Frenet_to_cylindrical_residual + 2*pi

    end function Frenet_to_cylindrical_residual

    ! -----------------------------------------------

    subroutine Frenet_to_cylindrical_1_point(phi0, total_R, total_z)
      ! Given a point on the axis with toroidal angle phi0, compute R and z for the associated point at r>0.

      implicit none

      real(dp), intent(in) :: phi0
      real(dp), intent(out) :: total_R, total_z
      real(dp) :: total_x, total_y, R0_at_phi0, z0_at_phi0, X_at_phi0, Y_at_phi0, Z_at_phi0
      real(dp) :: cosphi0, sinphi0
      real(dp) :: normal_R, normal_phi, normal_x, normal_y, normal_z
      real(dp) :: binormal_R, binormal_phi, binormal_x, binormal_y, binormal_z
      real(dp) :: tangent_R, tangent_phi, tangent_x, tangent_y, tangent_z
      
      sinphi0 = sin(phi0)
      cosphi0 = cos(phi0)
      R0_at_phi0   = periodic_splint(phi0,R0_spline)
      z0_at_phi0   = periodic_splint(phi0,z0_spline)
      X_at_phi0   = periodic_splint(phi0,X_spline)
      Y_at_phi0   = periodic_splint(phi0,Y_spline)
      normal_R     = periodic_splint(phi0,normal_R_spline)
      normal_phi   = periodic_splint(phi0,normal_phi_spline)
      normal_z     = periodic_splint(phi0,normal_z_spline)
      binormal_R   = periodic_splint(phi0,binormal_R_spline)
      binormal_phi = periodic_splint(phi0,binormal_phi_spline)
      binormal_z   = periodic_splint(phi0,binormal_z_spline)

      normal_x   =   normal_R * cosphi0 -   normal_phi * sinphi0
      normal_y   =   normal_R * sinphi0 +   normal_phi * cosphi0
      binormal_x = binormal_R * cosphi0 - binormal_phi * sinphi0
      binormal_y = binormal_R * sinphi0 + binormal_phi * cosphi0

      total_x = R0_at_phi0 * cosphi0 + X_at_phi0 * normal_x + Y_at_phi0 * binormal_x
      total_y = R0_at_phi0 * sinphi0 + X_at_phi0 * normal_y + Y_at_phi0 * binormal_y
      total_z = z0_at_phi0           + X_at_phi0 * normal_z + Y_at_phi0 * binormal_z

      if (order_r_squared) then
         Z_at_phi0    = periodic_splint(phi0,Z_spline)
         tangent_R    = periodic_splint(phi0,tangent_R_spline)
         tangent_phi  = periodic_splint(phi0,tangent_phi_spline)
         tangent_z    = periodic_splint(phi0,tangent_z_spline)

         tangent_x = tangent_R * cosphi0 - tangent_phi * sinphi0
         tangent_y = tangent_R * sinphi0 + tangent_phi * cosphi0

         total_x = total_x + Z_at_phi0 * tangent_x
         total_y = total_y + Z_at_phi0 * tangent_y
         total_z = total_z + Z_at_phi0 * tangent_z
      end if

      total_R = sqrt(total_x * total_x + total_y * total_y)

    end subroutine Frenet_to_cylindrical_1_point

  end subroutine quasisymmetry_Frenet_to_cylindrical_nonlinear

end module quasisymmetry_Frenet_to_cylindrical_mod
