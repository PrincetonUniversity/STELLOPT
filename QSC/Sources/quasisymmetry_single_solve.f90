subroutine quasisymmetry_single_solve

  use quasisymmetry_variables
  use vmec_input, only: ntor   !1/23/19,6/4/19.(9u6)

  implicit none

  real(dp) :: x
  integer :: iteration, new_N_phi
  real(dp), dimension(:), allocatable :: angle, sinangle, cosangle
  real :: solve_start_time, solve_end_time

  if (verbose) call cpu_time(solve_start_time)

  iota_tolerance_achieved = .false.
  elongation_tolerance_achieved = .false.
  iteration = 0
  N_phi = N_phi_original
  already_found_max_curvature = .false.
  skipped_solve = .false.

  if (sum(R0c) <= 1.0e-10) then
     if (verbose) print *,"R0 will be <= 0 at phi=0, so skipping solve."
     skipped_solve = .true.
     return
  end if

  if (sum(R0c(1::2)) - sum(R0c(2::2)) <= 1.0e-10) then
     if (verbose) print *,"R0 will be <= 0 at phi=pi/nfp, so skipping solve."
     skipped_solve = .true.
     return
  end if

  do 
     iteration = iteration + 1
     if (verbose) then
        print "(a)"," -------------------------------------------------------"
        print "(a,i4)"," Solving system using N_phi=",N_phi
     end if
     
     call quasisymmetry_init_phi()

     call quasisymmetry_init_axis()

     if (any(R0 < 1.0d-10)) then
        if (verbose) print *,"R0 is <= 0, so skipping solve."
        skipped_solve = .true.
        exit
     end if

     if (max_curvature > max_max_curvature_to_keep) then
        if (verbose) print *,"max_curvature > max_max_curvature_to_keep, so skipping solve."
        skipped_solve = .true.
        exit
     end if

     call quasisymmetry_solve()

     if (trim(resolution_option) == resolution_option_fixed) exit

     if (iteration > 1) then
        if (verbose) print "(a,es10.3)","                                      abs(iota - last_iota) =",abs(iota - last_iota)
        if (verbose) print "(a,es10.3)"," abs(max_elongation - last_max_elongation) / max_elongation =",abs(max_elongation - last_max_elongation) / max_elongation
        if (abs(iota - last_iota) <= iota_tolerance) then
           if (verbose) print *,"iota_tolerance (absolute) achieved."
           iota_tolerance_achieved = .true.
        end if
        !if (abs(max_elongation - last_max_elongation) <= elongation_tolerance) then  ! Absolute tolerance
        if (abs(max_elongation - last_max_elongation) / max_elongation <= elongation_tolerance) then  ! Relative tolerance
           if (verbose) print *,"elongation_tolerance (relative) achieved."
           elongation_tolerance_achieved = .true.
        end if
        if (max_elongation > max_precise_elongation) then
           if (verbose) print "(a)"," Ignoring elongation_tolerance and iota_tolerance since max_elongation > max_precise_elongation."
           elongation_tolerance_achieved = .true.
           iota_tolerance_achieved = .true.
        end if
        if (iota_tolerance_achieved .and. elongation_tolerance_achieved) exit
     end if

     last_iota = iota
     last_max_elongation = max_elongation

     new_N_phi = N_phi * 2 + 1
     write(0,*)'qs_single. new_N_phi,N_phi,ntor=',new_N_phi,N_phi,ntor  !6/4/19.(9u6)as (8k11g)
     if (new_N_phi > max_N_phi) then
        if (verbose) print *,"Stopping N_phi refinement since max_N_phi exceeded."
        exit
     end if
     N_phi = new_N_phi

  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now evaluate diagnostics that need only be evaluated at the final N_phi resolution
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (allocated(X1s_untwisted)) deallocate(X1s_untwisted)
  if (allocated(X1c_untwisted)) deallocate(X1c_untwisted)
  if (allocated(Y1s_untwisted)) deallocate(Y1s_untwisted)
  if (allocated(Y1c_untwisted)) deallocate(Y1c_untwisted)
  if (allocated(R1s)) deallocate(R1s)
  if (allocated(R1c)) deallocate(R1c)
  if (allocated(Z1s)) deallocate(Z1s)
  if (allocated(Z1c)) deallocate(Z1c)

  allocate(X1s_untwisted(N_phi))
  allocate(X1c_untwisted(N_phi))
  allocate(Y1s_untwisted(N_phi))
  allocate(Y1c_untwisted(N_phi))
  allocate(R1s(N_phi))
  allocate(R1c(N_phi))
  allocate(Z1s(N_phi))
  allocate(Z1c(N_phi))

  ! Derivatives of Y1c, Y1s, and X1c (without untwisting) are needed for computing the grad B tensor, and possibly for the higher-order-in-r terms.
  if (allocated(d_Y1c_d_zeta)) deallocate(d_Y1c_d_zeta)
  if (allocated(d_Y1s_d_zeta)) deallocate(d_Y1s_d_zeta)
  if (allocated(d_X1c_d_zeta)) deallocate(d_X1c_d_zeta)
  allocate(d_Y1c_d_zeta(N_phi))
  allocate(d_Y1s_d_zeta(N_phi))
  allocate(d_X1c_d_zeta(N_phi))
  d_Y1c_d_zeta = matmul(d_d_zeta,Y1c)
  d_Y1s_d_zeta = matmul(d_d_zeta,Y1s)
  d_X1c_d_zeta = matmul(d_d_zeta,X1c)

  if (trim(order_r_option) .ne. order_r_option_r1) call quasisymmetry_higher_order_in_r()

  if (verbose) then
     call cpu_time(solve_end_time)
     print *,"Time to solve equations, without diagnostics:",solve_end_time-solve_start_time
  end if

  iota_from_torsion = sum(torsion * d_l_d_phi) * d_phi * nfp / (2 * pi) + axis_helicity * nfp
  if (verbose) print *,"Iota expected from integrated torsion:",iota_from_torsion

  if (circular_cross_section_surface) then
     X1c = 1
     Y1s = 1
     Y1c = 0
  end if

  ! If helicity is nonzero, then the original X1s/X1c/Y1s/Y1c variables are defined with respect to a "poloidal" angle that
  ! is actually helical, with the theta=0 curve wrapping around the magnetic axis as you follow phi around toroidally. Therefore
  ! here we convert to an untwisted poloidal angle, such that the theta=0 curve does not wrap around the axis.
  if (axis_helicity == 0 .or. (.not. untwist)) then
     X1s_untwisted = 0
     X1c_untwisted = X1c
     Y1s_untwisted = Y1s
     Y1c_untwisted = Y1c
     if (trim(order_r_option) .ne. order_r_option_r1) then
        ! We have O(r^2) terms
        X20_untwisted = X20
        X2s_untwisted = X2s
        X2c_untwisted = X2c
        Y20_untwisted = Y20
        Y2s_untwisted = Y2s
        Y2c_untwisted = Y2c
        Z20_untwisted = Z20
        Z2s_untwisted = Z2s
        Z2c_untwisted = Z2c
     end if
     if (trim(order_r_option).ne.order_r_option_r1 .and. trim(order_r_option).ne.order_r_option_r2) then
        ! We have O(r^3) terms
        X3s1_untwisted = X3s1
        X3s3_untwisted = X3s3
        X3c1_untwisted = X3c1
        X3c3_untwisted = X3c3
        Y3s1_untwisted = Y3s1
        Y3s3_untwisted = Y3s3
        Y3c1_untwisted = Y3c1
        Y3c3_untwisted = Y3c3
        Z3s1_untwisted = Z3s1
        Z3s3_untwisted = Z3s3
        Z3c1_untwisted = Z3c1
        Z3c3_untwisted = Z3c3
     end if
  else
     allocate(angle(N_phi))
     allocate(sinangle(N_phi))
     allocate(cosangle(N_phi))
     angle = -axis_helicity * nfp * Boozer_toroidal_angle
     sinangle = sin(angle)
     cosangle = cos(angle)
     X1s_untwisted = X1s *   cosangle  + X1c * sinangle
     X1c_untwisted = X1s * (-sinangle) + X1c * cosangle
     Y1s_untwisted = Y1s *   cosangle  + Y1c * sinangle
     Y1c_untwisted = Y1s * (-sinangle) + Y1c * cosangle
     if (trim(order_r_option).ne.order_r_option_r1 .and. trim(order_r_option).ne.order_r_option_r2) then
        ! Then we have O(r^3) terms
        X3s1_untwisted = X3s1 *   cosangle  + X3c1 * sinangle
        X3c1_untwisted = X3s1 * (-sinangle) + X3c1 * cosangle
        Y3s1_untwisted = Y3s1 *   cosangle  + Y3c1 * sinangle
        Y3c1_untwisted = Y3s1 * (-sinangle) + Y3c1 * cosangle
        Z3s1_untwisted = Z3s1 *   cosangle  + Z3c1 * sinangle
        Z3c1_untwisted = Z3s1 * (-sinangle) + Z3c1 * cosangle
        sinangle = sin(3*angle)
        cosangle = cos(3*angle)
        X3s3_untwisted = X3s3 *   cosangle  + X3c3 * sinangle
        X3c3_untwisted = X3s3 * (-sinangle) + X3c3 * cosangle
        Y3s3_untwisted = Y3s3 *   cosangle  + Y3c3 * sinangle
        Y3c3_untwisted = Y3s3 * (-sinangle) + Y3c3 * cosangle
        Z3s3_untwisted = Z3s3 *   cosangle  + Z3c3 * sinangle
        Z3c3_untwisted = Z3s3 * (-sinangle) + Z3c3 * cosangle
     end if
     if (trim(order_r_option) .ne. order_r_option_r1) then
        ! Then we have O(r^2) terms
        X20_untwisted = X20
        Y20_untwisted = Y20
        Z20_untwisted = Z20
        sinangle = sin(2*angle)
        cosangle = cos(2*angle)
        X2s_untwisted = X2s *   cosangle  + X2c * sinangle
        X2c_untwisted = X2s * (-sinangle) + X2c * cosangle
        Y2s_untwisted = Y2s *   cosangle  + Y2c * sinangle
        Y2c_untwisted = Y2s * (-sinangle) + Y2c * cosangle
        Z2s_untwisted = Z2s *   cosangle  + Z2c * sinangle
        Z2c_untwisted = Z2s * (-sinangle) + Z2c * cosangle
     end if
     deallocate(sinangle,cosangle,angle)
  end if

  call quasisymmetry_Frenet_to_cylindrical_linear()
    
  call quasisymmetry_elongation_in_Rz_plane()

  call quasisymmetry_determine_B_helicity()

  call quasisymmetry_grad_B_tensor()

  if (verbose) then
     call cpu_time(solve_end_time)
     print *,"Time to solve equations and compute diagnostics:",solve_end_time-solve_start_time
  end if

end subroutine quasisymmetry_single_solve
