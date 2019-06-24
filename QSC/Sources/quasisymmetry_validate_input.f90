subroutine quasisymmetry_validate_input

  use quasisymmetry_variables

  implicit none

  integer :: j

  if (nfp<1) stop "Error! nfp must be positive."

  if (sign_G .ne. 1 .and. sign_G .ne. -1) stop "sign_G must be +1 or -1."

  if (N_phi < 0) stop "Error! N_phi must be positive"
  
  ! Ensure N_phi is always odd.
  if (mod(N_phi,2) ==0) N_phi = N_phi + 1

!  N_phi_original = N_phi   !6/6/19(9u6)as (7l23a), xferred to qs_rd_input.

  select case (trim(resolution_option))
  case (resolution_option_fixed)
  case (resolution_option_adaptive)
  case default
     print *,"Error! Invalid resolution_option:",resolution_option
     stop
  end select

  select case (trim(general_option))
  case (general_option_single)
  case (general_option_scan)
  case default
     print *,"Error! Invalid general_option:",general_option
     stop
  end select

  select case (trim(verbose_option))
  case (verbose_option_all)
  case (verbose_option_proc0)
  case (verbose_option_summary)
  case default
     print *,"Error! Invalid verbose_option:",verbose_option
     stop
  end select

  select case (trim(eta_bar_scan_option))
  case (eta_bar_scan_option_linear)
  case (eta_bar_scan_option_log)
  case (eta_bar_scan_option_2_sided_log)
  case default
     print *,"Error! Invalid eta_bar_scan_option:",eta_bar_scan_option
     stop
  end select

  select case (trim(sigma_initial_scan_option))
  case (sigma_initial_scan_option_linear)
  case (sigma_initial_scan_option_log)
  case (sigma_initial_scan_option_2_sided_log)
  case default
     print *,"Error! Invalid sigma_initial_scan_option:",sigma_initial_scan_option
     stop
  end select

  select case (trim(Fourier_scan_option))
  case (Fourier_scan_option_linear)
  case (Fourier_scan_option_2_sided_log)
  case (Fourier_scan_option_2_sided_log_except_Z0s1)
  case default
     print *,"Error! Invalid Fourier_scan_option:",Fourier_scan_option
     stop
  end select

  select case (trim(finite_r_option))
  case (finite_r_option_linear)
  case (finite_r_option_nonlinear)
  case default
     print *,"Error! Invalid finite_r_option:",finite_r_option
     stop
  end select

  select case (trim(order_r_option))
  case (order_r_option_r1)
  case (order_r_option_r2)
  case (order_r_option_r3_simplified)
  case (order_r_option_r3_simplified_with_Z3)
  case (order_r_option_r3_flux_constraint)
  case (order_r_option_r3_flux_constraint_const_B20)
  case (order_r_option_r3_B3)
  case (order_r_option_r3_X3s3_X3c3)
  !case (order_r_option_r3_X3s3_Y3s3)
  !case (order_r_option_r3_X3c3_Y3c3)
  !case (order_r_option_r3_Y3s3_Y3c3)
  case default
     print *,"Error! Invalid order_r_option:",order_r_option
     stop
  end select

!!$  if (order_r_squared .and. trim(finite_r_option)==finite_r_option_linear) then
!!$     print "(a)"," NOTE: Since order_r_squared==.true., finite_r_option is being set to 'nonlinear'."
!!$     finite_r_option = finite_r_option_nonlinear
!!$  end if

end subroutine quasisymmetry_validate_input
