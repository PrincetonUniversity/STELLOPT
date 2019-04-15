subroutine quasisymmetry_write_output

  use quasisymmetry_variables
  use ezcdf
  use vmec_input, only: lasym, ntor, RBC, RBS, ZBC, ZBS

  implicit none

  real :: time1, end_time
  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.

  ! Scalars:
  character(len=*), parameter :: &
       vn_general_option = "general_option", &
       vn_nfp = "nfp", &
       vn_sign_G = "sign_G", &
       vn_sign_psi = "sign_psi", &
       vn_resolution_option = "resolution_option", &
       vn_max_precise_elongation = "max_precise_elongation", &
       vn_sigma_initial_min = "sigma_initial_min", &
       vn_sigma_initial_max = "sigma_initial_max", &
       vn_sigma_initial_N_scan = "sigma_initial_N_scan", &
       vn_sigma_initial_scan_option = "sigma_initial_scan_option", &
       vn_eta_bar_min = "eta_bar_min", &
       vn_eta_bar_max = "eta_bar_max", &
       vn_eta_bar_N_scan = "eta_bar_N_scan", &
       vn_eta_bar_scan_option = "eta_bar_scan_option", &
       vn_Fourier_scan_option = "Fourier_scan_option", &
       vn_max_elongation_to_keep = "max_elongation_to_keep", &
       vn_max_max_curvature_to_keep = "max_max_curvature_to_keep", &
       vn_max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep = "max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep", &
       vn_min_iota_to_keep = "min_iota_to_keep", &
       vn_eta_bar = "eta_bar", &
       vn_sigma_initial = "sigma_initial", &
       vn_I2_over_B0 = "I2_over_B0", &
       vn_iota = "iota", &
       vn_max_elongation = "max_elongation", &
       vn_rms_curvature = "rms_curvature", &
       vn_max_curvature = "max_curvature", &
       vn_max_modBinv_sqrt_half_grad_B_colon_grad_B = "max_modBinv_sqrt_half_grad_B_colon_grad_B", &
       vn_axis_length = "axis_length", &
       vn_standard_deviation_of_R = "standard_deviation_of_R", &
       vn_standard_deviation_of_Z = "standard_deviation_of_Z", &
       vn_axis_helicity = "axis_helicity", &
       vn_B_helicity = "B_helicity", &
       vn_effective_nfp = "effective_nfp", &
       vn_Newton_tolerance_achieved = "Newton_tolerance_achieved", &
       vn_iota_tolerance_achieved = "iota_tolerance_achieved", &
       vn_elongation_tolerance_achieved = "elongation_tolerance_achieved", &
       vn_finite_r_option = "finite_r_option", &
       vn_N_phi = "N_phi", &
       vn_r = "r", &
       vn_B0 = "B0", &
       vn_mpol = "mpol", &
       vn_ntor = "ntor", &
       vn_lasym = "lasym", &
       vn_untwist = "untwist", &
       vn_order_r_squared = "order_r_squared", &
       vn_B2s = "B2s", &
       vn_B2c = "B2c", &
       vn_p2 = "p2", &
       vn_B20_mean = "B20_mean", &
       vn_B20_residual = "B20_residual"

  ! Arrays with dimension 1
  character(len=*), parameter :: &
       vn_iotas = "iotas", &
       vn_max_elongations = "max_elongations", &
       vn_rms_curvatures = "rms_curvatures", &
       vn_max_curvatures = "max_curvatures", &
       vn_max_modBinv_sqrt_half_grad_B_colon_grad_Bs = "max_modBinv_sqrt_half_grad_B_colon_grad_Bs", &
       vn_axis_lengths = "axis_lengths", &
       vn_standard_deviations_of_R = "standard_deviations_of_R", &
       vn_standard_deviations_of_Z = "standard_deviations_of_Z", &
       vn_axis_helicities = "axis_helicities", &
       vn_B_helicities = "B_helicities", &
       vn_effective_nfps = "effective_nfps", &
       vn_Newton_tolerance_achieveds = "Newton_tolerance_achieveds", &
       vn_iota_tolerance_achieveds = "iota_tolerance_achieveds", &
       vn_elongation_tolerance_achieveds = "elongation_tolerance_achieveds", &
       vn_R0s_min = "R0s_min", &
       vn_R0s_max = "R0s_max", &
       vn_R0s_N_scan = "R0s_N_scan", &
       vn_R0c_min = "R0c_min", &
       vn_R0c_max = "R0c_max", &
       vn_R0c_N_scan = "R0c_N_scan", &
       vn_Z0s_min = "Z0s_min", &
       vn_Z0s_max = "Z0s_max", &
       vn_Z0s_N_scan = "Z0s_N_scan", &
       vn_Z0c_min = "Z0c_min", &
       vn_Z0c_max = "Z0c_max", &
       vn_Z0c_N_scan = "Z0c_N_scan", &
       vn_scan_sigma_initial = "scan_sigma_initial", &
       vn_sigma_initial_values = "sigma_initial_values", &
       vn_scan_eta_bar = "scan_eta_bar", &
       vn_eta_bar_values = "eta_bar_values", &
       vn_R0c = "R0c", &
       vn_R0s = "R0s", &
       vn_Z0c = "Z0c", &
       vn_Z0s = "Z0s", &
       vn_phi = "phi", &
       vn_Boozer_toroidal_angle = "Boozer_toriodal_angle", &
       vn_R0 = "R0", &
       vn_z0 = "z0", &
       vn_curvature = "curvature", &
       vn_torsion = "torsion", &
       vn_sigma = "sigma", &
       vn_X1c = "X1c", &
       vn_Y1c = "Y1c", &
       vn_Y1s = "Y1s", &
       vn_X1c_untwisted = "X1c_untwisted", &
       vn_X1s_untwisted = "X1s_untwisted", &
       vn_Y1c_untwisted = "Y1c_untwisted", &
       vn_Y1s_untwisted = "Y1s_untwisted", &
       vn_R1c = "R1c", &
       vn_R1s = "R1s", &
       vn_z1c = "z1c", &
       vn_z1s = "z1s", &
       vn_elongation = "elongation", &
       vn_elongation_in_Rz_plane = "elongation_in_Rz_plane", &
       vn_d_l_d_phi = "d_l_d_phi", &
       vn_modBinv_sqrt_half_grad_B_colon_grad_B = "modBinv_sqrt_half_grad_B_colon_grad_B", &
       vn_X20 = "X20", &
       vn_X2s = "X2s", &
       vn_X2c = "X2c", &
       vn_Y20 = "Y20", &
       vn_Y2s = "Y2s", &
       vn_Y2c = "Y2c", &
       vn_Z20 = "Z20", &
       vn_Z2s = "Z2s", &
       vn_Z2c = "Z2c", &
       vn_X20_untwisted = "X20_untwisted", &
       vn_X2s_untwisted = "X2s_untwisted", &
       vn_X2c_untwisted = "X2c_untwisted", &
       vn_Y20_untwisted = "Y20_untwisted", &
       vn_Y2s_untwisted = "Y2s_untwisted", &
       vn_Y2c_untwisted = "Y2c_untwisted", &
       vn_Z20_untwisted = "Z20_untwisted", &
       vn_Z2s_untwisted = "Z2s_untwisted", &
       vn_Z2c_untwisted = "Z2c_untwisted", &
       vn_B20 = "B20", &
       vn_R20 = "R20", &
       vn_R2s = "R2s", &
       vn_R2c = "R2c", &
       vn_z20_cylindrical = "z20_cylindrical", &
       vn_z2s_cylindrical = "z2s_cylindrical", &
       vn_z2c_cylindrical = "z2c_cylindrical"

  ! Arrays with dimension 2
  character(len=*), parameter :: &
       vn_N_scan_array  = "N_scan_array", &
       vn_scan_R0c  = "scan_R0c", &
       vn_scan_R0s  = "scan_R0s", &
       vn_scan_Z0c  = "scan_Z0c", &
       vn_scan_Z0s  = "scan_Z0s", &
       vn_RBC = "RBC", &
       vn_RBS = "RBS", &
       vn_ZBC = "ZBC", &
       vn_ZBS = "ZBS"

!!$  ! Arrays with dimension 3
!!$  character(len=*), parameter :: &
!!$       vn_r_plasma  = "r_plasma"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now create variables that name the dimensions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       N_scan_dim = (/'N_scan'/), &
       axis_nmax_plus_1_dim = (/'axis_nmax_plus_1'/), &
       sigma_initial_N_scan_dim = (/'sigma_initial_N_scan'/), &
       eta_bar_N_scan_dim = (/'eta_bar_N_scan'/), &
       N_phi_dim = (/'N_phi'/)

  ! Arrays with dimension 2:
  ! The form of the array declarations here is inspired by
  ! http://stackoverflow.com/questions/21552430/gfortran-does-not-allow-character-arrays-with-varying-component-lengths
  character(len=*), parameter, dimension(2) :: &
       axis_nmax_plus_1_4_dim = (/ character(len=50) :: 'axis_nmax_plus_1','4'/), &
       N_scan_axis_nmax_plus_1_dim = (/ character(len=50) :: 'N_scan','axis_nmax_plus_1'/), &
       ntor_mpol_dim = (/ character(len=50) :: 'ntor_times_2_plus_1','mpol_plus_1' /)

!!$  ! Arrays with dimension 3:
!!$  character(len=*), parameter, dimension(3) :: &
!!$       xyz_ntheta_nzetal_plasma_dim = (/ character(len=50) :: 'xyz','ntheta_plasma','nzetal_plasma'/), &

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Only proc 0 writes.
  if (.not. proc0) return

  call cpu_time(time1)

  call cdf_open(ncid,output_filename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",output_filename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_general_option, general_option)

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_setatt(ncid, vn_nfp, 'Number of field periods, i.e. the number of identical toroidal segments, 5 for W7-X, 4 for HSX, etc. ' // &
       'Equivalent to the VMEC variable of the same name.')

  call cdf_define(ncid, vn_sign_G, sign_G)
  call cdf_define(ncid, vn_sign_psi, sign_psi)

  call cdf_define(ncid, vn_resolution_option, resolution_option)
  !call cdf_setatt(ncid, vn_resolution_option, 'Method used to define the geometry of the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_max_precise_elongation, max_precise_elongation)
  call cdf_define(ncid, vn_order_r_squared, order_r_squared)

  select case (trim(general_option))
  case (general_option_single)
     call cdf_define(ncid, vn_eta_bar, eta_bar)
     call cdf_define(ncid, vn_sigma_initial, sigma_initial)
     call cdf_define(ncid, vn_I2_over_B0, I2_over_B0)
     call cdf_define(ncid, vn_iota, iota)
     call cdf_define(ncid, vn_max_elongation, max_elongation)
     call cdf_define(ncid, vn_rms_curvature, rms_curvature)
     call cdf_define(ncid, vn_max_curvature, max_curvature)
     call cdf_define(ncid, vn_max_modBinv_sqrt_half_grad_B_colon_grad_B, max_modBinv_sqrt_half_grad_B_colon_grad_B)
     call cdf_define(ncid, vn_axis_length, axis_length)
     call cdf_define(ncid, vn_standard_deviation_of_R, standard_deviation_of_R)
     call cdf_define(ncid, vn_standard_deviation_of_Z, standard_deviation_of_Z)
     call cdf_define(ncid, vn_axis_helicity, axis_helicity)
     call cdf_define(ncid, vn_B_helicity, B_helicity)
     call cdf_define(ncid, vn_effective_nfp, effective_nfp)
     call cdf_define(ncid, vn_Newton_tolerance_achieved, Newton_tolerance_achieved)
     call cdf_define(ncid, vn_iota_tolerance_achieved, iota_tolerance_achieved)
     call cdf_define(ncid, vn_elongation_tolerance_achieved, elongation_tolerance_achieved)
     call cdf_define(ncid, vn_finite_r_option, finite_r_option)
     call cdf_define(ncid, vn_N_phi, N_phi)
     call cdf_define(ncid, vn_r, r)
     call cdf_define(ncid, vn_B0, B0)
     call cdf_define(ncid, vn_mpol, mpol_nonzero)
     call cdf_define(ncid, vn_ntor, ntor)
     call cdf_define(ncid, vn_lasym, lasym)
     call cdf_define(ncid, vn_untwist, untwist)
     if (order_r_squared) then
        call cdf_define(ncid, vn_B2s, B2s)
        call cdf_define(ncid, vn_B2c, B2c)
        call cdf_define(ncid, vn_p2, p2)
        call cdf_define(ncid, vn_B20_mean, B20_mean)
        call cdf_define(ncid, vn_B20_residual, B20_residual)
     end if
  case (general_option_scan)
     call cdf_define(ncid, vn_sigma_initial_min, sigma_initial_min)
     call cdf_define(ncid, vn_sigma_initial_max, sigma_initial_max)
     call cdf_define(ncid, vn_sigma_initial_N_scan, sigma_initial_N_scan)
     call cdf_define(ncid, vn_sigma_initial_scan_option, sigma_initial_scan_option)
     call cdf_define(ncid, vn_eta_bar_min, eta_bar_min)
     call cdf_define(ncid, vn_eta_bar_max, eta_bar_max)
     call cdf_define(ncid, vn_eta_bar_N_scan, eta_bar_N_scan)
     call cdf_define(ncid, vn_eta_bar_scan_option, eta_bar_scan_option)
     call cdf_define(ncid, vn_Fourier_scan_option, Fourier_scan_option)
     call cdf_define(ncid, vn_max_elongation_to_keep, max_elongation_to_keep)
     call cdf_define(ncid, vn_max_max_curvature_to_keep, max_max_curvature_to_keep)
     call cdf_define(ncid, vn_max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep, max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep)
     call cdf_define(ncid, vn_min_iota_to_keep, min_iota_to_keep)
  end select

  ! Arrays with dimension 1

  select case (trim(general_option))
  case (general_option_single)
     call cdf_define(ncid, vn_R0c, R0c(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0s, R0s(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c, Z0c(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s, Z0s(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_phi, phi, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Boozer_toroidal_angle, Boozer_toroidal_angle, dimname=N_phi_dim)
     call cdf_define(ncid, vn_R0, R0, dimname=N_phi_dim)
     call cdf_define(ncid, vn_z0, z0, dimname=N_phi_dim)
     call cdf_define(ncid, vn_curvature, curvature, dimname=N_phi_dim)
     call cdf_define(ncid, vn_torsion, torsion, dimname=N_phi_dim)
     call cdf_define(ncid, vn_sigma, sigma, dimname=N_phi_dim)
     call cdf_define(ncid, vn_X1c, X1c, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1c, Y1c, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1s, Y1s, dimname=N_phi_dim)
     call cdf_define(ncid, vn_X1c_untwisted, X1c_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_X1s_untwisted, X1s_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1c_untwisted, Y1c_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1s_untwisted, Y1s_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_R1c, R1c, dimname=N_phi_dim)
     call cdf_define(ncid, vn_R1s, R1s, dimname=N_phi_dim)
     call cdf_define(ncid, vn_z1c, z1c, dimname=N_phi_dim)
     call cdf_define(ncid, vn_z1s, z1s, dimname=N_phi_dim)
     call cdf_define(ncid, vn_elongation, elongation, dimname=N_phi_dim)
     call cdf_define(ncid, vn_elongation_in_Rz_plane, elongation_in_Rz_plane, dimname=N_phi_dim)
     call cdf_define(ncid, vn_d_l_d_phi, d_l_d_phi, dimname=N_phi_dim)
     call cdf_define(ncid, vn_modBinv_sqrt_half_grad_B_colon_grad_B, modBinv_sqrt_half_grad_B_colon_grad_B, dimname=N_phi_dim)
     if (order_r_squared) then
        call cdf_define(ncid, vn_X20, X20, dimname=N_phi_dim)
        call cdf_define(ncid, vn_X2s, X2s, dimname=N_phi_dim)
        call cdf_define(ncid, vn_X2c, X2c, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Y20, Y20, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Y2s, Y2s, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Y2c, Y2c, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Z20, Z20, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Z2s, Z2s, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Z2c, Z2c, dimname=N_phi_dim)
        call cdf_define(ncid, vn_X20_untwisted, X20_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_X2s_untwisted, X2s_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_X2c_untwisted, X2c_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Y20_untwisted, Y20_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Y2s_untwisted, Y2s_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Y2c_untwisted, Y2c_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Z20_untwisted, Z20_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Z2s_untwisted, Z2s_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_Z2c_untwisted, Z2c_untwisted, dimname=N_phi_dim)
        call cdf_define(ncid, vn_B20, B20, dimname=N_phi_dim)
        call cdf_define(ncid, vn_R20, R20, dimname=N_phi_dim)
        call cdf_define(ncid, vn_R2s, R2s, dimname=N_phi_dim)
        call cdf_define(ncid, vn_R2c, R2c, dimname=N_phi_dim)
        call cdf_define(ncid, vn_z20_cylindrical, z20_cylindrical, dimname=N_phi_dim)
        call cdf_define(ncid, vn_z2s_cylindrical, z2s_cylindrical, dimname=N_phi_dim)
        call cdf_define(ncid, vn_z2c_cylindrical, z2c_cylindrical, dimname=N_phi_dim)
     end if
  case (general_option_scan)
     call cdf_define(ncid, vn_iotas, iotas, dimname=N_scan_dim)
     call cdf_define(ncid, vn_max_elongations, max_elongations, dimname=N_scan_dim)
     call cdf_define(ncid, vn_rms_curvatures, rms_curvatures, dimname=N_scan_dim)
     call cdf_define(ncid, vn_max_curvatures, max_curvatures, dimname=N_scan_dim)
     call cdf_define(ncid, vn_max_modBinv_sqrt_half_grad_B_colon_grad_Bs, max_modBinv_sqrt_half_grad_B_colon_grad_Bs, dimname=N_scan_dim)
     call cdf_define(ncid, vn_axis_lengths, axis_lengths, dimname=N_scan_dim)
     call cdf_define(ncid, vn_standard_deviations_of_R, standard_deviations_of_R, dimname=N_scan_dim)
     call cdf_define(ncid, vn_standard_deviations_of_Z, standard_deviations_of_Z, dimname=N_scan_dim)
     call cdf_define(ncid, vn_axis_helicities, axis_helicities, dimname=N_scan_dim)
     call cdf_define(ncid, vn_B_helicities, B_helicities, dimname=N_scan_dim)
     call cdf_define(ncid, vn_effective_nfps, effective_nfps, dimname=N_scan_dim)
     call cdf_define(ncid, vn_Newton_tolerance_achieveds, Newton_tolerance_achieveds, dimname=N_scan_dim)
     call cdf_define(ncid, vn_iota_tolerance_achieveds, iota_tolerance_achieveds, dimname=N_scan_dim)
     call cdf_define(ncid, vn_elongation_tolerance_achieveds, elongation_tolerance_achieveds, dimname=N_scan_dim)
     call cdf_define(ncid, vn_R0s_min, R0s_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0s_max, R0s_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0s_N_scan, R0s_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0c_min, R0c_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0c_max, R0c_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0c_N_scan, R0c_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s_min, Z0s_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s_max, Z0s_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s_N_scan, Z0s_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c_min, Z0c_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c_max, Z0c_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c_N_scan, Z0c_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_sigma_initial, scan_sigma_initial, dimname=N_scan_dim)
     call cdf_define(ncid, vn_sigma_initial_values, sigma_initial_values, dimname=sigma_initial_N_scan_dim)
     call cdf_define(ncid, vn_scan_eta_bar, scan_eta_bar, dimname=N_scan_dim)
     call cdf_define(ncid, vn_eta_bar_values, eta_bar_values, dimname=eta_bar_N_scan_dim)
  end select

  ! Arrays with dimension 2

  select case (trim(general_option))
  case (general_option_single)
     call cdf_define(ncid, vn_RBC, RBC(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
     call cdf_define(ncid, vn_RBS, RBS(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
     call cdf_define(ncid, vn_ZBC, ZBC(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
     call cdf_define(ncid, vn_ZBS, ZBS(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
  case (general_option_scan)
     call cdf_define(ncid, vn_N_scan_array,  N_scan_array(1:axis_nmax+1,:), dimname=axis_nmax_plus_1_4_dim)
     call cdf_define(ncid, vn_scan_R0c,  scan_R0c, dimname=N_scan_axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_R0s,  scan_R0s, dimname=N_scan_axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_Z0c,  scan_Z0c, dimname=N_scan_axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_Z0s,  scan_Z0s, dimname=N_scan_axis_nmax_plus_1_dim)
  end select

  ! Arrays with dimension 3

  !call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_general_option, general_option)
  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_sign_G, sign_G)
  call cdf_write(ncid, vn_sign_psi, sign_psi)
  call cdf_write(ncid, vn_resolution_option, resolution_option)
  call cdf_write(ncid, vn_max_precise_elongation, max_precise_elongation)
  call cdf_write(ncid, vn_order_r_squared, order_r_squared)

  select case (trim(general_option))
  case (general_option_single)
     call cdf_write(ncid, vn_eta_bar, eta_bar)
     call cdf_write(ncid, vn_sigma_initial, sigma_initial)
     call cdf_write(ncid, vn_I2_over_B0, I2_over_B0)
     call cdf_write(ncid, vn_iota, iota)
     call cdf_write(ncid, vn_max_elongation, max_elongation)
     call cdf_write(ncid, vn_rms_curvature, rms_curvature)
     call cdf_write(ncid, vn_max_curvature, max_curvature)
     call cdf_write(ncid, vn_max_modBinv_sqrt_half_grad_B_colon_grad_B, max_modBinv_sqrt_half_grad_B_colon_grad_B)
     call cdf_write(ncid, vn_axis_length, axis_length)
     call cdf_write(ncid, vn_standard_deviation_of_R, standard_deviation_of_R)
     call cdf_write(ncid, vn_standard_deviation_of_Z, standard_deviation_of_Z)
     call cdf_write(ncid, vn_axis_helicity, axis_helicity)
     call cdf_write(ncid, vn_B_helicity, B_helicity)
     call cdf_write(ncid, vn_effective_nfp, effective_nfp)
     call cdf_write(ncid, vn_Newton_tolerance_achieved, Newton_tolerance_achieved)
     call cdf_write(ncid, vn_iota_tolerance_achieved, iota_tolerance_achieved)
     call cdf_write(ncid, vn_elongation_tolerance_achieved, elongation_tolerance_achieved)
     call cdf_write(ncid, vn_finite_r_option, finite_r_option)
     call cdf_write(ncid, vn_N_phi, N_phi)
     call cdf_write(ncid, vn_r, r)
     call cdf_write(ncid, vn_B0, B0)
     call cdf_write(ncid, vn_mpol, mpol_nonzero)
     call cdf_write(ncid, vn_ntor, ntor)
     call cdf_write(ncid, vn_lasym, lasym)
     call cdf_write(ncid, vn_untwist, untwist)
     if (order_r_squared) then
        call cdf_write(ncid, vn_B2s, B2s)
        call cdf_write(ncid, vn_B2c, B2c)
        call cdf_write(ncid, vn_p2, p2)
        call cdf_write(ncid, vn_B20_mean, B20_mean)
        call cdf_write(ncid, vn_B20_residual, B20_residual)
     end if
  case (general_option_scan)
     call cdf_write(ncid, vn_sigma_initial_min, sigma_initial_min)
     call cdf_write(ncid, vn_sigma_initial_max, sigma_initial_max)
     call cdf_write(ncid, vn_sigma_initial_N_scan, sigma_initial_N_scan)
     call cdf_write(ncid, vn_eta_bar_min, eta_bar_min)
     call cdf_write(ncid, vn_eta_bar_max, eta_bar_max)
     call cdf_write(ncid, vn_eta_bar_N_scan, eta_bar_N_scan)
     call cdf_write(ncid, vn_eta_bar_scan_option, eta_bar_scan_option)
     call cdf_write(ncid, vn_Fourier_scan_option, Fourier_scan_option)
     call cdf_write(ncid, vn_sigma_initial_scan_option, sigma_initial_scan_option)
     call cdf_write(ncid, vn_max_elongation_to_keep, max_elongation_to_keep)
     call cdf_write(ncid, vn_max_max_curvature_to_keep, max_max_curvature_to_keep)
     call cdf_write(ncid, vn_max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep, max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep)
     call cdf_write(ncid, vn_min_iota_to_keep, min_iota_to_keep)
  end select

  ! Arrays with dimension 1

  select case (trim(general_option))
  case (general_option_single)
     call cdf_write(ncid, vn_R0c, R0c(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0s, R0s(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c, Z0c(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s, Z0s(1:axis_nmax+1))
     call cdf_write(ncid, vn_phi, phi)
     call cdf_write(ncid, vn_Boozer_toroidal_angle, Boozer_toroidal_angle)
     call cdf_write(ncid, vn_R0, R0)
     call cdf_write(ncid, vn_z0, z0)
     call cdf_write(ncid, vn_curvature, curvature)
     call cdf_write(ncid, vn_torsion, torsion)
     call cdf_write(ncid, vn_sigma, sigma)
     call cdf_write(ncid, vn_X1c, X1c)
     call cdf_write(ncid, vn_Y1c, Y1c)
     call cdf_write(ncid, vn_Y1s, Y1s)
     call cdf_write(ncid, vn_X1c_untwisted, X1c_untwisted)
     call cdf_write(ncid, vn_X1s_untwisted, X1s_untwisted)
     call cdf_write(ncid, vn_Y1c_untwisted, Y1c_untwisted)
     call cdf_write(ncid, vn_Y1s_untwisted, Y1s_untwisted)
     call cdf_write(ncid, vn_R1c, R1c)
     call cdf_write(ncid, vn_R1s, R1s)
     call cdf_write(ncid, vn_z1c, z1c)
     call cdf_write(ncid, vn_z1s, z1s)
     call cdf_write(ncid, vn_elongation, elongation)
     call cdf_write(ncid, vn_elongation_in_Rz_plane, elongation_in_Rz_plane)
     call cdf_write(ncid, vn_d_l_d_phi, d_l_d_phi)
     call cdf_write(ncid, vn_modBinv_sqrt_half_grad_B_colon_grad_B, modBinv_sqrt_half_grad_B_colon_grad_B)
     if (order_r_squared) then
        call cdf_write(ncid, vn_X20, X20)
        call cdf_write(ncid, vn_X2s, X2s)
        call cdf_write(ncid, vn_X2c, X2c)
        call cdf_write(ncid, vn_Y20, Y20)
        call cdf_write(ncid, vn_Y2s, Y2s)
        call cdf_write(ncid, vn_Y2c, Y2c)
        call cdf_write(ncid, vn_Z20, Z20)
        call cdf_write(ncid, vn_Z2s, Z2s)
        call cdf_write(ncid, vn_Z2c, Z2c)
        call cdf_write(ncid, vn_X20_untwisted, X20_untwisted)
        call cdf_write(ncid, vn_X2s_untwisted, X2s_untwisted)
        call cdf_write(ncid, vn_X2c_untwisted, X2c_untwisted)
        call cdf_write(ncid, vn_Y20_untwisted, Y20_untwisted)
        call cdf_write(ncid, vn_Y2s_untwisted, Y2s_untwisted)
        call cdf_write(ncid, vn_Y2c_untwisted, Y2c_untwisted)
        call cdf_write(ncid, vn_Z20_untwisted, Z20_untwisted)
        call cdf_write(ncid, vn_Z2s_untwisted, Z2s_untwisted)
        call cdf_write(ncid, vn_Z2c_untwisted, Z2c_untwisted)
        call cdf_write(ncid, vn_B20, B20)
        call cdf_write(ncid, vn_R20, R20)
        call cdf_write(ncid, vn_R2s, R2s)
        call cdf_write(ncid, vn_R2c, R2c)
        call cdf_write(ncid, vn_z20_cylindrical, z20_cylindrical)
        call cdf_write(ncid, vn_z2s_cylindrical, z2s_cylindrical)
        call cdf_write(ncid, vn_z2c_cylindrical, z2c_cylindrical)
     end if
  case (general_option_scan)
     call cdf_write(ncid, vn_iotas, iotas)
     call cdf_write(ncid, vn_max_elongations, max_elongations)
     call cdf_write(ncid, vn_rms_curvatures, rms_curvatures)
     call cdf_write(ncid, vn_max_curvatures, max_curvatures)
     call cdf_write(ncid, vn_max_modBinv_sqrt_half_grad_B_colon_grad_Bs, max_modBinv_sqrt_half_grad_B_colon_grad_Bs)
     call cdf_write(ncid, vn_axis_lengths, axis_lengths)
     call cdf_write(ncid, vn_standard_deviations_of_R, standard_deviations_of_R)
     call cdf_write(ncid, vn_standard_deviations_of_Z, standard_deviations_of_Z)
     call cdf_write(ncid, vn_axis_helicities, axis_helicities)
     call cdf_write(ncid, vn_B_helicities, B_helicities)
     call cdf_write(ncid, vn_effective_nfps, effective_nfps)
     call cdf_write(ncid, vn_Newton_tolerance_achieveds, Newton_tolerance_achieveds)
     call cdf_write(ncid, vn_iota_tolerance_achieveds, iota_tolerance_achieveds)
     call cdf_write(ncid, vn_elongation_tolerance_achieveds, elongation_tolerance_achieveds)
     call cdf_write(ncid, vn_R0s_min, R0s_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0s_max, R0s_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0s_N_scan, R0s_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0c_min, R0c_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0c_max, R0c_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0c_N_scan, R0c_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s_min, Z0s_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s_max, Z0s_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s_N_scan, Z0s_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c_min, Z0c_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c_max, Z0c_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c_N_scan, Z0c_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_scan_sigma_initial, scan_sigma_initial)
     call cdf_write(ncid, vn_sigma_initial_values, sigma_initial_values)
     call cdf_write(ncid, vn_scan_eta_bar, scan_eta_bar)
     call cdf_write(ncid, vn_eta_bar_values, eta_bar_values)
  end select

  ! Arrays with dimension 2

  select case (trim(general_option))
  case (general_option_single)
     call cdf_write(ncid, vn_RBC, RBC(-ntor:ntor, 0:mpol_nonzero))
     call cdf_write(ncid, vn_RBS, RBS(-ntor:ntor, 0:mpol_nonzero))
     call cdf_write(ncid, vn_ZBC, ZBC(-ntor:ntor, 0:mpol_nonzero))
     call cdf_write(ncid, vn_ZBS, ZBS(-ntor:ntor, 0:mpol_nonzero))
  case (general_option_scan)
     call cdf_write(ncid, vn_N_scan_array,  N_scan_array(1:axis_nmax+1,:))
     call cdf_write(ncid, vn_scan_R0c,  scan_R0c)
     call cdf_write(ncid, vn_scan_R0s,  scan_R0s)
     call cdf_write(ncid, vn_scan_Z0c,  scan_Z0c)
     call cdf_write(ncid, vn_scan_Z0s,  scan_Z0s)
  end select

  ! Arrays with dimension 3

  !call cdf_write(ncid, vn_r_plasma, r_plasma)

  ! Finish up:
  call cdf_close(ncid)

  call cpu_time(end_time)
  print "(a,es10.3,a)"," Time to write output file:",end_time - time1," sec"

end subroutine quasisymmetry_write_output
