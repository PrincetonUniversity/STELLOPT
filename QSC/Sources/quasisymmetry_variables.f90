module quasisymmetry_variables

  use stel_kinds

  implicit none

  real(dp), parameter :: pi = 3.14159265358979d+0

  character(len=*), parameter :: &
       resolution_option_fixed = "fixed", &
       resolution_option_adaptive = "adaptive"
  character(len=50) :: resolution_option = resolution_option_fixed
  ! "fixed"    = Run using the specified N_phi.
  ! "adaptive" = Keep doubling N_phi (approximately, so N_phi remains odd) until iota_tolerance is achieved, or N_phi > max_N_phi.

  character(len=*), parameter :: &
       general_option_single = "single", &
       general_option_scan = "scan"
  character(len=50) :: general_option = general_option_single

  character(len=*), parameter :: &
       verbose_option_all = "all", &
       verbose_option_proc0 = "proc0", &
       verbose_option_summary = "summary"
  character(len=50) :: verbose_option = verbose_option_all

  character(len=*), parameter :: &
       eta_bar_scan_option_linear = "linear", &
       eta_bar_scan_option_log = "log", &
       eta_bar_scan_option_2_sided_log = "2_sided_log"
  character(len=50) :: eta_bar_scan_option = eta_bar_scan_option_linear

  character(len=*), parameter :: &
       sigma_initial_scan_option_linear = "linear", &
       sigma_initial_scan_option_log = "log", &
       sigma_initial_scan_option_2_sided_log = "2_sided_log"
  character(len=50) :: sigma_initial_scan_option = sigma_initial_scan_option_linear

  character(len=*), parameter :: &
       Fourier_scan_option_linear = "linear", &
       Fourier_scan_option_2_sided_log = "2_sided_log", &
       Fourier_scan_option_2_sided_log_except_Z0s1 = "2_sided_log_except_Z0s1"
  character(len=50) :: Fourier_scan_option = Fourier_scan_option_linear

  character(len=*), parameter :: &
       finite_r_option_linear = "linear", &
       finite_r_option_nonlinear = "nonlinear"
  character(len=50) :: finite_r_option = finite_r_option_linear

  character(len=*), parameter :: &
       order_r_option_r1 = "r1", &
       order_r_option_r2 = "r2", &
       order_r_option_r3_simplified = "r3_simplified", &
       order_r_option_r3_simplified_with_Z3 = "r3_simplified_with_Z3", &
       order_r_option_r3_flux_constraint = "r3_flux_constraint", &
       order_r_option_r3_flux_constraint_const_B20 = "r3_flux_constraint_const_B20", &
       order_r_option_r3_B3 = "r3_B3", &
       order_r_option_r3_X3s3_X3c3 = "r3_X3s3_X3c3", &
       order_r_option_r3_X3s3_Y3s3 = "r3_X3s3_Y3s3", &
       order_r_option_r3_X3c3_Y3c3 = "r3_X3c3_Y3c3", &
       order_r_option_r3_Y3s3_Y3c3 = "r3_Y3s3_Y3c3"
  character(len=50) :: order_r_option = order_r_option_r1

  real(dp) :: sigma_initial = 0

  integer :: nfp = 3

  integer :: sign_G = 1
  integer :: sign_psi = 1
  real(dp) :: I2_over_B0 = 0

  integer :: N_iterations = 20
  integer :: N_line_search = 10
  real(dp) :: Newton_tolerance = 1.0d-12
  real(dp) :: iota_tolerance = 1.0d-6
  real(dp) :: elongation_tolerance = 1.0d-2

  integer :: N_phi     !12/4/18.(7l23a) '= 15' out.
  integer :: N_phi_original = 15  !12/4/18.(7l23a)'= 15' in.
  integer :: max_N_phi = 100

!  integer, parameter :: max_axis_nmax = 1
!1/14/19.max_axis_nmax=1(orig:7m12a),10(8k8d),6(8k8e),1(8k9),6(8k11a),3(8x5).
  integer, parameter :: max_axis_nmax = 3
  integer :: axis_nmax = 1
  real(dp), dimension(max_axis_nmax + 1) :: R0s, R0c, Z0s, Z0c ! Fourier coefficients for the magnetic axis  
  real(dp) :: eta_bar   !3/1/19.QSC2 has B1s_over_B0, B1c_over_B0 here.
  integer :: max_n   !1/29/19.(8o7)xferred fra qs_rd_input.

  real(dp) :: max_max_curvature_to_keep = 5.0d+0
  real(dp) :: min_iota_to_keep = 0.05d+0
  real(dp) :: max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep = 5.0d+0

  integer :: matrix_size, axis_helicity, B_helicity, effective_nfp
  real(dp) :: last_iota, last_max_elongation, d_phi
  real(dp), dimension(:,:), allocatable :: d_d_phi, d_d_zeta
  real(dp), dimension(:), allocatable :: phi_extended, R0_extended, Z0_extended
  real(dp), dimension(:), allocatable :: phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp, Boozer_toroidal_angle
  real(dp), dimension(:), allocatable :: d_l_d_phi, curvature, torsion, B1Squared_over_curvatureSquared
  real(dp), dimension(:,:), allocatable :: tangent_cylindrical, normal_cylindrical, binormal_cylindrical
  real(dp), dimension(:,:), allocatable :: tangent_Cartesian, normal_Cartesian, binormal_Cartesian
  real(dp), dimension(:), allocatable :: sigma, X1s, X1c, Y1s, Y1c, R1s, R1c, Z1s, Z1c, elongation, elongation_in_Rz_plane
  real(dp), dimension(:), allocatable :: X1s_untwisted, X1c_untwisted, Y1s_untwisted, Y1c_untwisted
  real(dp) :: B0_over_abs_G0, abs_G0_over_B0, iota, max_elongation, rms_curvature, max_curvature, axis_length
  real(dp) :: max_precise_elongation = 10 ! Above this value, we won't do a precise solve, just take maxval over the phi grid.
  real(dp) :: max_elongation_to_keep = 10 ! Discard solutions with max(elongation) higher than this value. Set to e.g. 1.0e200 to keep all solutions.
  real(dp), dimension(:,:), allocatable :: Jacobian
  real(dp), dimension(:), allocatable :: residual, step_direction
  logical :: already_found_max_curvature, skipped_solve
  logical :: consider_only_nfp = .false.
  integer :: dimension_Fourier = 0
  real(dp), dimension(:,:), allocatable :: sin_n_phi, cos_n_phi

  character(len=200) :: output_filename
  character(len=200) :: vmec_template_filename = ''
  character(len=200) :: new_vmec_filename

  real(dp) :: start_time, total_time

  real(dp), dimension(max_axis_nmax+1) :: R0s_min, R0s_max, R0c_min, R0c_max, Z0s_min, Z0s_max, Z0c_min, Z0c_max
  real(dp) :: eta_bar_min = 1, eta_bar_max = 1, sigma_initial_min = 0, sigma_initial_max = 0
  integer, dimension(max_axis_nmax+1) :: R0s_N_scan=0, R0c_N_scan=0, Z0s_N_scan=0, Z0c_N_scan=0
  integer :: eta_bar_N_scan=0, sigma_initial_N_scan=0
  integer*8 :: N_scan
  real(dp), dimension(:), allocatable :: iotas, max_elongations, rms_curvatures, max_curvatures, axis_lengths, eta_bar_values, sigma_initial_values
  real(dp), dimension(:), allocatable :: max_modBinv_sqrt_half_grad_B_colon_grad_Bs
  real(dp), dimension(:), allocatable :: standard_deviations_of_R, standard_deviations_of_Z
  real(dp) :: standard_deviation_of_R, standard_deviation_of_Z
  integer, dimension(:), allocatable :: axis_helicities, B_helicities, effective_nfps
  logical, dimension(:), allocatable :: iota_tolerance_achieveds, elongation_tolerance_achieveds, Newton_tolerance_achieveds
  logical :: iota_tolerance_achieved, elongation_tolerance_achieved, Newton_tolerance_achieved
  integer, dimension(max_axis_nmax+1, 4) :: N_scan_array
  logical :: untwist = .true.
  real(dp), dimension(:), allocatable :: modBinv_sqrt_half_grad_B_colon_grad_B
  real(dp) :: max_modBinv_sqrt_half_grad_B_colon_grad_B
  real(dp), dimension(:), allocatable :: B0_order_a_squared_to_cancel

  real(dp), dimension(:), allocatable :: scan_eta_bar, scan_sigma_initial
  real(dp), dimension(:,:), allocatable :: scan_R0c, scan_R0s, scan_Z0c, scan_Z0s
  real(dp) :: r = 0.1d+0
  integer :: mpol_nonzero

  real(dp), dimension(:), allocatable :: X20, X2s, X2c, Y20, Y2s, Y2c, Z20, Z2s, Z2c, B20
  real(dp), dimension(:), allocatable :: X20_untwisted, X2s_untwisted, X2c_untwisted
  real(dp), dimension(:), allocatable :: Y20_untwisted, Y2s_untwisted, Y2c_untwisted, Z20_untwisted, Z2s_untwisted, Z2c_untwisted
  real(dp), dimension(:), allocatable :: R20, R2s, R2c, z20_cylindrical, z2s_cylindrical, z2c_cylindrical
  real(dp) :: B2s = 0.0d+0
  real(dp) :: B2c = 0.0d+0
  real(dp) :: B3s3_input = 0.0d+0
  real(dp) :: B3c3_input = 0.0d+0
  real(dp) :: B0 = 1.0d+0
  real(dp) :: p2 = 0.0d+0
  real(dp), parameter :: mu0 = 1.25663706143592d-6
  real(dp) :: B20_mean, B20_residual, iota2
  real(dp) :: Y3c1_initial

  real(dp), dimension(:), allocatable :: d_X1c_d_zeta, d_Y1c_d_zeta, d_Y1s_d_zeta
  real(dp), dimension(:), allocatable :: X3s1, X3s3, X3c1, X3c3, Y3s1, Y3s3, Y3c1, Y3c3, Z3s1, Z3s3, Z3c1, Z3c3
  real(dp), dimension(:), allocatable :: B3s1, B3s3, B3c1, B3c3
  real(dp), dimension(:), allocatable :: R3s1, R3s3, R3c1, R3c3, z3s1_cylindrical, z3s3_cylindrical, z3c1_cylindrical, z3c3_cylindrical
  real(dp), dimension(:), allocatable :: X3s1_untwisted, X3s3_untwisted, X3c1_untwisted, X3c3_untwisted
  real(dp), dimension(:), allocatable :: Y3s1_untwisted, Y3s3_untwisted, Y3c1_untwisted, Y3c3_untwisted
  real(dp), dimension(:), allocatable :: Z3s1_untwisted, Z3s3_untwisted, Z3c1_untwisted, Z3c3_untwisted
  real(dp) :: iota_from_torsion
  logical :: circular_cross_section_surface = .false.
  integer :: finite_r_nonlinear_N_theta = 20

  integer :: N_procs, mpi_rank
  logical :: proc0, verbose = .true.

  character(len=120) :: xtqsc  !hm-10/14/18,6/4/19.

  namelist / quasisymmetry / resolution_option, general_option, verbose_option, nfp, sign_G, sign_psi, I2_over_B0, vmec_template_filename, r, &
       N_iterations, N_line_search, Newton_tolerance, iota_tolerance, elongation_tolerance, max_precise_elongation, max_elongation_to_keep, N_phi, max_N_phi, &
       R0s, R0c, Z0s, Z0c, eta_bar, sigma_initial, eta_bar_scan_option, sigma_initial_scan_option, Fourier_scan_option, consider_only_nfp, &
       R0s_min, R0s_max, R0s_N_scan, R0c_min, R0c_max, R0c_N_scan, Z0s_min, Z0s_max, Z0s_N_scan, Z0c_min, Z0c_max, Z0c_N_scan, &
       eta_bar_min, eta_bar_max, eta_bar_N_scan, sigma_initial_min, sigma_initial_max, sigma_initial_N_scan, max_max_curvature_to_keep, min_iota_to_keep, &
       finite_r_option, order_r_option, B0, B2s, B2c, p2, untwist, max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep, B3s3_input, B3c3_input, Y3c1_initial, &
       circular_cross_section_surface, finite_r_nonlinear_N_theta

end module quasisymmetry_variables

