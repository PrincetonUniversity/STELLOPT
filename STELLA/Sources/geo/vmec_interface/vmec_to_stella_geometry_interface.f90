! vmec_to_stella_geometry_interface.f90
! Initial code written by Matt Landreman, University of Maryland in August 2017.
! Skip down ~25 lines for detailed description of the input and output parameters.
! Modified 2018-2019 by Michael Barnes

module vmec_to_stella_geometry_interface_mod

  implicit none

  private

  public :: vmec_to_stella_geometry_interface
  public :: get_nominal_vmec_zeta_grid
  public :: read_vmec_equilibrium

  real :: theta_pest_target, zeta0
  real, dimension(2) :: vmec_radial_weight_full, vmec_radial_weight_half
  integer, dimension(2) :: vmec_radial_index_full, vmec_radial_index_half

  logical :: lasym
  integer :: nfp, isigng
  integer :: ns, mnmax, mnmax_nyq
  integer :: mpol, ntor
  real :: Aminor
  real, dimension (:), allocatable :: xm, xn
  real, dimension (:), allocatable :: xm_nyq, xn_nyq
  real, dimension (:,:), allocatable :: rmnc, rmns
  real, dimension (:,:), allocatable :: lmnc, lmns
  real, dimension (:,:), allocatable :: zmnc, zmns
  real, dimension (:,:), allocatable :: bmnc, bmns
  real, dimension (:,:), allocatable :: gmnc, gmns
  real, dimension (:,:), allocatable :: bsupumnc, bsupumns
  real, dimension (:,:), allocatable :: bsupvmnc, bsupvmns
  real, dimension (:,:), allocatable :: bsubumnc, bsubumns
  real, dimension (:,:), allocatable :: bsubvmnc, bsubvmns
  real, dimension (:,:), allocatable :: bsubsmnc, bsubsmns
  real, dimension (:), allocatable :: phi, phip, iotas, iotaf, presf

contains

  subroutine read_vmec_equilibrium (vmec_filename)

    use read_wout_mod, only: read_wout_file, read_wout_deallocate
    use read_wout_mod, only: nfp_vmec => nfp
    use read_wout_mod, only: lasym_vmec => lasym
    use read_wout_mod, only: isigng_vmec => isigng
    use read_wout_mod, only: Aminor_vmec => Aminor
    use read_wout_mod, only: ns_vmec => ns
    use read_wout_mod, only: mnmax_nyq_vmec => mnmax_nyq
    use read_wout_mod, only: mnmax_vmec => mnmax
    use read_wout_mod, only: mpol_vmec => mpol
    use read_wout_mod, only: ntor_vmec => ntor
    use read_wout_mod, only: xm_vmec => xm
    use read_wout_mod, only: xn_vmec => xn
    use read_wout_mod, only: xm_nyq_vmec => xm_nyq
    use read_wout_mod, only: xn_nyq_vmec => xn_nyq
    use read_wout_mod, only: phi_vmec => phi
    use read_wout_mod, only: phip_vmec => phip
    use read_wout_mod, only: lmnc_vmec => lmnc
    use read_wout_mod, only: lmns_vmec => lmns
    use read_wout_mod, only: rmnc_vmec => rmnc
    use read_wout_mod, only: rmns_vmec => rmns
    use read_wout_mod, only: zmnc_vmec => zmnc
    use read_wout_mod, only: zmns_vmec => zmns
    use read_wout_mod, only: bmnc_vmec => bmnc
    use read_wout_mod, only: bmns_vmec => bmns
    use read_wout_mod, only: gmnc_vmec => gmnc
    use read_wout_mod, only: gmns_vmec => gmns
    use read_wout_mod, only: bsupumnc_vmec => bsupumnc
    use read_wout_mod, only: bsupvmnc_vmec => bsupvmnc
    use read_wout_mod, only: bsupumns_vmec => bsupumns
    use read_wout_mod, only: bsupvmns_vmec => bsupvmns
    use read_wout_mod, only: bsubumnc_vmec => bsubumnc
    use read_wout_mod, only: bsubvmnc_vmec => bsubvmnc
    use read_wout_mod, only: bsubumns_vmec => bsubumns
    use read_wout_mod, only: bsubvmns_vmec => bsubvmns
    use read_wout_mod, only: bsubsmnc_vmec => bsubsmnc
    use read_wout_mod, only: bsubsmns_vmec => bsubsmns
    use read_wout_mod, only: iotas_vmec => iotas
    use read_wout_mod, only: iotaf_vmec => iotaf
    use read_wout_mod, only: presf_vmec => presf

    implicit none

    ! vmec_filename is the vmec wout_* file that will be read.
    character (*), intent (in) :: vmec_filename

    integer :: ierr, iopen

    !*********************************************************************
    ! Read in everything from the vmec wout file using libstell.
    !*********************************************************************

    write (*,'(A)') "############################################################"
    write (*,'(A)') "                       MAGNETIC FIELD"
    write (*,'(A)') "############################################################"
    write (*,*) "About to read VMEC wout file: '",trim(vmec_filename),"'."
    call read_wout_file(vmec_filename, ierr, iopen)
    if (iopen .ne. 0) stop 'error opening wout file'
    if (ierr .ne. 0) stop 'error reading wout file'
    write (*,*) "Successfully read VMEC data from '",trim(vmec_filename),"'."

    nfp = nfp_vmec
    lasym = lasym_vmec
    isigng = isigng_vmec
    aminor = aminor_vmec
    ns = ns_vmec
    mnmax = mnmax_vmec
    mnmax_nyq = mnmax_nyq_vmec
    mpol = mpol_vmec
    ntor = ntor_vmec

    write (*,*) " "
    write (*,*) "  Characteristics of the magnetic field:"
    write (*,'(A44, I1)') "      Number of field periods (nfp):"//REPEAT(' ',50),nfp
    write (*,'(A44, L1)') "      Stellarator-asymmetric? (lasym):"//REPEAT(' ',50),lasym

    if (.not.allocated(rmnc)) then
       allocate (xm(mnmax)) ; xm = xm_vmec
       allocate (xn(mnmax)) ; xn = xn_vmec
       allocate (xm_nyq(mnmax_nyq)) ; xm_nyq = xm_nyq_vmec
       allocate (xn_nyq(mnmax_nyq)) ; xn_nyq = xn_nyq_vmec
       allocate (rmnc(mnmax,ns)) ; rmnc = rmnc_vmec
       allocate (lmns(mnmax,ns)) ; lmns = lmns_vmec
       allocate (zmns(mnmax,ns)) ; zmns = zmns_vmec
       allocate (bmnc(mnmax_nyq,ns)) ; bmnc = bmnc_vmec
       allocate (gmnc(mnmax_nyq,ns)) ; gmnc = gmnc_vmec
       allocate (bsupumnc(mnmax_nyq,ns)) ; bsupumnc = bsupumnc_vmec
       allocate (bsupvmnc(mnmax_nyq,ns)) ; bsupvmnc = bsupvmnc_vmec
       allocate (bsubumnc(mnmax_nyq,ns)) ; bsubumnc = bsubumnc_vmec
       allocate (bsubvmnc(mnmax_nyq,ns)) ; bsubvmnc = bsubvmnc_vmec
       allocate (bsubsmns(mnmax_nyq,ns)) ; bsubsmns = bsubsmns_vmec
       allocate (phi(ns)) ; phi = phi_vmec
       allocate (phip(ns)) ; phip = phip_vmec
       allocate (iotas(ns)) ; iotas = iotas_vmec
       allocate (iotaf(ns)) ; iotaf = iotaf_vmec
       allocate (presf(ns)) ; presf = presf_vmec
       if (lasym) then
          allocate (rmns(mnmax,ns)) ; rmns = rmns_vmec
          allocate (lmnc(mnmax,ns)) ; lmnc = lmnc_vmec
          allocate (zmnc(mnmax,ns)) ; zmnc = zmnc_vmec
          allocate (bmns(mnmax_nyq,ns)) ; bmns = bmns_vmec
          allocate (gmns(mnmax_nyq,ns)) ; gmns = gmns_vmec
          allocate (bsupumns(mnmax_nyq,ns)) ; bsupumns = bsupumns_vmec
          allocate (bsupvmns(mnmax_nyq,ns)) ; bsupvmns = bsupvmns_vmec
          allocate (bsubumns(mnmax_nyq,ns)) ; bsubumns = bsubumns_vmec
          allocate (bsubvmns(mnmax_nyq,ns)) ; bsubvmns = bsubvmns_vmec
          allocate (bsubsmnc(mnmax_nyq,ns)) ; bsubsmnc = bsubsmnc_vmec
       end if
    end if

    ! deallocate all arrays opened externally in read_wout_mod
    call read_wout_deallocate
    
  end subroutine read_vmec_equilibrium

  subroutine get_nominal_vmec_zeta_grid (nzgrid, zeta_center, number_of_field_periods_stella, &
       number_of_field_periods_device, zeta)

    implicit none

    ! 2*nzgrid+1 is the number of zeta grid points for the nominal zeta grid
    integer, intent (in) :: nzgrid
    ! The zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
    ! yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
    real, intent (in) :: zeta_center
    ! number_of_field_periods_device is the number of field periods sampled by stella
    real, intent (in out) :: number_of_field_periods_stella
    ! number_of_field_periods_device is the number of field periods for the device
    real, intent (out) :: number_of_field_periods_device
    ! On exit, zeta holds the nominal grid points in the toroidal angle zeta
    real, dimension (-nzgrid:), intent (out) :: zeta

    real, parameter :: pi = 3.1415926535897932d+0
    integer :: j

    number_of_field_periods_device = nfp

    if (number_of_field_periods_stella < 0.0) &
         number_of_field_periods_stella = number_of_field_periods_device

    zeta = [( zeta_center + (pi*j*number_of_field_periods_stella)/(nfp*nzgrid), j=-nzgrid, nzgrid )]

  end subroutine get_nominal_vmec_zeta_grid

  subroutine vmec_to_stella_geometry_interface(nalpha, alpha0, nzgrid, &
       zeta_center, number_of_field_periods_to_include, &
       desired_normalized_toroidal_flux, vmec_surface_option, verbose, &
       normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, nfp_out, &
       sign_toroidal_flux, &
       alpha, zeta, bmag, gradpar_zeta, grad_alpha_grad_alpha, &
       grad_alpha_grad_psi, grad_psi_grad_psi, gds23, gds24, gds25, gds26, &
       gbdrift_alpha, gbdrift0_psi, cvdrift_alpha, cvdrift0_psi, &
       theta_vmec, B_sub_zeta, B_sub_theta_vmec)

    use fzero_mod, only: fzero

    implicit none

    !*********************************************************************
    ! Input parameters
    !*********************************************************************

    ! nalpha is the number of grid points in the alpha coordinate:
    integer, intent(in) :: nalpha

    ! alpha0 is the first alpha value to include in the alpha grid
    real, intent(in) :: alpha0

    ! The zeta grid has nzgrid*2+1 points, including the "repeated" point at index -nzgrid and +nzgrid.
    integer, intent(in) :: nzgrid

    ! The zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
    ! yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
    real, intent(in) :: zeta_center

    ! If number_of_field_periods_to_include is > 0, then this parameter does what you think:
    ! the extent of the toroidal in zeta will be 2*pi*number_of_field_periods_to_include/nfp.
    ! If number_of_field_periods_to_include is <= 0, the entire 2*pi toroidal domain will be included.
    real, intent(in) :: number_of_field_periods_to_include

    ! The parameter desired_normalized_toroidal_flux determines which flux surface from the VMEC file will be used
    ! for the computation. This parameter should lie in the interval [0,1].
    real, intent(in) :: desired_normalized_toroidal_flux

    ! If vmec_surface_option = 0, the magnetic surface specified by desired_normalized_toroidal_flux will be used,
    ! by interpolating between the surfaces available in the vmec file.
    ! If vmec_surface_option = 1, the magnetic surface on vmec's HALF radial mesh will be used that is closest to desired_normalized_toroidal_flux.
    ! If vmec_surface_option = 2, the magnetic surface on vmec's FULL radial mesh will be used that is closest to desired_normalized_toroidal_flux.    
    ! Other values of vmec_surface_option will cause the program to abort with an error.
    integer, intent(in) :: vmec_surface_option

    ! If verbose is .true. in the vmec_parameters namelist, lots of diagnostic information is printed.
    logical, intent(in) :: verbose

    !*********************************************************************
    ! Output quantities
    !*********************************************************************

    ! On exit, normalized_toroidal_flux_used holds the flux surface that was actually used for the geometry,
    ! as measured by psi_toroidal / psi_{toroidal,edge}
    real, intent(out) :: normalized_toroidal_flux_used

    ! Safety factor q = 1/iota
    real, intent(out) :: safety_factor_q

    ! Magnetic shear shat = (x/q) * (d q / d x) where x = Aminor_p * sqrt(psi_toroidal / psi_{toroidal,edge})
    ! and Aminor_p is the minor radius calculated by VMEC.
    real, intent(out) :: shat

    ! L_reference is the reference length used for stella normalization, in meters.
    real, intent(out) :: L_reference

    ! B_reference is the reference magnetic field strength used for stella normalization, in Tesla.
    real, intent(out) :: B_reference

    ! nfp is the number of field periods given by VMEC
    real, intent(out) :: nfp_out

    integer, intent (out) :: sign_toroidal_flux

    ! On exit, alpha holds the grid points in alpha = theta_p - iota * zeta, where theta_p is the PEST toroidal angle
    real, dimension (:), intent(out) :: alpha

    ! On exit, zeta holds the grid points in the toroidal angle zeta
    real, dimension (-nzgrid:), intent(out) :: zeta

    real, dimension (:,-nzgrid:), intent(out) :: theta_vmec

    real, dimension (:,-nzgrid:), intent (out) :: bmag
    ! gradpar_zeta = b . grad zeta, with zeta the physical toroidal angle
    ! taken to increase in the counter-clockwise direction
    real, dimension (:,-nzgrid:), intent (out) :: gradpar_zeta
    ! grad alpha . grad alpha in units of 1/L_ref^2, with alpha = theta_pest - iota * zeta
    real, dimension (:,-nzgrid:), intent (out) :: grad_alpha_grad_alpha
    ! grad alpha . grad psi_t in units of B_ref, with alpha = theta_pest - iota * zeta
    real, dimension (:,-nzgrid:), intent (out) :: grad_alpha_grad_psi
    ! grad psi_t . grad psi_t in units of (a*B_ref)^2, with alpha = theta_pest - iota * zeta
    real, dimension (:,-nzgrid:), intent (out) :: grad_psi_grad_psi
    ! 2 * bhat/B x (grad B / B) . grad alpha * B_ref * L_ref^2
    real, dimension (:,-nzgrid:), intent (out) :: gbdrift_alpha
    ! 2 * bhat/B x (bhat . grad bhat) . grad alpha * B_ref * L_ref^2
    real, dimension (:,-nzgrid:), intent (out) :: cvdrift_alpha

    real, dimension (:,-nzgrid:), intent (out) :: gds23, gds24, gds25, gds26
!    real, dimension (:,-nzgrid:), intent (out) :: gbdrift, gbdrift0, cvdrift, cvdrift0
    real, dimension (:,-nzgrid:), intent (out) :: gbdrift0_psi, cvdrift0_psi

    real, dimension (:,-nzgrid:), intent (out) :: B_sub_theta_vmec, B_sub_zeta

    !*********************************************************************
    ! Variables used internally by this subroutine
    !*********************************************************************

    real, parameter :: pi = 3.1415926535897932d+0
    real, parameter :: zero = 0.0d+0
    real, parameter :: one = 1.0d+0
    real, parameter :: mu_0 = 4*pi*(1.0d-7)

    integer :: j, index, izeta, ialpha, isurf, m, n, imn, imn_nyq
    real :: angle, sin_angle, cos_angle, temp, edge_toroidal_flux_over_2pi
!    integer :: fzero_flag
    real :: number_of_field_periods_to_include_final
    real :: dphi, iota, min_dr2, ds, d_pressure_d_s, d_iota_d_s, scale_factor
    real :: theta_vmec_min, theta_vmec_max, sqrt_s
    real, dimension(:), allocatable :: dr2, normalized_toroidal_flux_full_grid, normalized_toroidal_flux_half_grid
    real, dimension(:), allocatable :: d_pressure_d_s_on_half_grid, d_iota_d_s_on_half_grid
!    real :: root_solve_absolute_tolerance, root_solve_relative_tolerance
    logical :: non_Nyquist_mode_available, found_imn
    real, dimension(:,:), allocatable :: B, sqrt_g, R, B_dot_grad_theta_pest_over_B_dot_grad_zeta, temp2D
    real, dimension(:,:), allocatable :: d_B_d_theta_vmec, d_B_d_zeta, d_B_d_s
    real, dimension(:,:), allocatable :: d_R_d_theta_vmec, d_R_d_zeta, d_R_d_s
    real, dimension(:,:), allocatable :: d_Z_d_theta_vmec, d_Z_d_zeta, d_Z_d_s
    real, dimension(:,:), allocatable :: d_X_d_theta_vmec, d_X_d_zeta, d_X_d_s
    real, dimension(:,:), allocatable :: d_Y_d_theta_vmec, d_Y_d_zeta, d_Y_d_s
    real, dimension(:,:), allocatable :: d_Lambda_d_theta_vmec, d_Lambda_d_zeta, d_Lambda_d_s
    !real, dimension(:,:), allocatable :: B_sub_s, B_sub_theta_vmec, B_sub_zeta
    real, dimension(:,:), allocatable :: B_sub_s
    real, dimension(:,:), allocatable :: B_sup_theta_vmec, B_sup_zeta
    real, dimension(:), allocatable :: d_B_d_s_mnc, d_B_d_s_mns
    real, dimension(:), allocatable :: d_R_d_s_mnc, d_R_d_s_mns
    real, dimension(:), allocatable :: d_Z_d_s_mnc, d_Z_d_s_mns
    real, dimension(:), allocatable :: d_Lambda_d_s_mnc, d_Lambda_d_s_mns
    real, dimension(:,:), allocatable :: grad_s_X, grad_s_Y, grad_s_Z
    real, dimension(:,:), allocatable :: grad_theta_vmec_X, grad_theta_vmec_Y, grad_theta_vmec_Z
    real, dimension(:,:), allocatable :: grad_theta_pest_X, grad_theta_pest_Y, grad_theta_pest_Z
    real, dimension(:,:), allocatable :: grad_zeta_X, grad_zeta_Y, grad_zeta_Z
    real, dimension(:,:), allocatable :: grad_psi_X, grad_psi_Y, grad_psi_Z
    real, dimension(:,:), allocatable :: grad_alpha_X, grad_alpha_Y, grad_alpha_Z
    real, dimension(:,:), allocatable :: B_cross_grad_B_dot_grad_alpha, B_cross_grad_B_dot_grad_alpha_alternate
    real, dimension(:,:), allocatable :: B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_alternate
    real, dimension(:,:), allocatable :: grad_B_X, grad_B_Y, grad_B_Z
    real, dimension(:,:), allocatable :: B_X, B_Y, B_Z
    real, dimension (:,:), allocatable :: gradzeta_grady, gradzeta_gradx
    real, dimension (:,:), allocatable :: gradtheta_grady, gradtheta_gradx

    logical :: theta_converged

    !*********************************************************************
    ! VMEC variables of interest:
    ! ns = number of flux surfaces used by VMEC
    ! nfp = number of field periods, e.g. 5 for W7-X, 4 for HSX
    ! iotas = rotational transform (1/q) on the half grid.
    ! iotaf = rotational transform on the full grid.
    ! presf = pressure on the full grid.
    !
    ! All VMEC quantities (B, pressure, etc) are in SI units.
    ! 
    ! In VMEC, quantities on the half grid have the same number of array elements (ns) as quantities on the full grid,
    ! but the first array element is 0.
    !
    !*********************************************************************

    !*********************************************************************
    ! Beginning of executable statements.
    !*********************************************************************

!    if (verbose) print *,"Entering subroutine vmec_to_stella_geometry_interface."

    !*********************************************************************
    ! Do some validation.
    !*********************************************************************

    if (nalpha<1) then
       print *,"Error! nalpha must be >= 1. Instead it is",nalpha
       stop
    end if

    if (nzgrid<1) then
       print *,"Error! nzgrid must be >= 1. Instead it is",nzgrid
       stop
    end if

    if (desired_normalized_toroidal_flux <= 0) then
       print *,"Error! desired_normalized_toroidal_flux must be >0. Instead it is",desired_normalized_toroidal_flux
       stop
    end if

    if (desired_normalized_toroidal_flux > 1) then
       print *,"Error! desired_normalized_toroidal_flux must be <= 1. Instead it is",desired_normalized_toroidal_flux
       stop
    end if

    nfp_out = nfp

    ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq arrays are sometimes
    ! not populated. The next few lines here provide a workaround:
    if (maxval(abs(xm_nyq)) < 1 .and. maxval(abs(xn_nyq)) < 1) then
       if (mnmax_nyq == mnmax) then
          if (verbose) print *,"xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
          xm_nyq = xm
          xn_nyq = xn
       else
          print *,"Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
          stop
       end if
    end if

    edge_toroidal_flux_over_2pi = phi(ns) / (2*pi) * isigng ! isigns is called signgs in the wout*.nc file. Why is this signgs here?

    ! this gives the sign of the edge toroidal flux
    sign_toroidal_flux = int(sign(1.1,edge_toroidal_flux_over_2pi))

    if (verbose) write (*,'(A43, I2)') "      Sign of the toroidal flux from VMEC:"//REPEAT(' ',50), sign_toroidal_flux

    ! Set reference length and magnetic field for stella's normalization, 
    ! using the choices made by Pavlos Xanthopoulos in GIST:
    L_reference = Aminor ! Note that 'Aminor' in read_wout_mod is called 'Aminor_p' in the wout*.nc file.
!    B_reference = 2 * edge_toroidal_flux_over_2pi / (L_reference * L_reference)
    B_reference = 2 * abs(edge_toroidal_flux_over_2pi) / (L_reference * L_reference)

    if (verbose) then
       write (*,*) "  "
       write (*,*) "  Reference values for the stella normalization:"
       write (*,'(A42, F15.12, A7)') "      Reference length (minor radius a):"//REPEAT(' ',50),L_reference," meters"
       write (*,'(A42, F15.12, A6)') "      Reference magnetic field strength:"//REPEAT(' ',50),B_reference," Tesla"
    end if

    ! --------------------------------------------------------------------------------
    ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
    ! --------------------------------------------------------------------------------

    ! 'phi' is vmec's array of the toroidal flux (not divided by 2pi!) on vmec's radial grid.
    if (abs(phi(1)) > 1d-14) then
       print *,"Error! VMEC phi array does not begin with 0."
       print *,"phi:",phi
       stop
    end if

    dphi = phi(2) - phi(1)
    do j=3,ns
       if (abs(phi(j)-phi(j-1)-dphi) > 1d-11) then
          print *,"Error! VMEC phi array is not uniformly spaced."
          print *,"phi:",phi
          stop
       end if
    end do

    ! The variable called 'phips' in the wout file is called just 'phip' in read_wout_mod.F.
    ! phips is on the half-mesh, so skip first point.
    do j=2,ns
       if (abs(phip(j)+phi(ns)/(2*pi)) > 1d-11) then
          print *,"Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
          print *,"phip(s):",phip
          stop
       end if
    end do

    ! The first mode in the m and n arrays should be m=n=0:
    if (xm(1) .ne. 0) stop "First element of xm in the wout file should be 0."
    if (xn(1) .ne. 0) stop "First element of xn in the wout file should be 0."
    if (xm_nyq(1) .ne. 0) stop "First element of xm_nyq in the wout file should be 0."
    if (xn_nyq(1) .ne. 0) stop "First element of xn_nyq in the wout file should be 0."

    ! Lambda should be on the half mesh, so its value at radial index 1 should be 0 for all (m,n)
    if (maxval(abs(lmns(:,1))) > 0) then
       print *,"Error! Expected lmns to be on the half mesh, but its value at radial index 1 is nonzero."
       print *,"Here comes lmns(:,1):", lmns(:,1)
       stop
    end if
    if (lasym) then
       if (maxval(abs(lmnc(:,1))) > 0) then
          print *,"Error! Expected lmnc to be on the half mesh, but its value at radial index 1 is nonzero."
          print *,"Here comes lmnc(:,1):", lmnc(:,1)
          stop
       end if
    end if

    ! --------------------------------------------------------------------------------
    ! End of sanity checks.
    ! --------------------------------------------------------------------------------

    allocate(normalized_toroidal_flux_full_grid(ns))
    normalized_toroidal_flux_full_grid = [( real(j-1)/(ns-1), j=1,ns )]

    ! Build an array of the half grid points:
    allocate(normalized_toroidal_flux_half_grid(ns-1))
    do j = 1,ns-1
       normalized_toroidal_flux_half_grid(j) = (normalized_toroidal_flux_full_grid(j) + normalized_toroidal_flux_full_grid(j+1))*(0.5d+0)
    end do

    !*********************************************************************
    ! Determine which flux surface to use, based on 
    ! desired_normalized_toroidal_flux and vmec_surface_option.
    !*********************************************************************

    ! Possible values of vmec_surface_option:
    ! 0 = Use the exact radius requested.
    ! 1 = Use the nearest value of the VMEC half grid.
    ! 2 = Use the nearest value of the VMEC full grid.

    select case (vmec_surface_option)
    case (0)
       ! Use exact radius requested.
       normalized_toroidal_flux_used = desired_normalized_toroidal_flux

    case (1)
       ! Use nearest value of the VMEC half grid

       ! Compute differences
       allocate(dr2(ns-1))
       dr2 = (normalized_toroidal_flux_half_grid - desired_normalized_toroidal_flux) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns-1
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       normalized_toroidal_flux_used = normalized_toroidal_flux_half_grid(index)
       deallocate(dr2)

    case (2)
       ! Use nearest value of the VMEC full grid

       ! Compute differences
       allocate(dr2(ns))
       dr2 = (normalized_toroidal_flux_full_grid - desired_normalized_toroidal_flux) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       normalized_toroidal_flux_used = normalized_toroidal_flux_full_grid(index)
       deallocate(dr2)

    case default
       print *,"Error! vmec_surface_option must be 0, 1, or 2. It is instead ",vmec_surface_option
       stop
    end select

    ! --------------------------------------------------------------------------------
    ! Done choosing the actual radius to use.
    ! --------------------------------------------------------------------------------

    ! In general, we get quantities for stella by linear interpolation, taking a weighted average of the quantity from
    ! 2 surfaces in the VMEC file. Sometimes the weights are 0 and 1, i.e. no interpolation is needed.

    ! For any VMEC quantity Q on the full grid, the value used in stella will be
    !  Q_stella = Q(vmec_radial_index_full(1))*vmec_radial_weight_full(1) + Q(vmec_radial_index_full(2))*vmec_radial_weight_full(2)

    ! For any VMEC quantity Q on the half grid, the value used in stella will be
    !  Q_stella = Q(vmec_radial_index_half(1))*vmec_radial_weight_half(1) + Q(vmec_radial_index_half(2))*vmec_radial_weight_half(2)


    ! Handle quantities for the full grid
    if (normalized_toroidal_flux_used>1) then
       stop "Error! normalized_toroidal_flux_used cannot be >1"
    elseif (normalized_toroidal_flux_used<0) then
       stop "Error! normalized_toroidal_flux_used cannot be <0"
    elseif (normalized_toroidal_flux_used==1) then
       vmec_radial_index_full(1) = ns-1
       vmec_radial_index_full(2) = ns
       vmec_radial_weight_full(1) = zero
    else
       ! normalized_toroidal_flux_used is >= 0 and <1
       ! This is the most common case.
       vmec_radial_index_full(1) = floor(normalized_toroidal_flux_used*(ns-1))+1
       vmec_radial_index_full(2) = vmec_radial_index_full(1) + 1
       vmec_radial_weight_full(1) = vmec_radial_index_full(1) - normalized_toroidal_flux_used*(ns-one)
    end if
    vmec_radial_weight_full(2) = one - vmec_radial_weight_full(1)

    ! Handle quantities for the half grid
    if (normalized_toroidal_flux_used < normalized_toroidal_flux_half_grid(1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the magnetic axis.) Results are likely to be inaccurate."

       ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
       vmec_radial_index_half(1) = 2
       vmec_radial_index_half(2) = 3
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_used) / (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_half_grid(1))

    elseif (normalized_toroidal_flux_used > normalized_toroidal_flux_half_grid(ns-1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the last closed flux surface.) Results may be inaccurate."
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_used) &
            / (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_half_grid(ns-2))

    elseif (normalized_toroidal_flux_used == normalized_toroidal_flux_half_grid(ns-1)) then
       ! We are exactly at the last point of the half grid
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = zero
    else
       ! normalized_toroidal_flux_used is inside the half grid.
       ! This is the most common case.
       vmec_radial_index_half(1) = floor(normalized_toroidal_flux_used*(ns-1) + 0.5d+0)+1
       if (vmec_radial_index_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmec_radial_index_half(1) = 2
       end if
       vmec_radial_index_half(2) = vmec_radial_index_half(1) + 1
       vmec_radial_weight_half(1) = vmec_radial_index_half(1) - normalized_toroidal_flux_used*(ns-one) - (0.5d+0)
    end if
    vmec_radial_weight_half(2) = one-vmec_radial_weight_half(1)

!     if (verbose) then
!        if (abs(vmec_radial_weight_half(1)) < 1e-14) then
!           print "(a,i4,a,i4,a)","   Using radial index ",vmec_radial_index_half(2)," of ",ns," from vmec's half mesh."
!        elseif (abs(vmec_radial_weight_half(2)) < 1e-14) then
!           print "(a,i4,a,i4,a)","   Using radial index ",vmec_radial_index_half(1)," of ",ns," from vmec's half mesh."
!        else
!           print "(a,i4,a,i4,a,i4,a)", "   Interpolating using radial indices ",vmec_radial_index_half(1)," and ",vmec_radial_index_half(2),&
!                " of ",ns," from vmec's half mesh."
!           print "(a,f17.14,a,f17.14)", "   Weights for half mesh = ",vmec_radial_weight_half(1)," and ",vmec_radial_weight_half(2)
!           print "(a,i4,a,i4,a,i4,a)", "   Interpolating using radial indices ",vmec_radial_index_full(1)," and ",vmec_radial_index_full(2),&
!                " of ",ns," from vmec's full mesh."
!           print "(a,f17.14,a,f17.14)", "   Weights for full mesh = ",vmec_radial_weight_full(1)," and ",vmec_radial_weight_full(2)
!        end if
!     end if



    !*********************************************************************
    ! Evaluate several radial-profile functions at the flux surface
    ! we ended up choosing.
    !*********************************************************************
    if (verbose) write (*,*) " "
    if (verbose) write (*,*) "  Radial-profile functions at the chosen flux surface:"
    iota = iotas(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + iotas(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    if (verbose) write (*,'(A21, F15.12)') "      iota:"//REPEAT(' ',50),iota
    safety_factor_q = 1/iota

    allocate(d_iota_d_s_on_half_grid(ns))
    d_iota_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    if (verbose) write (*,'(A21, ES20.12E3)') "      ds:"//REPEAT(' ',50),ds
    d_iota_d_s_on_half_grid(2:ns) = (iotaf(2:ns) - iotaf(1:ns-1)) / ds
    d_iota_d_s =  &
         d_iota_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_iota_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    deallocate(d_iota_d_s_on_half_grid)
    if (verbose) write (*,'(A21, ES20.12E3)') "      diota/ds:"//REPEAT(' ',50),d_iota_d_s
    ! shat = (r/q)(dq/dr) where r = a sqrt(s).
    !      = - (r/iota) (d iota / d r) = -2 (s/iota) (d iota / d s)
    shat = (-2 * normalized_toroidal_flux_used / iota) * d_iota_d_s

    allocate(d_pressure_d_s_on_half_grid(ns))
    d_pressure_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    d_pressure_d_s_on_half_grid(2:ns) = (presf(2:ns) - presf(1:ns-1)) / ds
    d_pressure_d_s =  &
         d_pressure_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_pressure_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    deallocate(d_pressure_d_s_on_half_grid)
    if (verbose) write (*,'(A21, ES20.12E3)') "      dpressure/ds:"//REPEAT(' ',50),d_pressure_d_s
    if (verbose) write (*,*) " "

    !*********************************************************************
    ! Set up the coordinate grids.
    !*********************************************************************

    alpha = [( alpha0 + ((j-1)*2*pi) / nalpha, j=1, nalpha )]

!!$    if (number_of_field_periods_to_include > nfp) then
!!$       print *,"Error! number_of_field_periods_to_include > nfp"
!!$       print *,"  number_of_field_periods_to_include =",number_of_field_periods_to_include
!!$       print *,"  nfp =",nfp
!!$       stop
!!$    end if
    number_of_field_periods_to_include_final = number_of_field_periods_to_include
    if (number_of_field_periods_to_include <= 0) then
       number_of_field_periods_to_include_final = nfp
       if (verbose) print *,"   Since number_of_field_periods_to_include was <= 0, it is being reset to nfp =",nfp
    end if

    zeta = [( zeta_center + (pi*j*number_of_field_periods_to_include_final)/(nfp*nzgrid), j=-nzgrid, nzgrid )]

    !*********************************************************************
    ! We know theta_pest = alpha + iota * zeta, but we need to determine
    ! theta_vmec = theta_pest - Lambda.
    !*********************************************************************

    if (verbose) print *,"Beginning root solves to determine theta_vmec."
    do izeta = -nzgrid, nzgrid
       zeta0 = zeta(izeta)
       do ialpha = 1,nalpha
          theta_pest_target = alpha(ialpha) + iota * zeta0
          ! Guess that theta_vmec will be within 0.3 radians of theta_pest:
          theta_vmec_min = theta_pest_target - 0.3
          theta_vmec_max = theta_pest_target + 0.3
          
          call get_root (theta_vmec_min, theta_vmec_max, theta_vmec(ialpha,izeta), theta_converged)
          ! In the 4th argument, we are telling the root-finder (fzero) to use theta_pest as the initial guess for theta_vmec.
!          call fzero(fzero_residual, theta_vmec_min, theta_vmec_max, theta_pest_target, &
!               root_solve_relative_tolerance, root_solve_absolute_tolerance, fzero_flag)
          ! Note: fzero returns its answer in theta_vmec_min.
!          theta_vmec(ialpha,izeta) = theta_vmec_min
!          if (fzero_flag == 4) then
!             stop "ERROR: fzero returned error 4: no sign change in residual"
!          else if (fzero_flag > 2) then
!             print *,"WARNING: fzero returned an error code:",fzero_flag
!          end if
          if (.not.theta_converged) then
             write (*,*) "ERROR: could not find root needed to compute theta_vmec. aborting"
             stop
          end if
       end do
    end do

!    if (verbose) then
!       do izeta = -nzgrid, nzgrid
!          do ialpha = 1, nalpha
!             write (*,*) 'theta_vmec', alpha(ialpha), zeta(izeta), theta_vmec(ialpha,izeta)
!          end do
!          write (*,*) 
!       end do
!    end if

    !*********************************************************************
    ! Initialize geometry arrays
    !*********************************************************************

    bmag = 0
    gradpar_zeta = 0
    grad_alpha_grad_alpha = 0.0
    grad_alpha_grad_psi = 0.0
    grad_psi_grad_psi = 0.0
    gds23 = 0.0
    gds24 = 0.0
    gds25 = 0.0
    gds26 = 0.0
    gbdrift_alpha = 0
    gbdrift0_psi = 0
    cvdrift_alpha = 0
    cvdrift0_psi = 0
    B_sub_theta_vmec = 0
    B_sub_zeta = 0

    allocate(B(nalpha,-nzgrid:nzgrid))
    allocate(temp2D(nalpha,-nzgrid:nzgrid))
    allocate(sqrt_g(nalpha,-nzgrid:nzgrid))
    allocate(R(nalpha,-nzgrid:nzgrid))
    allocate(d_B_d_theta_vmec(nalpha,-nzgrid:nzgrid))
    allocate(d_B_d_zeta(nalpha,-nzgrid:nzgrid))
    allocate(d_B_d_s(nalpha,-nzgrid:nzgrid))
    allocate(d_R_d_theta_vmec(nalpha,-nzgrid:nzgrid))
    allocate(d_R_d_zeta(nalpha,-nzgrid:nzgrid))
    allocate(d_R_d_s(nalpha,-nzgrid:nzgrid))
    allocate(d_Z_d_theta_vmec(nalpha,-nzgrid:nzgrid))
    allocate(d_Z_d_zeta(nalpha,-nzgrid:nzgrid))
    allocate(d_Z_d_s(nalpha,-nzgrid:nzgrid))
    allocate(d_Lambda_d_theta_vmec(nalpha,-nzgrid:nzgrid))
    allocate(d_Lambda_d_zeta(nalpha,-nzgrid:nzgrid))
    allocate(d_Lambda_d_s(nalpha,-nzgrid:nzgrid))
    allocate(B_sub_s(nalpha,-nzgrid:nzgrid))
    !allocate(B_sub_theta_vmec(nalpha,-nzgrid:nzgrid))
    !allocate(B_sub_zeta(nalpha,-nzgrid:nzgrid))
    allocate(B_sup_theta_vmec(nalpha,-nzgrid:nzgrid))
    allocate(B_sup_zeta(nalpha,-nzgrid:nzgrid))

    allocate(d_B_d_s_mnc(ns))
    allocate(d_B_d_s_mns(ns))
    allocate(d_R_d_s_mnc(ns))
    allocate(d_R_d_s_mns(ns))
    allocate(d_Z_d_s_mnc(ns))
    allocate(d_Z_d_s_mns(ns))
    allocate(d_Lambda_d_s_mnc(ns))
    allocate(d_Lambda_d_s_mns(ns))

    allocate(d_X_d_s(nalpha, -nzgrid:nzgrid))
    allocate(d_X_d_theta_vmec(nalpha, -nzgrid:nzgrid))
    allocate(d_X_d_zeta(nalpha, -nzgrid:nzgrid))
    allocate(d_Y_d_s(nalpha, -nzgrid:nzgrid))
    allocate(d_Y_d_theta_vmec(nalpha, -nzgrid:nzgrid))
    allocate(d_Y_d_zeta(nalpha, -nzgrid:nzgrid))

    allocate(grad_s_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_s_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_s_Z(nalpha, -nzgrid:nzgrid))
    allocate(grad_theta_vmec_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_theta_vmec_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_theta_vmec_Z(nalpha, -nzgrid:nzgrid))
    allocate(grad_theta_pest_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_theta_pest_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_theta_pest_Z(nalpha, -nzgrid:nzgrid))
    allocate(grad_zeta_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_zeta_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_zeta_Z(nalpha, -nzgrid:nzgrid))
    allocate(grad_psi_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_psi_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_psi_Z(nalpha, -nzgrid:nzgrid))
    allocate(grad_alpha_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_alpha_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_alpha_Z(nalpha, -nzgrid:nzgrid))
    
    allocate(B_X(nalpha, -nzgrid:nzgrid))
    allocate(B_Y(nalpha, -nzgrid:nzgrid))
    allocate(B_Z(nalpha, -nzgrid:nzgrid))
    allocate(grad_B_X(nalpha, -nzgrid:nzgrid))
    allocate(grad_B_Y(nalpha, -nzgrid:nzgrid))
    allocate(grad_B_Z(nalpha, -nzgrid:nzgrid))
    allocate(B_cross_grad_B_dot_grad_alpha(nalpha, -nzgrid:nzgrid))
    allocate(B_cross_grad_B_dot_grad_alpha_alternate(nalpha, -nzgrid:nzgrid))
    allocate(B_cross_grad_s_dot_grad_alpha(nalpha, -nzgrid:nzgrid))
    allocate(B_cross_grad_s_dot_grad_alpha_alternate(nalpha, -nzgrid:nzgrid))

    allocate (gradzeta_grady(nalpha,-nzgrid:nzgrid))
    allocate (gradzeta_gradx(nalpha,-nzgrid:nzgrid))
    allocate (gradtheta_grady(nalpha,-nzgrid:nzgrid))
    allocate (gradtheta_gradx(nalpha,-nzgrid:nzgrid))

    B = 0
    sqrt_g = 0
    R = 0
    d_B_d_theta_vmec = 0
    d_B_d_zeta = 0
    d_B_d_s = 0
    d_R_d_theta_vmec = 0
    d_R_d_zeta = 0
    d_R_d_s = 0
    d_Z_d_theta_vmec = 0
    d_Z_d_zeta = 0
    d_Z_d_s = 0
    d_Lambda_d_theta_vmec = 0
    d_Lambda_d_zeta = 0
    d_Lambda_d_s = 0
    B_sub_s = 0
    B_sub_theta_vmec = 0
    B_sub_zeta = 0
    B_sup_theta_vmec = 0
    B_sup_zeta = 0

    !*********************************************************************
    ! Now that we know the grid points in theta_vmec, we can evaluate
    ! all the geometric quantities on the grid points.
    !*********************************************************************

    do imn_nyq = 1, mnmax_nyq ! All the quantities we need except R, Z, and Lambda use the _nyq mode numbers.

       m = int(xm_nyq(imn_nyq))
       n = int(xn_nyq(imn_nyq)/nfp)

       if (abs(m) >= mpol .or. abs(n) > ntor) then
          non_Nyquist_mode_available = .false.
       else
          non_Nyquist_mode_available = .true.
          ! Find the imn in the non-Nyquist arrays that corresponds to the same m and n.
          found_imn = .false.
          do imn = 1,mnmax
             if (xm(imn)==m .and. xn(imn)==n*nfp) then
                found_imn = .true.
                exit
             end if
          end do
          if ((xm(imn) .ne. m) .or. (xn(imn) .ne. n*nfp)) stop "Something went wrong!"
          if (.not. found_imn) stop "Error! imn could not be found matching the given imn_nyq."
       end if

       ! All quantities are multiplied by a variable scale_factor which can in principle depend on m and n.
       ! For now we just set scale_factor = 1. In the future, scale_factor could be used to lower the
       ! symmetry-breaking Fourier components, or filter out certain Fourier components in some way.
       scale_factor = 1

       ! -----------------------------------------------------
       ! First, consider just the stellarator-symmetric terms:
       ! -----------------------------------------------------

       ! Evaluate the radial derivatives we will need:

       ! B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
       ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

       d_B_d_s_mnc(2:ns-1) = (bmnc(imn_nyq,3:ns) - bmnc(imn_nyq,2:ns-1)) / ds
       ! Simplistic extrapolation at the endpoints:
       d_B_d_s_mnc(1) = d_B_d_s_mnc(2)
       d_B_d_s_mnc(ns) = d_B_d_s_mnc(ns-1)

       if (non_Nyquist_mode_available) then
          ! R is on the full mesh:
          d_R_d_s_mnc(2:ns) = (rmnc(imn,2:ns) - rmnc(imn,1:ns-1)) / ds
          d_R_d_s_mnc(1) = 0

          ! Z is on the full mesh:
          d_Z_d_s_mns(2:ns) = (zmns(imn,2:ns) - zmns(imn,1:ns-1)) / ds
          d_Z_d_s_mns(1) = 0

          ! Lambda is on the half mesh:
          d_Lambda_d_s_mns(2:ns-1) = (lmns(imn,3:ns) - lmns(imn,2:ns-1)) / ds
          ! Simplistic extrapolation at the endpoints:
          d_Lambda_d_s_mns(1) = d_Lambda_d_s_mns(2)
          d_Lambda_d_s_mns(ns) = d_Lambda_d_s_mns(ns-1)
       else
          d_R_d_s_mnc = 0
          d_Z_d_s_mns = 0
          d_Lambda_d_s_mns = 0
       end if

       ! End of evaluating radial derivatives.

       do ialpha = 1,nalpha
          do izeta = -nzgrid, nzgrid
             angle = m * theta_vmec(ialpha,izeta) - n * nfp * zeta(izeta)
             cos_angle = cos(angle)
             sin_angle = sin(angle)

             do isurf = 1,2

                ! Handle |B|:
                temp = bmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B(ialpha,izeta) = B(ialpha,izeta) + temp * cos_angle
                d_B_d_theta_vmec(ialpha,izeta) = d_B_d_theta_vmec(ialpha,izeta) - m * temp * sin_angle
                d_B_d_zeta(ialpha,izeta) = d_B_d_zeta(ialpha,izeta) + n * nfp * temp * sin_angle       

                ! Handle Jacobian:
                temp = gmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                sqrt_g(ialpha,izeta) = sqrt_g(ialpha,izeta) + temp * cos_angle

                ! Handle B sup theta:
                temp = bsupumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sup_theta_vmec(ialpha,izeta) = B_sup_theta_vmec(ialpha,izeta) + temp * cos_angle

                ! Handle B sup zeta:
                temp = bsupvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sup_zeta(ialpha,izeta) = B_sup_zeta(ialpha,izeta) + temp * cos_angle

                ! Handle B sub theta:
                temp = bsubumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sub_theta_vmec(ialpha,izeta) = B_sub_theta_vmec(ialpha,izeta) + temp * cos_angle

                ! Handle B sub zeta:
                temp = bsubvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sub_zeta(ialpha,izeta) = B_sub_zeta(ialpha,izeta) + temp * cos_angle

                ! Handle B sub psi.
                ! Unlike the other components of B, this one is on the full mesh.
                temp = bsubsmns(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                temp = temp*scale_factor
                B_sub_s(ialpha,izeta) = B_sub_s(ialpha,izeta) + temp * sin_angle

                ! Handle d B / d s
                ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                temp = d_B_d_s_mnc(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                temp = temp*scale_factor
                d_B_d_s(ialpha,izeta) = d_B_d_s(ialpha,izeta) + temp * cos_angle

                ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                if (non_Nyquist_mode_available) then

                   ! Handle R, which is on the full mesh
                   temp = rmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   R(ialpha,izeta) = R(ialpha,izeta) + temp * cos_angle
                   d_R_d_theta_vmec(ialpha,izeta) = d_R_d_theta_vmec(ialpha,izeta) - temp * m * sin_angle
                   d_R_d_zeta(ialpha,izeta)  = d_R_d_zeta(ialpha,izeta)  + temp * n * nfp * sin_angle

                   ! Handle Z, which is on the full mesh
                   temp = zmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   !Z(ialpha,izeta) = Z(ialpha,izeta) + temp * sin_angle  ! We don't actually need Z itself, only derivatives of Z.
                   d_Z_d_theta_vmec(ialpha,izeta) = d_Z_d_theta_vmec(ialpha,izeta) + temp * m * cos_angle
                   d_Z_d_zeta(ialpha,izeta)  = d_Z_d_zeta(ialpha,izeta)  - temp * n * nfp * cos_angle

                   ! Handle Lambda:
                   temp = lmns(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   ! We don't need Lambda itself, just its derivatives.
                   d_Lambda_d_theta_vmec(ialpha,izeta) = d_Lambda_d_theta_vmec(ialpha,izeta) + m * temp * cos_angle
                   d_Lambda_d_zeta(ialpha,izeta) = d_Lambda_d_zeta(ialpha,izeta) - n * nfp * temp * cos_angle       

                   ! Handle d R / d s
                   ! Since R is on the full mesh, its radial derivative is on the half mesh.
                   temp = d_R_d_s_mnc(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   d_R_d_s(ialpha,izeta) = d_R_d_s(ialpha,izeta) + temp * cos_angle

                   ! Handle d Z / d s
                   ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                   temp = d_Z_d_s_mns(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   d_Z_d_s(ialpha,izeta) = d_Z_d_s(ialpha,izeta) + temp * sin_angle

                   ! Handle d Lambda / d s
                   ! Since Lambda is on the half mesh, its radial derivative is on the full mesh.
                   temp = d_Lambda_d_s_mns(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   d_Lambda_d_s(ialpha,izeta) = d_Lambda_d_s(ialpha,izeta) + temp * sin_angle

                end if
             end do
          end do
       end do

       ! -----------------------------------------------------
       ! Now consider the stellarator-asymmetric terms.
       ! -----------------------------------------------------

       if (lasym) then

          ! Evaluate the radial derivatives we will need:

          ! B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
          ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

          d_B_d_s_mns(2:ns-1) = (bmns(imn_nyq,3:ns) - bmns(imn_nyq,2:ns-1)) / ds
          ! Simplistic extrapolation at the endpoints:
          d_B_d_s_mns(1) = d_B_d_s_mns(2)
          d_B_d_s_mns(ns) = d_B_d_s_mns(ns-1)

          if (non_Nyquist_mode_available) then
             ! R is on the full mesh:
             d_R_d_s_mns(2:ns) = (rmns(imn,2:ns) - rmns(imn,1:ns-1)) / ds
             d_R_d_s_mns(1) = 0

             ! Z is on the full mesh:
             d_Z_d_s_mnc(2:ns) = (zmnc(imn,2:ns) - zmnc(imn,1:ns-1)) / ds
             d_Z_d_s_mnc(1) = 0

             ! Lambda is on the half mesh:
             d_Lambda_d_s_mnc(2:ns-1) = (lmnc(imn_nyq,3:ns) - lmnc(imn_nyq,2:ns-1)) / ds
             ! Simplistic extrapolation at the endpoints:
             d_Lambda_d_s_mnc(1) = d_Lambda_d_s_mnc(2)
             d_Lambda_d_s_mnc(ns) = d_Lambda_d_s_mnc(ns-1)
          else
             d_R_d_s_mns = 0
             d_Z_d_s_mnc = 0
             d_Lambda_d_s_mnc = 0
          end if

          ! End of evaluating radial derivatives.

          do ialpha = 1,nalpha
             do izeta = -nzgrid, nzgrid
                angle = m * theta_vmec(ialpha,izeta) - n * nfp * zeta(izeta)
                cos_angle = cos(angle)
                sin_angle = sin(angle)

                do isurf = 1,2

                   ! Handle |B|:
                   temp = bmns(imn_nyq,vmec_radial_index_half(1)) * vmec_radial_weight_half(1)
                   temp = temp * scale_factor
                   B(ialpha,izeta) = B(ialpha,izeta) + temp * sin_angle
                   d_B_d_theta_vmec(ialpha,izeta) = d_B_d_theta_vmec(ialpha,izeta) + m * temp * cos_angle
                   d_B_d_zeta(ialpha,izeta) = d_B_d_zeta(ialpha,izeta) - n * nfp * temp * cos_angle

                   ! Handle Jacobian:
                   temp = gmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   sqrt_g(ialpha,izeta) = sqrt_g(ialpha,izeta) + temp * sin_angle

                   ! Handle B sup theta:
                   temp = bsupumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sup_theta_vmec(ialpha,izeta) = B_sup_theta_vmec(ialpha,izeta) + temp * sin_angle

                   ! Handle B sup zeta:
                   temp = bsupvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sup_zeta(ialpha,izeta) = B_sup_zeta(ialpha,izeta) + temp * sin_angle

                   ! Handle B sub theta:
                   temp = bsubumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sub_theta_vmec(ialpha,izeta) = B_sub_theta_vmec(ialpha,izeta) + temp * sin_angle

                   ! Handle B sub zeta:
                   temp = bsubvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sub_zeta(ialpha,izeta) = B_sub_zeta(ialpha,izeta) + temp * sin_angle

                   ! Handle B sub psi.
                   ! Unlike the other components of B, this one is on the full mesh.
                   temp = bsubsmnc(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   B_sub_s(ialpha,izeta) = B_sub_s(ialpha,izeta) + temp * cos_angle

                   ! Handle d B / d s.
                   ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                   temp = d_B_d_s_mns(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   d_B_d_s(ialpha,izeta) = d_B_d_s(ialpha,izeta) + temp * sin_angle

                   ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                   if (non_Nyquist_mode_available) then

                      ! Handle R, which is on the full mesh
                      temp = rmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      R(ialpha,izeta) = R(ialpha,izeta) + temp * sin_angle
                      d_R_d_theta_vmec(ialpha,izeta) = d_R_d_theta_vmec(ialpha,izeta) + temp * m * cos_angle
                      d_R_d_zeta(ialpha,izeta)  = d_R_d_zeta(ialpha,izeta)  - temp * n * nfp * cos_angle

                      ! Handle Z, which is on the full mesh
                      temp = zmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      ! Z(ialpha,izeta) = Z(ialpha,izeta) + temp * cos_angle   ! We don't actually need Z itself, only derivatives of Z.
                      d_Z_d_theta_vmec(ialpha,izeta) = d_Z_d_theta_vmec(ialpha,izeta) - temp * m * sin_angle
                      d_Z_d_zeta(ialpha,izeta)  = d_Z_d_zeta(ialpha,izeta)  + temp * n * nfp * sin_angle

                      ! Handle Lambda, which is on the half mesh
                      temp = lmnc(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      ! We don't actually need Lambda itself, only derivatives of Lambda.
                      d_Lambda_d_theta_vmec(ialpha,izeta) = d_Lambda_d_theta_vmec(ialpha,izeta) - temp * m * sin_angle
                      d_Lambda_d_zeta(ialpha,izeta)  = d_Lambda_d_zeta(ialpha,izeta)  + temp * n * nfp * sin_angle

                      ! Handle d R / d s.
                      ! Since R is on the full mesh, its radial derivative is on the half mesh.
                      temp = d_R_d_s_mns(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      d_R_d_s(ialpha,izeta) = d_R_d_s(ialpha,izeta) + temp * sin_angle

                      ! Handle d Z / d s.
                      ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                      temp = d_Z_d_s_mnc(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      d_Z_d_s(ialpha,izeta) = d_Z_d_s(ialpha,izeta) + temp * cos_angle

                      ! Handle d Lambda / d s.
                      ! Since Lambda is on the half mesh, its radial derivative is on the full mesh.
                      temp = d_Lambda_d_s_mnc(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      d_Lambda_d_s(ialpha,izeta) = d_Lambda_d_s(ialpha,izeta) + temp * cos_angle
                   end if
                end do
             end do
          end do
       end if
    end do

    !*********************************************************************
    ! Sanity check: If the conversion to theta_pest has been done 
    ! correctly, we should find that 
    ! (B dot grad theta_pest) / (B dot grad zeta) = iota.
    ! Let's verify this:
    !*********************************************************************

    allocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta(nalpha, -nzgrid:nzgrid))
    ! Compute (B dot grad theta_pest) / (B dot grad zeta):
    B_dot_grad_theta_pest_over_B_dot_grad_zeta = (B_sup_theta_vmec * (1 + d_Lambda_d_theta_vmec) + B_sup_zeta * d_Lambda_d_zeta) / B_sup_zeta 
    temp2D = iota
    call test_arrays(B_dot_grad_theta_pest_over_B_dot_grad_zeta, temp2D, .false., 0.01, 'iota')
    deallocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta)

    !*********************************************************************
    ! Using R(theta,zeta) and Z(theta,zeta), compute the Cartesian
    ! components of the gradient basis vectors using the dual relations:
    !*********************************************************************

    do izeta = -nzgrid, nzgrid
       cos_angle = cos(zeta(izeta))
       sin_angle = sin(zeta(izeta))

       ! X = R * cos(zeta)
       d_X_d_theta_vmec(:,izeta) = d_R_d_theta_vmec(:,izeta) * cos_angle
       d_X_d_zeta(:,izeta) = d_R_d_zeta(:,izeta) * cos_angle - R(:,izeta) * sin_angle
       d_X_d_s(:,izeta) = d_R_d_s(:,izeta) * cos_angle

       ! Y = R * sin(zeta)
       d_Y_d_theta_vmec(:,izeta) = d_R_d_theta_vmec(:,izeta) * sin_angle
       d_Y_d_zeta(:,izeta) = d_R_d_zeta(:,izeta) * sin_angle + R(:,izeta) * cos_angle
       d_Y_d_s(:,izeta) = d_R_d_s(:,izeta) * sin_angle

!!$       ! Y = -R * sin(zeta)
!!$       d_Y_d_theta_vmec(:,izeta) = -d_R_d_theta_vmec(:,izeta) * sin_angle
!!$       d_Y_d_zeta(:,izeta) = -(d_R_d_zeta(:,izeta) * sin_angle + R(:,izeta) * cos_angle)
!!$       d_Y_d_s(:,izeta) = -d_R_d_s(:,izeta) * sin_angle

    end do

    ! Use the dual relations to get the Cartesian components of grad s, grad theta_vmec, and grad zeta:
    grad_s_X = (d_Y_d_theta_vmec * d_Z_d_zeta - d_Z_d_theta_vmec * d_Y_d_zeta) / sqrt_g
    grad_s_Y = (d_Z_d_theta_vmec * d_X_d_zeta - d_X_d_theta_vmec * d_Z_d_zeta) / sqrt_g
    grad_s_Z = (d_X_d_theta_vmec * d_Y_d_zeta - d_Y_d_theta_vmec * d_X_d_zeta) / sqrt_g

    grad_theta_vmec_X = (d_Y_d_zeta * d_Z_d_s - d_Z_d_zeta * d_Y_d_s) / sqrt_g
    grad_theta_vmec_Y = (d_Z_d_zeta * d_X_d_s - d_X_d_zeta * d_Z_d_s) / sqrt_g
    grad_theta_vmec_Z = (d_X_d_zeta * d_Y_d_s - d_Y_d_zeta * d_X_d_s) / sqrt_g

    grad_zeta_X = (d_Y_d_s * d_Z_d_theta_vmec - d_Z_d_s * d_Y_d_theta_vmec) / sqrt_g
    grad_zeta_Y = (d_Z_d_s * d_X_d_theta_vmec - d_X_d_s * d_Z_d_theta_vmec) / sqrt_g
    grad_zeta_Z = (d_X_d_s * d_Y_d_theta_vmec - d_Y_d_s * d_X_d_theta_vmec) / sqrt_g
    ! End of the dual relations.

    ! next get grad_theta_pest = grad (theta_vmec + Lambda)
    grad_theta_pest_X = (1.0 + d_Lambda_d_theta_vmec)*grad_theta_vmec_X &
         + d_Lambda_d_zeta*grad_zeta_X + d_Lambda_d_s*grad_s_X
    grad_theta_pest_Y = (1.0 + d_Lambda_d_theta_vmec)*grad_theta_vmec_Y &
         + d_Lambda_d_zeta*grad_zeta_Y + d_Lambda_d_s*grad_s_Y
    grad_theta_pest_Z = (1.0 + d_Lambda_d_theta_vmec)*grad_theta_vmec_Z &
         + d_Lambda_d_zeta*grad_zeta_Z + d_Lambda_d_s*grad_s_Z

    ! Sanity check: grad_zeta_X should be -sin(zeta) / R:
    do izeta = -nzgrid,nzgrid
       temp2D(:,izeta) = -sin(zeta(izeta)) / R(:,izeta)
    end do
    
    call test_arrays(grad_zeta_X, temp2D, .false., 1.0e-2, 'grad_zeta_X')
    grad_zeta_X = temp2D ! We might as well use the exact value, which is in temp2D.

    ! Sanity check: grad_zeta_Y should be cos(zeta) / R:
    do izeta = -nzgrid,nzgrid
       temp2D(:,izeta) = cos(zeta(izeta)) / R(:,izeta)
    end do
    call test_arrays(grad_zeta_Y, temp2D, .false., 1.0e-2, 'grad_zeta_Y')
    grad_zeta_Y = temp2D ! We might as well use the exact value, which is in temp2D.

    ! grad_zeta_Z should be 0:
    call test_arrays(grad_zeta_Z, temp2D, .true., 1.0e-14, 'grad_zeta_Z')
    grad_zeta_Z = 0

    !*********************************************************************
    ! Compute the Cartesian components of other quantities we need:
    !*********************************************************************

    grad_psi_X = grad_s_X * edge_toroidal_flux_over_2pi
    grad_psi_Y = grad_s_Y * edge_toroidal_flux_over_2pi
    grad_psi_Z = grad_s_Z * edge_toroidal_flux_over_2pi

    ! Form grad alpha = grad (theta_vmec + Lambda - iota * zeta)
    do izeta = -nzgrid,nzgrid
       grad_alpha_X(:,izeta) = (d_Lambda_d_s(:,izeta) - zeta(izeta) * d_iota_d_s) * grad_s_X(:,izeta)
       grad_alpha_Y(:,izeta) = (d_Lambda_d_s(:,izeta) - zeta(izeta) * d_iota_d_s) * grad_s_Y(:,izeta)
       grad_alpha_Z(:,izeta) = (d_Lambda_d_s(:,izeta) - zeta(izeta) * d_iota_d_s) * grad_s_Z(:,izeta)
    end do
    grad_alpha_X = grad_alpha_X + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_X + (-iota + d_Lambda_d_zeta) * grad_zeta_X
    grad_alpha_Y = grad_alpha_Y + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Y + (-iota + d_Lambda_d_zeta) * grad_zeta_Y
    grad_alpha_Z = grad_alpha_Z + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Z + (-iota + d_Lambda_d_zeta) * grad_zeta_Z

    grad_B_X = d_B_d_s * grad_s_X + d_B_d_theta_vmec * grad_theta_vmec_X + d_B_d_zeta * grad_zeta_X
    grad_B_Y = d_B_d_s * grad_s_Y + d_B_d_theta_vmec * grad_theta_vmec_Y + d_B_d_zeta * grad_zeta_Y
    grad_B_Z = d_B_d_s * grad_s_Z + d_B_d_theta_vmec * grad_theta_vmec_Z + d_B_d_zeta * grad_zeta_Z

    !temp2D = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_R_d_zeta + (iota - d_Lambda_d_zeta) * d_R_d_theta_vmec) / sqrt_g
!!$    do izeta = -nzgrid,nzgrid
!!$       !B_X(:,izeta) = temp2D(:,izeta) * cos(zeta(izeta))
!!$       !B_Y(:,izeta) = temp2D(:,izeta) * sin(zeta(izeta))
!!$       sin_angle = sin(zeta(izeta))
!!$       cos_angle = cos(zeta(izeta))
!!$       B_X(:,izeta) = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec(:,izeta)) * (d_R_d_zeta(:,izeta)*cos_angle - R(:,izeta)*sin_angle) &
!!$            + (iota - d_Lambda_d_zeta(:,izeta)) * d_R_d_theta_vmec(:,izeta)*cos_angle) / sqrt_g(:,izeta)
!!$       B_Y(:,izeta) = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec(:,izeta)) * (d_R_d_zeta(:,izeta)*sin_angle + R(:,izeta)*cos_angle) &
!!$            + (iota - d_Lambda_d_zeta(:,izeta)) * d_R_d_theta_vmec(:,izeta)*sin_angle) / sqrt_g(:,izeta)
!!$    end do
    B_X = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_X_d_zeta + (iota - d_Lambda_d_zeta) * d_X_d_theta_vmec) / sqrt_g
    B_Y = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Y_d_zeta + (iota - d_Lambda_d_zeta) * d_Y_d_theta_vmec) / sqrt_g
    B_Z = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Z_d_zeta + (iota - d_Lambda_d_zeta) * d_Z_d_theta_vmec) / sqrt_g

    sqrt_s = sqrt(normalized_toroidal_flux_used)

    !*********************************************************************
    ! Sanity tests: Verify that the Jacobian equals the appropriate
    ! cross product of the basis vectors.
    !*********************************************************************

    temp2D = 0 &
         + d_X_d_s * d_Y_d_theta_vmec * d_Z_d_zeta &
         + d_Y_d_s * d_Z_d_theta_vmec * d_X_d_zeta &
         + d_Z_d_s * d_X_d_theta_vmec * d_Y_d_zeta &
         - d_Z_d_s * d_Y_d_theta_vmec * d_X_d_zeta &
         - d_X_d_s * d_Z_d_theta_vmec * d_Y_d_zeta &
         - d_Y_d_s * d_X_d_theta_vmec * d_Z_d_zeta
    call test_arrays(sqrt_g, temp2D, .false., 3.0e-3, 'sqrt_g')

    temp2D = 0 &
         + grad_s_X * grad_theta_vmec_Y * grad_zeta_Z &
         + grad_s_Y * grad_theta_vmec_Z * grad_zeta_X &
         + grad_s_Z * grad_theta_vmec_X * grad_zeta_Y &
         - grad_s_Z * grad_theta_vmec_Y * grad_zeta_X &
         - grad_s_X * grad_theta_vmec_Z * grad_zeta_Y &
         - grad_s_Y * grad_theta_vmec_X * grad_zeta_Z
    call test_arrays(1/sqrt_g, temp2D, .false., 1.0e-2, '1/sqrt_g')
    
    !*********************************************************************
    ! Sanity tests: Verify that 
    ! \vec{B} dot (each of the covariant and contravariant basis vectors)
    ! matches the corresponding term from VMEC.
    !*********************************************************************

    call test_arrays(B_X * d_X_d_theta_vmec + B_Y * d_Y_d_theta_vmec + B_Z * d_Z_d_theta_vmec, B_sub_theta_vmec, .false., 1.0e-2, 'B_sub_theta_vmec')
    call test_arrays(B_X * d_X_d_s          + B_Y * d_Y_d_s          + B_Z * d_Z_d_s,          B_sub_s,          .false., 1.0e-2, 'B_sub_s')
    call test_arrays(B_X * d_X_d_zeta       + B_Y * d_Y_d_zeta       + B_Z * d_Z_d_zeta,       B_sub_zeta,       .false., 1.0e-2, 'B_sub_zeta')

    call test_arrays(B_X *          grad_s_X + B_Y *          grad_s_Y + B_Z *          grad_s_Z,           temp2D,  .true., 1.0e-2, 'B_sup_s')
    call test_arrays(B_X *       grad_zeta_X + B_Y *       grad_zeta_Y + B_Z *       grad_zeta_Z,       B_sup_zeta, .false., 1.0e-2, 'B_sup_zeta')
    call test_arrays(B_X * grad_theta_vmec_X + B_Y * grad_theta_vmec_Y + B_Z * grad_theta_vmec_Z, B_sup_theta_vmec, .false., 1.0e-2, 'B_sup_theta_vmec')

    !*********************************************************************
    ! For gbdrift, we need \vect{B} cross grad |B| dot grad alpha.
    ! For cvdrift, we also need \vect{B} cross grad s dot grad alpha.
    ! Let us compute both of these quantities 2 ways, and make sure the two
    ! approaches give the same answer (within some tolerance).
    !*********************************************************************

    B_cross_grad_s_dot_grad_alpha = (B_sub_zeta * (1 + d_Lambda_d_theta_vmec) &
         - B_sub_theta_vmec * (d_Lambda_d_zeta - iota) ) / sqrt_g

    B_cross_grad_s_dot_grad_alpha_alternate = 0 &
         + B_X * grad_s_Y * grad_alpha_Z &
         + B_Y * grad_s_Z * grad_alpha_X &
         + B_Z * grad_s_X * grad_alpha_Y &
         - B_Z * grad_s_Y * grad_alpha_X &
         - B_X * grad_s_Z * grad_alpha_Y &
         - B_Y * grad_s_X * grad_alpha_Z 

    call test_arrays(B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_alternate, &
         .false., 1.0e-2, 'B_cross_grad_s_dot_grad_alpha')

    do izeta = -nzgrid,nzgrid
       B_cross_grad_B_dot_grad_alpha(:,izeta) = 0 &
            + (B_sub_s(:,izeta) * d_B_d_theta_vmec(:,izeta) * (d_Lambda_d_zeta(:,izeta) - iota) &
            + B_sub_theta_vmec(:,izeta) * d_B_d_zeta(:,izeta) * (d_Lambda_d_s(:,izeta) - zeta(izeta) * d_iota_d_s) &
            + B_sub_zeta(:,izeta) * d_B_d_s(:,izeta) * (1 + d_Lambda_d_theta_vmec(:,izeta)) &
            - B_sub_zeta(:,izeta) * d_B_d_theta_vmec(:,izeta) * (d_Lambda_d_s(:,izeta) - zeta(izeta) * d_iota_d_s) &
            - B_sub_theta_vmec(:,izeta) * d_B_d_s(:,izeta) * (d_Lambda_d_zeta(:,izeta) - iota) &
            - B_sub_s(:,izeta) * d_B_d_zeta(:,izeta) * (1 + d_Lambda_d_theta_vmec(:,izeta))) / sqrt_g(:,izeta)
    end do

    B_cross_grad_B_dot_grad_alpha_alternate = 0 &
         + B_X * grad_B_Y * grad_alpha_Z &
         + B_Y * grad_B_Z * grad_alpha_X &
         + B_Z * grad_B_X * grad_alpha_Y &
         - B_Z * grad_B_Y * grad_alpha_X &
         - B_X * grad_B_Z * grad_alpha_Y &
         - B_Y * grad_B_X * grad_alpha_Z 

    call test_arrays(B_cross_grad_B_dot_grad_alpha, B_cross_grad_B_dot_grad_alpha_alternate, &
         .false., 1.0e-2, 'B_cross_grad_B_dot_grad_alpha')

    !*********************************************************************
    ! Finally, assemble the quantities needed for stella.
    !*********************************************************************

    ! See the latex note gs2_full_surface_stellarator_geometry in the "doc" directory for a derivation of the formulae that follow.

    bmag = B / B_reference

    gradpar_zeta = L_reference * B_sup_zeta / B

    ! grad alpha . grad alpha in units of 1/L_ref^2, with alpha = theta_pest - iota * zeta
    grad_alpha_grad_alpha = L_reference*L_reference*(grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z)

    ! this is grad alpha . grad psi_t in units of B_reference
    grad_alpha_grad_psi = (grad_alpha_X*grad_psi_X + grad_alpha_Y*grad_psi_Y + grad_alpha_Z*grad_psi_Z) / B_reference

    ! this is grad psi_t . grad psi_t in units of (B_reference*L_reference)^2
    grad_psi_grad_psi = (grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z) &
         / (L_reference * L_reference * B_reference * B_reference)

    ! dx/dpsitor = sign_torflux/rhotor/Lref/Bref
    ! dy/dalpha = Lref*rhotor

    ! this is shat * grad x . grad y = shat * (grad psi_tor . grad alpha) * dx/dpsi_tor * dy/dalpha
    ! = sign_torflux * shat * (grad psi_tor . grad alpha) / Bref
!    gds21 = (grad_alpha_X * grad_psi_X + grad_alpha_Y * grad_psi_Y + grad_alpha_Z * grad_psi_Z) &
!         * sign_toroidal_flux * shat / B_reference

    ! this is shat^2 * | grad x | ^2 = shat^2 * |grad psi_tor|^2 * (dx/dpsitor)^2
    ! = shat^2 * |grad psi_tor|^2 * / rhotor^2 / Lref^2 / Bref^2
!    gds22 = (grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z) &
!         * shat * shat / (L_reference * L_reference * B_reference * B_reference * normalized_toroidal_flux_used)

    ! this is (grad zeta . grad x_stella) / bmag^2 = (grad zeta . grad psitor) * dx/dpsitor / bmag^2
    ! = (grad zeta . grad psitor) * sign_torflux/rhotor/Lref/Bref/ bmag^2
    gradzeta_gradx = sign_toroidal_flux &
         * (grad_zeta_X * grad_psi_X + grad_zeta_Y * grad_psi_Y + grad_zeta_Z * grad_psi_Z) &
         / (L_reference * B_reference * sqrt(normalized_toroidal_flux_used) * bmag**2)

    ! this is (grad zeta . grad y_stella) / bmag^2 = (grad zeta . grad alpha) * dy/dalpha / bmag^2
    ! = (grad zeta . grad alpha) * Lref * rhotor / bmag^2
    gradzeta_grady = (grad_zeta_X * grad_alpha_X + grad_zeta_Y * grad_alpha_Y + grad_zeta_Z * grad_alpha_Z) &
         * L_reference * sqrt(normalized_toroidal_flux_used) / bmag**2

    ! this is (grad theta_pest . grad x_stella) / bmag^2
    gradtheta_gradx = (grad_theta_pest_X * grad_psi_X + grad_theta_pest_Y * grad_psi_Y &
         + grad_theta_pest_Z * grad_psi_Z) &
         / (L_reference * B_reference * sqrt(normalized_toroidal_flux_used) * bmag**2)

    ! this is (grad theta_pest . grad y_stella) / bmag^2
    gradtheta_grady = (grad_theta_pest_X * grad_alpha_X + grad_theta_pest_Y * grad_alpha_Y &
         + grad_theta_pest_Z * grad_alpha_Z) &
         * L_reference * sqrt(normalized_toroidal_flux_used) / bmag**2

    ! this is psitor/|psitor|*((grad y_stella . grad zeta)*(grad x_stella . grad y_stella) 
    ! - (grad x_stella . grad zeta)*|grad y_stella|^2) / (B/Bref)^2
    gds23 = sign_toroidal_flux*(gradzeta_grady * sign_toroidal_flux * grad_alpha_grad_psi &
         - gradzeta_gradx * grad_alpha_grad_alpha * normalized_toroidal_flux_used)
!    gds23 = sign_toroidal_flux*(gradzeta_grady * gds21/shat - gradzeta_gradx * gds2)

    ! this is psitor/|psitor| * ((grad y_stella . grad zeta) * |grad x_stella|^2 
    ! - (grad x_stella . grad zeta)*(grad x_stella . grad y_stella)) / (B/Bref)^2
    gds24 = (gradzeta_grady * grad_psi_grad_psi/normalized_toroidal_flux_used &
         - gradzeta_gradx * sign_toroidal_flux * grad_alpha_grad_psi) * sign_toroidal_flux
!    gds24 = (gradzeta_grady * gds22/shat**2 - gradzeta_gradx * gds21/shat) * sign_toroidal_flux

    ! this is ((grad y_stella . grad theta_pest)*(grad x_stella . grad y_stella)
    ! - (grad x_stella . grad theta_pest)*|grad y_stella|^2) / (B/Bref)^2
    gds25 = gradtheta_grady * sign_toroidal_flux * grad_alpha_grad_psi &
         - gradtheta_gradx * grad_alpha_grad_alpha * normalized_toroidal_flux_used
!    gds25 = gradtheta_grady * gds21/shat - gradtheta_gradx * gds2

    ! this is ((grad y_stella . grad theta_pest) * |grad x_stella|^2 
    ! - (grad x_stella . grad theta_pest)*(grad x_stella . grad y_stella)) / 2 / (B/Bref)^2
    gds26 = 0.5*(gradtheta_grady * grad_psi_grad_psi/normalized_toroidal_flux_used &
         - gradtheta_gradx * sign_toroidal_flux * grad_alpha_grad_psi)
!    gds26 = 0.5*(gradtheta_grady * gds22/shat**2 - gradtheta_gradx * gds21/shat)

    gbdrift_alpha = 2 * B_reference * L_reference * L_reference * B_cross_grad_B_dot_grad_alpha &
         / (B * B * B)
!    gbdrift = 2 * B_reference * L_reference * L_reference * sqrt_s * B_cross_grad_B_dot_grad_alpha &

    gbdrift0_psi = -edge_toroidal_flux_over_2pi &
         * (B_sub_theta_vmec * d_B_d_zeta - B_sub_zeta * d_B_d_theta_vmec) / sqrt_g &
         * 2 * shat / (B * B * B)
!    gbdrift0 = (B_sub_theta_vmec * d_B_d_zeta - B_sub_zeta * d_B_d_theta_vmec) / sqrt_g * edge_toroidal_flux_over_2pi &
!    gbdrift0 = abs(edge_toroidal_flux_over_2pi) &
!         * (B_sub_theta_vmec * d_B_d_zeta - B_sub_zeta * d_B_d_theta_vmec) / sqrt_g &
!         * 2 * shat / (B * B * B * sqrt_s)
    ! In the above 2-line expression for gbdrift0, the first line is \vec{B} \times \nabla B \cdot \nabla \psi.

    cvdrift_alpha = gbdrift_alpha + 2 * B_reference * L_reference * L_reference * mu_0 * d_pressure_d_s &
         * B_cross_grad_s_dot_grad_alpha / (B * B * B * B)
!    cvdrift = gbdrift + 2 * B_reference * L_reference * L_reference * sqrt_s * mu_0 * d_pressure_d_s &
!         * B_cross_grad_s_dot_grad_alpha / (B * B * B * B)

    cvdrift0_psi = gbdrift0_psi

    !*********************************************************************
    ! Free all arrays that were allocated.
    !*********************************************************************

    deallocate(B)
    deallocate(temp2D)
    deallocate(sqrt_g)
    deallocate(R)
    deallocate(d_B_d_theta_vmec)
    deallocate(d_B_d_zeta)
    deallocate(d_B_d_s)
    deallocate(d_R_d_theta_vmec)
    deallocate(d_R_d_zeta)
    deallocate(d_R_d_s)
    deallocate(d_Z_d_theta_vmec)
    deallocate(d_Z_d_zeta)
    deallocate(d_Z_d_s)
    deallocate(d_Lambda_d_theta_vmec)
    deallocate(d_Lambda_d_zeta)
    deallocate(d_Lambda_d_s)
    deallocate(B_sub_s)
    !deallocate(B_sub_theta_vmec)
    !deallocate(B_sub_zeta)
    deallocate(B_sup_theta_vmec)
    deallocate(B_sup_zeta)

    deallocate(d_B_d_s_mnc)
    deallocate(d_B_d_s_mns)
    deallocate(d_R_d_s_mnc)
    deallocate(d_R_d_s_mns)
    deallocate(d_Z_d_s_mnc)
    deallocate(d_Z_d_s_mns)
    deallocate(d_Lambda_d_s_mnc)
    deallocate(d_Lambda_d_s_mns)

    deallocate(d_X_d_s)
    deallocate(d_X_d_theta_vmec)
    deallocate(d_X_d_zeta)
    deallocate(d_Y_d_s)
    deallocate(d_Y_d_theta_vmec)
    deallocate(d_Y_d_zeta)

    deallocate(grad_s_X)
    deallocate(grad_s_Y)
    deallocate(grad_s_Z)
    deallocate(grad_theta_vmec_X)
    deallocate(grad_theta_vmec_Y)
    deallocate(grad_theta_vmec_Z)
    deallocate(grad_theta_pest_X)
    deallocate(grad_theta_pest_Y)
    deallocate(grad_theta_pest_Z)
    deallocate(grad_zeta_X)
    deallocate(grad_zeta_Y)
    deallocate(grad_zeta_Z)
    deallocate(grad_psi_X)
    deallocate(grad_psi_Y)
    deallocate(grad_psi_Z)
    deallocate(grad_alpha_X)
    deallocate(grad_alpha_Y)
    deallocate(grad_alpha_Z)

    deallocate(B_X)
    deallocate(B_Y)
    deallocate(B_Z)
    deallocate(grad_B_X)
    deallocate(grad_B_Y)
    deallocate(grad_B_Z)
    deallocate(B_cross_grad_B_dot_grad_alpha)
    deallocate(B_cross_grad_B_dot_grad_alpha_alternate)
    deallocate(B_cross_grad_s_dot_grad_alpha)
    deallocate(B_cross_grad_s_dot_grad_alpha_alternate)

    deallocate (gradzeta_grady)
    deallocate (gradzeta_gradx)
    deallocate (gradtheta_grady)
    deallocate (gradtheta_gradx)

    deallocate(normalized_toroidal_flux_full_grid)
    deallocate(normalized_toroidal_flux_half_grid)

    if (verbose) then
       write (*,*) "Leaving vmec_to_stella_geometry_interface."
       write (*,*)
    end if

    if (allocated(rmnc)) then
       deallocate (xm, xn)
       deallocate (xm_nyq, xn_nyq)
       deallocate (rmnc, lmns, zmns, bmnc, gmnc)
       deallocate (bsupumnc, bsupvmnc, bsubumnc, bsubvmnc, bsubsmns)
       deallocate (phi, phip, iotas, iotaf, presf)
    end if
    if (allocated(rmns)) then
       deallocate (rmns, lmnc, zmnc, bmns, gmns)
       deallocate (bsupumns, bsupvmns, bsubumns, bsubvmns, bsubsmnc)
    end if

  contains

    subroutine test_arrays(array1, array2, should_be_0, tolerance, name)
      ! This subroutine is used for verifying the geometry arrays.
      ! When should_be_0 = .true., the subroutine verifies that |array1| = 0 to within 
      !     an absolute tolerance specified by 'tolerance'. array2 is ignored in this case.
      ! When should_be_0 = .false., the subroutine verifies that array1 = array2 
      !     to within a relative tolerance specified by 'tolerance'.

      implicit none

      real, dimension(nalpha,-nzgrid:nzgrid) :: array1, array2
      real :: tolerance
      character(len=*) :: name
      logical :: should_be_0
      real :: max_value, max_difference

      if (should_be_0) then
         max_value = maxval(abs(array1))
!         if (verbose) print *,"  maxval(abs(",trim(name),")):",max_value,"(should be << 1.)"
         if (max_value > tolerance) then
            print *,"Error! ",trim(name)," should be 0, but instead it is:"
            do ialpha = 1,nalpha
               print *,array1(ialpha,:)
            end do
            stop
         end if
      else
         max_difference = maxval(abs(array1 - array2)) / maxval(abs(array1) + abs(array2))
!         if (verbose) print *,"  Relative difference between two methods for computing ",trim(name),":",max_difference,"(should be << 1.)"
         if (max_difference > tolerance) then
            print *,"Error! Two methods for computing ",trim(name)," disagree. Here comes method 1:"
            do ialpha = 1,nalpha
               print *,array1(ialpha,:)
            end do
            print *,"Here comes method 2:"
            do ialpha = 1,nalpha
               print *,array2(ialpha,:)
            end do
            print *,"Here comes the difference:"
            do ialpha = 1,nalpha
               print *,array1(ialpha,:) - array2(ialpha,:)
            end do
            stop
         end if
      end if

    end subroutine test_arrays

  end subroutine vmec_to_stella_geometry_interface

  subroutine get_root (a0, b0, root, converged)

    implicit none
    
    real, intent (in) :: a0, b0
    real, intent (out) :: root
    logical, intent (out) :: converged

    integer, parameter :: itmax_bracket = 10
    integer, parameter :: itmax_root = 10
    real, parameter :: tol = 1.0e-10
    integer :: it
    real :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,eps
    
    a = a0
    b = b0
    fa = fzero_residual(a)
    fb = fzero_residual(b)
    do it = 1, itmax_bracket
       eps = epsilon(a)
       if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
          write (*,*)
          write (*,*) 'in vmec_to_stella_geometry_interface, theta_min=', a, ' and theta_max=', b, ' do not bracket root.'
          write (*,*) 'f(theta_min)=', fa, 'and f(theta_max)=', fb, '.'
          a = a-0.3
          b = b+0.3
          write (*,*) 'Trying again with values ', a, ' and ', b, ' .'
          fa = fzero_residual(a)
          fb = fzero_residual(b)
       else
          exit
       end if
    end do

    c=b
    fc=fb
    do it = 1, itmax_root
       if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
          c=a
          fc=fa
          d=b-a
          e=d
       end if
       if (abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       end if
       tol1=2.0*eps*abs(b)+0.5*tol
       xm=0.5*(c-b)
       if (abs(xm) <= tol1 .or. fb == 0.0) then
          root = b
          converged = .true.
          exit
       end if
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa
          if (a==c) then
             p=2.0*xm*s
             q=1.0-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
             q=(q-1.0)*(r-1.0)*(s-1.0)
          end if
          if (p > 0.0) q=-q
          p=abs(p)
          if (2.0*p < min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          end if
       else
          d=xm
          e=d
       end if
       a=b
       fa=fb
       b=b+merge(d,sign(tol1,xm), abs(d) > tol1)
       fb = fzero_residual(b)
    end do

  end subroutine get_root

  function fzero_residual(theta_vmec_try)

    ! Note that lmns and lmnc use the non-Nyquist xm, xn, and mnmax.
    ! Also note that lmns and lmnc are on the HALF grid.

    implicit none
    
    real :: theta_vmec_try, fzero_residual
    real :: angle, sinangle, cosangle
    integer :: imn, which_surface

    ! residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
    fzero_residual = theta_vmec_try - theta_pest_target

    do imn = 1, mnmax
       angle = xm(imn)*theta_vmec_try - xn(imn)*zeta0
       sinangle = sin(angle)
       cosangle = cos(angle)
       do which_surface = 1,2
          fzero_residual = fzero_residual + vmec_radial_weight_half(which_surface) * lmns(imn,vmec_radial_index_half(which_surface)) * sinangle
          if (lasym) then
             fzero_residual = fzero_residual + vmec_radial_weight_half(which_surface) * lmnc(imn,vmec_radial_index_half(which_surface)) * cosangle
          end if
       end do
    end do

  end function fzero_residual

end module vmec_to_stella_geometry_interface_mod
