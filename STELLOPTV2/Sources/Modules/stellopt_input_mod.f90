!-----------------------------------------------------------------------
!     Module:        stellopt_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/24/2012
!     Description:   This module contains the STELLOPT input namelist and
!                    subroutine which initializes and reads the
!                    STELLOPT input namelist.
!-----------------------------------------------------------------------
      MODULE stellopt_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stellopt_runtime
      USE stellopt_vars
      USE stellopt_targets
      USE safe_open_mod, ONLY: safe_open
      USE diagno_runtime, ONLY: DIAGNO_VERSION
      USE vmec0, ONLY: version_
      USE vmec_input, ONLY: lasym_local => lasym
      USE vmec_params, ONLY: version_vmec=> version_
      USE mpi_params                                                    ! MPI
!DEC$ IF DEFINED (GENE)
      !USE par_other, ONLY: svn_gene => svn_rev, release_gene => release !OLD SVN Version
      USE par_other, ONLY: svn_gene => git_master, release_gene => git_branch
!DEC$ ENDIF        
!DEC$ IF DEFINED (BEAMS3D_OPT)
      USE beams3d_runtime, ONLY: BEAMS3D_VERSION
!DEC$ ENDIF        
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF        
!-----------------------------------------------------------------------
!     Input Namelists
!         &optimum
!            nfunc_max          Maximum number of function evaluations
!                               1: Single iteration
!                               2: Single iteration with Jacobian eval (if equil_type=='LMDIF')
!            equil_type         Equilibrium Code:
!                                  'VMEC2000' (default)
!                                  'SPEC'
!            opt_type           Optimizer Type
!                                  'LMDIF'    (default)
!                                  'GADE'
!            ftol               Determines tollerance in sum of squares (LMDIF)
!            xtol               Determines relative errror in approximate solution (LMDIF)
!            gtol               Determines orthagonality of solution (LMDIF)
!            epsfcn             Determines jacobian step length (LMDIF)
!                               Determines crossover factor (GADE)
!                               Determines local attractor scaling, global = 1 (PSO)
!            factor             Determines initial step size (LMDIF)
!                               Determines mutation scaling factor (GADE)
!                               Determines Velocity scaling factor (PSO)
!            mode               Determines if scaling is automatic (1) or user(0) (LMDIF)
!                               Determines strategy (GADE) 
!                               Detremines number of divisions if > numprocs (MAP)
!            cr_strategy        Crossover strategy (GADE, 0:exponential, 1: binomial)
!            npopulation        Size of population (defaults to nproc if -1 or not set)
!            lkeep_mins         Keep minimum files.
!            lphiedge_opt       Logical to control PHIEDGE variation
!            lcurtor_opt        Logical to control CURTOR variation
!            lbcrit_opt         Logical to control BCRIT variation
!            lpscale_opt        Logical to control PRES_SCALE variation
!            lextcur_opt        Logical array to control EXTCUR varaition
!            laphi_opt          Logical array to control APHI variation
!            lam_opt            Logical array to control AM variation
!            lac_opt            Logical array to control AC variation
!            lai_opt            Logical array to control AI variation
!            lah_opt            Logical array to control AH variation
!            lat_opt            Logical array to control AT variation
!            lam_s_opt          Logical array to control AM_AUX_S variation
!            lam_f_opt          Logical array to control AM_AUX_F variation
!            lac_s_opt          Logical array to control AC_AUX_S variation
!            lac_f_opt          Logical array to control AC_AUX_F variation
!            lai_s_opt          Logical array to control AI_AUX_S variation
!            lai_f_opt          Logical array to control AI_AUX_F variation                         
!            lphi_s_opt         Logical array to control PHI_AUX_S variation (Electrostatic potential) 
!            lphi_f_opt         Logical array to control PHI_AUX_F variation (Electrostatic potential)
!            lah_f_opt          Logical array to control AH_AUX_F variation
!            lat_f_opt          Logical array to control AT_AUX_F variation
!            lbound_opt         Logical array to control Boudnary variation
!            lrho_opt           Logical array to control HB Boudnary variation
!            rho_exp            Integer controling value of HB Boundary exponent (default 2)
!            dphiedge_opt       Scale factor for PHIEDGE variation
!            dcurtor_opt        Scale factor for CURTOR variation
!            dbcrit_opt         Scale factor for BCRIT variation
!            dextcur_opt        Scale factor array for EXTCUR variation
!            daphi_opt          Scale factor array for APHI variation
!            dam_opt            Scale factor array for AM variation
!            dac_opt            Scale factor array for AC variation
!            dai_opt            Scale factor array for AI variation
!            dah_opt            Scale factor array for AH variation
!            dat_opt            Scale factor array for AT variation
!            dam_s_opt          Scale factor array for AM_AUX_S variation
!            dam_f_opt          Scale factor array for AM_AUX_F variation
!            dac_s_opt          Scale factor array for AC_AUX_S variation
!            dac_f_opt          Scale factor array for AC_AUX_F variation
!            dai_s_opt          Scale factor array for AI_AUX_S variation
!            dai_f_opt          Scale factor array for AI_AUX_F variation
!            dphi_s_opt         Scale factor array for PHI_AUX_S variation
!            dphi_f_opt         Scale factor array for PHI_AUX_F variation
!            dbound_opt         Scale factor array for Boundarys variation
!            mboz               Boozer Tranfrom Poloidal Fourier Content
!            nboz               Boozer Tranform Toroidal Fourier Content
!            target_phiedge     Target for PHIEDGE
!            target_RBtor       Target for RBtor
!            target_R0          Target for R axis (phi=0 plane)
!            target_Z0          Target for R axis (phi=0 plane)
!            target_curtor      Target for CURTOR
!            target_curtor_max  Target for maximum toroidal current
!            target_volume      Target for volume
!            target_beta        Target for Beta
!            target_betat       Target for Beta Toroidal
!            target_betap       Target for Beta Poloidal
!            target_wp          Target for stored energy
!            target_apsect_max  Target upper limit on aspect ratio
!            width_apsec_max    Width of TANH function which defines wall in parameter space
!            target_aspect      Target for aspect ratio
!            target_gradp_max   Target upper limit on grad(p)
!            width_gradp_max    Width of TANH function which defines wall in parameter space
!            target_pmin        Target lower limit on min(p)
!            width_pmin         Width of TANH function which defines wall in parameter space
!            target_curvature   Target Edge Curvature Kertosis
!            target_press       Target array for pressure
!            target_ne_line     Target array for line integrated electron density
!            r0_ne_line         R starting point for line integrated measurement (R1: endpoint) [m]
!            phi0_ne_line       PHI starting point for line integrated measurement (PHI1: endpoint) [rad]
!            z0_ne_line         Z starting point for line integrated measurement (Z1: endpoint) [m]
!            target_faraday     Target array for line integrated Faraday rotation.
!            r0_faraday         R starting point for line integrated measurement (R1: endpoint) [m]
!            phi0_faraday       PHI starting point for line integrated measurement (PHI1: endpoint) [rad]
!            z0_faraday         Z starting point for line integrated measurement (Z1: endpoint) [m]
!            target_sxr         Target array for soft X-Ray Diagnostic.
!            r0_sxr             R starting point for soft X-Ray Chord (R1: endpoint) [m]
!            phi0_sxr           PHI starting point for soft X-Ray Chord (PHI1: endpoint) [rad]
!            z0_sxr             Z starting point for soft X-Ray Chord (Z1: endpoint) [m]
!            target_te          Target array for electron temperature
!            r_te               R electron temperature location array
!            z_te               Z electron temperature location array
!            phi_te             PHI electron temperature location array
!            s_te               s electron temperature location array
!            target_ne          Target array for electron density
!            r_ne               R electron density location array
!            z_ne               Z electron density location array
!            phi_ne             PHI electron density location array
!            s_ne               s electron denstiy location array
!            target_ti          Target array for ion temperature
!            r_ti               R ion temperature location array
!            z_ti               Z ion temperature location array
!            phi_ti             PHI ion temperature location array
!            s_ti               s ion temperature location array
!            target_vphi        Target array for toroidal rotation
!            r_vphi             R toroidal rotation location array
!            z_vphi             Z toroidal rotation location array
!            phi_vphi           PHI toroidal rotation location array
!            s_vphi             s toroidal rotation location array
!            r_iota             R Rotational Transform location array
!            z_iota             Z Rotational Transform location array
!            phi_iota           PHI Rotational Transform location array
!            s_iota             s Rotational Transform location array
!            target_mse         Target arry of MSE gammas
!            r_mse              R MSE datapoint array
!            phi_mse            PHI MSE datapoint array
!            z_mse              Z MSE datapoint array
!            s_mse              s MSE datapoint array
!            a1_mse             A1 MSE coefficient (A1-A6 defined this way)
!            vac_mse            Vacuum MSE signal
!            lmse_extcur        Logical Controlling Which currents are factored out
!            target_bprobe      Magnetic Field Measurement
!            target_fluxloop    Fluxloop Measurement
!            target_segrog      Rogowski Coil Measurement
!            magdiag_coil       Coils file for magnetic diagnostics
!            target_vessel      Minimum distance between equilibrium and vessel
!            vessel_string      Vessel file
!            target_balloon     Target array for ballooning stability calculation
!            balloon_theta      Array indicating values of theta for ballooning calculation
!            balloon_zeta       Array indicating values of zeta for ballooning calculation
!            target_bootstrap   Target array for boostrap current <J*B> (difference between VMEC and bootstrap>
!            target_neo         Target array for neoclassical transport (should be zero)
!            target_Jstar       Target array for dJ*/du (should = 0)
!            NumJstar           Division of trapped mu's
!            target_helicity    Target arry for Boozer Spectrum Helicity (=0)
!            helicity           Complex number specifying target helicity
!            target_resjac      Target Resonant Jacobian modes in Boozer Spectrum (iota_res = nfp*n/m)
!            xm_resjac          Array of corresponding m values to target
!            xn_resjac          Array of corresponding n values to target
!            target_separatrix  Target array minimum distance between plasma and separatrix
!            r_separatrix       Array of (ntheta,nzeta) radial separatrix values [m]
!            z_separatrix       Array of (ntheta,nzeta) vertical separatrix values [m] 
!            phi_separatrix     Array of (ntheta,nzeta) toroidal angle separatrix values [rad]
!            target_txport      Array of target values for turbulent transport optimization
!            s_txport           Array of normalized toroidal fluxes to calculated turbulent transport
!            target_dkes        Array of target values for drift kinetic optimization
!            nu_dkes            Array of nu values for DKES optimization [m] (Efield from phi array)
!            target_limter      Target array minimum distance between plasma and limiter surface
!            r_limiter          Array of (ntheta,nzeta) radial limiter values [m]
!            z_limiter          Array of (ntheta,nzeta) vertical limiter values [m] 
!            phi_limiter        Array of (ntheta,nzeta) toroidal angle limiter values [rad]
!            txport_proxy       String of proxy function name.
!
!   NOTE:  All varaibles which start with target have an similar
!          varaible starting with sigma which defines the error bars.
!-----------------------------------------------------------------------
      NAMELIST /optimum/ nfunc_max, equil_type, opt_type,&
                         ftol, xtol, gtol, epsfcn, factor, refit_param, &
                         cr_strategy, mode, lkeep_mins, lrefit,&
                         npopulation, noptimizers, &
                         lphiedge_opt, lcurtor_opt, lbcrit_opt, &
                         lpscale_opt, lmix_ece_opt,&
                         lextcur_opt, laphi_opt, lam_opt, lac_opt, &
                         lai_opt, lah_opt, lat_opt, lam_s_opt, &
                         lam_f_opt, lac_s_opt, lac_f_opt, lai_s_opt, &
                         lai_f_opt, lne_f_opt, lte_f_opt, lti_f_opt,&
                         lbeamj_f_opt, lbootj_f_opt, lzeff_f_opt, &
                         lth_f_opt, lphi_s_opt, lphi_f_opt, &
                         lrho_opt, ldeltamn_opt, lbound_opt, laxis_opt, lmode_opt, &
                         lne_opt, lte_opt, lti_opt, lth_opt, lzeff_opt, &
                         lah_f_opt, lat_f_opt, lcoil_spline, &
                         dphiedge_opt, dcurtor_opt, dbcrit_opt, &
                         dpscale_opt, dmix_ece_opt,&
                         dextcur_opt, daphi_opt, dam_opt, dac_opt, &
                         dai_opt, dah_opt, dat_opt, dam_s_opt, &
                         dam_f_opt, dac_s_opt, dac_f_opt, dai_s_opt, &
                         dai_f_opt, dne_f_opt, dte_f_opt, dti_f_opt,&
                         dth_f_opt, dphi_s_opt, dphi_f_opt,&
                         dbeamj_f_opt, dbootj_f_opt, dzeff_f_opt, dbound_opt, &
                         drho_opt, ddeltamn_opt, &
                         dne_opt, dte_opt, dti_opt, dth_opt, dzeff_opt, &
                         dah_f_opt, dat_f_opt, daxis_opt, &
                         dmix_ece_opt, dcoil_spline, &
                         ne_aux_s, te_aux_s, ti_aux_s, th_aux_s, phi_aux_s,&
                         beamj_aux_s, bootj_aux_s, zeff_aux_s, &
                         ne_aux_f, te_aux_f, ti_aux_f, th_aux_f, phi_aux_f,&
                         beamj_aux_f, bootj_aux_f, zeff_aux_f, &
                         ne_opt, te_opt, ti_opt, th_opt, zeff_opt, &
                         ne_type, te_type, ti_type, th_type, &
                         beamj_type, bootj_type, zeff_type, &
                         ne_min, te_min, ti_min, th_min, beamj_f_min, &
                         bootj_f_min, zeff_min, zeff_f_min, &
                         ne_max, te_max, ti_max, th_max, beamj_f_max, &
                         bootj_f_max, zeff_max, zeff_f_max, &
                         ah_f_min, at_f_min, &
                         ah_f_max, at_f_max, &
                         raxis_min, raxis_max, &
                         zaxis_min, zaxis_max, &
                         rbc_min, rbc_max, zbs_min, zbs_max, &
                         rbs_min, rbs_max, zbc_min, zbc_max, &
                         mboz, nboz, rho_exp, &
                         coil_splinesx,coil_splinesy,coil_splinesz,&
                         coil_splinefx,coil_splinefy,coil_splinefz,&
                         coil_splinefx_min,coil_splinefy_min,coil_splinefz_min,&
                         coil_splinefx_max,coil_splinefy_max,coil_splinefz_max,&
                         target_phiedge, sigma_phiedge, &
                         target_rbtor, sigma_rbtor, &
                         target_r0, sigma_r0, target_z0, sigma_z0, &
                         target_curtor, sigma_curtor, &
                         target_curtor_max, sigma_curtor_max, &
                         target_volume, sigma_volume, &
                         target_beta, sigma_beta, &
                         target_betapol, sigma_betapol, &
                         target_betator, sigma_betator, &
                         target_wp, sigma_wp, &
                         target_aspect, sigma_aspect, &
                         target_extcur, sigma_extcur, &
                         target_aspect_max, sigma_aspect_max, width_aspect_max, &
                         target_gradp_max, sigma_gradp_max, width_gradp_max, &
                         target_pmin, sigma_pmin, width_pmin, &
                         target_curvature, sigma_curvature, &
                         target_kappa, sigma_kappa, phi_kappa, &
                         target_kappa_box, sigma_kappa_box, phi_kappa_box, &
                         target_kappa_avg, sigma_kappa_avg, &
                         target_magwell, sigma_magwell, &
                         target_press, sigma_press, r_press, z_press, phi_press, s_press,&
                         target_te, sigma_te, r_te, z_te, phi_te, s_te,&
                         target_ne, sigma_ne, r_ne, z_ne, phi_ne, s_ne,&
                         target_ne_line, sigma_ne_line, r0_ne_line, phi0_ne_line, z0_ne_line,&
                         r1_ne_line, phi1_ne_line, z1_ne_line, &
                         target_te_line, sigma_te_line, r0_te_line, phi0_te_line, z0_te_line,&
                         r1_te_line, phi1_te_line, z1_te_line, cutoff_te_line, &
                         target_ti_line, sigma_ti_line, r0_ti_line, phi0_ti_line, z0_ti_line,&
                         r1_ti_line, phi1_ti_line, z1_ti_line, &
                         target_ti, sigma_ti, r_ti, z_ti, phi_ti, s_ti,&
                         target_vphi, sigma_vphi, r_vphi, z_vphi, phi_vphi, s_vphi, qm_ratio,&
                         target_iota, sigma_iota, r_iota, z_iota, phi_iota, s_iota,&
                         target_vaciota, sigma_vaciota, r_vaciota, z_vaciota, phi_vaciota, s_vaciota,&
                         target_mse, sigma_mse, r_mse, z_mse, phi_mse, &
                         a1_mse, a2_mse, a3_mse, a4_mse, a5_mse, a6_mse, &
                         a7_mse, vac_mse, lmse_extcur, &
                         target_faraday, sigma_faraday, r0_faraday, phi0_faraday, z0_faraday,&
                         r1_faraday, phi1_faraday, z1_faraday, &
                         target_sxr, sigma_sxr, r0_sxr, phi0_sxr, z0_sxr,&
                         r1_sxr, phi1_sxr, z1_sxr, &
                         target_bprobe, target_fluxloop, target_segrog,&
                         sigma_bprobe, sigma_fluxloop, sigma_segrog, magdiag_coil,&
                         target_vessel, sigma_vessel, vessel_string, &
                         phiedge_min, phiedge_max, curtor_min, curtor_max, &
                         bcrit_min, bcrit_max, pscale_min, pscale_max, &
                         mix_ece_min, mix_ece_max, &
                         extcur_min, extcur_max, &
                         am_min, am_max, ac_min, ac_max, ai_min, ai_max, &
                         ah_min, ah_max, at_min, at_max, aphi_min, aphi_max, &
                         am_f_min, ac_f_min, ai_f_min, phi_f_min, ne_f_min, &
                         te_f_min, ti_f_min, th_f_min, &
                         am_f_max, ac_f_max, ai_f_max, phi_f_max, ne_f_max, &
                         te_f_max, ti_f_max, th_f_max, bound_min, bound_max, &
                         delta_min, delta_max, &
                         target_balloon, sigma_balloon, balloon_theta, balloon_zeta,&
                         target_bootstrap,sigma_bootstrap, target_neo, sigma_neo,&
                         target_Jstar, sigma_Jstar, NumJstar,&
                         target_helicity, sigma_helicity, helicity,&
                         target_helicity_old, sigma_helicity_old, &
                         target_resjac, sigma_resjac, xm_resjac, xn_resjac,&
                         target_separatrix, sigma_separatrix, &
                         r_separatrix, z_separatrix, phi_separatrix, &
                         target_limiter, sigma_limiter, &
                         r_limiter, z_limiter, phi_limiter, &
                         lglobal_txport, nz_txport, nalpha_txport, alpha_start_txport, alpha_end_txport, &
                         target_txport, sigma_txport, s_txport, txport_proxy,&
                         target_dkes, sigma_dkes, nu_dkes, &
                         target_jdotb,sigma_jdotb,target_bmin,sigma_bmin,&
                         target_bmax,sigma_bmax,target_jcurv,sigma_jcurv,&
                         target_orbit,sigma_orbit,nu_orbit,nv_orbit,&
                         mass_orbit,Z_orbit,vperp_orbit,&
                         np_orbit,vll_orbit,mu_orbit, target_coil_bnorm,&
                         sigma_coil_bnorm, nu_bnorm, nv_bnorm,&
                         target_ece,sigma_ece,freq_ece, mix_ece, vessel_ece, mirror_ece, &
                         antennaposition_ece, targetposition_ece, rbeam_ece, rfocus_ece, &
                         targettype_ece, antennatype_ece, nra_ece, nphi_ece, &
                         target_kink, sigma_kink,mlmnb_kink,mlmns_kink,ivac_kink,&
                         nj_kink, nk_kink, lssl_kink, lssd_kink, mmaxdf_kink, nmaxdf_kink, &
                         lregcoil_winding_surface_separation_opt, &
                         dregcoil_winding_surface_separation_opt, &
                         lregcoil_current_density_opt, &
                         dregcoil_current_density_opt, &
                         target_regcoil_winding_surface_separation, &
                         sigma_regcoil_winding_surface_separation, &
                         target_regcoil_bnorm, sigma_regcoil_bnorm, &
                         target_regcoil_chi2_b, sigma_regcoil_chi2_b, &
                         target_regcoil_current_density, sigma_regcoil_current_density, &
                         regcoil_winding_surface_separation, &
                         regcoil_current_density
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_stellopt_input:   Reads optimum namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_stellopt_input(filename, istat, ithread)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      INTEGER, INTENT(in) :: ithread
      LOGICAL :: lexist
      INTEGER :: i, iunit, local_master
      ! Initializations
      nfunc_max       = 5000
      opt_type        = 'LMDIF'
      equil_type      = 'VMEC2000'
      ftol            = 1.0D-06
      xtol            = 1.0D-06
      gtol            = 0.0
      epsfcn          = 1.0D-06
      mode            = 1       ! Default in case user forgets
      factor          = 100.
      cr_strategy     = 0
      npopulation     = -1
      noptimizers     = -1
      refit_param     = 0.75
      rho_exp         = 4
      lkeep_mins      = .FALSE.
      lrefit          = .FALSE.
      lphiedge_opt    = .FALSE.
      lcurtor_opt     = .FALSE.
      lpscale_opt     = .FALSE.
      lbcrit_opt      = .FALSE.
      lmix_ece_opt    = .FALSE.
      lextcur_opt(:)  = .FALSE.
      laphi_opt(:)    = .FALSE.
      lam_opt(:)      = .FALSE.
      lac_opt(:)      = .FALSE.
      lai_opt(:)      = .FALSE.
      lah_opt(:)      = .FALSE.
      lat_opt(:)      = .FALSE.
      lne_opt(:)      = .FALSE.
      lzeff_opt(:)    = .FALSE.
      lte_opt(:)      = .FALSE.
      lti_opt(:)      = .FALSE.
      lth_opt(:)      = .FALSE.
      lam_s_opt(:)    = .FALSE.
      lam_f_opt(:)    = .FALSE.
      lac_s_opt(:)    = .FALSE.
      lac_f_opt(:)    = .FALSE.
      lai_s_opt(:)    = .FALSE.
      lai_f_opt(:)    = .FALSE.
      lphi_s_opt(:)   = .FALSE.
      lphi_f_opt(:)   = .FALSE.
      lne_f_opt(:)    = .FALSE.
      lzeff_f_opt(:)  = .FALSE.
      lte_f_opt(:)    = .FALSE.
      lti_f_opt(:)    = .FALSE.
      lth_f_opt(:)    = .FALSE.
      lbeamj_f_opt(:) = .FALSE.
      lbootj_f_opt(:) = .FALSE.
      lah_f_opt(:)    = .FALSE.
      lat_f_opt(:)    = .FALSE.
      lbound_opt(:,:) = .FALSE.
      lrho_opt(:,:)   = .FALSE.
      ldeltamn_opt(:,:) = .FALSE.
      lmode_opt(:,:)  = .FALSE.
      laxis_opt(:)    = .FALSE.
      lcoil_spline(:,:) = .FALSE.
      dphiedge_opt    = -1.0
      dcurtor_opt     = -1.0
      dpscale_opt     = -1.0
      dbcrit_opt      = -1.0
      dmix_ece_opt    = -1.0
      dextcur_opt(:)  = -1.0
      daphi_opt(:)    = -1.0
      dam_opt(:)      = -1.0
      dac_opt(:)      = -1.0
      dai_opt(:)      = -1.0
      dah_opt(:)      = -1.0
      dat_opt(:)      = -1.0
      dne_opt(:)      = -1.0
      dzeff_opt(:)    = -1.0
      dte_opt(:)      = -1.0
      dti_opt(:)      = -1.0
      dth_opt(:)      = -1.0
      dam_s_opt(:)    = -1.0
      dam_f_opt(:)    = -1.0
      dac_s_opt(:)    = -1.0
      dac_f_opt(:)    = -1.0
      dai_s_opt(:)    = -1.0
      dai_f_opt(:)    = -1.0
      dphi_s_opt(:)   = -1.0
      dphi_f_opt(:)   = -1.0
      dne_f_opt(:)    = -1.0
      dzeff_f_opt(:)  = -1.0
      dte_f_opt(:)    = -1.0
      dti_f_opt(:)    = -1.0
      dth_f_opt(:)    = -1.0
      dbeamj_f_opt(:) = -1.0
      dbootj_f_opt(:) = -1.0
      dah_f_opt(:)    = -1.0
      dat_f_opt(:)    = -1.0
      daxis_opt(:)    = -1.0
      dbound_opt(:,:) = -1.0
      drho_opt(:,:) = -1.0
      ddeltamn_opt(:,:) = -1.0
      dcoil_spline(:,:) = -1.0
      lregcoil_winding_surface_separation_opt    = .FALSE.
      dregcoil_winding_surface_separation_opt    = -1.0
      lregcoil_current_density_opt    = .FALSE.
      dregcoil_current_density_opt    = -1.0
      IF (.not.ltriangulate) THEN  ! This is done because values may be set by trinagulate
         phiedge_min     = -bigno;  phiedge_max     = bigno
         curtor_min      = -bigno;  curtor_max      = bigno
         bcrit_min       = -bigno;  bcrit_max       = bigno
         pscale_min      = 0.0   ;  pscale_max      = bigno
         extcur_min      = -bigno;  extcur_max      = bigno
         aphi_min        = -bigno;  aphi_max        = bigno
         am_min          = -bigno;  am_max          = bigno
         ac_min          = -bigno;  ac_max          = bigno
         ai_min          = -bigno;  ai_max          = bigno
         ah_min          = -bigno;  ah_max          = bigno
         at_min          = -bigno;  at_max          = bigno
         am_f_min        = 0.0;     am_f_max        = bigno
         ac_f_min        = -bigno;  ac_f_max        = bigno
         ai_f_min        = -bigno;  ai_f_max        = bigno
         phi_f_min       = -bigno;  phi_f_max       = bigno
         ah_f_min        = -bigno;  ah_f_max        = bigno
         at_f_min        = -bigno;  at_f_max        = bigno
         raxis_min        = -bigno; raxis_max        = bigno
         zaxis_min        = -bigno; zaxis_max        = bigno
         rbc_min         = -bigno;  rbc_max         = bigno
         zbs_min         = -bigno;  zbs_max         = bigno
         bound_min       = -bigno;  bound_max       = bigno
         delta_min       = -bigno;  delta_max       = bigno
      END IF
      mix_ece_min     = 0.0   ;  mix_ece_max     = 1.0
      ne_min          = -bigno;  ne_max          = bigno
      zeff_min        = -bigno;  zeff_max        = bigno
      te_min          = -bigno;  te_max          = bigno
      ti_min          = -bigno;  ti_max          = bigno
      th_min          = -bigno;  th_max          = bigno
      ne_f_min        = 0.0;     ne_f_max        = bigno_ne
      zeff_f_min      = 0.0;     zeff_f_max      = bigno
      te_f_min        = 0.0;     te_f_max        = bigno
      ti_f_min        = 0.0;     ti_f_max        = bigno
      th_f_min        = 0.0;     th_f_max        = bigno
      beamj_f_min     = -bigno;  beamj_f_max     = bigno
      bootj_f_min     = -bigno;  bootj_f_max     = bigno
      coil_splinefx_min       = -bigno;  coil_splinefx_max       = bigno
      coil_splinefy_min       = -bigno;  coil_splinefy_max       = bigno
      coil_splinefz_min       = -bigno;  coil_splinefz_max       = bigno
      target_regcoil_winding_surface_separation = 0.0
      sigma_regcoil_winding_surface_separation = bigno
      regcoil_winding_surface_separation = 1.0
      regcoil_winding_surface_separation_min = 1.0e-3
      regcoil_winding_surface_separation_max = 2.0
      target_regcoil_current_density = 0.0
      sigma_regcoil_current_density = bigno
      regcoil_current_density = 8.0e6
      regcoil_current_density_min = 0.0
      regcoil_current_density_max = bigno
      ne_type         = 'akima_spline'
      zeff_type       = 'akima_spline'
      te_type         = 'akima_spline'
      ti_type         = 'akima_spline'
      th_type         = 'akima_spline'
      beamj_type      = 'power_series'
      bootj_type      = 'power_series'
      ne_opt(0:20)       = 0.0
      zeff_opt(0:20)     = 0.0
      te_opt(0:20)       = 0.0
      ti_opt(0:20)       = 0.0
      th_opt(0:20)       = 0.0
      ne_aux_s(:)     = -1.0
      ne_aux_s(1:5)   = (/0.0,0.25,0.50,0.75,1.0/)
      ne_aux_f(:)     = 0.0 ! Do this so profile norm doesn't get screwy
      ne_aux_f(1:5)   = 1.0 ! Do this so we can scale T to P
      zeff_aux_s(:)   = -1.0
      zeff_aux_s(1:5) = (/0.0,0.25,0.50,0.75,1.0/)
      zeff_aux_f(:)   = 0.0 
      zeff_aux_f(1:5) = 1.0 
      te_aux_s(:)     = -1.0
      te_aux_f(:)     = 0.0
      ti_aux_s(:)     = -1.0
      ti_aux_f(:)     = 0.0
      th_aux_s(:)     = -1.0
      th_aux_f(:)     = 0.0 ! Probably need to recast th as ph later
      phi_aux_s(:)    = -1.0
      phi_aux_f(:)    = 0.0
      beamj_aux_s(:)  = -1.0
      beamj_aux_s(1:5) = (/0.0,0.25,0.50,0.75,1.0/)
      beamj_aux_f(:)  = 0.0
      bootj_aux_s(:)  = -1.0
      bootj_aux_s(1:5) = (/0.0,0.25,0.50,0.75,1.0/)
      bootj_aux_f(:)  = 0.0
      coil_splinesx(:,:) = -1
      coil_splinesy(:,:) = -1
      coil_splinesz(:,:) = -1
      coil_splinefx(:,:) = 0
      coil_splinefy(:,:) = 0
      coil_splinefz(:,:) = 0
      mboz            = 64
      nboz            = 64
      target_phiedge  = 0.0
      sigma_phiedge   = bigno
      target_rbtor    = 0.0
      sigma_rbtor     = bigno
      target_r0       = 0.0
      sigma_r0        = bigno
      target_z0       = 0.0
      sigma_z0        = bigno
      target_curtor   = 0.0
      sigma_curtor    = bigno
      target_curtor_max = 0.0
      sigma_curtor_max  = bigno
      target_volume   = 0.0
      sigma_volume    = bigno
      target_beta     = 0.0
      sigma_beta      = bigno
      target_betapol  = 0.0
      sigma_betapol   = bigno
      target_betator   = 0.0
      sigma_betator    = bigno
      target_wp        = 0.0
      sigma_wp         = bigno
      target_aspect    = 0.0
      sigma_aspect     = bigno
      target_aspect_max= 0.0
      sigma_aspect_max = bigno
      width_aspect_max = 0.5 ! Ideally we'd want it to be something like N*EPSFCN
      target_gradp_max = 0.0
      sigma_gradp_max  = bigno
      width_gradp_max  = 0.5
      target_pmin      = 0.0
      sigma_pmin       = bigno
      width_pmin       = 0.5
      target_curvature = 0.0
      sigma_curvature  = bigno
      target_magwell   = 0.0
      sigma_magwell    = bigno
      target_kappa     = 0.0
      sigma_kappa      = bigno
      phi_kappa        = 0.0
      target_kappa_box = 0.0
      sigma_kappa_box  = bigno
      phi_kappa_box    = 0.0
      target_kappa_avg = 0.0
      sigma_kappa_avg  = bigno
      target_kink(:)  = 0.0
      sigma_kink(:)   = bigno
      mlmnb_kink      = 264
      mlmns_kink(:)   = 91
      ivac_kink       = 24
      nj_kink         = 512
      nk_kink         = 512
      mmaxdf_kink     = 127
      nmaxdf_kink     = 31
      lssl_kink(:)    = 4096
      lssd_kink(:)    = 4096
      target_vessel   = 0.0
      sigma_vessel    = bigno
      vessel_string   = ''
      target_te(:)    = 0.0
      sigma_te(:)     = bigno
      target_extcur   = 0.0
      sigma_extcur    = bigno
      norm_press      = 1.0
      r_press(:)         = 0.0
      z_press(:)         = 0.0
      phi_press(:)       = 0.0
      s_press(:)         = -1.0
      target_press(:)    = 0.0
      sigma_press(:)     = bigno
      r_te(:)         = 0.0
      z_te(:)         = 0.0
      phi_te(:)       = 0.0
      s_te(:)         = -1.0
      target_ne(:)    = 0.0
      sigma_ne(:)     = bigno_ne
      r_ne(:)         = 0.0
      z_ne(:)         = 0.0
      phi_ne(:)       = 0.0
      s_ne(:)         = -1.0
      target_ne_line(:) = 0.0
      sigma_ne_line(:) = bigno_ne
      r0_ne_line(:)   = 0.0
      phi0_ne_line(:) = 0.0
      z0_ne_line(:)   = 0.0
      r1_ne_line(:)   = 0.0
      phi1_ne_line(:) = 0.0
      z1_ne_line(:)   = 0.0
      target_te_line(:) = 0.0
      sigma_te_line(:) = bigno
      r0_te_line(:)   = 0.0
      phi0_te_line(:) = 0.0
      z0_te_line(:)   = 0.0
      r1_te_line(:)   = 0.0
      phi1_te_line(:) = 0.0
      z1_te_line(:)   = 0.0
      cutoff_te_line  = 3000.0
      target_ti_line(:) = 0.0
      sigma_ti_line(:) = bigno
      r0_ti_line(:)   = 0.0
      phi0_ti_line(:) = 0.0
      z0_ti_line(:)   = 0.0
      r1_ti_line(:)   = 0.0
      phi1_ti_line(:) = 0.0
      z1_ti_line(:)   = 0.0
      target_ti(:)    = 0.0
      sigma_ti(:)     = bigno
      r_ti(:)         = 0.0
      z_ti(:)         = 0.0
      phi_ti(:)       = 0.0
      s_ti(:)         = -1.0
      target_vphi(:)  = 0.0
      sigma_vphi(:)   = bigno
      r_vphi(:)       = 0.0
      z_vphi(:)       = 0.0
      phi_vphi(:)     = 0.0
      s_vphi(:)       = -1.0
      qm_ratio        = 9.57883058054E+07
      target_iota(:)  = 0.0
      sigma_iota(:)   = bigno
      r_iota(:)       = 0.0
      z_iota(:)       = 0.0
      phi_iota(:)     = 0.0
      s_iota(:)       = -1.0
      target_vaciota(:)  = 0.0
      sigma_vaciota(:)   = bigno
      r_vaciota(:)       = 0.0
      z_vaciota(:)       = 0.0
      phi_vaciota(:)     = 0.0
      s_vaciota(:)       = -1.0
      target_mse(:)   = 0.0
      sigma_mse(:)    = bigno
      r_mse(:)        = 0.0
      z_mse(:)        = 0.0
      phi_mse(:)      = 0.0
      s_mse(:)        = -1.0
      a1_mse(:)       = 0.0
      a2_mse(:)       = 0.0
      a3_mse(:)       = 0.0
      a4_mse(:)       = 0.0
      a5_mse(:)       = 0.0
      a6_mse(:)       = 0.0
      a7_mse(:)       = 0.0
      vac_mse(:)      = 0.0
      lmse_extcur(:)  = .false.
      target_faraday(:) = 0.0
      sigma_faraday(:) = bigno_ne
      r0_faraday(:)   = 0.0
      phi0_faraday(:) = 0.0
      z0_faraday(:)   = 0.0
      r1_faraday(:)   = 0.0
      phi1_faraday(:) = 0.0
      z1_faraday(:)   = 0.0
      target_sxr(:) = 0.0
      sigma_sxr(:) = bigno_ne
      r0_sxr(:)   = 0.0
      phi0_sxr(:) = 0.0
      z0_sxr(:)   = 0.0
      r1_sxr(:)   = 0.0
      phi1_sxr(:) = 0.0
      z1_sxr(:)   = 0.0
      target_ece(:,:) = 0
      sigma_ece(:,:) = bigno
      mix_ece        = 0.8
      vessel_ece     = ''
      mirror_ece     = ''
      targettype_ece  = 'cyl'
      antennatype_ece = 'cyl'
      antennaposition_ece = 0
      targetposition_ece = 0
      rbeam_ece = 0
      rfocus_ece = 0
      nra_ece = 0
      nphi_ece = 8
      target_bprobe(:)= 0.0
      sigma_bprobe(:) = bigno
      target_fluxloop(:) = 0.0
      sigma_fluxloop(:) = bigno
      target_segrog(:)= 0.0
      sigma_segrog(:) = bigno
      magdiag_coil = ''
      target_balloon(:)= 0.0
      sigma_balloon(:) = bigno
      balloon_theta(:)= -1.0
      balloon_zeta(:) = -1.0
      target_bootstrap(:) = 0.0
      sigma_bootstrap(:) = bigno
      target_neo(:)   = 0.0
      sigma_neo(:)    = bigno
      target_Jstar(:) = 0.0
      sigma_Jstar(:)  = bigno
      NumJstar        = 4
      target_helicity(:) = 0.0
      sigma_helicity(:)  = bigno
      helicity           = CMPLX(0.0,0.0)
      target_helicity_old(:) = 0.0
      sigma_helicity_old(:)  = bigno
      target_resjac(:)  = 0.0
      sigma_resjac(:)   = bigno
      xn_resjac(:)      = 0
      xm_resjac(:)      = -99
      target_separatrix = 0.0
      sigma_separatrix  = bigno
      r_separatrix      = 0.0
      z_separatrix      = 0.0
      phi_separatrix    = 0.0
      lglobal_txport    = .FALSE.
      nz_txport         = 128
      nalpha_txport     = 1
      alpha_start_txport= 0.0
      alpha_end_txport  = pi2/2
      target_txport     = 0.0
      sigma_txport      = bigno
      s_txport          = -1.0
      txport_proxy      = 'prox1f'
      target_orbit      = 0.0
      sigma_orbit       = bigno
      mass_orbit        = 6.64465675E-27 ! Default to He4
      Z_orbit           = 2.0
      nu_orbit          = 0
      nv_orbit          = 0
      np_orbit          = 0
      vll_orbit         = 0
      mu_orbit          = 0
      vperp_orbit       = 0
      target_dkes       = 0.0
      sigma_dkes        = bigno
      nu_dkes           = 0.01
      target_jdotb      = 0.0
      sigma_jdotb       = bigno
      target_jcurv      = 0.0
      sigma_jcurv       = bigno
      target_bmin       = 0.0
      sigma_bmin        = bigno
      target_bmax       = 0.0
      sigma_bmax        = bigno
      target_limiter    = 0.0
      sigma_limiter     = bigno
      r_limiter         = 0.0
      z_limiter         = 0.0
      phi_limiter       = 0.0
      target_coil_bnorm = 0.0
      sigma_coil_bnorm  = bigno
      nu_bnorm          = 256
      nv_bnorm          = 64
      target_regcoil_bnorm = 0.0
      sigma_regcoil_bnorm  = bigno
      target_regcoil_chi2_b = 0.0
      sigma_regcoil_chi2_b  = bigno
      target_regcoil_current_density = 8.0e6
      sigma_regcoil_current_density  = bigno
      ! Read name list
      lexist            = .false.
      istat=0
      iunit=12
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) CALL handle_err(FILE_EXIST_ERR,TRIM(filename),istat)
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat /= 0) CALL handle_err(FILE_OPEN_ERR,TRIM(filename),istat)
      READ(iunit,NML=optimum,IOSTAT=istat)
      IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'OPTIMUM in: '//TRIM(filename),istat)
      CALL FLUSH(iunit)
      CLOSE(iunit)
      ! Fix String vars
      equil_type=TRIM(equil_type)
      equil_type=ADJUSTL(equil_type)
      opt_type  =TRIM(opt_type)
      opt_type  = ADJUSTL(opt_type)
      magdiag_coil = TRIM(magdiag_coil)
      magdiag_coil = ADJUSTL(magdiag_coil)
      ne_type = ADJUSTL(ne_type)
      te_type = ADJUSTL(te_type)
      ti_type = ADJUSTL(ti_type)
      th_type = ADJUSTL(th_type)
      beamj_type = ADJUSTL(beamj_type)
      bootj_type = ADJUSTL(bootj_type)
      ! Coil Optimization
      IF (ANY(ANY(lcoil_spline,2),1)) lcoil_geom = .true.
      ! If fixed boundary optimization or mapping turn off restart
      IF (ANY(ANY(lbound_opt,2),1) .or. opt_type=='map') lno_restart = .true.
      ! Test for proper normalization on ne profile
      IF (ANY(ne_opt > 1.0E+15) .or. ANY(ne_aux_f > 1.0E+15)) THEN
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Electron density profile is set too large.  The'
            WRITE(6,*) '  NE_OPT and NE_AUX_F arrays should be normalized to'
            WRITE(6,*) '  1.0E+18 [m^-3].'
            WRITE(6,*) '  '
            WRITE(6,*) '  This DOES NOT apply to the TARGET_NE or SIGMA_NE'
            WRITE(6,*) '  arrays which are normalized to 1.0E+00 [m^-3].'
            WRITE(6,*) '  '
            WRITE(6,*) '  STELLOPT will now renormalize.'
         END IF
         ne_opt = ne_opt / ne_norm
         ne_aux_f = ne_aux_f / ne_norm
      END IF
      ! Print code messages
      CALL tolower(equil_type)
      IF ((myid == master) .and. (TRIM(equil_type(1:4)) == 'vmec') ) THEN
         WRITE(6,*)        " Equilibrium calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========   Parallel Variational Moments Equilibrium Code (v "//TRIM(version_vmec)//")      ========="
         WRITE(6,"(2X,A)") "=========       Variational Moments Equilibrium Code (v "//TRIM(version_vmec)//")           ========="
         WRITE(6,"(2X,A)") "=========                (S. Hirshman, J. Whitson)                      ========="
         WRITE(6,"(2X,A)") "=========         http://vmecwiki.pppl.wikispaces.net/VMEC              ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
      !IF ((myid == master) .and. (TRIM(equil_type(1:8)) == 'paravmec') ) THEN
      !   WRITE(6,*)        " Equilibrium calculation provided by: "
      !   WRITE(6,"(2X,A)") "================================================================================="
      !   WRITE(6,"(2X,A)") "=========   Parallel Variational Moments Equilibrium Code (v "//TRIM(version_vmec)//")      ========="
      !   WRITE(6,"(2X,A)") "=========                (S. Hirshman, J. Whitson)                      ========="
      !   WRITE(6,"(2X,A)") "=========         http://vmecwiki.pppl.wikispaces.net/VMEC              ========="
      !   WRITE(6,"(2X,A)") "================================================================================="
      !   WRITE(6,*)        "    "
      !END IF
!DEC$ IF DEFINED (BEAMS3D_OPT)
      IF (myid == master .and. ANY(sigma_orbit < bigno) ) THEN
         WRITE(6,*)               " Energetic Particle calculation provided by: "
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                      BEAMS3D (v",BEAMS3D_VERSION,")                         ========="
         WRITE(6,"(2X,A)")        "=========                  (M. McMillan, S. Lazerson)                   ========="
         WRITE(6,"(2X,A)")        "=========                       lazerson@pppl.gov                       ========="
         WRITE(6,"(2X,A)")        "=========          http://vmecwiki.pppl.wikispaces.net/BEAMS3D          ========="
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_orbit < bigno)) THEN
         sigma_orbit = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the BEAMS3D code.   '
            WRITE(6,*) '  Optimization of particle orbits not possible.'
            WRITE(6,*) '  Disabling energetic particle targets.'
         END IF
      END IF
!DEC$ ENDIF
      IF (myid == master .and. ANY(sigma_bootstrap < bigno)) THEN
         WRITE(6,*)        " Bootstrap calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,A,A)") "=========                    BOOTSJ (v",version_,")                             ========="
         WRITE(6,"(2X,A)") "=========            (J. Tolliver, K. Shaing, P. Moroz)                 ========="
         WRITE(6,"(2X,A)") "=========        http://vmecwiki.pppl.wikispaces.net/BOOTSJ             ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
      IF (myid == master .and. ANY(sigma_balloon < bigno)) THEN
         WRITE(6,*)        " Ballooning stability calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                    COBRAVMEC (v",4.10,")                         ========="
         WRITE(6,"(2X,A)") "=========                 (R. Sanchez, S. Hirshman)                     ========="
         WRITE(6,"(2X,A)") "=========                   raul.sanchez@uc3m.es                        ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ IF DEFINED (TERPSICHORE)
      IF (myid == master .and. ANY(sigma_kink < bigno)) THEN
         WRITE(6,*)        " Kink stability calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                    TERPSICHORE (v2016)                        ========="
         WRITE(6,"(2X,A)") "=========    (D. V. ANDERSON, W. A. COOPER, R. GRUBER AND U. SCHWENN)   ========="
         WRITE(6,"(2X,A)") "=========                   wilfred.cooper@epfl.ch                      ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_kink < bigno)) THEN
         sigma_kink(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the TERPSICHORE code.   '
            WRITE(6,*) '  Optimization of kink stability not possible.'
            WRITE(6,*) '  Disabling kink stability targets.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (TRAVIS)
      IF (myid == master .and. ANY(sigma_ece < bigno)) THEN
         WRITE(6,*)        " ECE Radiation calculation provided by: "
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                            TRAVIS                             ========="
         CALL printversion_f77
         WRITE(6,"(2X,A)")        "=========                    (N. Marushchenko)                          ========="
         WRITE(6,"(2X,A)")        "=========              nikolai.marushchenko@ipp.mpg.de                  ========="
         WRITE(6,"(2X,A)")        "================================================================================="
         WRITE(6,*)               "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_ece < bigno)) THEN
         sigma_ece = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the TRAVIS code.   '
            WRITE(6,*) '  Optimization of ECE Radiation not possible.'
            WRITE(6,*) '  Disabling ECE Radiation targets.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (COILOPTPP)
      IF (myid == master .and. (sigma_coil_bnorm < bigno)) THEN
         WRITE(6,*)        " Stellarator Coil Optimization provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========                            COILOPT++                          ========="
         WRITE(6,"(2X,A)") "=========                    (J. Breslau, S. Lazerson)                  ========="
         WRITE(6,"(2X,A)") "=========                        jbreslau@pppl.gov                      ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (sigma_coil_bnorm < bigno) THEN
         sigma_coil_bnorm = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Coil optimization with the COILOPT++'
            WRITE(6,*) '  code has been disabled.  Coil optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (REGCOIL)
      IF (myid == master .and. ((sigma_regcoil_bnorm < bigno) .or. &
                                (sigma_regcoil_chi2_b < bigno) .or. &
                                (sigma_regcoil_current_density < bigno) )) THEN
         WRITE(6,*)        " Stellarator REGCOIL Optimization provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========                            REGCOIL                            ========="
         WRITE(6,"(2X,A)") "=========                        (M. Landreman)                         ========="
         WRITE(6,"(2X,A)") "=========               Matt dot Landreman at gmail dot com             ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (myid == master .and. ((sigma_regcoil_bnorm < bigno) .or. &
                                (sigma_regcoil_chi2_b < bigno) .or &
                                (sigma_regcoil_current_density < bigno) ) THEN
         sigma_regcoil_bnorm = bigno
         sigma_regcoil_chi2_b = bigno
         sigma_regcoil_current_density = bigno
         WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(6,*) '  Coil optimization with the REGCOIL
         WRITE(6,*) '  code has been disabled.  Coil optimziation'
         WRITE(6,*) '  has been turned off.  Contact your vendor for'
         WRITE(6,*) '  further information.'
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (DKES_OPT)
      IF (myid == master .and. ANY(sigma_dkes < bigno)) THEN
         WRITE(6,*)        " Drift-Kinetic Equation Solver (DKES) provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========           Drift Kinetic Equation Solver, Variational          ========="
         WRITE(6,"(2X,A)") "=========                    (S.Hirshman, W. van Rij)                   ========="
         WRITE(6,"(2X,A)") "=========                      hirshmansp@ornl.gov                      ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_dkes < bigno)) THEN
         sigma_dkes(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Drift-kinetic optimization with the DKES'
            WRITE(6,*) '  code has been disabled.  Drift-kinetic optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
      IF (myid == master .and. (ANY(sigma_fluxloop < bigno) .or. ANY(sigma_bprobe < bigno) .or. ANY(sigma_segrog < bigno) )) THEN
         WRITE(6,*)        " Magnetic Diagnostic calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A,F5.2,A)") "=========                    DIAGNO (v",DIAGNO_VERSION,")                             ========="
         WRITE(6,"(2X,A)") "=========            (S.Lazerson, H Gardner, J. Geiger)                 ========="
         WRITE(6,"(2X,A)") "=========                   lazerson@pppl.gov                           ========="
         WRITE(6,"(2X,A)") "=========       http://vmecwiki.pppl.wikispaces.net/DIAGNO              ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ IF DEFINED (NEO_OPT)
      IF (myid == master .and. ANY(sigma_neo < bigno)) THEN
         WRITE(6,*)        " Neoclassical Transport calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========                      NEO (v3.02)                              ========="
         WRITE(6,"(2X,A)") "=========               (W.Kernbichler and S.Kasilov)                   ========="
         WRITE(6,"(2X,A)") "=========               kernbichler@itp.tu-graz.ac.at                   ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_neo < bigno)) THEN
         sigma_neo(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Neoclassical transport optimization with the NEO'
            WRITE(6,*) '  code has been disabled.  Neoclassical optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (TXPORT_OPT)
      IF (myid == master .and. ANY(sigma_txport < bigno)) THEN
         WRITE(6,*)        " Geometry Interface to Turbulent Transport provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)")       "=========        Geometry Interface for Stellarators and Tokamaks       ========="
         WRITE(6,"(2X,A)")       "=========          (P.Xanthopoulos, W.A.Cooper, and Yu.Turkin)          ========="
         WRITE(6,"(2X,A)")       "=========          pax@ipp.mpg.de  http://www.ipp.mpg.de/~pax/          ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
         WRITE(6,"(2X,A)") "     NOTICE: New TXPORT variables now used to control execution COORDINATES,"
         WRITE(6,"(2X,A)") "             IN_OUT, and SETUP namelists are ignored.  LGLOBAL_TXPORT,"
         WRITE(6,"(2X,A)") "             NZ_TXPORT, NALPHA_TXPORT, ALPHA0_TXPORT have been added to the"
         WRITE(6,"(2X,A)") "             OPTIMUM namelist."
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      IF (ANY(sigma_txport < bigno)) THEN
         sigma_txport(:) = bigno
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  Turbulent transport optimization with the GIST/TXPORT'
            WRITE(6,*) '  code has been disabled.  Turbulent optimziation'
            WRITE(6,*) '  has been turned off.  Contact your vendor for'
            WRITE(6,*) '  further information.'
         END IF
      END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (GENE)
      CALL tolower(txport_proxy)
      IF (myid == master .and. ANY(sigma_txport < bigno) .and. (TRIM(txport_proxy(1:4)) == 'gene') ) THEN
         WRITE(6,*)        " Turbulent Transport calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========       Gyrokinetic Electromagnetic Numerical Experiment        ========="
         WRITE(6,"(2X,A)") "=========                  (F. Jenko, P. Xanthopoulos)                  ========="
         WRITE(6,"(2X,A)") "=========              http://www2.ipp.mpg.de/~fsj/gene/                ========="
         !WRITE(6,"(2X,A)") "=========              GENE11  (release "//release_gene//")                  ========="
         !WRITE(6,"(2X,A)") "=========                      (svn-revision: "//svn_gene//")     ========="
         WRITE(6,"(2X,A)") "=========              GENE11  (git-branch "//release_gene//")                  ========="
         WRITE(6,"(2X,A)") "=========                      (git-master: "//svn_gene//")     ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
!DEC$ ELSE
      CALL tolower(txport_proxy)
      IF (ANY(sigma_txport < bigno) .and. (TRIM(txport_proxy(1:4)) == 'gene')) THEN
         txport_proxy = 'prox1d'
         IF (myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  STELLOPT has not been linked to the GENE code.      '
            WRITE(6,*) '  Optimization with linear GENE for turblent'
            WRITE(6,*) '  transport not possible.  Defaulting to proxy function'
            WRITE(6,*) '        txport_proxy = prox1d'
         END IF
      END IF
!DEC$ ENDIF
      ! Force some behavior
      lbooz(1) = .FALSE.
      target_balloon(1)   = 0.0;  sigma_balloon(1)   = bigno
      target_bootstrap(1) = 0.0;  sigma_bootstrap(1) = bigno
      target_neo(1)       = 0.0;  sigma_neo(1)       = bigno
      target_dkes(1)      = 0.0;  sigma_dkes(1)      = bigno
      target_dkes(2)      = 0.0;  sigma_dkes(2)      = bigno
      target_helicity(1)  = 0.0;  sigma_helicity(1)  = bigno
      END SUBROUTINE read_stellopt_input
      
      SUBROUTINE write_optimum_namelist(iunit,istat)
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(in) :: istat
      INTEGER :: ik, n, m, u, v
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,L1,2(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3))"
      WRITE(iunit,'(A)') '&OPTIMUM'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(A)') '!       Optimizer Run Control Parameters'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,outint) 'NFUNC_MAX',nfunc_max
      WRITE(iunit,outstr) 'EQUIL_TYPE',TRIM(equil_type)
      WRITE(iunit,outstr) 'OPT_TYPE',TRIM(opt_type)
      WRITE(iunit,outflt) 'FTOL',ftol
      WRITE(iunit,outflt) 'XTOL',xtol
      WRITE(iunit,outflt) 'GTOL',gtol
      WRITE(iunit,outflt) 'EPSFCN',epsfcn
      WRITE(iunit,outflt) 'FACTOR',factor
      WRITE(iunit,outint) 'MODE',mode
      WRITE(iunit,outint) 'CR_STRATEGY',cr_strategy
      WRITE(iunit,outint) 'NPOPULATION',npopulation
      WRITE(iunit,outint) 'NOPTIMIZERS',noptimizers
      WRITE(iunit,outboo) 'LKEEP_MINS',lkeep_mins
      !WRITE(iunit,outboo) 'LREFIT',lrefit
      !WRITE(iunit,outflt) 'REFIT_PARAM',refit_param
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(A)') '!       Optimized Quantities'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      IF (lregcoil_winding_surface_separation_opt) THEN
         WRITE(iunit,'(2X,A,E22.14)') &
                'REGCOIL_WINDING_SURFACE_SEPARATION = ', &
                regcoil_winding_surface_separation
         WRITE(iunit,onevar) 'LREGCOIL_WINDING_SURFACE_SEPARATION', & 
                lregcoil_winding_surface_separation_opt, &
                'REGCOIL_WINDING_SURFACE_SEPARATION_MIN', &
                regcoil_winding_surface_separation_min, &
                'REGCOIL_WINDING_SURFACE_SEPARATION_MAX', &
               regcoil_winding_surface_separation_max
         IF (dregcoil_winding_surface_separation_opt > 0) &
                 WRITE(iunit,outflt) 'DREGCOIL_WINDING_SURFACE_SEPARATION', &
                 dregcoil_winding_surface_separation_opt
      END IF
      IF (lregcoil_current_density_opt) THEN
         WRITE(iunit,'(2X,A,E22.14)') &
                'REGCOIL_CURRENT_DENSITY = ', &
                regcoil_current_density
         WRITE(iunit,onevar) 'LREGCOIL_CURRENT_DENSITY', & 
                lregcoil_current_density_opt, &
                'REGCOIL_CURRENT_DENSITY_MIN', &
                regcoil_current_density_min, &
                'REGCOIL_CURRENT_DENSITY_MAX', &
               regcoil_current_density_max
         IF (dregcoil_current_density_opt > 0) &
                 WRITE(iunit,outflt) 'DREGCOIL_CURRENT_DENSITY', &
                 dregcoil_current_density_opt
      END IF
      IF (lphiedge_opt) THEN
         WRITE(iunit,onevar) 'LPHIEDGE_OPT',lphiedge_opt,'PHIEDGE_MIN',phiedge_min,'PHIEDGE_MAX',phiedge_max
         IF (dphiedge_opt > 0) WRITE(iunit,outflt) 'DPHIEDGE_OPT',dphiedge_opt
      END IF
      IF (lcurtor_opt) THEN
         WRITE(iunit,onevar) 'LCURTOR_OPT',lcurtor_opt,'CURTOR_MIN',curtor_min,'CURTOR_MAX',curtor_max
         IF (dcurtor_opt > 0)  WRITE(iunit,outflt) 'DCURTOR_OPT',dcurtor_opt
      END IF
      IF (lpscale_opt) THEN
         WRITE(iunit,onevar) 'LPSCALE_OPT',lpscale_opt,'PSCALE_MIN',pscale_min,'PSCALE_MAX',pscale_max
         IF (dpscale_opt > 0)  WRITE(iunit,outflt) 'DPSCALE_OPT',dpscale_opt
      END IF
      IF (lbcrit_opt) THEN
         WRITE(iunit,onevar) 'LBCRIT_OPT',lbcrit_opt,'BCRIT_MIN',bcrit_min,'BCRIT_MAX',bcrit_max
         IF (dbcrit_opt > 0)   WRITE(iunit,outflt) 'DBCRIT_OPT',dbcrit_opt
      END IF
      IF (lmix_ece_opt) THEN
         WRITE(iunit,onevar) 'LMIX_ECE_OPT',lmix_ece_opt,'MIX_ECE_MIN',mix_ece_min,'MIX_ECE_MAX',mix_ece_max
         IF (dmix_ece_opt > 0)   WRITE(iunit,outflt) 'DMIX_ECE_OPT',dmix_ece_opt
      END IF
      IF (ANY(lextcur_opt)) THEN
        n=0
        DO ik = LBOUND(lextcur_opt,DIM=1), UBOUND(lextcur_opt,DIM=1)
           IF(lextcur_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LEXTCUR_OPT',ik,lextcur_opt(ik),'EXTCUR_MIN',ik,extcur_min(ik),'EXTCUR_MAX',ik,extcur_max(ik)
        END DO
        IF (ANY(dextcur_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DEXTCUR_OPT',(dextcur_opt(ik), ik = 1, n)
      END IF
      IF (ANY(laphi_opt)) THEN
        n=0
        DO ik = LBOUND(laphi_opt,DIM=1), UBOUND(laphi_opt,DIM=1)
           IF(laphi_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LAPHI_OPT',ik,laphi_opt(ik),'APHI_MIN',ik,aphi_min(ik),'APHI_MAX',ik,aphi_max(ik)
        END DO
        IF (ANY(daphi_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAPHI_OPT',(daphi_opt(ik), ik = 1, n)
      END IF
      IF (ANY(lam_opt)) THEN
        n=0
        DO ik = LBOUND(lam_opt,DIM=1), UBOUND(lam_opt,DIM=1)
           IF(lam_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LAM_OPT',ik,lam_opt(ik),'AM_MIN',ik,am_min(ik),'AM_MAX',ik,am_max(ik)
        END DO
        IF (ANY(dam_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAM_OPT',(dam_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lac_opt)) THEN
        n=0
        DO ik = LBOUND(lac_opt,DIM=1),UBOUND(lac_opt,DIM=1)
           IF(lac_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LAC_OPT',ik,lac_opt(ik),'AC_MIN',ik,ac_min(ik),'AC_MAX',ik,ac_max(ik)
        END DO
        IF (ANY(dac_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAC_OPT',(dac_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lai_opt)) THEN
        n=0
        DO ik = LBOUND(lai_opt,DIM=1),UBOUND(lai_opt,DIM=1)
           IF(lai_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LAI_OPT',ik,lai_opt(ik),'AI_MIN',ik,ai_min(ik),'AI_MAX',ik,ai_max(ik)
        END DO
        IF (ANY(dai_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAI_OPT',(dai_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lah_opt)) THEN
        n=0
        DO ik = LBOUND(lah_opt,DIM=1), UBOUND(lah_opt,DIM=1)
           IF(lah_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LAH_OPT',ik,lah_opt(ik),'AH_MIN',ik,ah_min(ik),'AH_MAX',ik,ah_max(ik)
        END DO
        IF (ANY(dah_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAH_OPT',(dah_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lat_opt)) THEN
        n=0
        DO ik = LBOUND(lat_opt,DIM=1),UBOUND(lat_opt,DIM=1)
           IF(lat_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LAT_OPT',ik,lat_opt(ik),'AT_MIN',ik,at_min(ik),'AT_MAX',ik,at_max(ik)
        END DO
        IF (ANY(dat_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAT_OPT',(dat_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lne_opt)) THEN
        n=0
        DO ik = LBOUND(lne_opt,DIM=1), UBOUND(lne_opt,DIM=1)
           IF(lne_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LNE_OPT',ik,lne_opt(ik),'NE_MIN',ik,ne_min(ik),'NE_MAX',ik,ne_max(ik)
        END DO
        IF (ANY(dne_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DNE_OPT',(dne_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lzeff_opt)) THEN
        n=0
        DO ik = LBOUND(lzeff_opt,DIM=1),UBOUND(lzeff_opt,DIM=1)
           IF(lzeff_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LZEFF_OPT',ik,lzeff_opt(ik),'ZEFF_MIN',ik,zeff_min(ik),'ZEFF_MAX',ik,zeff_max(ik)
        END DO
        IF (ANY(dne_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DZEFF_OPT',(dne_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lte_opt)) THEN
        n=0
        DO ik = LBOUND(lte_opt,DIM=1),UBOUND(lte_opt,DIM=1)
           IF(lte_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LTE_OPT',ik,lte_opt(ik),'TE_MIN',ik,te_min(ik),'TE_MAX',ik,te_max(ik)
        END DO
        IF (ANY(dte_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DTE_OPT',(dte_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lti_opt)) THEN
        n=0
        DO ik = LBOUND(lti_opt,DIM=1),UBOUND(lti_opt,DIM=1)
           IF(lti_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LTI_OPT',ik,lti_opt(ik),'TI_MIN',ik,ti_min(ik),'TI_MAX',ik,ti_max(ik)
        END DO
        IF (ANY(dti_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DTI_OPT',(dti_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lth_opt)) THEN
        n=0
        DO ik = LBOUND(lth_opt,DIM=1),UBOUND(lth_opt,DIM=1)
           IF(lth_opt(ik)) n=ik
        END DO
        DO ik = 0, n
           WRITE(iunit,vecvar) 'LTH_OPT',ik,lth_opt(ik),'TH_MIN',ik,th_min(ik),'TH_MAX',ik,th_max(ik)
        END DO
        IF (ANY(dth_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DTH_OPT',(dth_opt(ik), ik = 0, n)
      END IF
      IF (ANY(lam_s_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lam_s_opt,DIM=1)
           IF(lam_s_opt(ik)) n=ik
        END DO
        WRITE(iunit,"(2X,A,1X,'=',10(2X,L1))") 'LAM_S_OPT',(lam_s_opt(ik), ik = 1, n)
        IF (ANY(dam_s_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAM_S_OPT',(dam_s_opt(ik), ik = 1, n)
      END IF
      IF (ANY(lam_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lam_f_opt,DIM=1)
           IF(lam_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LAM_F_OPT',ik,lam_f_opt(ik),'AM_F_MIN',ik,am_f_min(ik),'AM_F_MAX',ik,am_f_max(ik)
        END DO
        IF (ANY(dam_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAM_F_OPT',(dam_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lac_s_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lac_s_opt,DIM=1)
           IF(lac_s_opt(ik)) n=ik
        END DO
        WRITE(iunit,"(2X,A,1X,'=',10(2X,L1))") 'LAC_S_OPT',(lac_s_opt(ik), ik = 1, n)
        IF (ANY(dac_s_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAC_S_OPT',(dac_s_opt(ik), ik = 1, n)
      END IF
      IF (ANY(lac_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lac_f_opt,DIM=1)
           IF(lac_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LAC_F_OPT',ik,lac_f_opt(ik),'AC_F_MIN',ik,ac_f_min(ik),'AC_F_MAX',ik,ac_f_max(ik)
        END DO
        IF (ANY(dac_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAC_F_OPT',(dac_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lai_s_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lai_s_opt,DIM=1)
           IF(lai_s_opt(ik)) n=ik
        END DO
        WRITE(iunit,"(2X,A,1X,'=',10(2X,L1))") 'LAI_S_OPT',(lai_s_opt(ik), ik = 1, n)
        IF (ANY(dai_s_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAI_S_OPT',(dai_s_opt(ik), ik = 1, n)
      END IF
      IF (ANY(lai_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lai_f_opt,DIM=1)
           IF(lai_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LAI_F_OPT',ik,lai_f_opt(ik),'AI_F_MIN',ik,ai_f_min(ik),'AI_F_MAX',ik,ai_f_max(ik)
        END DO
        IF (ANY(dai_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAI_F_OPT',(dai_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lphi_s_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lphi_s_opt,DIM=1)
           IF(lphi_s_opt(ik)) n=ik
        END DO
        WRITE(iunit,"(2X,A,1X,'=',10(2X,L1))") 'LPHI_S_OPT',(lphi_s_opt(ik), ik = 1, n)
        IF (ANY(dphi_s_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DPHI_S_OPT',(dphi_s_opt(ik), ik = 1, n)
      END IF
      IF (ANY(lphi_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lphi_f_opt,DIM=1)
           IF(lphi_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LPHI_F_OPT',ik,lphi_f_opt(ik),'PHI_F_MIN',ik,phi_f_min(ik),'PHI_F_MAX',ik,phi_f_max(ik)
        END DO
        IF (ANY(dphi_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DPHI_F_OPT',(dphi_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lne_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lne_f_opt,DIM=1)
           IF(lne_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LNE_F_OPT',ik,lne_f_opt(ik),'NE_F_MIN',ik,ne_f_min(ik),'NE_F_MAX',ik,ne_f_max(ik)
        END DO
        IF (ANY(dne_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DNE_F_OPT',(dne_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lzeff_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lzeff_f_opt,DIM=1)
           IF(lzeff_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LZEFF_F_OPT',ik,lzeff_f_opt(ik),'ZEFF_F_MIN',ik,zeff_f_min(ik),'ZEFF_F_MAX',ik,zeff_f_max(ik)
        END DO
        IF (ANY(dne_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DZEFF_F_OPT',(dzeff_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lte_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lte_f_opt,DIM=1)
           IF(lte_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LTE_F_OPT',ik,lte_f_opt(ik),'TE_F_MIN',ik,te_f_min(ik),'TE_F_MAX',ik,te_f_max(ik)
        END DO
        IF (ANY(dte_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DTE_F_OPT',(dte_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lti_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lti_f_opt,DIM=1)
           IF(lti_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LTI_F_OPT',ik,lti_f_opt(ik),'TI_F_MIN',ik,ti_f_min(ik),'TI_F_MAX',ik,ti_f_max(ik)
        END DO
        IF (ANY(dti_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'dti_f_opt',(dti_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lth_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lth_f_opt,DIM=1)
           IF(lth_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LTH_F_OPT',ik,lth_f_opt(ik),'TH_F_MIN',ik,th_f_min(ik),'TH_F_MAX',ik,th_f_max(ik)
        END DO
        IF (ANY(dth_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DTH_F_OPT',(dth_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lah_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lah_f_opt,DIM=1)
           IF(lah_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LAH_F_OPT',ik,lah_f_opt(ik),'AH_F_MIN',ik,ah_f_min(ik),'AH_F_MAX',ik,ah_f_max(ik)
        END DO
        IF (ANY(dah_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAH_F_OPT',(dah_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lat_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lat_f_opt,DIM=1)
           IF(lat_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LAT_F_OPT',ik,lat_f_opt(ik),'AT_F_MIN',ik,at_f_min(ik),'AT_F_MAX',ik,at_f_max(ik)
        END DO
        IF (ANY(dat_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DAT_F_OPT',(dat_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lbeamj_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lbeamj_f_opt,DIM=1)
           IF(lbeamj_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LBEAMJ_F_OPT',ik,lbeamj_f_opt(ik),'BEAMJ_F_MIN',ik,beamj_f_min(ik),'BEAMJ_F_MAX',ik,beamj_f_max(ik)
        END DO
        IF (ANY(dbeamj_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DBEAMJ_F_OPT',(dbeamj_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(lbootj_f_opt)) THEN
        n=0
        DO ik = 1,UBOUND(lbootj_f_opt,DIM=1)
           IF(lbootj_f_opt(ik)) n=ik
        END DO
        DO ik = 1, n
           WRITE(iunit,vecvar) 'LBOOTJ_F_OPT',ik,lbootj_f_opt(ik),'BOOTJ_F_MIN',ik,bootj_f_min(ik),'BOOTJ_F_MAX',ik,bootj_f_max(ik)
        END DO
        IF (ANY(dbootj_f_opt > 0)) WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'DBOOTJ_F_OPT',(dbootj_f_opt(ik), ik = 1, n)
      END IF
      
      IF (ANY(laxis_opt)) THEN
         DO n = LBOUND(laxis_opt,DIM=1), UBOUND(laxis_opt,DIM=1)
            IF (laxis_opt(n) .and. (raxis_min(n)>-bigno .or. raxis_max(n)<bigno .or. zaxis_min(n)>-bigno .or. zaxis_max(n)<bigno)) THEN
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',1X,L1,5(2X,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LAXIS_OPT(',n,')',laxis_opt(n),&
                 'RAXIS_MIN(',n,')',raxis_min(n),&
                 'RAXIS_MAX(',n,')',raxis_max(n),&
                 'ZAXIS_MIN(',n,')',zaxis_min(n),&
                 'ZAXIS_MAX(',n,')',zaxis_max(n),&
                 'DAXIS_OPT(',n,')',daxis_opt(n)
            ELSEIF (laxis_opt(n)) THEN
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',1X,L1,1(2X,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LAXIS_OPT(',n,')',laxis_opt(n),&
                 'DAXIS_OPT(',n,')',daxis_opt(n)
            END IF
         END DO
      END IF
      IF (ANY(lrho_opt)) THEN
         DO m = LBOUND(lrho_opt,DIM=2), UBOUND(lrho_opt,DIM=2)
            DO n = LBOUND(lrho_opt,DIM=1), UBOUND(lrho_opt,DIM=1)
               IF(lrho_opt(n,m) .and. (bound_min(n,m)>-bigno .or. bound_max(n,m)<bigno)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,4(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LRHO_OPT(',n,',',m,')',lrho_opt(n,m),&
                 'BOUND_MIN(',n,',',m,')',bound_min(n,m),&
                 'BOUND_MAX(',n,',',m,')',bound_max(n,m),&
                 'DRHO_OPT(',n,',',m,')',drho_opt(n,m)
               ELSEIF (lrho_opt(n,m)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,1(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LRHO_OPT(',n,',',m,')',lrho_opt(n,m),&
                 'DRHO_OPT(',n,',',m,')',drho_opt(n,m)
               END IF
            END DO
         END DO
         WRITE(iunit,outint) 'RHO_EXP',rho_exp
      END IF
      IF (ANY(ldeltamn_opt)) THEN
         DO m = LBOUND(ldeltamn_opt,DIM=2), UBOUND(ldeltamn_opt,DIM=2)
            DO n = LBOUND(ldeltamn_opt,DIM=1), UBOUND(ldeltamn_opt,DIM=1)
               IF(ldeltamn_opt(n,m) .and. (delta_min(n,m)>-bigno .or. delta_max(n,m)<bigno)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,4(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LDELTAMN_OPT(',n,',',m,')',ldeltamn_opt(n,m),&
                 'DELTA_MIN(',n,',',m,')',delta_min(n,m),&
                 'DELTA_MAX(',n,',',m,')',delta_max(n,m),&
                 'DDELTAMN_OPT(',n,',',m,')',ddeltamn_opt(n,m)
               ELSEIF (ldeltamn_opt(n,m)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,1(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LDELTAMN_OPT(',n,',',m,')',ldeltamn_opt(n,m),&
                 'DDELTAMN_OPT(',n,',',m,')',ddeltamn_opt(n,m)
               END IF
            END DO
         END DO
      END IF
      IF (ANY(lmode_opt)) THEN
         DO m = LBOUND(lmode_opt,DIM=2), UBOUND(lmode_opt,DIM=2)
           DO n = LBOUND(lmode_opt,DIM=1), UBOUND(lmode_opt,DIM=1)
               IF(lmode_opt(n,m) .and. (bound_min(n,m)>-bigno .or. bound_max(n,m)<bigno)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,4(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LMODE_OPT(',n,',',m,')',lmode_opt(n,m),&
                 'BOUND_MIN(',n,',',m,')',bound_min(n,m),&
                 'BOUND_MAX(',n,',',m,')',bound_max(n,m),&
                 'DBOUND_OPT(',n,',',m,')',dbound_opt(n,m)
               ELSEIF (lrho_opt(n,m)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,1(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LMODE_OPT(',n,',',m,')',lmode_opt(n,m),&
                 'DBOUND_OPT(',n,',',m,')',dbound_opt(n,m)
               END IF
           END DO
        END DO
      END IF

      
      IF (ANY(lbound_opt)) THEN
         DO m = LBOUND(lbound_opt,DIM=2), UBOUND(lbound_opt,DIM=2)
           DO n = LBOUND(lbound_opt,DIM=1), UBOUND(lbound_opt,DIM=1)
              IF(lbound_opt(n,m)) THEN
                 WRITE(iunit,"(2X,A,I4.3,A,I4.3,A,1X,'=',1X,L1,5(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                 'LBOUND_OPT(',n,',',m,')',lbound_opt(n,m),&
                 'RBC_MIN(',n,',',m,')',rbc_min(n,m),&
                 'RBC_MAX(',n,',',m,')',rbc_max(n,m),&
                 'ZBS_MIN(',n,',',m,')',zbs_min(n,m),&
                 'ZBS_MAX(',n,',',m,')',zbs_max(n,m),&
                 'DBOUND_OPT(',n,',',m,')',dbound_opt(n,m)
                 IF (lasym_local) THEN
                    WRITE(iunit,"(4(2X,A,I4.3,A,I4.3,A,1X,'=',1X,E22.14))")&
                    'RBS_MIN(',n,',',m,')',rbs_min(n,m),&
                    'RBS_MAX(',n,',',m,')',rbs_max(n,m),&
                    'ZBC_MIN(',n,',',m,')',zbc_min(n,m),&
                    'ZBC_MAX(',n,',',m,')',zbc_max(n,m)
                 END IF
              END IF
           END DO
        END DO
      END IF

      IF (ANY(lcoil_spline)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!       Coil Splines'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         ! For now assumes sx,sy, and sz are the same size.
         DO n = LBOUND(lcoil_spline,DIM=1), UBOUND(lcoil_spline,DIM=1)
            IF (ANY(coil_splinesx(n,:)>-1)) THEN
               WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
               WRITE(iunit,'(A,I4.3)') '!       Coil Number ',n
               ik = MINLOC(coil_splinesx(n,:),DIM=1)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,L1))") 'COIL_SPLINESX(',n,',:)',(lcoil_spline(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'DCOIL_SPLINE(',n,',:)',(dcoil_spline(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINESX(',n,',:)',(coil_splinesx(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFX(',n,',:)',(coil_splinefx(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINESY(',n,',:)',(coil_splinesy(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFY(',n,',:)',(coil_splinefy(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINESZ(',n,',:)',(coil_splinesz(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFZ(',n,',:)',(coil_splinefz(n,m), m = 1, ik)
               ! Min/Max
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFX_MIN(',n,',:)',(coil_splinefx_min(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFX_MAX(',n,',:)',(coil_splinefx_max(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFY_MIN(',n,',:)',(coil_splinefy_min(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFY_MAX(',n,',:)',(coil_splinefy_max(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFZ_MIN(',n,',:)',(coil_splinefz_min(n,m), m = 1, ik)
               WRITE(iunit,"(2X,A,I4.3,A,1X,'=',5(2X,E22.14))") 'COIL_SPLINEFZ_MAX(',n,',:)',(coil_splinefz_max(n,m), m = 1, ik)
            END IF
         END DO
      END IF
      
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(A)') '!       Profile Functions'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      ! NE
      n = MINLOC(ne_opt,DIM=1)
      m = MINLOC(ne_aux_s(2:),DIM=1)
      IF (n > 1 .or. m > 4)  WRITE(iunit,outstr) 'NE_TYPE',TRIM(ne_type)
      IF (n > 1) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(2X,E22.14))") 'NE_OPT',(ne_opt(ik), ik = 0, 20)
      END IF
      IF (m > 5) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'NE_AUX_S',(ne_aux_s(ik), ik=1,m)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'NE_AUX_F',(ne_aux_f(ik), ik=1,m)
      END IF
      ! ZEFF
      n = MINLOC(zeff_opt,DIM=1)
      m = MINLOC(zeff_aux_s(2:),DIM=1)
      IF (n > 1 .or. m > 4)  WRITE(iunit,outstr) 'ZEFF_TYPE',TRIM(zeff_type)
      IF (n > 1) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(2X,E22.14))") 'ZEFF_OPT',(zeff_opt(ik), ik = 0, 20)
      END IF
      IF (m > 4) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'ZEFF_AUX_S',(zeff_aux_s(ik), ik=1,m)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'ZEFF_AUX_F',(zeff_aux_f(ik), ik=1,m)
      END IF
      ! TE
      n = MINLOC(te_opt,DIM=1)
      m = MINLOC(te_aux_s(2:),DIM=1)
      IF (n > 1 .or. m > 4)  WRITE(iunit,outstr) 'TE_TYPE',TRIM(te_type)
      IF (n > 1) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(2X,E22.14))") 'TE_OPT',(te_opt(ik), ik = 0, 20)
      END IF
      IF (m > 4) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'TE_AUX_S',(te_aux_s(ik), ik=1,m)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'TE_AUX_F',(te_aux_f(ik), ik=1,m)
      END IF
      ! TI
      n = MINLOC(ti_opt,DIM=1)
      m = MINLOC(ti_aux_s(2:),DIM=1)
      IF (n > 1 .or. m > 4)  WRITE(iunit,outstr) 'TI_TYPE',TRIM(ti_type)
      IF (n > 1) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(2X,E22.14))") 'TI_OPT',(ti_opt(ik), ik = 0, 20)
      END IF
      IF (m > 4) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'TI_AUX_S',(ti_aux_s(ik), ik=1,m)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'TI_AUX_F',(ti_aux_f(ik), ik=1,m)
      END IF
      ! Currents
      ik = MINLOC(beamj_aux_s(2:),DIM=1)
      IF (ik > 2) THEN
         WRITE(iunit,outstr) 'BEAMJ_TYPE',TRIM(beamj_type)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'BEAMJ_AUX_S',(beamj_aux_s(n), n=1,ik)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'BEAMJ_AUX_F',(beamj_aux_f(n), n=1,ik)
      END IF
      ik = MINLOC(bootj_aux_s(2:),DIM=1)
      IF (ik > 2) THEN
         WRITE(iunit,outstr) 'BOOTJ_TYPE',TRIM(bootj_type)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'BOOTJ_AUX_S',(bootj_aux_s(n), n=1,ik)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'BOOTJ_AUX_F',(bootj_aux_f(n), n=1,ik)
      END IF
      ! E-static potential
      ik = MINLOC(phi_aux_s(2:),DIM=1)
      IF (ik > 4) THEN
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'PHI_AUX_S',(phi_aux_s(n), n=1,ik)
         WRITE(iunit,"(2X,A,1X,'=',5(1X,E22.14))") 'PHI_AUX_F',(phi_aux_f(n), n=1,ik)
      END IF
      IF (ANY(lbooz)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          BOOZER COORDINATE TRANSFORMATION'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outint) 'MBOZ',mboz 
         WRITE(iunit,outint) 'NBOZ',nboz
      END IF
      IF (ANY(sigma_helicity < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          BOOZER COORDINATE HELICITY'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outcmp) 'HELICITY',NINT(REAL(helicity)),NINT(AIMAG(helicity))
         n=0
         DO ik = 1,UBOUND(sigma_helicity,DIM=1)
            IF(sigma_helicity(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_helicity(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_HELICITY(',ik,') = ',target_helicity(ik), &
                          'SIGMA_HELICITY(',ik,') = ',sigma_helicity(ik)
         END DO
      END IF
      IF (ANY(sigma_helicity_old < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          BOOZER COORDINATE HELICITY (OLD)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outcmp) 'HELICITY',NINT(REAL(helicity)),NINT(AIMAG(helicity))
         n=0
         DO ik = 1,UBOUND(sigma_helicity_old,DIM=1)
            IF(sigma_helicity_old(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_helicity_old(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_HELICITY_OLD(',ik,') = ',target_helicity_old(ik), &
                          'SIGMA_HELICITY_OLD(',ik,') = ',sigma_helicity_old(ik)
         END DO
      END IF
      IF (ANY(sigma_resjac < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          BOOZER Resonant Modes'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_resjac,DIM=1)
            IF(sigma_resjac(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_resjac(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_RESJAC(',ik,') = ',target_resjac(ik), &
                          'SIGMA_RESJAC(',ik,') = ',sigma_resjac(ik), &
                          'XM_RESJAC(',ik,') = ',xm_resjac(ik), &
                          'XN_RESJAC(',ik,') = ',xn_resjac(ik)
         END DO
      END IF
      IF (ANY(sigma_balloon < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          BALLOONING CALCULATION'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=COUNT(balloon_theta >= 0.0)
         WRITE(iunit,"(2X,A,1X,'=',10(2X,E22.14))") 'BALLOON_THETA',(balloon_theta(ik), ik = 1, n)
         n=COUNT(balloon_zeta >= 0.0)
         WRITE(iunit,"(2X,A,1X,'=',10(2X,E22.14))") 'BALLOON_ZETA',(balloon_zeta(ik), ik = 1, n)
         n=0
         DO ik = 1,UBOUND(sigma_balloon,DIM=1)
            IF(sigma_balloon(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_balloon(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_BALLOON(',ik,') = ',target_balloon(ik), &
                          'SIGMA_BALLOON(',ik,') = ',sigma_balloon(ik)
         END DO
      END IF
      IF (ANY(sigma_bootstrap < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          BOOTSTRAP CALCULATION'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_bootstrap,DIM=1)
            IF(sigma_bootstrap(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
           IF (sigma_bootstrap(ik) < bigno)  WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_BOOTSTRAP(',ik,') = ',target_bootstrap(ik), &
                          'SIGMA_BOOTSTRAP(',ik,') = ',sigma_bootstrap(ik)
         END DO
      END IF
      IF (ANY(sigma_neo < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          NEOCLASSICAL TRANSPORT (NEO)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_neo,DIM=1)
            IF(sigma_neo(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_neo(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_NEO(',ik,') = ',target_neo(ik), &
                          'SIGMA_NEO(',ik,') = ',sigma_neo(ik)
         END DO
      END IF
      IF (ANY(sigma_kink < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          KINK STABILITY (TERPSICHORE)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,"(2(2X,A,1X,'=',1X,I5))") 'MLMNB_KINK',mlmnb_kink,'IVAC_KINK',ivac_kink
         WRITE(iunit,"(2(2X,A,1X,'=',1X,I5))") 'MMAXDF_KINK',mmaxdf_kink,'NMAXDF_KINK',nmaxdf_kink
         DO ik = 1, UBOUND(sigma_kink,DIM=1)
            IF (sigma_kink(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14),5(2X,A,I3.3,A,I6))") &
                          'TARGET_KINK(',ik,') = ',target_kink(ik), &
                          'SIGMA_KINK(',ik,') = ',sigma_kink(ik),&
                          'MLMNS_KINK(',ik,') = ',mlmns_kink(ik),&
                          'NJ_KINK(',ik,') = ',nj_kink(ik),&
                          'NK_KINK(',ik,') = ',nk_kink(ik),&
                          'LSSL_KINK(',ik,') = ',lssl_kink(ik),&
                          'LSSD_KINK(',ik,') = ',lssd_kink(ik)
         END DO
      END IF
      IF (ANY(sigma_dkes < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          DRIFT-KINETICS (DKES)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_dkes,DIM=1)
            IF(sigma_dkes(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_dkes(ik) < bigno) WRITE(iunit,"(3(2X,A,I3.3,A,E22.14))") &
                          'TARGET_DKES(',ik,') = ',target_dkes(ik), &
                          'SIGMA_DKES(',ik,') = ',sigma_dkes(ik), &
                          'NU_DKES(',ik,') = ',nu_dkes(ik)
         END DO
      END IF
      IF (ANY(sigma_jdotb < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Parllel Current (<J.B>)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_jdotb,DIM=1)
            IF(sigma_jdotb(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
           IF (sigma_jdotb(ik) < bigno)  WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_JDOTB(',ik,') = ',target_jdotb(ik), &
                          'SIGMA_JDOTB(',ik,') = ',sigma_jdotb(ik)
         END DO
      END IF
      IF (ANY(sigma_magwell < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Magnetic Well (W>0 Stable)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_magwell,DIM=1)
            IF(sigma_magwell(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
           IF (sigma_magwell(ik) < bigno)  WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_MAGWELL(',ik,') = ',target_magwell(ik), &
                          'SIGMA_MAGWELL(',ik,') = ',sigma_magwell(ik)
         END DO
      END IF
      IF (ANY(sigma_jcurv < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Toroidal Current (<JCURV>)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_jcurv,DIM=1)
            IF(sigma_jcurv(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_jcurv(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_JCURV(',ik,') = ',target_jcurv(ik), &
                          'SIGMA_JCURV(',ik,') = ',sigma_jcurv(ik)
         END DO
      END IF
      IF (ANY(sigma_bmin < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          |B|_min'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_bmin,DIM=1)
            IF(sigma_bmin(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_bmin(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_BMIN(',ik,') = ',target_bmin(ik), &
                          'SIGMA_BMIN(',ik,') = ',sigma_bmin(ik)
         END DO
      END IF
      IF (ANY(sigma_bmax < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          |B|_max'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         n=0
         DO ik = 1,UBOUND(sigma_bmax,DIM=1)
            IF(sigma_bmax(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_bmax(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_BMAX(',ik,') = ',target_bmax(ik), &
                          'SIGMA_BMAX(',ik,') = ',sigma_bmax(ik)
         END DO
      END IF
      IF (ANY(sigma_Jstar < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          TRAPPED PARTICLE CONFINEMENT (J*)'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,"(2X,A,1X,'=',1(2X,I4.4))") 'NumJstar',NumJstar
         n=0
         DO ik = 1,UBOUND(sigma_Jstar,DIM=1)
            IF(sigma_Jstar(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_Jstar(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_JSTAR(',ik,') = ',target_Jstar(ik), &
                          'SIGMA_JSTAR(',ik,') = ',sigma_Jstar(ik)
         END DO
      END IF
      IF (ANY(sigma_txport < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          TURBULENT TRANSPORT'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outstr) 'TXPORT_PROXY',TRIM(txport_proxy)
         WRITE(iunit,outboo) 'LGLOBAL_TXPORT',lglobal_txport
         WRITE(iunit,outint) 'NZ_TXPORT',nz_txport
         WRITE(iunit,outint) 'NALPHA_TXPORT',nalpha_txport
         WRITE(iunit,outflt) 'ALPHA_START_TXPORT',alpha_start_txport
         WRITE(iunit,outflt) 'ALPHA_END_TXPORT',alpha_end_txport
         n=0
         DO ik = 1,UBOUND(sigma_txport,DIM=1)
            IF(sigma_txport(ik) < bigno) n=ik
         END DO
         DO ik = 1, n
            IF (sigma_txport(ik) < bigno) WRITE(iunit,"(3(2X,A,I3.3,A,E22.14))") &
                          'S_TXPORT(',ik,') = ',s_txport(ik), &
                          'TARGET_TXPORT(',ik,') = ',target_txport(ik), &
                          'SIGMA_TXPORT(',ik,') = ',sigma_txport(ik)
         END DO
      END IF
      IF (ANY(sigma_orbit < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          ORBIT OPTIMIZATION'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outflt) 'MASS_ORBIT',mass_orbit
         WRITE(iunit,outflt) 'Z_ORBIT',Z_orbit
         WRITE(iunit,outint) 'NU_ORBIT',nu_orbit 
         WRITE(iunit,outint) 'NV_ORBIT',nv_orbit
         n=0
         DO ik = 1,UBOUND(sigma_orbit,DIM=1)
            IF(sigma_orbit(ik) < bigno) WRITE(iunit,"(2(2X,A,I3.3,A,E22.14))") &
                          'TARGET_ORBIT(',ik,') = ',target_orbit(ik), &
                          'SIGMA_ORBIT(',ik,') = ',sigma_orbit(ik)
         END DO
         WRITE(iunit,outint) 'NP_ORBIT',np_orbit
         DO ik = 1, np_orbit
            WRITE(iunit,"(3(2X,A,I3.3,A,E22.14))") &
                          'VLL_ORBIT(',ik,') = ',VLL_orbit(ik), &
                          'MU_ORBIT(',ik,') = ',MU_orbit(ik),&
                          'VPERP_ORBIT(',ik,') = ',VPERP_orbit(ik)
         END DO
      END IF
      IF (sigma_coil_bnorm < bigno) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          COIL OPTIMIZATION'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outint) 'NU_BNORM',nu_bnorm 
         WRITE(iunit,outint) 'NV_BNORM',nv_bnorm
         WRITE(iunit,outflt) 'TARGET_COIL_BNORM',target_coil_bnorm
         WRITE(iunit,outflt) 'SIGMA_COIL_BNORM',sigma_coil_bnorm
      END IF
      IF ((sigma_regcoil_bnorm < bigno) .or. &
          (sigma_regcoil_chi2_b < bigno) .or.  &
          (sigma_regcoil_current_density < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          REGCOIL BNORM OPTIMIZATION'  
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outint) 'NU_BNORM',nu_bnorm 
         WRITE(iunit,outint) 'NV_BNORM',nv_bnorm
         WRITE(iunit,outflt) 'TARGET_COIL_BNORM',target_regcoil_bnorm
         WRITE(iunit,outflt) 'SIGMA_COIL_BNORM',sigma_regcoil_bnorm
         WRITE(iunit,outflt) 'TARGET_COIL_CHI2_B',target_regcoil_chi2_b
         WRITE(iunit,outflt) 'SIGMA_COIL_CHI2_B',sigma_regcoil_chi2_b
         WRITE(iunit,outflt) 'TARGET_CURRENT_DENSITY',target_regcoil_current_density
         WRITE(iunit,outflt) 'SIGMA_CURRENT_DENSITY',sigma_regcoil_current_density
      END IF
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(A)') '!         EQUILIBRIUM/GEOMETRY OPTIMIZATION PARAMETERS' 
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      IF (sigma_phiedge < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_PHIEDGE',target_phiedge
         WRITE(iunit,outflt) 'SIGMA_PHIEDGE',sigma_phiedge
      END IF 
      IF (sigma_curtor < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_CURTOR',target_curtor
         WRITE(iunit,outflt) 'SIGMA_CURTOR',sigma_curtor
      END IF 
      IF (sigma_curtor_max < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_CURTOR_MAX',target_curtor_max
         WRITE(iunit,outflt) 'SIGMA_CURTOR_MAX',sigma_curtor_max
      END IF 
      IF (sigma_rbtor < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_RBTOR',target_rbtor
         WRITE(iunit,outflt) 'SIGMA_RBTOR',sigma_rbtor
      END IF 
      IF (sigma_r0 < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_R0',target_r0
         WRITE(iunit,outflt) 'SIGMA_R0',sigma_r0
      END IF 
      IF (sigma_z0 < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_Z0',target_z0
         WRITE(iunit,outflt) 'SIGMA_Z0',sigma_z0
      END IF 
      IF (sigma_volume < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_VOLUME',target_volume
         WRITE(iunit,outflt) 'SIGMA_VOLUME',sigma_volume
      END IF 
      IF (sigma_beta < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_BETA',target_beta
         WRITE(iunit,outflt) 'SIGMA_BETA',sigma_beta
      END IF 
      IF (sigma_betapol < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_BETAPOL',target_betapol
         WRITE(iunit,outflt) 'SIGMA_BETAPOL',sigma_betapol
      END IF 
      IF (sigma_betator < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_BETATOR',target_betator
         WRITE(iunit,outflt) 'SIGMA_BETATOR',sigma_betator
      END IF 
      IF (sigma_wp < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_WP',target_wp
         WRITE(iunit,outflt) 'SIGMA_WP',sigma_wp
      END IF 
      IF (sigma_aspect < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_ASPECT',target_aspect
         WRITE(iunit,outflt) 'SIGMA_ASPECT',sigma_aspect
      END IF 
      IF (sigma_curvature < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_CURVATURE',target_curvature
         WRITE(iunit,outflt) 'SIGMA_CURVATURE',sigma_curvature
      END IF 
      IF (sigma_kappa < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_KAPPA',target_kappa
         WRITE(iunit,outflt) 'SIGMA_KAPPA',sigma_kappa
         WRITE(iunit,outflt) 'PHI_KAPPA',phi_kappa
      END IF 
      IF (sigma_kappa_box < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_KAPPA_BOX',target_kappa_box
         WRITE(iunit,outflt) 'SIGMA_KAPPA_BOX',sigma_kappa_box
         WRITE(iunit,outflt) 'PHI_KAPPA_BOX',phi_kappa_box
      END IF 
      IF (sigma_kappa_avg < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_KAPPA_AVG',target_kappa_avg
         WRITE(iunit,outflt) 'SIGMA_KAPPA_AVG',sigma_kappa_avg
      END IF 
      IF (sigma_aspect_max < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_ASPECT_MAX',target_aspect_max
         WRITE(iunit,outflt) 'SIGMA_ASPECT_MAX',sigma_aspect_max
         WRITE(iunit,outflt) 'WIDTH_ASPECT_MAX',width_aspect_max
      END IF          
      IF (sigma_gradp_max < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_GRADP_MAX',target_gradp_max
         WRITE(iunit,outflt) 'SIGMA_GRADP_MAX',sigma_gradp_max
         WRITE(iunit,outflt) 'WIDTH_GRADP_MAX',width_gradp_max
      END IF          
      IF (sigma_pmin < bigno) THEN
         WRITE(iunit,outflt) 'TARGET_PMIN',target_pmin
         WRITE(iunit,outflt) 'SIGMA_PMIN',sigma_pmin
         WRITE(iunit,outflt) 'WIDTH_PMIN',width_pmin
      END IF
      IF (ANY(sigma_extcur < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Coil Current Optimization'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_extcur,DIM=1)
            IF (sigma_extcur(ik) < bigno) THEN
               WRITE(iunit,"(2(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'TARGET_EXTCUR(',ik,')',target_extcur(ik),& 
                  'SIGMA_EXTCUR(',ik,')',sigma_extcur(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_press < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Plasma Pressure OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outflt) 'NORM_PRESS',norm_press
         DO ik = 1, UBOUND(sigma_press,DIM=1)
            IF (sigma_press(ik) < bigno .and. s_press(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_PRESS(',ik,')',r_press(ik),&
                  'PHI_PRESS(',ik,')',phi_press(ik),& 
                  'Z_PRESS(',ik,')',z_press(ik),&
                  'TARGET_PRESS(',ik,')',target_press(ik),& 
                  'SIGMA_PRESS(',ik,')',sigma_press(ik)
            ELSE IF (sigma_press(ik) < bigno .and. s_press(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_PRESS(',ik,')',s_press(ik),&
                  'TARGET_PRESS(',ik,')',target_press(ik),& 
                  'SIGMA_PRESS(',ik,')',sigma_press(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_ne < bigno_ne)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          ELECTRON DENSITY OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_ne,DIM=1)
            IF (sigma_ne(ik) < bigno_ne .and. s_ne(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_NE(',ik,')',r_ne(ik),&
                  'PHI_NE(',ik,')',phi_ne(ik),& 
                  'Z_NE(',ik,')',z_ne(ik),&
                  'TARGET_NE(',ik,')',target_ne(ik),& 
                  'SIGMA_NE(',ik,')',sigma_ne(ik)
            ELSE IF (sigma_ne(ik) < bigno_ne .and. s_ne(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_NE(',ik,')',s_ne(ik),&
                  'TARGET_NE(',ik,')',target_ne(ik),& 
                  'SIGMA_NE(',ik,')',sigma_ne(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_ne_line < bigno_ne)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          LINE INTEGRATED ELECTRON DENSITY OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_ne_line,DIM=1)
            IF (sigma_ne_line(ik) < bigno_ne) THEN
               WRITE(iunit,"(8(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R0_NE_LINE(',ik,')',r0_ne_line(ik),&
                  'PHI0_NE_LINE(',ik,')',phi0_ne_line(ik),&
                  'Z0_NE_LINE(',ik,')',z0_ne_line(ik),&
                  'R1_NE_LINE(',ik,')',r1_ne_line(ik),&
                  'PHI1_NE_LINE(',ik,')',phi1_ne_line(ik),&
                  'Z1_NE_LINE(',ik,')',z1_ne_line(ik),&
                  'TARGET_NE_LINE(',ik,')',target_ne_line(ik),&
                  'SIGMA_NE_LINE(',ik,')',sigma_ne_line(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_te_line < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          LINE INTEGRATED ELECTRON TEMPERATURE OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,outflt) 'CUTOFF_TE_LINE',cutoff_te_line
         DO ik = 1, UBOUND(sigma_te_line,DIM=1)
            IF (sigma_te_line(ik) < bigno) THEN
               WRITE(iunit,"(8(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R0_TE_LINE(',ik,')',r0_te_line(ik),&
                  'PHI0_TE_LINE(',ik,')',phi0_te_line(ik),&
                  'Z0_TE_LINE(',ik,')',z0_te_line(ik),&
                  'R1_TE_LINE(',ik,')',r1_te_line(ik),&
                  'PHI1_TE_LINE(',ik,')',phi1_te_line(ik),&
                  'Z1_TE_LINE(',ik,')',z1_te_line(ik),&
                  'TARGET_TE_LINE(',ik,')',target_te_line(ik),&
                  'SIGMA_TE_LINE(',ik,')',sigma_te_line(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_ti_line < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          LINE INTEGRATED ION TEMPERATURE OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_ti_line,DIM=1)
            IF (sigma_ti_line(ik) < bigno) THEN
               WRITE(iunit,"(8(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R0_TI_LINE(',ik,')',r0_ti_line(ik),&
                  'PHI0_TI_LINE(',ik,')',phi0_ti_line(ik),&
                  'Z0_TI_LINE(',ik,')',z0_ti_line(ik),&
                  'R1_TI_LINE(',ik,')',r1_ti_line(ik),&
                  'PHI1_TI_LINE(',ik,')',phi1_ti_line(ik),&
                  'Z1_TI_LINE(',ik,')',z1_ti_line(ik),&
                  'TARGET_TI_LINE(',ik,')',target_ti_line(ik),&
                  'SIGMA_TI_LINE(',ik,')',sigma_ti_line(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_te < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          ELECTRON TEMPERATURE OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_te,DIM=1)
            IF (sigma_te(ik) < bigno .and. s_te(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_TE(',ik,')',r_te(ik),&
                  'PHI_TE(',ik,')',phi_te(ik),& 
                  'Z_TE(',ik,')',z_te(ik),&
                  'TARGET_TE(',ik,')',target_te(ik),& 
                  'SIGMA_TE(',ik,')',sigma_te(ik)
            ELSE IF (sigma_te(ik) < bigno .and. s_te(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_TE(',ik,')',s_te(ik),&
                  'TARGET_TE(',ik,')',target_te(ik),& 
                  'SIGMA_TE(',ik,')',sigma_te(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_ti < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          ION TEMPERATURE OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_ti,DIM=1)
            IF (sigma_ti(ik) < bigno .and. s_ti(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_TI(',ik,')',r_ti(ik),&
                  'PHI_TI(',ik,')',phi_ti(ik),& 
                  'Z_TI(',ik,')',z_ti(ik),&
                  'TARGET_TI(',ik,')',target_ti(ik),& 
                  'SIGMA_TI(',ik,')',sigma_ti(ik)
            ELSE IF (sigma_ti(ik) < bigno .and. s_ti(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_TI(',ik,')',s_ti(ik),&
                  'TARGET_TI(',ik,')',target_ti(ik),& 
                  'SIGMA_TI(',ik,')',sigma_ti(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_vphi < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          TOROIDAL ROTATION OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(2X,A,E22.14)') 'QM_RATIO = ',qm_ratio
         DO ik = 1, UBOUND(sigma_vphi,DIM=1)
            IF (sigma_vphi(ik) < bigno .and. s_vphi(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_VPHI(',ik,')',r_vphi(ik),&
                  'PHI_VPHI(',ik,')',phi_vphi(ik),& 
                  'Z_VPHI(',ik,')',z_vphi(ik),&
                  'TARGET_VPHI(',ik,')',target_vphi(ik),& 
                  'SIGMA_VPHI(',ik,')',sigma_vphi(ik)
            ELSE IF (sigma_vphi(ik) < bigno .and. s_vphi(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_VPHI(',ik,')',s_vphi(ik),&
                  'TARGET_VPHI(',ik,')',target_vphi(ik),& 
                  'SIGMA_VPHI(',ik,')',sigma_vphi(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_ece < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          ECE Reflectometry OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(2X,A,E22.14)') 'MIX_ECE = ',mix_ece
         WRITE(iunit,'(2X,A,I3.3)') 'NRA_ECE = ',nra_ece
         WRITE(iunit,'(2X,A,I3.3)') 'NPHI_ECE = ',nphi_ece
         IF (LEN_TRIM(vessel_ece) > 1) WRITE(iunit,outstr) 'VESSEL_ECE',TRIM(vessel_ece)
         IF (LEN_TRIM(mirror_ece) > 1) WRITE(iunit,outstr) 'MIRROR_ECE',TRIM(mirror_ece)
         IF (LEN_TRIM(targettype_ece) > 1) WRITE(iunit,outstr) 'TARGETTYPE_ECE',TRIM(targettype_ece)
         IF (LEN_TRIM(antennatype_ece) > 1) WRITE(iunit,outstr) 'ANTENNATYPE_ECE',TRIM(antennatype_ece)
         DO u = 1, UBOUND(sigma_ece,DIM=1)
               IF (ALL(sigma_ece(u,:) >= bigno)) CYCLE
               WRITE(iunit,"(2X,A,I3.3,A,1X,'=',1X,3E22.14)")'ANTENNAPOSITION_ECE(',u,',1:3)',antennaposition_ece(u,1:3)
               WRITE(iunit,"(2X,A,I3.3,A,1X,'=',1X,3E22.14)")'TARGETPOSITION_ECE(',u,',1:3)',targetposition_ece(u,1:3)
               WRITE(iunit,"(2X,A,I3.3,A,1X,'=',1X,3E22.14)")'RBEAM_ECE(',u,',1:3)',rbeam_ece(u,1:3)
               WRITE(iunit,"(2X,A,I3.3,A,1X,'=',1X,3E22.14)")'RFOCUS_ECE(',u,',1:3)',rfocus_ece(u,1:3)
               DO v = 1, UBOUND(sigma_ece,DIM=2)
                  IF (sigma_ece(u,v) >= bigno) CYCLE
                  WRITE(iunit,"(3(5X,A,I3.3,A,I3.3,A,1X,'=',1X,E22.14))") &
                        'TARGET_ECE(',u,',',v,')',target_ece(u,v),&
                        'SIGMA_ECE(',u,',',v,')',sigma_ece(u,v),& 
                        'FREQ_ECE(',u,',',v,')',freq_ece(u,v)
               END DO
         END DO
      END IF
      IF (ANY(sigma_iota < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Rotational Transform OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_iota,DIM=1)
            IF (sigma_iota(ik) < bigno .and. s_iota(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_IOTA(',ik,')',r_iota(ik),&
                  'PHI_IOTA(',ik,')',phi_iota(ik),& 
                  'Z_IOTA(',ik,')',z_iota(ik),&
                  'TARGET_IOTA(',ik,')',target_iota(ik),& 
                  'SIGMA_IOTA(',ik,')',sigma_iota(ik)
            ELSE IF (sigma_iota(ik) < bigno .and. s_iota(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_IOTA(',ik,')',s_iota(ik),&
                  'TARGET_IOTA(',ik,')',target_iota(ik),& 
                  'SIGMA_IOTA(',ik,')',sigma_iota(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_vaciota < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          Vacuum Rotational Transform OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_vaciota,DIM=1)
            IF (sigma_vaciota(ik) < bigno .and. s_vaciota(ik) < 0) THEN
               WRITE(iunit,"(5(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_VACIOTA(',ik,')',r_vaciota(ik),&
                  'PHI_VACIOTA(',ik,')',phi_vaciota(ik),& 
                  'Z_VACIOTA(',ik,')',z_vaciota(ik),&
                  'TARGET_VACIOTA(',ik,')',target_vaciota(ik),& 
                  'SIGMA_VACIOTA(',ik,')',sigma_vaciota(ik)
            ELSE IF (sigma_vaciota(ik) < bigno .and. s_vaciota(ik) >= 0) THEN
               WRITE(iunit,"(3(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_VACIOTA(',ik,')',s_vaciota(ik),&
                  'TARGET_VACIOTA(',ik,')',target_vaciota(ik),& 
                  'SIGMA_VACIOTA(',ik,')',sigma_vaciota(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_faraday < bigno_ne)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          FARADAY ROTATION OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(sigma_faraday,DIM=1)
            IF (sigma_faraday(ik) < bigno_ne) THEN
               WRITE(iunit,"(8(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R0_FARADAY(',ik,')',r0_faraday(ik),&
                  'PHI0_FARADAY(',ik,')',phi0_faraday(ik),&
                  'Z0_FARADAY(',ik,')',z0_faraday(ik),&
                  'R1_FARADAY(',ik,')',r1_faraday(ik),&
                  'PHI1_FARADAY(',ik,')',phi1_faraday(ik),&
                  'Z1_FARADAY(',ik,')',z1_faraday(ik),&
                  'TARGET_FARADAY(',ik,')',target_faraday(ik),&
                  'SIGMA_FARADAY(',ik,')',sigma_faraday(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_mse < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          MOTIONAL STARK EFFECT OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO ik = 1, UBOUND(lmse_extcur,DIM=1)
            IF (lmse_extcur(ik)) WRITE(iunit,"(2X,A,I3.3,A,1X,'=',1X,L1)") 'LMSE_EXTCUR(',ik,')',lmse_extcur(ik)
         END DO
         DO ik = 1, UBOUND(sigma_mse,DIM=1)
            IF (sigma_mse(ik) < bigno .and. s_mse(ik) < 0) THEN
               WRITE(iunit,"(13(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'R_MSE(',ik,')',r_mse(ik),&
                  'PHI_MSE(',ik,')',phi_mse(ik),& 
                  'Z_MSE(',ik,')',z_mse(ik),&
                  'A1_MSE(',ik,')',a1_mse(ik),&
                  'A2_MSE(',ik,')',a2_mse(ik),&
                  'A3_MSE(',ik,')',a3_mse(ik),&
                  'A4_MSE(',ik,')',a4_mse(ik),&
                  'A5_MSE(',ik,')',a5_mse(ik),&
                  'A6_MSE(',ik,')',a6_mse(ik),&
                  'A7_MSE(',ik,')',a7_mse(ik),&
                  'TARGET_MSE(',ik,')',target_mse(ik),& 
                  'SIGMA_MSE(',ik,')',sigma_mse(ik),& 
                  'VAC_MSE(',ik,')',vac_mse(ik)
            ELSE IF (sigma_mse(ik) < bigno .and. s_mse(ik) >= 0) THEN
               WRITE(iunit,"(11(2X,A,I3.3,A,1X,'=',1X,E22.14))") &
                  'S_MSE(',ik,')',z_mse(ik),&
                  'A1_MSE(',ik,')',a1_mse(ik),&
                  'A2_MSE(',ik,')',a2_mse(ik),&
                  'A3_MSE(',ik,')',a3_mse(ik),&
                  'A4_MSE(',ik,')',a4_mse(ik),&
                  'A5_MSE(',ik,')',a5_mse(ik),&
                  'A6_MSE(',ik,')',a6_mse(ik),&
                  'A7_MSE(',ik,')',a7_mse(ik),&
                  'TARGET_MSE(',ik,')',target_mse(ik),& 
                  'SIGMA_MSE(',ik,')',sigma_mse(ik),& 
                  'VAC_MSE(',ik,')',vac_mse(ik)
            END IF
         END DO
      END IF
      IF (ANY(sigma_bprobe < bigno) .or. ANY(sigma_fluxloop < bigno) .or. ANY(sigma_segrog < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          MAGNETIC DIAGNOSTIC OPTIMIZATION'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         IF (LEN_TRIM(magdiag_coil) > 1) WRITE(iunit,outstr) 'MAGDIAG_COIL',TRIM(magdiag_coil)
         DO ik = 1, UBOUND(sigma_bprobe,DIM=1)
            IF (target_bprobe(ik) /= 0.0) THEN
               WRITE(iunit,"(2(2X,A,I3.3,A,1X,'=',1X,E22.14))")&
                  'TARGET_BPROBE(',ik,')',target_bprobe(ik),&
                  'SIGMA_BPROBE(',ik,')',sigma_bprobe(ik)
            END IF
         END DO
         DO ik = 1, UBOUND(sigma_fluxloop,DIM=1)
            IF (target_fluxloop(ik) /= 0.0) THEN
               WRITE(iunit,"(2(2X,A,I3.3,A,1X,'=',1X,E22.14))")&
                  'TARGET_FLUXLOOP(',ik,')',target_fluxloop(ik),&
                  'SIGMA_FLUXLOOP(',ik,')',sigma_fluxloop(ik)
            END IF
         END DO
         DO ik = 1, UBOUND(sigma_segrog,DIM=1)
            IF (target_segrog(ik) /= 0.0) THEN
               WRITE(iunit,"(2(2X,A,I3.3,A,1X,'=',1X,E22.14))")&
                  'TARGET_SEGROG(',ik,')',target_segrog(ik),&
                  'SIGMA_SEGROG(',ik,')',sigma_segrog(ik)
            END IF
         END DO
      END IF
      IF (sigma_vessel < bigno) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          VACCUM VESSEL LIMITER'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         IF (LEN_TRIM(vessel_string) > 1) WRITE(iunit,outstr) 'VESSEL_STRING',TRIM(vessel_string)
         WRITE(iunit,outflt) 'TARGET_VESSEL',target_vessel
         WRITE(iunit,outflt) 'SIGMA_VESSEL',sigma_vessel
      END IF
      IF (ANY(sigma_separatrix < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          TARGET SEPARATRIX'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO u = 1, nu_max
            DO v = 1, nv_max
               IF (sigma_separatrix(u,v) < bigno) &
                  WRITE(iunit,"(5(2X,A,I3.3,',',I3.3,A,E22.14))") &
                  'R_SEPARATRIX(',u,v,') = ',r_separatrix(u,v),&
                  'PHI_SEPARATRIX(',u,v,') = ',phi_separatrix(u,v),&
                  'Z_SEPARATRIX(',u,v,') = ',z_separatrix(u,v),&
                  'TARGET_SEPARATRIX(',u,v,') = ',target_separatrix(u,v),&
                  'SIGMA_SEPARATRIX(',u,v,') = ',sigma_separatrix(u,v)
            END DO
         END DO
      END IF
      IF (ANY(sigma_limiter < bigno)) THEN
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         WRITE(iunit,'(A)') '!          TARGET LIMITER'
         WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
         DO u = 1, nu_max
            DO v = 1, nv_max
               IF (sigma_limiter(u,v) < bigno) &
                  WRITE(iunit,"(5(2X,A,I3.3,',',I3.3,A,E22.14))") &
                  'R_LIMITER(',u,v,') = ',r_limiter(u,v),&
                  'PHI_LIMITER(',u,v,') = ',phi_limiter(u,v),&
                  'Z_LIMITER(',u,v,') = ',z_limiter(u,v),&
                  'TARGET_LIMITER(',u,v,') = ',target_limiter(u,v),&
                  'SIGMA_LIMITER(',u,v,') = ',sigma_limiter(u,v)
            END DO
         END DO
      END IF
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_optimum_namelist
      
      END MODULE stellopt_input_mod
