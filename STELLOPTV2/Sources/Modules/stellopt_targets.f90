!-----------------------------------------------------------------------
!     Module:        stellopt_targets
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/25/2012
!     Description:   This module contains the STELLOPT target variables.
!-----------------------------------------------------------------------
      MODULE stellopt_targets       
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vparams, ONLY: nsd

!-----------------------------------------------------------------------
!     Module Variables
!            mboz               Poloidal Fourier Content for Boozer
!            nboz               Toroidal Fourier Content for Boozer
!            target_phiedge     Target for PHIEDGE
!            target_curtor      Target for CURTOR      
!            target_curtor_max  Target for upper limit on CURTOR      
!            target_volume	     Target for volume
!            target_beta        Target for Beta     
!            target_rbtor	     Target for R*Btor
!            target_r0   	     Target for Raxis (phi=0)
!            target_z0   	     Target for Zaxis (phi=0)
!            target_wp          Target for stored energy
!            target_aspect	     Target for aspect ratio
!            target_aspect_max  Target for upper limit on aspect ratio
!            target_bprobe      Target B-Field probe array
!            target_linene      Target line_integrated density array
!            target_press	     Target array for pressure
!            target_te          Target array for electron temperature
!            target_ne          Target array for electron density
!            target_ti          Target array for ion temperature
!            r_te               R electron temperature location array
!            z_te               Z electron temperature location array
!            phi_te             PHI electron temperature location array
!            s_te               s electron temperature location array
!            r_ne               R electron density location array
!            z_ne               Z electron density location array
!            phi_ne             PHI electron density location array
!            s_ne               s electron denstiy location array
!            r_ti               R ion temperature location array
!            z_ti               Z ion temperature location array
!            phi_ti             PHI ion temperature location array
!            s_ti               s ion temperature location array
!            r_iota             R Rotational Transform location array
!            z_iota             Z Rotational Transform location array
!            phi_iota           PHI Rotational Transform location array
!            s_iota             s Rotational Transform location array
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL     ::  lneed_magdiag
      LOGICAL, DIMENSION(nsd)  :: lbooz
      INTEGER     ::  mboz, nboz, NumJstar
      INTEGER, PARAMETER :: nprof = 512
      INTEGER, PARAMETER :: nprobes = 2048
      INTEGER, PARAMETER :: nu_max = 361
      INTEGER, PARAMETER :: nv_max = 361
      INTEGER, PARAMETER :: nsys   = 16
      INTEGER, PARAMETER :: npart_max = 16384
      REAL(rprec) ::  target_phiedge, sigma_phiedge
      REAL(rprec) ::  target_curtor, sigma_curtor
      REAL(rprec) ::  target_curtor_max, sigma_curtor_max
      REAL(rprec) ::  target_volume, sigma_volume
      REAL(rprec) ::  target_beta, sigma_beta
      REAL(rprec) ::  target_betator, sigma_betator
      REAL(rprec) ::  target_betapol, sigma_betapol
      REAL(rprec) ::  target_wp, sigma_wp
      REAL(rprec) ::  target_aspect, sigma_aspect
      REAL(rprec) ::  target_rbtor, sigma_rbtor
      REAL(rprec) ::  target_r0, sigma_r0
      REAL(rprec) ::  target_z0, sigma_z0
      REAL(rprec) ::  target_aspect_max, sigma_aspect_max, width_aspect_max
      REAL(rprec) ::  target_gradp_max, sigma_gradp_max, width_gradp_max
      REAL(rprec) ::  target_pmin, sigma_pmin, width_pmin
      REAL(rprec) ::  target_curvature, sigma_curvature
      REAL(rprec) ::  target_kappa, sigma_kappa
      REAL(rprec), PARAMETER ::  bigno_ne = 1.0E27
      REAL(rprec) ::  norm_press
      REAL(rprec) ::  qm_ratio
      REAL(rprec) ::  cutoff_te_line
      REAL(rprec), DIMENSION(nprof) ::  target_press, sigma_press, &
                                        r_press, z_press, phi_press, s_press  
      REAL(rprec), DIMENSION(nprof) ::  target_te, sigma_te, &
                                        r_te, z_te, phi_te, s_te  
      REAL(rprec), DIMENSION(nprof) ::  target_ne, sigma_ne, &
                                        r_ne, z_ne, phi_ne, s_ne  
      REAL(rprec), DIMENSION(nprof) ::  target_ne_line,sigma_ne_line, &
                                        r0_ne_line, phi0_ne_line, z0_ne_line, &
                                        r1_ne_line, phi1_ne_line, z1_ne_line
      REAL(rprec), DIMENSION(nprof) ::  target_te_line,sigma_te_line, &
                                        r0_te_line, phi0_te_line, z0_te_line, &
                                        r1_te_line, phi1_te_line, z1_te_line
      REAL(rprec), DIMENSION(nprof) ::  target_ti_line,sigma_ti_line, &
                                        r0_ti_line, phi0_ti_line, z0_ti_line, &
                                        r1_ti_line, phi1_ti_line, z1_ti_line
      REAL(rprec), DIMENSION(nprof) ::  target_faraday,sigma_faraday, &
                                        r0_faraday, phi0_faraday, z0_faraday, &
                                        r1_faraday, phi1_faraday, z1_faraday
      REAL(rprec), DIMENSION(nprof) ::  target_sxr,sigma_sxr, &
                                        r0_sxr, phi0_sxr, z0_sxr, &
                                        r1_sxr, phi1_sxr, z1_sxr
      REAL(rprec), DIMENSION(nprof) ::  target_ti, sigma_ti, &
                                        r_ti, z_ti, phi_ti, s_ti   
      REAL(rprec), DIMENSION(nprof) ::  target_vphi, sigma_vphi, &
                                        r_vphi, z_vphi, phi_vphi, s_vphi   
      REAL(rprec), DIMENSION(nprof) ::  target_iota, sigma_iota, &
                                        r_iota, z_iota, phi_iota, s_iota   
      REAL(rprec), DIMENSION(nprof) ::  target_vaciota, sigma_vaciota, &
                                        r_vaciota, z_vaciota, phi_vaciota, s_vaciota   
      REAL(rprec), DIMENSION(nprof) ::  target_mse, sigma_mse, &
                                        r_mse, z_mse, phi_mse, s_mse,&
                                        a1_mse, a2_mse, a3_mse, a4_mse,&
                                        a5_mse, a6_mse, a7_mse, vac_mse
      LOGICAL,     DIMENSION(nprof) ::  lmse_extcur
      REAL(rprec), DIMENSION(nprobes)  ::  target_bprobe, sigma_bprobe       ! Note this number is hardcoded in chisq_brobes SAL 2/10/14
      REAL(rprec), DIMENSION(nprof) ::  target_segrog, sigma_segrog, &
                                        target_fluxloop, sigma_fluxloop
      CHARACTER(256)                ::  magdiag_coil
      REAL(rprec), DIMENSION(nprof) ::  target_extcur, sigma_extcur 
      CHARACTER(256)                ::  vessel_string
      REAL(rprec)                   ::  target_vessel, sigma_vessel
      REAL(rprec), DIMENSION(nprof) ::  balloon_theta, balloon_zeta
      REAL(rprec), DIMENSION(nsd)   ::  target_bmin, sigma_bmin
      REAL(rprec), DIMENSION(nsd)   ::  target_bmax, sigma_bmax
      REAL(rprec), DIMENSION(nsd)   ::  target_jcurv, sigma_jcurv
      REAL(rprec), DIMENSION(nsd)   ::  target_jdotb, sigma_jdotb
      REAL(rprec), DIMENSION(nsd)   ::  target_balloon, sigma_balloon
      REAL(rprec), DIMENSION(nsd)   ::  target_bootstrap, sigma_bootstrap
      REAL(rprec), DIMENSION(nsd)   ::  target_neo, sigma_neo
      REAL(rprec), DIMENSION(nsd)   ::  target_Jstar, sigma_Jstar
      REAL(rprec), DIMENSION(nsd)   ::  target_helicity, sigma_helicity
      COMPLEX                       ::  helicity
      REAL(rprec), DIMENSION(nsd)   ::  target_resjac, sigma_resjac, &
                                        xm_resjac, xn_resjac
      LOGICAL                       ::  lglobal_txport
      INTEGER                       ::  nz_txport, nalpha_txport
      REAL(rprec)                   ::  alpha0_txport
      REAL(rprec), DIMENSION(nsd)   ::  target_txport, sigma_txport, &
                                        s_txport
      CHARACTER(256)                ::  txport_proxy
      REAL(rprec), DIMENSION(nsd)   ::  target_DKES, sigma_DKES, nu_DKES
      REAL(rprec), DIMENSION(nu_max,nv_max) ::  target_separatrix, sigma_separatrix, &
                                                r_separatrix, z_separatrix, phi_separatrix
      REAL(rprec), DIMENSION(nu_max,nv_max) ::  target_limiter, sigma_limiter, &
                                                r_limiter, z_limiter, phi_limiter

      INTEGER                            :: nu_orbit, nv_orbit, np_orbit
      REAL(rprec)                        :: mass_orbit, Z_orbit
      REAL(rprec), DIMENSION(nsd)        ::  target_orbit, sigma_orbit
      REAL(rprec), DIMENSION(npart_max)  :: vll_orbit, mu_orbit, vperp_orbit
      REAL(rprec), DIMENSION(nsys) ::  target_kink, sigma_kink
      INTEGER, DIMENSION(nsys) :: mlmns_kink, lssl_kink, lssd_kink,nj_kink,nk_kink
      INTEGER     ::  mlmnb_kink, ivac_kink, mmaxdf_kink, nmaxdf_kink
      REAL(rprec), DIMENSION(nsys,nprof) :: target_ece, sigma_ece,&
                                            freq_ece
      INTEGER                            :: nra_ece, nphi_ece
      REAL(rprec), DIMENSION(nsys,3)     :: antennaposition_ece, targetposition_ece,rbeam_ece,rfocus_ece
      CHARACTER(256)                     :: vessel_ece,mirror_ece,targettype_ece,antennatype_ece
      
      INTEGER     ::  numws
      REAL(rprec) ::  target_coil_bnorm, sigma_coil_bnorm
      INTEGER     ::  nu_bnorm,nv_bnorm
                                        
      INTEGER, PARAMETER :: jtarget_aspect     = 100
      INTEGER, PARAMETER :: jtarget_rbtor      = 1001
      INTEGER, PARAMETER :: jtarget_r0         = 1002
      INTEGER, PARAMETER :: jtarget_z0         = 1003
      INTEGER, PARAMETER :: jtarget_curvature  = 1004
      INTEGER, PARAMETER :: jtarget_kappa      = 1005
      INTEGER, PARAMETER :: jtarget_beta       = 101
      INTEGER, PARAMETER :: jtarget_betapol    = 1011
      INTEGER, PARAMETER :: jtarget_betator    = 1012
      INTEGER, PARAMETER :: jtarget_curtor     = 102
      INTEGER, PARAMETER :: jtarget_curtor_max = 1021
      INTEGER, PARAMETER :: jtarget_phiedge    = 103
      INTEGER, PARAMETER :: jtarget_volume     = 104
      INTEGER, PARAMETER :: jtarget_wp         = 105
      INTEGER, PARAMETER :: jtarget_aspect_max = 106
      INTEGER, PARAMETER :: jtarget_gradp_max  = 107
      INTEGER, PARAMETER :: jtarget_pmin       = 108
      INTEGER, PARAMETER :: jtarget_extcur     = 109
      INTEGER, PARAMETER :: jtarget_vessel     = 110
      INTEGER, PARAMETER :: jtarget_separatrix = 111
      INTEGER, PARAMETER :: jtarget_limiter    = 112
      INTEGER, PARAMETER :: jtarget_ne         = 200
      INTEGER, PARAMETER :: jtarget_line_ne    = 2001
      INTEGER, PARAMETER :: jtarget_te         = 201
      INTEGER, PARAMETER :: jtarget_line_te    = 2011
      INTEGER, PARAMETER :: jtarget_ti         = 202
      INTEGER, PARAMETER :: jtarget_line_ti    = 2021
      INTEGER, PARAMETER :: jtarget_vphi       = 2022
      INTEGER, PARAMETER :: jtarget_press      = 203
      INTEGER, PARAMETER :: jtarget_iota       = 300  
      INTEGER, PARAMETER :: jtarget_iprime     = 301
      INTEGER, PARAMETER :: jtarget_vaciota    = 302  
      INTEGER, PARAMETER :: jtarget_mse        = 401
      INTEGER, PARAMETER :: jtarget_faraday    = 402
      INTEGER, PARAMETER :: jtarget_sxr        = 403
      INTEGER, PARAMETER :: jtarget_ece        = 404
      INTEGER, PARAMETER :: jtarget_bprobe     = 501
      INTEGER, PARAMETER :: jtarget_segrog     = 502
      INTEGER, PARAMETER :: jtarget_fluxloop   = 503
      INTEGER, PARAMETER :: jtarget_balloon    = 601
      INTEGER, PARAMETER :: jtarget_kink       = 6011
      INTEGER, PARAMETER :: jtarget_bootstrap  = 602
      INTEGER, PARAMETER :: jtarget_neo        = 603
      INTEGER, PARAMETER :: jtarget_Jstar      = 604
      INTEGER, PARAMETER :: jtarget_helicity   = 605
      INTEGER, PARAMETER :: jtarget_resjac     = 606
      INTEGER, PARAMETER :: jtarget_txport     = 607
      INTEGER, PARAMETER :: jtarget_dkes       = 608
      INTEGER, PARAMETER :: jtarget_jdotb      = 609
      INTEGER, PARAMETER :: jtarget_jcurv      = 6091
      INTEGER, PARAMETER :: jtarget_bmin       = 610
      INTEGER, PARAMETER :: jtarget_bmax       = 611
      INTEGER, PARAMETER :: jtarget_orbit      = 612
      INTEGER, PARAMETER :: jtarget_coil_bnorm = 613
      
      
      CONTAINS
      
      SUBROUTINE write_targets(iunit,var_num)
      INTEGER, INTENT(in) :: iunit, var_num
      CHARACTER*(*), PARAMETER ::  out_format = '(5X,A)'
      CHARACTER*(*), PARAMETER ::  out_format_1D = '(5X,A,I3.3,A)'
      CHARACTER*(*), PARAMETER ::  out_format_2D = '(5X,A,I3.3,A,I3.3,A)'
      SELECT CASE(var_num)
         CASE(jtarget_aspect)
            WRITE(iunit, out_format) 'Aspect Ratio'
         CASE(jtarget_aspect_max)
            WRITE(iunit, out_format) 'Max Aspect Ratio (upper limit)'
         CASE(jtarget_beta)
            WRITE(iunit, out_format) 'Plasma Beta'
         CASE(jtarget_betapol)
            WRITE(iunit, out_format) 'Plasma Beta (Poloidal)'
         CASE(jtarget_betator)
            WRITE(iunit, out_format) 'Plasma Beta (Toroidal)'
         CASE(jtarget_curvature)
            WRITE(iunit, out_format) 'Boundary Curvature (kertosis)'
         CASE(jtarget_kappa)
            WRITE(iunit, out_format) 'Boundary Elongation'
         CASE(jtarget_curtor)
            WRITE(iunit, out_format) 'Net Toroidal Current'
         CASE(jtarget_curtor_max)
            WRITE(iunit, out_format) 'Net Toroidal Current (max)'
         CASE(jtarget_phiedge)
            WRITE(iunit, out_format) 'Total Enclosed Toroidal Flux'
         CASE(jtarget_volume)
            WRITE(iunit, out_format) 'Plasma Volume'
         CASE(jtarget_wp)
            WRITE(iunit, out_format) 'Plasma Stored Energy'
         CASE(jtarget_gradp_max)
            WRITE(iunit, out_format) 'Max Pressure Gradient (upper limit)'
         CASE(jtarget_pmin)
            WRITE(iunit, out_format) 'Min Pressure'
         CASE(jtarget_rbtor)
            WRITE(iunit, out_format) 'R*Btor'
         CASE(jtarget_r0)
            WRITE(iunit, out_format) 'R0 (phi=0)'
         CASE(jtarget_z0)
            WRITE(iunit, out_format) 'Z0 (phi=0)'
         CASE(jtarget_extcur)
            WRITE(iunit, out_format) 'External currents'
         CASE(jtarget_press)
            WRITE(iunit, out_format) 'Plasma Pressure'
         CASE(jtarget_ne)
            WRITE(iunit, out_format) 'Electron Density'
         CASE(jtarget_line_ne)
            WRITE(iunit, out_format) 'Line Integrated Electron Density'
         CASE(jtarget_line_te)
            WRITE(iunit, out_format) 'Line Integrated Electron Temperature'
         CASE(jtarget_line_ti)
            WRITE(iunit, out_format) 'Line Integrated Ion Temperature'
         CASE(jtarget_te)
            WRITE(iunit, out_format) 'Electron Temperature'
         CASE(jtarget_ti)
            WRITE(iunit, out_format) 'Ion Temperature'
         CASE(jtarget_vphi)
            WRITE(iunit, out_format) 'Toroidal Rotation'
         CASE(jtarget_iota)
            WRITE(iunit, out_format) 'Rotational Transform'
         CASE(jtarget_vaciota)
            WRITE(iunit, out_format) 'Vacuum Rotational Transform'
         CASE(jtarget_iprime)
            WRITE(iunit, out_format) 'Current Profile (I'')'
         CASE(jtarget_mse)
            WRITE(iunit, out_format) 'Motional Stark Effect Diagnostic'
         CASE(jtarget_faraday)
            WRITE(iunit, out_format) 'Faraday Rotation'
         CASE(jtarget_sxr)
            WRITE(iunit, out_format) 'Soft X-Ray'
         CASE(jtarget_ece)
            WRITE(iunit, out_format) 'ECE Reflectometry Diagnostic'
         CASE(jtarget_bprobe)
            WRITE(iunit, out_format) 'Magnetic Field Probe'
         CASE(jtarget_fluxloop)
            WRITE(iunit, out_format) 'Fluxloop'
         CASE(jtarget_segrog)
            WRITE(iunit, out_format) 'Rogowski Coil'
         CASE(jtarget_balloon)
            WRITE(iunit, out_format) 'Ballooning Stability'
         CASE(jtarget_kink)
            WRITE(iunit, out_format) 'Kink Stability'
         CASE(jtarget_bootstrap)
            WRITE(iunit, out_format) 'Bootstrap Current'
         CASE(jtarget_neo)
            WRITE(iunit, out_format) 'Neoclassical Transport'
         CASE(jtarget_Jstar)
            WRITE(iunit, out_format) 'Trapped Particle J*'
         CASE(jtarget_helicity)
            WRITE(iunit, out_format) 'Boozer Spectrum Helicity'
         CASE(jtarget_txport)
            WRITE(iunit, out_format) 'Turbulent Transport'
         CASE(jtarget_orbit)
            WRITE(iunit, out_format) 'Particle Orbits (BEAMS3D)'
         CASE(jtarget_dkes)
            WRITE(iunit, out_format) 'Drift-Kinetics (DKES)'
         CASE(jtarget_jdotb)
            WRITE(iunit, out_format) '<J.B>'
         CASE(jtarget_jcurv)
            WRITE(iunit, out_format) 'Toroidal Curent'
         CASE(jtarget_bmin)
            WRITE(iunit, out_format) '|B| Minimum'
         CASE(jtarget_bmax)
            WRITE(iunit, out_format) '|B| Maximum'
         CASE(jtarget_resjac)
            WRITE(iunit, out_format) 'Boozer Resonant Modes'
         CASE(jtarget_separatrix)
            WRITE(iunit, out_format) 'Separatrix'
         CASE(jtarget_limiter)
            WRITE(iunit, out_format) 'Limiter'
         CASE(jtarget_coil_bnorm)
            WRITE(iunit, out_format) 'COILOPT++ Normal Field'
      END SELECT
      END SUBROUTINE write_targets
      
      
      
      END MODULE stellopt_targets
