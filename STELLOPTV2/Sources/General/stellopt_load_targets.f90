!-----------------------------------------------------------------------
!     Subroutine:    stellopt_load_targets
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine calls the various chi-squared
!                    functions.  It is important to note that on the
!                    first call to a chi-squared function only the
!                    number of targets should be returned.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_load_targets(m, fvec, iflag, ncnt)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_targets
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        m       Number of function dimensions
!        fvec    Output array of function values
!        iflag   Processor number
!        ncnt    Current function evaluation
!----------------------------------------------------------------------
      INTEGER, INTENT(in)      :: ncnt
      INTEGER, INTENT(inout)   :: m,iflag
      REAL(rprec), INTENT(out) :: fvec(m)
      
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, iunit,m_sav

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      mtargets=0
      !------------- SCALAR TARGETS ----------------------------
      ! PHIEDGE
      IF (sigma_phiedge < bigno)  &
         CALL chisq_phiedge(target_phiedge,sigma_phiedge,ncnt,iflag)
      ! CURTOR
      IF (sigma_curtor < bigno)  &
         CALL chisq_curtor(target_curtor,sigma_curtor,ncnt,iflag)
      ! CURTOR (MAX)
      IF (sigma_curtor_max < bigno)  &
         CALL chisq_curtor_max(target_curtor_max,sigma_curtor_max,ncnt,iflag)
      ! RBtor
      IF (sigma_rbtor < bigno)  &
         CALL chisq_rbtor(target_rbtor,sigma_rbtor,ncnt,iflag)
      ! R0
      IF (sigma_r0 < bigno)  &
         CALL chisq_r0(target_r0,sigma_r0,ncnt,iflag)
      ! Z0
      IF (sigma_z0 < bigno)  &
         CALL chisq_z0(target_z0,sigma_z0,ncnt,iflag)
      ! B0
      IF (sigma_b0 < bigno)  &
         CALL chisq_b0(target_b0,sigma_b0,ncnt,iflag)
      ! VOLUME
      IF (sigma_volume < bigno)  &
         CALL chisq_volume(target_volume,sigma_volume,ncnt,iflag)
      ! BETA (total)
      IF (sigma_beta < bigno)  &
         CALL chisq_beta(target_beta,sigma_beta,ncnt,iflag)
      ! BETA (poloidal)
      IF (sigma_betapol < bigno)  &
         CALL chisq_betapol(target_betapol,sigma_betapol,ncnt,iflag)
      ! BETA (toroidal)
      IF (sigma_betator < bigno)  &
         CALL chisq_betator(target_betator,sigma_betator,ncnt,iflag)
      ! STORED ENERGY
      IF (sigma_wp < bigno)  &
         CALL chisq_wp(target_wp,sigma_wp,ncnt,iflag)
      ! ASPECT RATIO
      IF (sigma_aspect < bigno)  &
         CALL chisq_aspect(target_aspect,sigma_aspect,ncnt,iflag)
      ! CURVATURE (kertosis)
      IF (sigma_curvature < bigno)  &
         CALL chisq_curvature(target_curvature,sigma_curvature,ncnt,iflag)
      ! KAPPA (ellipticity)
      IF (sigma_kappa < bigno)  &
         CALL chisq_kappa(target_kappa,sigma_kappa,ncnt,iflag)
      ! KAPPA (ellipticity box)
      IF (sigma_kappa_box < bigno)  &
         CALL chisq_kappa_box(target_kappa_box,sigma_kappa_box,ncnt,iflag)
      ! KAPPA (ellipticity avg)
      IF (sigma_kappa_avg < bigno)  &
         CALL chisq_kappa_avg(target_kappa_avg,sigma_kappa_avg,ncnt,iflag)
      ! ASPECT RATIO (MAX)
      IF (sigma_aspect_max < bigno)  &
         CALL chisq_aspect(target_aspect_max,sigma_aspect_max,ncnt,iflag)
      ! PRESSURE (MIN)
      IF (sigma_pmin < bigno)  &
         CALL chisq_pmin(target_pmin,sigma_pmin,ncnt,iflag)
         
      !------------- ARRAY TARGETS ----------------------------
      ! EXTERNAL CURRENTS
      IF (ANY(sigma_extcur < bigno))  &
         CALL chisq_extcur(target_extcur,sigma_extcur,ncnt,iflag)
         
      !------------- LINE INTEGRATED TARGETS -------------------
      ! LINE ELECTRON DENSITY
      IF (ANY(sigma_ne_line < bigno_ne)) &
         CALL chisq_line_ne(target_ne_line, sigma_ne_line, ncnt,iflag)
      ! LINE ELECTRON TEMPERATURE
      IF (ANY(sigma_te_line < bigno)) &
         CALL chisq_line_te(target_te_line, sigma_te_line, ncnt,iflag)
      ! LINE ION TEMPERATURE
      IF (ANY(sigma_ti_line < bigno)) &
         CALL chisq_line_ti(target_ti_line, sigma_ti_line, ncnt,iflag)
      ! XICS Brightness
      IF (ANY(sigma_xics_bright < bigno)) &
         CALL chisq_xics_bright(target_xics_bright, sigma_xics_bright, ncnt,iflag)
      ! XICS
      IF (ANY(sigma_xics < bigno)) &
         CALL chisq_xics(target_xics, sigma_xics, ncnt,iflag)
      ! SOFT X-RAYS
      IF (ANY(sigma_sxr < bigno)) &
         CALL chisq_sxr(target_sxr, sigma_sxr, ncnt,iflag)
      ! FARADAY ROTATION
      IF (ANY(sigma_faraday < bigno_ne)) &
         CALL chisq_faraday(target_faraday, sigma_faraday, ncnt,iflag)
         
      !------------- OTHER TARGETS -------------------
      !  ECE Reflectometry
      IF (ANY(sigma_ece < bigno)) &
         CALL chisq_ece(target_ece, sigma_ece, ncnt,iflag)
         
      !------------- PROFILE TARGETS ---------------------------
      ! PRESSURE PROFILE
      IF (ANY(sigma_press < bigno)) &
         CALL chisq_press(target_press, sigma_press, ncnt,iflag)
      ! ELECTRON DENSITY
      IF (ANY(sigma_ne < bigno_ne)) &
         CALL chisq_ne(target_ne, sigma_ne, ncnt,iflag)
      ! ELECTRON TEMPERATURE
      IF (ANY(sigma_te < bigno)) &
         CALL chisq_te(target_te, sigma_te, ncnt,iflag)
      ! ION TEMPERATURE
      IF (ANY(sigma_ti < bigno)) &
         CALL chisq_ti(target_ti, sigma_ti, ncnt,iflag)
      ! ION TEMPERATURE
      IF (ANY(sigma_vphi < bigno)) &
         CALL chisq_vphi(target_vphi, sigma_vphi, ncnt,iflag)
      ! ROTATIONAL TRANSFORM
      IF (ANY(sigma_iota < bigno)) &
         CALL chisq_iota(target_iota, sigma_iota, ncnt,iflag)
      ! Vacuum ROTATIONAL TRANSFORM
      IF (ANY(sigma_vaciota < bigno)) &
         CALL chisq_vaciota(target_vaciota, sigma_vaciota, ncnt,iflag)
      ! Parallel Current <J.B>
      IF (ANY(sigma_jdotb < bigno)) &
         CALL chisq_jdotb(target_jdotb, sigma_jdotb, ncnt,iflag)
      ! Magnetic Well
      IF (ANY(sigma_magwell < bigno)) &
         CALL chisq_magwell(target_magwell, sigma_magwell, ncnt,iflag)
      ! B_min
      IF (ANY(sigma_bmin < bigno)) &
         CALL chisq_bmin(target_bmin, sigma_bmin, ncnt,iflag)
      ! B_max
      IF (ANY(sigma_bmax < bigno)) &
         CALL chisq_bmax(target_bmax, sigma_bmax, ncnt,iflag)
      ! MSE Diagnostic
      IF (ANY(sigma_mse < bigno)) &
         CALL chisq_mse(target_mse, sigma_mse, ncnt,iflag)
      ! CURRENT DENSITY <JCURV>
      IF (ANY(sigma_jcurv < bigno)) &
         CALL chisq_jcurv(target_jcurv, sigma_jcurv, ncnt,iflag)
         
      !------------- COIL GEOMETRY TARGETS ---------------------
      ! Coil lengths
      IF (ANY(sigma_coillen < bigno)) &
         CALL chisq_coillen(target_coillen, sigma_coillen, ncnt, iflag)
      IF (sigma_coilsep < bigno) &
         CALL chisq_coilsep(target_coilsep, sigma_coilsep, ncnt, iflag)
      IF (ANY(sigma_coilcrv < bigno)) &
         CALL chisq_coilcrv(target_coilcrv, sigma_coilcrv, ncnt, iflag)
      IF (ANY(sigma_coilself < bigno)) &
         CALL chisq_coilself(target_coilself, sigma_coilself, ncnt, iflag)

      !------------- EXTERNAL TARGETS --------------------------
      !  This section of the code relys upon external libraries
      !  for calculation of target parameters
      ! B-Probes
      IF (ANY(sigma_bprobe < bigno)) &
         CALL chisq_bprobes(target_bprobe, sigma_bprobe, ncnt,iflag)
      ! Fluxloops
      IF (ANY(sigma_fluxloop < bigno)) &
         CALL chisq_fluxloops(target_fluxloop, sigma_fluxloop, ncnt,iflag)
      ! Rogowski
      IF (ANY(sigma_segrog < bigno)) &
         CALL chisq_segrog(target_segrog, sigma_segrog, ncnt,iflag)
      ! Vessel (limiting)
      IF (sigma_vessel < bigno) &
         CALL chisq_vessel(target_vessel, sigma_vessel, ncnt,iflag)
      ! Separatrix
      IF (ANY(sigma_separatrix < bigno)) &
         CALL chisq_separatrix(target_separatrix, sigma_separatrix, ncnt,iflag)
      ! Limiter
      IF (ANY(sigma_limiter < bigno)) &
         CALL chisq_limiter(target_limiter, sigma_limiter, ncnt,iflag)
      ! Ballooning
      IF (ANY(sigma_balloon < bigno)) &
         CALL chisq_balloon(target_balloon, sigma_balloon, ncnt,iflag)
      ! Boostrap
      IF (ANY(sigma_bootstrap < bigno)) &
         CALL chisq_bootstrap(target_bootstrap, sigma_bootstrap, ncnt,iflag)
      ! NEO
      IF (ANY(sigma_neo < bigno)) &
         CALL chisq_neo(target_neo, sigma_neo, ncnt,iflag)
      ! DKES
      IF (ANY(sigma_dkes < bigno)) &
         CALL chisq_dkes(target_dkes, sigma_dkes, ncnt,iflag)
      ! TXPORT
      IF (ANY(sigma_txport < bigno)) &
         CALL chisq_txport(target_txport, sigma_txport, ncnt,iflag)
      ! Orbit
      IF (ANY(sigma_orbit < bigno)) &
         CALL chisq_orbit(target_orbit, sigma_orbit, ncnt,iflag)
      ! |Bmn| Helicity
      IF (ANY(sigma_helicity < bigno)) &
         CALL chisq_helicity(target_helicity, sigma_helicity, ncnt,iflag)
      ! |Bmn| Helicity (OLD)
      IF (ANY(sigma_helicity_old < bigno)) &
         CALL chisq_helicity_ornl(target_helicity_old, sigma_helicity_old, ncnt,iflag)
      ! J*
      IF (ANY(sigma_Jstar < bigno)) &
         CALL chisq_jstar(target_Jstar, sigma_Jstar, ncnt,iflag)
      ! Resonant Jacobian
      IF (ANY(sigma_resjac < bigno)) &
         CALL chisq_resjac(target_resjac, sigma_resjac, ncnt,iflag)
      ! Coil Optimization
      IF (sigma_coil_bnorm < bigno) &
         CALL chisq_coil_bnorm(target_coil_bnorm, sigma_coil_bnorm, ncnt,iflag)
      ! REGCOIL Coil Optimization (CHI2_B targets)
      IF (ANY(sigma_regcoil_chi2_b < bigno)) THEN
         CALL chisq_regcoil_chi2_b(target_regcoil_chi2_b, sigma_regcoil_chi2_b, ncnt,iflag)
      END IF
      ! Kink
      IF (ANY(sigma_kink < bigno)) &
         CALL chisq_kink(target_kink, sigma_kink, ncnt,iflag)

      ! Return if an initialization call
      IF (ncnt < 0) RETURN
      
      ! Check some stuff
      IF (mtargets .ne. m) THEN; iflag=-2; RETURN; END IF
      
      ! Calculate fvec
      !PRINT *,m,mtargets,fvec
      fvec(1:m) = (vals(1:m)-targets(1:m))/ABS(sigmas(1:m))
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_load_targets
