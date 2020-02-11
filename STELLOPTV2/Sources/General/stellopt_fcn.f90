!-----------------------------------------------------------------------
!     Subroutine:    stellopt_fcn
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine calculates the function which is
!                    minimized by STELLOPT.  Originally developed for
!                    the lmdif function.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_fcn(m, n, x, fvec, iflag, ncnt)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE stellopt_targets
      USE equil_utils, ONLY: eval_prof_spline
      USE vmec_input
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
!                             animec_flag, flow_flag
      USE vmec_main, ONLY:  multi_ns_grid
      USE read_wout_mod, ONLY: read_wout_file, write_wout_file, read_wout_deallocate
      USE mpi_params                                                    ! MPI
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        m       Number of function dimensions
!        n       Number of function variables
!        x       Vector of function variables
!        fvec    Output array of function values
!        iflag   Processor number
!        ncnt    Current function evaluation
!----------------------------------------------------------------------
      INTEGER, INTENT(in)      ::  m, n, ncnt
      INTEGER, INTENT(inout)   :: iflag
      REAL(rprec), INTENT(inout)  :: x(n)
      REAL(rprec), INTENT(out) :: fvec(m)
      
!-----------------------------------------------------------------------
!     Local Variables
!        iflag       Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL ::  lscreen
      INTEGER ::  nvar_in, dex, dex2, ik, istat, iunit, pass, mf,nf
      INTEGER ::  vctrl_array(5)
      REAL(rprec) :: norm_aphi, norm_am, norm_ac, norm_ai, norm_ah,&
                     norm_at, norm_ne, norm_te, norm_ti, norm_th, &
                     norm_phi, norm_zeff, norm_emis_xics, &
                     norm_beamj, norm_bootj, temp
      INTEGER, PARAMETER     :: max_refit = 2
      REAL(rprec), PARAMETER :: ec  = 1.60217653D-19
      CHARACTER(len = 16)     :: temp_str
      CHARACTER(len = 128)    :: reset_string
      CHARACTER(len = 256)    :: ctemp_str
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Load variables first
      norm_aphi = 1; norm_am = 1; norm_ac = 1; norm_ai = 1
      norm_ah   = 1; norm_at = 1; norm_phi = 1; norm_zeff = 1
      norm_ne   = 1; norm_te = 1; norm_ti  = 1; norm_th = 1
      norm_beamj = 1; norm_bootj = 1; norm_emis_xics = 1

      ! Save variables
      DO nvar_in = 1, n
         IF (var_dex(nvar_in) == iaphi .and. arr_dex(nvar_in,2) == norm_dex) norm_aphi = x(nvar_in)
         IF (var_dex(nvar_in) == iam .and. arr_dex(nvar_in,2) == norm_dex) norm_am = x(nvar_in)
         IF (var_dex(nvar_in) == iac .and. arr_dex(nvar_in,2) == norm_dex) norm_ac = x(nvar_in)
         IF (var_dex(nvar_in) == iai .and. arr_dex(nvar_in,2) == norm_dex) norm_ai = x(nvar_in)
         IF (var_dex(nvar_in) == iah .and. arr_dex(nvar_in,2) == norm_dex) norm_ah = x(nvar_in)
         IF (var_dex(nvar_in) == iat .and. arr_dex(nvar_in,2) == norm_dex) norm_at = x(nvar_in)
         IF (var_dex(nvar_in) == ine .and. arr_dex(nvar_in,2) == norm_dex) norm_ne = x(nvar_in)
         IF (var_dex(nvar_in) == izeff .and. arr_dex(nvar_in,2) == norm_dex) norm_zeff = x(nvar_in)
         IF (var_dex(nvar_in) == ite .and. arr_dex(nvar_in,2) == norm_dex) norm_te = x(nvar_in)
         IF (var_dex(nvar_in) == iti .and. arr_dex(nvar_in,2) == norm_dex) norm_ti = x(nvar_in)
         IF (var_dex(nvar_in) == ith .and. arr_dex(nvar_in,2) == norm_dex) norm_th = x(nvar_in)
         IF (var_dex(nvar_in) == iam_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_am = x(nvar_in)
         IF (var_dex(nvar_in) == iai_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_ai = x(nvar_in)
         IF (var_dex(nvar_in) == iphi_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_phi = x(nvar_in)
         IF (var_dex(nvar_in) == iac_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_ac = x(nvar_in)
         IF (var_dex(nvar_in) == ine_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_ne = x(nvar_in)
         IF (var_dex(nvar_in) == ite_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_te = x(nvar_in)
         IF (var_dex(nvar_in) == iti_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_ti = x(nvar_in)
         IF (var_dex(nvar_in) == ith_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_th = x(nvar_in)
         IF (var_dex(nvar_in) == iah_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_ah = x(nvar_in)
         IF (var_dex(nvar_in) == iat_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_at = x(nvar_in)
         IF (var_dex(nvar_in) == izeff_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_zeff = x(nvar_in)
         IF (var_dex(nvar_in) == ibeamj_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_beamj = x(nvar_in)
         IF (var_dex(nvar_in) == ibootj_aux_f .and. arr_dex(nvar_in,2) == norm_dex) norm_bootj = x(nvar_in)
         IF (var_dex(nvar_in) == iemis_xics_f .and. arr_dex(nvar_in,2) == norm_dex) norm_emis_xics = x(nvar_in)
      END DO

      ! Unpack array (minus RBC/ZBS/RBS/ZBC)
      DO nvar_in = 1, n
         IF (arr_dex(nvar_in,2) == norm_dex) cycle
         IF (var_dex(nvar_in) == ixval) xval = x(nvar_in)
         IF (var_dex(nvar_in) == iyval) yval = x(nvar_in)
         IF (var_dex(nvar_in) == iphiedge) phiedge = x(nvar_in)
         IF (var_dex(nvar_in) == icurtor) curtor = x(nvar_in)
         IF (var_dex(nvar_in) == ipscale) pres_scale = x(nvar_in)
         IF (var_dex(nvar_in) == imixece) mix_ece = x(nvar_in)
         IF (var_dex(nvar_in) == ixics_v0) xics_v0 = x(nvar_in)
         !IF (var_dex(nvar_in) == ifocus_winding_surface_separation) &
         !       focus_winding_surface_separation = x(nvar_in)
         IF (var_dex(nvar_in) == iregcoil_winding_surface_separation) &
                regcoil_winding_surface_separation = x(nvar_in)
         IF (var_dex(nvar_in) == ibcrit) bcrit = x(nvar_in)
         IF (var_dex(nvar_in) == iextcur) extcur(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iaphi) aphi(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iam) am(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iac) ac(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iai) ai(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iah) ah(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iat) at(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ine) ne_opt(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == izeff) zeff_opt(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ite) te_opt(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iti) ti_opt(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ith) th_opt(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iam_aux_s) am_aux_s(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iam_aux_f) am_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iac_aux_s) ac_aux_s(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iac_aux_f) ac_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ibeamj_aux_f) beamj_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ibootj_aux_f) bootj_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iai_aux_s) ai_aux_s(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iai_aux_f) ai_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iphi_aux_f) phi_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ine_aux_f) ne_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == izeff_aux_f) zeff_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ite_aux_f) te_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iti_aux_f) ti_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == ith_aux_f) th_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iah_aux_f) ah_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iat_aux_f) at_aux_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iemis_xics_f) emis_xics_f(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iraxis_cc) raxis_cc(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == izaxis_cs) zaxis_cs(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == iraxis_cs) raxis_cs(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == izaxis_cc) zaxis_cc(arr_dex(nvar_in,1)) = x(nvar_in)
         IF (var_dex(nvar_in) == irhobc)     rhobc(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ideltamn)   deltamn(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == icoil_splinefx)   coil_splinefx(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == icoil_splinefy)   coil_splinefy(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == icoil_splinefz)   coil_splinefz(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ifocus_ds_rbound_c) focus_ds_rbound_c(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ifocus_ds_rbound_s) focus_ds_rbound_s(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ifocus_ds_zbound_c) focus_ds_zbound_c(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ifocus_ds_zbound_s) focus_ds_zbound_s(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == iregcoil_rcws_rbound_c) regcoil_rcws_rbound_c(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == iregcoil_rcws_rbound_s) regcoil_rcws_rbound_s(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == iregcoil_rcws_zbound_c) regcoil_rcws_zbound_c(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == iregcoil_rcws_zbound_s) regcoil_rcws_zbound_s(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == iRosenbrock_X) Rosenbrock_X(arr_dex(nvar_in,1)) = x(nvar_in)
      END DO

      ! Adust Boundary Representation
      IF (ANY(var_dex == irhobc)) THEN
         CALL unique_boundary(rbc,zbs,rhobc,mpol1d,ntord,mpol-1,ntor,mpol-1,rho_exp)
      END IF
      IF (ANY(var_dex == ideltamn)) THEN
         CALL unique_boundary_PG(rbc,zbs,deltamn,ntord,mpol1d,mpol-1,ntor)
      END IF

      ! Unpack RBC/ZBS/RBS/ZBC
      DO nvar_in = 1, n
         IF (var_dex(nvar_in) == ibound_rbc) rbc(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ibound_rbs) rbs(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ibound_zbc) zbc(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == ibound_zbs) zbs(arr_dex(nvar_in,1),arr_dex(nvar_in,2)) = x(nvar_in)
         IF (var_dex(nvar_in) == imodemn) THEN
            nf = arr_dex(nvar_in,1)
            mf = arr_dex(nvar_in,2)
            rbc(nf,mf) = x(nvar_in)
            zbs(nf,mf) = x(nvar_in)
            IF (mf == 0) THEN
               raxis_cc(nf) = x(nvar_in)
               zaxis_cs(nf) = x(nvar_in)
            END IF
         END IF
      END DO

      ! Apply normalization
      aphi = aphi * norm_aphi
      am   = am   * norm_am
      ac   = ac   * norm_ac
      ai   = ai   * norm_ai
      ne_opt = ne_opt * norm_ne
      zeff_opt = zeff_opt * norm_zeff
      te_opt = te_opt * norm_te
      ti_opt = ti_opt * norm_ti
      th_opt = th_opt * norm_th
      am_aux_f = am_aux_f * norm_am
      ac_aux_f = ac_aux_f * norm_ac
      ai_aux_f = ai_aux_f * norm_ai
      phi_aux_f = phi_aux_f * norm_phi
      ne_aux_f = ne_aux_f * norm_ne
      zeff_aux_f = zeff_aux_f * norm_zeff
      te_aux_f = te_aux_f * norm_te
      ti_aux_f = ti_aux_f * norm_ti
      th_aux_f = th_aux_f * norm_th
      ah_aux_f = ah_aux_f * norm_ah
      at_aux_f = at_aux_f * norm_at
      beamj_aux_f = beamj_aux_f * norm_beamj
      bootj_aux_f = bootj_aux_f * norm_bootj
      emis_xics_f = emis_xics_f * norm_emis_xics

      ! Handle cleanup
      IF (iflag < -2) THEN
         CALL stellopt_clean_up(ncnt,iflag)
         iflag = 0
         ! Now normalize arrays otherwise we'll be multiplying by normalizations on next iteration for non-varied quantities
         aphi = aphi / norm_aphi
         am   = am   / norm_am
         ac   = ac   / norm_ac
         ai   = ai   / norm_ai
         ne_opt = ne_opt / norm_ne
         zeff_opt = zeff_opt / norm_zeff
         te_opt = te_opt / norm_te
         ti_opt = ti_opt / norm_ti
         th_opt = th_opt / norm_th
         am_aux_f = am_aux_f / norm_am
         ac_aux_f = ac_aux_f / norm_ac
         ai_aux_f = ai_aux_f / norm_ai
         phi_aux_f = phi_aux_f / norm_phi
         ne_aux_f = ne_aux_f / norm_ne
         zeff_aux_f = zeff_aux_f / norm_zeff
         te_aux_f = te_aux_f / norm_te
         ti_aux_f = ti_aux_f / norm_ti
         th_aux_f = th_aux_f / norm_th
         ah_aux_f = ah_aux_f / norm_ah
         at_aux_f = at_aux_f / norm_at
         bootj_aux_f = bootj_aux_f / norm_bootj
         beamj_aux_f = beamj_aux_f / norm_beamj
         emis_xics_f = emis_xics_f / norm_emis_xics
         RETURN
      END IF

      ! Handle lscreen
      lscreen = .false.
      if (iflag < 0) lscreen=.true.
      istat = iflag
      !PRINT *,myid,iflag,iflag+myid,MOD(iflag+myid,4)
      ! Generate Random errors
!      IF (ncnt > n) THEN
!         IF (MOD(iflag+myid,4) < 1) THEN
!            fvec(1:m) = 10*SQRT(bigno/m)
!            iflag = 0
!            aphi = aphi / norm_aphi
!            am   = am   / norm_am
!            ac   = ac   / norm_ac
!            ai   = ai   / norm_ai
!            ne_opt = ne_opt / norm_ne
!            zeff_opt = zeff_opt / norm_zeff
!            te_opt = te_opt / norm_te
!            ti_opt = ti_opt / norm_ti
!            th_opt = th_opt / norm_th
!            am_aux_f = am_aux_f / norm_am
!            ac_aux_f = ac_aux_f / norm_ac
!            ai_aux_f = ai_aux_f / norm_ai
!            phi_aux_f = phi_aux_f / norm_phi
!            ne_aux_f = ne_aux_f / norm_ne
!            zeff_aux_f = zeff_aux_f / norm_zeff
!            te_aux_f = te_aux_f / norm_te
!            ti_aux_f = ti_aux_f / norm_ti
!            th_aux_f = th_aux_f / norm_th
!            RETURN
!         END IF
!      END IF

      ! Handle making a temporary string
      IF (iflag .eq. -1) istat = 0
      WRITE(temp_str,'(i5)') istat
      proc_string = TRIM(TRIM(id_string) // '_opt' // TRIM(ADJUSTL(temp_str)))

      ! Handle coil geometry variations
      IF (lcoil_geom) THEN
         CALL stellopt_spline_to_coil(lscreen)
         ctemp_str = 'write_mgrid'
         CALL stellopt_paraexe(ctemp_str,proc_string,lscreen)
      END IF

      IF (iflag .eq. -1) THEN 
         IF (lverb) WRITE(6,*) '---------------------------  EQUILIBRIUM CALCULATION  ------------------------'
      END IF

      ! Assume we've already read the stellopt input namelist and any input files.
      CALL tolower(equil_type)
         SELECT CASE (TRIM(equil_type))
            CASE('vmec2000_old','animec','flow','satire')
            CASE('paravmec','parvmec','vmec2000')
               iflag = 0
               CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
               iflag = ier_paraexe
               IF (lscreen .and. lverb) WRITE(6,*)  '-------------------------  PARAVMEC CALCULATION DONE  -----------------------'
            CASE('vboot')
               iflag = 0
               CALL stellopt_vboot(lscreen,iflag)
            CASE('vmec2000_oneeq')
               IF (iflag .eq. -1) THEN
                  iflag = 0
                  CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
                  iflag = ier_paraexe
                  IF (lscreen .and. lverb) WRITE(6,*)  '-------------------------  PARAVMEC CALCULATION DONE  -----------------------'
               ELSE
                  CALL read_wout_deallocate
                  CALL read_wout_file(TRIM(id_string)//'_opt0',iflag)
                  CALL write_wout_file('wout_'//TRIM(proc_string)//'.nc',iflag)
                  CALL stellopt_prof_to_vmec(proc_string,iflag)
                  iflag = 0
               END IF
            CASE('spec')
            CASE('test')
               !Do Nothing
         END SELECT
         ! Check profiles for negative values of pressure
         dex = MINLOC(am_aux_s(2:),DIM=1)
         IF (dex > 2) THEN
            IF (ANY(am_aux_f(1:dex) < 0)) iflag = -55
            IF (ALL(am_aux_f(1:dex) == 0)) iflag = -55
         END IF
         IF (pres_scale < 0) iflag = -55
         ! Now call any functions necessary to read or load the
         ! equilibrium output.  Things like conversion to other
         ! coordinate systems should be put here.  Note that it should be
         ! a function call which is handles every equil_type.  Note these
         ! functions should handle iflag by returning immediately if
         ! iflag is set to a negative number upon entry.
         CALL stellopt_load_equil(lscreen,iflag)

         ! Calls to secondary codes
         proc_string_old = proc_string ! So we can find the DIAGNO files
         IF (ANY(sigma_balloon < bigno)) CALL stellopt_balloon(lscreen,iflag)
         ctemp_str = 'booz_xform'
         IF (ANY(lbooz) .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
         ctemp_str = 'bootsj'
         IF (ANY(sigma_bootstrap < bigno) .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
         ctemp_str = 'diagno'
         IF (lneed_magdiag .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
         ctemp_str = 'neo'
         IF (ANY(sigma_neo < bigno) .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
!DEC$ IF DEFINED (TERPSICHORE)
         ctemp_str = 'terpsichore'
         IF (ANY(sigma_kink < bigno) .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
!DEC$ ENDIF
!DEC$ IF DEFINED (TRAVIS)
         ctemp_str = 'travis'
         IF (ANY(sigma_ece < bigno) .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
!DEC$ ENDIF
!DEC$ IF DEFINED (DKES_OPT)
         IF (ANY(sigma_dkes < bigno)) CALL stellopt_dkes(lscreen,iflag)
!DEC$ ENDIF

         ! NOTE ALL parallel secondary codes go here
!DEC$ IF DEFINED (TXPORT_OPT)
         IF (ANY(sigma_txport < bigno)) CALL stellopt_txport(lscreen,iflag)
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
         IF (ANY(sigma_orbit < bigno)) CALL stellopt_orbits(lscreen,iflag)
!DEC$ ENDIF
!DEC$ IF DEFINED (COILOPTPP)
         ctemp_str = 'coilopt++'
         IF (sigma_coil_bnorm < bigno .and. (iflag>=0)) CALL stellopt_paraexe(ctemp_str,proc_string,lscreen); iflag = ier_paraexe
!DEC$ ENDIF
!DEC$ IF DEFINED (FOCUS)
         ctemp_str = 'focus'
         write(6,*),'<----stellopt_fcn checking focus request'
         IF ( ( ANY(sigma_focus_bn < bigno) ) .and. (iflag >=0)) then
            write(6,*),'<----stellopt_fcn calling paraexe for ',ctemp_str
            CALL stellopt_paraexe(ctemp_str,proc_string,lscreen)
         END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (REGCOIL)
         IF ( ( ANY(sigma_regcoil_chi2_b < bigno) .or.   &
                ANY(sigma_regcoil_lambda < bigno) .or.    &
                ANY(sigma_regcoil_max_K < bigno) .or.    &
                ANY(sigma_regcoil_current_potential < bigno) .or.    &
                ANY(sigma_regcoil_max_current_potential < bigno) .or.    &
	        ANY(sigma_regcoil_rms_K < bigno) .or.    &
	        ANY(sigma_regcoil_chi2_k < bigno) .or.    &
	        ANY(sigma_regcoil_max_bnormal < bigno) .or.    &
	        ANY(sigma_regcoil_area_coil < bigno) .or.    &
	        ANY(sigma_regcoil_area_plasma < bigno) .or.    &
	        ANY(sigma_regcoil_area_diff < bigno) .or.    &
	        ANY(sigma_regcoil_volume_coil < bigno) .or.    &
	        ANY(sigma_regcoil_volume_plasma < bigno) .or.    &
	        ANY(sigma_regcoil_volume_diff < bigno) .or.    &
                ANY(sigma_regcoil_c2p_dist_min < bigno) .or.    &
	        ANY(sigma_regcoil_bnormal_total < bigno) .or.  &
	        ANY(sigma_regcoil_K2 < bigno)   ) .and. (iflag >=0)) then
           !Skipping parallelization for now - instead, call regcoi_driver directly.
           CALL stellopt_regcoil_driver(proc_string, lscreen, iflag)
           ! Unstested parallel call is commented out here.
           ! ctemp_str = 'regcoil'
           ! CALL stellopt_paraexe(ctemp_str,proc_string,lscreen)
         end if
!DEC$ ENDIF

         ! Now we load target values if an error was found then
         ! exagerate the fvec values so that those directions are not
         ! searched this levenberg step
         IF (iflag == 0) THEN
            CALL stellopt_load_targets(m,fvec,iflag,ncnt)
            WHERE(ABS(fvec) > bigno) fvec = bigno
            ier_paraexe = 0
         ELSE
            IF (lscreen) RETURN ! Make sure we can do at least the initial integration
            fvec(1:m) = 10*SQRT(bigno/m)
            iflag = 0 ! Because we wish to continue
            ier_paraexe = 0
         END IF
      ! Now normalize arrays otherwise we'll be multiplying by normalizations on next iteration for non-varied quantities
      aphi = aphi / norm_aphi
      am   = am   / norm_am
      ac   = ac   / norm_ac
      ai   = ai   / norm_ai
      ne_opt = ne_opt / norm_ne
      zeff_opt = zeff_opt / norm_zeff
      te_opt = te_opt / norm_te
      ti_opt = ti_opt / norm_ti
      th_opt = th_opt / norm_th
      am_aux_f = am_aux_f / norm_am
      ac_aux_f = ac_aux_f / norm_ac
      ai_aux_f = ai_aux_f / norm_ai
      phi_aux_f = phi_aux_f / norm_phi
      ne_aux_f = ne_aux_f / norm_ne
      zeff_aux_f = zeff_aux_f / norm_zeff
      te_aux_f = te_aux_f / norm_te
      ti_aux_f = ti_aux_f / norm_ti
      th_aux_f = th_aux_f / norm_th
      ah_aux_f = ah_aux_f / norm_ah
      at_aux_f = at_aux_f / norm_at
      bootj_aux_f = bootj_aux_f / norm_bootj
      beamj_aux_f = beamj_aux_f / norm_beamj
      emis_xics_f = emis_xics_f / norm_emis_xics
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_fcn
