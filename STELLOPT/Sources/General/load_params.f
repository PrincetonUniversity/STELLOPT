      SUBROUTINE load_params(xc_opt, nvar, iflag, lscreen)
      USE optim
      USE bootsj_input
      USE optim_params, ONLY: optimum
      USE safe_open_mod
      USE legendre_params
      USE vmec_input
      USE vparams, ONLY: zero
      USE coilsnamin, ONLY: lmodular, lsaddle, lbcoil, lvf, lsurfv,
     1    lsadsfv, ltfc, ltfcv, lmodcur, lsadcur
      USE mpi_params                                                     !MPI
      USE vacfield_mod, ONLY: angles, shifts, nstell_coils, vac_extcur
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: iflag
      INTEGER, INTENT(in) :: nvar
      REAL(rprec), DIMENSION(nvar), INTENT(in) :: xc_opt
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p96 = 0.96_dp
      REAL(rprec), PARAMETER :: ec  = 1.60217653e-19
!DEC$ IF DEFINED (NETCDF)
      CHARACTER(LEN=*), PARAMETER :: mgrid_ext = '.nc'
!DEC$ ELSE
      CHARACTER(LEN=*), PARAMETER :: mgrid_ext = '.bin'
!DEC$ ENDIF
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      LOGICAL :: lneteti
      INTEGER :: nb, ik, nvar_in, ireset, j1, iunit, istat
      INTEGER :: nvariables, nzero,ntemp
      INTEGER :: num_coefs, start_coefs
      LOGICAL :: firstp=.true., firstj=.true.,
     1  first=.true., second=.false.
      REAL(rprec) :: sum0, sum1, newmass, hs
      real(rprec) :: am0_9_opt                                           !COBRA
      real(rprec) :: tm0_9_opt                                           !LEGENDRE
      real(rprec) :: integral_pressure, integral_current, integral_iota
      CHARACTER(len=LEN_TRIM(home_dir)+20) :: version
      CHARACTER(len=LEN_TRIM(min_wout_file)+20) :: temp
      CHARACTER(LEN=8) :: screen
      LOGICAL :: lreset0, ex, lneed_mgrid
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL unique_boundary, unique_boundary_PG
      REAL(rprec), EXTERNAL :: Diode
      REAL(rprec) , EXTERNAL :: eval_prof
C-----------------------------------------------

!     SAL - 07/06/11 - Read the previous miminum file for am if fitting
!                      Do this here before we modify anything
!     VMEC 8.47 contains this info in wout file so we'll use the wout
!     values.
!      if (lpres_prof_fit) then
!         call read_input_am(trim(min_input_file),istat,am)
!      end if   
      
      lneteti   = .FALSE.
      nvar_in = 0
      second = (.NOT. first)

      lneed_mgrid = lcoil_geom .or. lvac_opt

!
!     SET UP RBC, ZBS BOUNDARY ARRAYS FROM ARRAY XC_OPT USED INTERNALLY
!     BY OPTIMIZATION CODE (THIS IS NOT THE VMEC XC ARRAY!)
!     NOTE: IF NITER_OPT = 1 AND THIS COMES FROM A FREE-BDY RUN, THEN
!     THE BOUNDARY MAY NOT BE IN THE EXPECTED FORM, SO KEEP RBC, ZBS
!     UNCHANGED
!

      IF (.not.lfreeb .and. .not.lvac_opt) THEN
         IF (nopt_boundary .le. 1) THEN
            DO ik = 1, irm0_bdy
               j1 = ik + nvar_in
               rbc(nrz0_opt(ik),0) = xc_opt(j1)
            END DO
            DO ik = 1+irm0_bdy, irm0_bdy+izm0_bdy
               j1 = ik + nvar_in
               zbs(nrz0_opt(ik),0) = xc_opt(j1)
            END DO
            nvar_in = nvar_in + irm0_bdy + izm0_bdy
            DO ik = 1, irho_bdy
               j1 = ik + nvar_in
               rhobc(nbrho_opt(ik),mbrho_opt(ik)) = xc_opt(j1)
            END DO
         ELSE IF (nopt_boundary .eq. 2) THEN
            DO ik = 1, irho_bdy
              j1 = ik + nvar_in
              delta_mn(nbrho_opt(ik),mbrho_opt(ik)) = xc_opt(j1)
            END DO
         ENDIF

         nvar_in = nvar_in + irho_bdy

         IF (.not.lniter1 .and. nopt_boundary.le.1)
     1       CALL unique_boundary(rbc, zbs, rhobc, mpol1d, ntord, 
     2                            mpol1_opt, ntor_opt, mrho1_opt)
         IF (.not.lniter1 .and. nopt_boundary.eq.2)
     1       CALL unique_boundary_PG(rbc, zbs, delta_mn, ntord, mpol1d,
     2                            mpol1_opt, ntor_opt)

      ELSE                                                               !Free-boundary optimization
         IF (lvac_opt) THEN
!        EXTRACT TILT, ANGLE OF COILS FROM xc_opt
            nb = 0
            DO ik = 1, nstell_coils
               nb = nb+1
               angles(ik) % as_array(1) = xc_opt(nb)    !theta
               nb = nb+1
               angles(ik) % as_array(2) = xc_opt(nb)    !phi
               nb = nb+1
               angles(ik) % as_array(3) = xc_opt(nb)    !rot-angle
               nb = nb+1
               shifts(ik) % as_array(1) = xc_opt(nb)    !x-shift
               nb = nb+1
               shifts(ik) % as_array(2) = xc_opt(nb)    !y-shift
               nb = nb+1
               shifts(ik) % as_array(3) = xc_opt(nb)    !z-shift
            END DO
            nvar_in = nvar_in+nb
            nb = 0
            DO ik = 1, SIZE(lextcur)
               IF (lextcur(ik)) THEN
                  nb = nb + 1
                  vac_extcur(ik) = xc_opt(nvar_in+nb)
               END IF
            END DO
            nvar_in = nvar_in+nb
         END IF

         IF( .not. lcoil_geom) THEN
            nb = 0
            DO ik = 1, SIZE(lextcur)
               IF (lextcur(ik)) THEN
                  nb = nb + 1
                  extcur(ik) = xc_opt(nvar_in+nb)
               END IF
            END DO
            IF (nb .ne. nextcur_opt)
     1         STOP 'Error counting EXTERNAL coils!'
            nvar_in = nvar_in + nextcur_opt

         ELSE                                                            !lcoil_geom=T
!************************************************************************
!     Load and count modular coil variables
!************************************************************************
            IF (lmodular) THEN
               CALL load_modular_coils (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
               IF (lmodcur) THEN
                  CALL load_modular_currents
     1                 (nvariables, xc_opt(nvar_in+1))
                  nvar_in = nvar_in + nvariables
               END IF
            END IF
!************************************************************************
!     Load and count saddle coil variables (IF lsaddle = true)
!************************************************************************
            IF (lsaddle) THEN
               CALL load_saddle_coils (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
               IF (lsadcur) THEN
                  CALL load_saddle_currents
     1                 (nvariables, xc_opt(nvar_in+1))
                  nvar_in = nvar_in + nvariables
               END IF
            END IF
!************************************************************************
!     Load and count background coil variables (IF lbcoil = true)
!************************************************************************
            IF (lbcoil) THEN
               CALL load_bg_currents (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            END IF
!************************************************************************
!     Load and count vf coil variables (IF lvf = true)
!************************************************************************
            IF (lvf) THEN
               CALL load_vf_currents (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            END IF

!************************************************************************
!     Load tf coil currents (IF ltfc = true)
!************************************************************************
            IF (ltfc .and. ltfcv) THEN
               CALL load_tf_coils (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            END IF
!************************************************************************
!     Load modular winding surface variables (IF lsurfv = true)
!************************************************************************
            IF (lsurfv) THEN
               CALL load_modular_wsurf (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            END IF
!************************************************************************
!     Load saddle winding surface variables (IF lsadsfv = true)
!************************************************************************
            IF (lsadsfv) THEN
               CALL load_saddle_wsurf (nvariables, xc_opt(nvar_in+1))
               nvar_in = nvar_in + nvariables
            END IF

         END IF               ! LCOIL_GEOM

      END IF      ! LFREEB


!************************************************************************
!     LOAD AI or AC PROFILE COEFFICIENTS (IF LIOTA_PROF_OPT = TRUE or LCUR_PROF_OPT = TRUE)
!     PICK LAST COEFFICIENT SO THAT FOR NCURR=1, <JZETA(s=1)> = 0, TO AVOID
!     SURFACE CURRENT FORMATION
!************************************************************************

      IF (num_ai .gt. 11) STOP 'num_ai>11 in STELLOPT load_params'

      IF (lcur_prof_opt) THEN
         CALL tolower(pcurr_type)
         SELECT CASE(TRIM(pcurr_type))
         CASE('gauss_trunc')
            STOP 'PCURR_TYPE = gauss_trunc not implemented in STELLOPT'
         CASE('two_power')
            STOP 'PCURR_TYPE = two_power not implemented in STELLOPT'
         CASE('power_series_i')
            STOP 
     1        'PCURR_TYPE = power_series_i not implemented in STELLOPT'
         CASE('akima_spline_ip', 'akima_spline_i',
     1        'cubic_spline_ip', 'cubic_spline_i')
            nvar_in = nvar_in + 1
            curtor = xc_opt(nvar_in)
            num_coefs=minloc(ac_aux_s(2:),dim=1)                        ! Find zeros
            IF (lcur_opt_edge0) THEN
               ac_aux_f(num_coefs) = 0
               num_coefs = num_coefs -1                                 ! Don't vary edge coefficient
            END IF
            IF (num_coefs > 4) THEN
               DO ik = 1, num_coefs
                  nvar_in = nvar_in + 1
                  ac_aux_f(ik) = xc_opt(nvar_in)
               END DO
            ELSE
               STOP 'CHECK AC_AUX_S for number of values (>4)'
            END IF
         CASE DEFAULT
            nzero = 10  ! num_ai    !  ac() entry to use to force j(a)=0.
            IF(l_legendre) THEN                                             !LEGENDRE
               tc(0:num_ai-1) = xc_opt(1+nvar_in:num_ai+nvar_in)
               IF (lj1zero) tc(num_ai) = -SUM(tc(0:num_ai-1))
               CALL legendre_to_power(num_ai-1, a_leg_inv, b_leg_inv,
     1            tc, ac(0))
            ELSE
               ntemp = 1
               num_ai = 11
               do ik = 0, 10
                  if (ac_mask(ik) .ne. 0) then
                     ac(ik) = xc_opt(nvar_in+ntemp)
                     ntemp=ntemp+1                                           ! Keep track of moments turned on
                  end if
               end do
               ! first element is curtor
               sum0 = zero; sum1 = zero
               DO j1 = 1, num_ai-1
                  sum0 = sum0 + ac(j1)
                  sum1 = sum1 + ac(j1)/(j1+1)
               END DO
               ac(0) = ac(0) - sum1                                          ! make sure the integral is correct
               IF (lcur_opt_edge0) THEN
                  num_ai = 9
                  ac(num_ai) = (num_ai+1)*(ac(0)+sum0) / (10-num_ai)
                  ac(nzero) = -( ac(0) + sum0 + ac(num_ai))                  ! force edge to zero
               END IF
            ENDIF
            sum1 = ac(0)
            DO j1 = 1, 10
               sum1 = sum1 + ac(j1)/(j1+1)
            END DO
            curtor = sum1                                                   !Varies curtor this way
            nvar_in = nvar_in + ntemp - 1                                   !Only some moments turned on
         END SELECT
      ELSE IF (liota_prof_opt) THEN
         CALL tolower(piota_type)
         SELECT CASE(TRIM(piota_type))
         CASE('sum_atan')
            STOP 'PIOTA_TYPE = sum_atan not implemented in STELLOPT'
         CASE('akima_spline','cubic_spline')
            num_coefs=minloc(ai_aux_s(2:),dim=1)                        ! Find zeros
            IF (num_coefs > 4) THEN
               DO ik = 2, num_coefs-1
                  nvar_in = nvar_in + 1
                  ai_aux_f(ik) = xc_opt(nvar_in) * ai_aux_f(1)
               END DO
               ! Calculate and apply integral Iota
!               num_coefs=minloc(ai_aux_s(2:),dim=1)                        ! Find zeros for integral
!               nvar_in = nvar_in + 1
!               integral_iota = 0
!               DO ik = 1, num_coefs-1
!                  integral_iota = integral_iota
!     1                      + ai_aux_f(ik) * 
!     2                           ( ai_aux_s(ik+1) - ai_aux_s(ik) )
!               END DO
!               ai_aux_f(1:num_coefs) = ai_aux_f(1:num_coefs) *
!     1                              xc_opt(nvar_in) / integral_iota
            ELSE
               STOP 'CHECK AI_AUX_S for number of values (>4)'
            END IF
         CASE('pedestal')
            STOP 'PIOTA_TYPE = pedestal not implemented in STELLOPT'
         CASE('rational')
            STOP 'PIOTA_TYPE = rational not implemented in STELLOPT'
         CASE('nice_quadratic')
            STOP 'PIOTA_TYPE=nice_quadratic not implemented in STELLOPT'
         CASE DEFAULT                                                    ! Power Series
            IF(l_legendre) THEN                                             !LEGENDRE
               ti(0:num_ai-1) = xc_opt(1+nvar_in:num_ai+nvar_in)
               CALL legendre_to_power(num_ai-1, a_leg_inv, b_leg_inv,
     1            ti, ai(0))
               nvar_in = nvar_in + num_ai
            ELSE
               DO j1 = 2, num_ai
                  nvar_in = nvar_in + 1
                  ai(j1-1) = xc_opt(nvar_in)*ai(0)
               END DO
               nvar_in = nvar_in+1
               ai(0:10) = ai(0:10) * xc_opt(nvar_in) /
     1            sum( ai(0:10) /
     2                (/ 1,2,3,4,5,6,7,8,9,10,11 /) )
            ENDIF
         END SELECT
      END IF

!     Alternative way to vary curtor
      IF (ncurr.eq.1 .and. (.not.lcur_prof_opt) .and.
     1   (ABS(sigma_maxcurrent).lt.bigno .or. ledge_current .or.
     2      ABS(sigma_curtor) < bigno) ) THEN
         nvar_in = nvar_in + 1
         curtor = xc_opt(nvar_in)
      END IF

!*****************************************************************************
!       VARIATION OF PHI-EDGE, USING LIMITED (PHIEDGE_MIN,MAX) AND "DIODE" MODEL
      IF (lphiedge) THEN
         nvar_in = nvar_in + 1
         IF (phiedge_diode) THEN
            phiedge = .5_dp*((phiedge_max + phiedge_min) +
     1             (phiedge_max - phiedge_min)*Diode(xc_opt(nvar_in)))
         ELSE
            phiedge = xc_opt(nvar_in)
         END IF
      END IF
      
!************************************************************************
!       CODE added by S. Lazerson (05/14/12)
!       MSE Radial Electric Field coefficients
!************************************************************************
      IF (lmse_er) THEN
         num_coefs=minloc(er_aux_s(2:),dim=1)                        ! Find zeros
         IF (num_coefs > 4) THEN
            DO ik = 1, num_coefs
               nvar_in = nvar_in + 1
               er_aux_f(ik) = xc_opt(nvar_in)
            END DO
         ELSE
            STOP 'CHECK ER_AUX_S for number of values (>4)'
         END IF
         num_coefs=minloc(ez_aux_s(2:),dim=1)                        ! Find zeros
         IF (num_coefs > 4) THEN
            DO ik = 1, num_coefs
               nvar_in = nvar_in + 1
               ez_aux_f(ik) = xc_opt(nvar_in)
            END DO
         ELSE
            STOP 'CHECK EZ_AUX_S for number of values (>4)'
         END IF
      END IF
      
!************************************************************************
!       CODE added by S. Lazerson (06/04/12)
!       NE, TE, TI coefficients
!************************************************************************
      num_coefs = minloc(ne_aux_s(2:), dim=1)
      IF (num_coefs > 4) THEN
!        NE may not go to zero where Te goes to zero
!         IF (lpres_opt_edge0) THEN
!            ne_aux_f(num_coefs) = 0.0
!            num_coefs = num_coefs-1
!         END IF
         nvar_in = nvar_in + 1
         sum1 = xc_opt(nvar_in)
         DO ik = 1, num_coefs
            nvar_in = nvar_in + 1
            ne_aux_f(ik) = xc_opt(nvar_in)*sum1
         END DO
         lneteti = .true.
      END IF
      num_coefs = minloc(te_aux_s(2:), dim=1)
      IF (num_coefs > 4) THEN
         IF (lpres_opt_edge0) THEN
            te_aux_f(num_coefs) = 0.0
            num_coefs = num_coefs-1
         END IF
         nvar_in = nvar_in + 1
         sum1 = xc_opt(nvar_in)
         DO ik = 1, num_coefs
            nvar_in = nvar_in + 1
            te_aux_f(ik) = xc_opt(nvar_in)*sum1
         END DO
         lneteti = .true.
      END IF
      num_coefs = minloc(ti_aux_s(2:), dim=1)
      IF (num_coefs > 4) THEN
         IF (lpres_opt_edge0) THEN
            ti_aux_f(num_coefs) = 0.0
            num_coefs = num_coefs-1
         END IF
         nvar_in = nvar_in + 1
         sum1 = xc_opt(nvar_in)
         DO ik = 1, num_coefs
            nvar_in = nvar_in + 1
            ti_aux_f(ik) = xc_opt(nvar_in)*sum1
         END DO
         lneteti = .true.
      END IF
      
      IF (lneteti) THEN
         nvar_in = nvar_in + 1
         factor_p_prof = xc_opt(nvar_in)
         num_coefs = minloc(am_aux_s(2:), dim=1)
         DO ik = 1, num_coefs
            am_aux_f(ik) = factor_p_prof * ne_aux_f(ik) * ec
     1                     * (te_aux_f(ik)+ti_aux_f(ik))
         END DO
      END IF
      
!************************************************************************
!       CODE added by S. Lazerson (08/16/11)
!       ANIMEC: Vary Anisotropic coefficients
!************************************************************************
      IF (lanimec) THEN
         IF (lani_bcrit) THEN                                             ! Critical Field Strength
            nvar_in = nvar_in + 1
            bcrit = xc_opt(nvar_in)
         END IF
         IF (lani_tperp) THEN                                             ! T_perp/T_|| for hot part.
            DO ik = 0, 10
               IF (at_mask(ik) .ne. 0) THEN
                  nvar_in = nvar_in + 1
                  at(ik) = xc_opt(nvar_in)
               END IF
            END DO
         END IF
         IF (lani_phot) THEN                                              ! P_hot/P_iso ratio
            DO ik = 0, 10
               IF (ah_mask(ik) .ne. 0) THEN
                  nvar_in = nvar_in + 1
                  ah(ik) = xc_opt(nvar_in)
               END IF
            END DO
         END IF
      END IF
   
!************************************************************************
!       CODE added by R.SANCHEZ (01/19/99).and E. Lazarus (05/25/05)
!       COBRA: Vary mass coefficients to taylor pressure profile to ballooning modes.
!       New profile options added by Sam Lazerson (04/20/12).
!************************************************************************
      IF (lpres_prof_opt) THEN
         CALL tolower(pmass_type)
         SELECT CASE(TRIM(pmass_type))
         CASE('gauss_trunc')
            STOP 'PMASS_TYPE = gauss_trunc not implemented in STELLOPT'
         CASE('two_power')
            STOP 'PMASS_TYPE = two_power not implemented in STELLOPT'
         CASE('two_lorentz')
            STOP 'PMASS_TYPE = two_lorentz not implemented in STELLOPT'
         CASE('akima_spline','cubic_spline')
            num_coefs=minloc(am_aux_s(2:),dim=1)                        ! Find zeros
            IF (lpres_opt_edge0) THEN
               am_aux_f(num_coefs) = 0                                  ! Set edge to zero
            END IF
            IF (lpres_opt_edgegr0) THEN
               am_aux_f(num_coefs-1) = am_aux_f(num_coefs)              ! Set gradient to zero
            END IF
            IF (num_coefs > 4) THEN
               ! Get Values
               DO ik = 2, num_coefs-2
                  nvar_in = nvar_in + 1
                  am_aux_f(ik) = xc_opt(nvar_in)*am_aux_f(1)            ! Vary ratios to ease workload 
               END DO
               ! Get Edge Gradient Values
               IF (.not. lpres_opt_edgegr0) THEN
                  nvar_in = nvar_in + 1
                  am_aux_f(num_coefs-1) = xc_opt(nvar_in)*am_aux_f(1)
               END IF   
               ! Get Edge Value
               IF (.not. lpres_opt_edge0) THEN
                  nvar_in = nvar_in + 1
                  am_aux_f(num_coefs) = xc_opt(nvar_in)*am_aux_f(1)
               END IF
            ELSE
               STOP 'CHECK AM_AUX_S for number of values (>4)'
            END IF
            ! Calculate and apply integral pressure
            num_coefs=minloc(am_aux_s(2:),dim=1)                        ! Find zeros for integral
            nvar_in = nvar_in + 1
            integral_pressure = 0
            DO ik = 1, num_coefs-1
               integral_pressure = integral_pressure
     1                      + am_aux_f(ik) * 
     2                           ( am_aux_s(ik+1) - am_aux_s(ik) )
            END DO
            am_aux_f(1:num_coefs) = am_aux_f(1:num_coefs) *
     1                              xc_opt(nvar_in) / integral_pressure
            ! Match Beta
 !           IF (     abs(sigma_beta)    < bigno 
 !    1          .or. abs(sigma_eplasma) < bigno) THEN
 !              nvar_in = nvar_in + 1
 !              factor_p_prof = xc_opt(nvar_in)/am_aux_f(1)
 !           END IF
            ! Test for positive pressure profile
            IF (ANY(am_aux_f(1:num_coefs) < 0)) THEN
               ierr_vmec = 1
               iflag = 0
               RETURN       ! can't run vmec, so get out of here
            END IF
         CASE('pedestal')
            nvar_in = nvar_in + 1
            am(0) = xc_opt(nvar_in)
            DO j1 = 1, MIN(15, pres_opt_nmax)
               nvar_in = nvar_in + 1
               am(j1) = xc_opt(nvar_in)*am(0)
            END DO
            ! Note am(16) isn't used
            ! Note am(20) is automatically set
            DO j1 = 17, 19
               nvar_in = nvar_in + 1
               am(j1) = xc_opt(nvar_in)
            END DO
            ! set floor for am(19) so code will run.
            IF (am(19) .le. 0) am(19) = 1.0E-4
            !STOP 'PMASS_TYPE = pedestal not implemented in STELLOPT'
         CASE('rational')
            STOP 'PMASS_TYPE = rational not implemented in STELLOPT'
         CASE DEFAULT                                                   ! power_series
            IF(.not. l_legendre) THEN                                        !LEGENDRE
               am0_9_opt = zero
               am(0) = 1
               am(1:10) = 0
               ! Get AM(1) through AM(8)
               DO j1 = 1, MIN(8, pres_opt_nmax)
                  nvar_in = nvar_in + 1
                  am(j1) = xc_opt(nvar_in)
               END DO
               IF( lpres_opt_edge0 .and. lpres_opt_edgegr0) THEN
                  am(9) = - sum( am(1:8)*
     1                    (/ 9,8,7,6,5,4,3,2 /) )                        ! use am(9 ) to guarantee p_edge gradient=0.
                  am0_9_opt = sum(am(0:9))
                        am(10) = - am0_9_opt                                      ! use am(10) to guarantee p_edge=0.
               ELSE IF (lpres_opt_edge0) THEN
                  IF (pres_opt_nmax >= 9) THEN
                     nvar_in = nvar_in + 1
                     am(9) = xc_opt(nvar_in)
                  END IF
                  am0_9_opt = sum(am(0:9))
                  am(10) = - am0_9_opt                                      !use am(10) to guarantee p_edge=0.
               ELSE IF( lpres_opt_edgegr0 ) THEN
                  IF (pres_opt_nmax >= 9) THEN
                     nvar_in = nvar_in + 1
                     am(9) = xc_opt(nvar_in)
                  END IF
                  am(10) = - sum( am(1:9)*
     1                       (/ 1,2,3,4,5,6,7,8,9 /) ) / 10      !use am(10) to guarantee p_edge gradient=0.
               ELSE IF(pres_opt_nmax >= 9) THEN
                  nvar_in = nvar_in + 1
                  am(9) = xc_opt(nvar_in)
                  IF (pres_opt_nmax >= 10) THEN
                     nvar_in = nvar_in + 1
                     am(10) = xc_opt(nvar_in)
                  ENDIF
               END IF
!              normalize to s-integral of pressure profile
               nvar_in = nvar_in+1
               am(0:10) = am(0:10) * xc_opt(nvar_in) /
     1            sum( am(0:10) /
     2                (/ 1,2,3,4,5,6,7,8,9,10,11 /) )
            ELSE                                                             !LEGENDRE
               DO j1 = 0, num_am-1
                  nvar_in = nvar_in + 1
                  tm(j1) = xc_opt(nvar_in)
               ENDDO
               IF (lp1zero) tm(num_am) = -SUM(tm(0:num_am-1))                !Since tm(1) = 1 for all m
               CALL legendre_to_power(num_am, a_leg_inv, b_leg_inv, tm,
     1                                am(0))
            ENDIF
            ! Scaling factor for Beta Matching
            IF (     abs(sigma_beta)    < bigno 
     1          .or. abs(sigma_eplasma) < bigno) THEN
               nvar_in = nvar_in + 1
               factor_p_prof = xc_opt(nvar_in) / am(0)
            END IF
            ! Test for positive Pressure profile
            hs = 1._rprec / (nrad-1)
            DO j1 = 2, nrad
               IF( eval_prof(am, (j1-1.5_rprec)*hs ) < 0) THEN
                  ierr_vmec = 1
                  iflag = 0
                  RETURN       ! can't run vmec, so get out of here
               END IF
            END DO
         END SELECT
      ENDIF
      
!************************************************************************


      IF (nvar_in .ne. nvar) THEN
         PRINT *,'MYID: ',myid,' nvar_in = ', nvar_in,' nvar = ', nvar
         STOP 'nvar_in != nvar in STELLOPT load_params'
      END IF

      lreset0 = lreset_opt
      IF (iflag .eq. -1) lreset0 = .true.
      ireset = 0

!************************************************************************
!     WRITE OUT DATA TO INPUT FILE (WILL BE READ BY XVMEC AND XCOILGEOM)
!************************************************************************
      lrecon = .false.
      IF (lreset0) THEN
         raxis_cc(:) = raxis_old(:,1);  raxis_cs(:) = raxis_old(:,2)
         zaxis_cs(:) = zaxis_old(:,1);  zaxis_cc(:) = zaxis_old(:,2)
      END IF

      IF (lneed_mgrid) mgrid_file='mgrid_' // TRIM(opt_ext) // mgrid_ext

      CALL write_indata (input_file, iflag)
      
      IF (iflag .ne. 0) PRINT*,"write_indata (input_file, iflag)",
     1 " iflag=",iflag," input_file ",TRIM(input_file)
      IF (iflag .ne. 0) RETURN

!***********************************************************************
!     GENERATE MGRID FILE AND PARSE EXTERNAL CURRENTS ARRAY.
!     INPUT_FILE MUST BE CLOSED (ABOVE) TO BE READ IN GENERATE MGRID CALL
!     MUST RE-WRITE INPUT FILE TO REFLECT UPDATED CURRENTS.
!***********************************************************************
      IF (lneed_mgrid) THEN
         CALL generate_mgrid (opt_ext, (lscreen .and. (myid.eq.master)),
     1                        iflag)
         IF (iflag .ne. 0) RETURN

!***********************************************************************
!        RE-WRITE DATA TO INPUT FILE (WILL BE READ BY XVMEC), NOW THAT
!        WE HAVE EXTCUR VALUES. ALSO RECOMPUTES NEXTCUR_VMEC VALUE.
!***********************************************************************
         CALL write_indata (input_file, iflag)
         IF (iflag .ne. 0) RETURN
      END IF

!********************************************************************
!     EXECUTE VMEC CODE AND LOAD VMEC_INPUT MODULE VALUES
!***********************************************************************
      IF (lanimec) THEN
         version = TRIM(home_dir) // 'xanimec' // TRIM(exe_suffix)
      ELSE
         version = TRIM(home_dir) // 'xvmec2000' // TRIM(exe_suffix)
      END IF
      INQUIRE(file=TRIM(version), exist=ex, iostat=istat)
      IF (ex) THEN
         screen = 'noscreen'
         IF (lscreen) screen = 'screen'
         WRITE (temp,'(1x,a,1x,l1)') screen, lreset0
         IF (.not.lreset0) WRITE                                 
     1       (temp,'(a,1x,a)') TRIM(temp), TRIM(min_wout_file)
!DEC$ IF .NOT.DEFINED (WIN32)
         IF (.not.lscreen) temp = TRIM(temp) // ' > /dev/null'
!DEC$ ENDIF
         iunit = 99
         CALL load_physics_codes (version, 'input', temp, 'wout',
     1        opt_ext, iunit, iflag)
      ELSE
         iflag = -6
         RETURN
      ENDIF

      IF (lneed_mgrid) CALL system(remove // TRIM(mgrid_file))
      FIRST = .false.
      CALL flush(6)
      END SUBROUTINE load_params


      FUNCTION Diode(x)
      USE stel_constants
      REAL(rprec) :: x, Diode

      Diode = 2*ATAN(x)/pi

      END FUNCTION Diode
