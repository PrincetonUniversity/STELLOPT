      SUBROUTINE load_target(xc_opt, fvec, extension, nopt, ncnt,
     1                       iflag, lscreen)
!      USE read_boozer_mod, ONLY: mnboz_b, mboz_b, nboz_b, 
!     1                           read_boozer_deallocate          
      USE read_boozer_mod
      USE optim
      USE boozer_params
      USE chisq_mod
      USE safe_open_mod
      USE vmec_input, ONLY: lfreeb, extcur, phiedge
      USE vparams, ONLY: zero, one, twopi, mu0
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nopt, iflag
      INTEGER, INTENT(in) :: ncnt
      REAL(rprec), DIMENSION(*) :: fvec, xc_opt
      CHARACTER(LEN=*), INTENT(in)  :: extension
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
      LOGICAL :: lprint_bmin = .false.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: icontrol, jcount, no_elements, nsval,
     1   n, k, i, nmodB, istat, iunit, nweight(ntargets),
     2   mboz_surf, mskip, iunit_err
      REAL(rprec) :: chisq_tot, epsm
      REAL(rprec), DIMENSION(ntargets) :: chisq, mean_weight
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmod_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bmin_b, bmax_b
      REAL(rprec) :: hs, dmax_j, bmin_0, bmin_pi2, bmin_pi,
     2  bmax_0, bmax_pi2, bmax_pi, dtheta, sqrt_nsurf
      CHARACTER(len=LEN_TRIM(home_dir)+20) :: version
      LOGICAL :: ex, lscreen_master
C-----------------------------------------------
!
!     DESCRIPTION OF LOCAL INPUT AND CONTROL VARIABLES
!     (SEE LIBSTELL module OPTIM_PARAMS NAMELIST for description of
!      variable in that namelist)
!
!   Transport optimization related
!---------------------------------
!       lprint_bmn (local var)  controls if boozer spectra is printed to a file (=t)
!                               only true if niter_opt = 1 and lbmn = t or
!                               pseudo-sym is on.
!
!       sigma_jstar          2D Array (nrad by 4) of sigmas for Jstar
!                            at 4 values of velocity pitch-angle eps/mu on each surface
!
!       ldkes_opt            =T optimize diffusion coefficients from DKES
!       sigma_dkes           =weights for DKES diffusion coefficients
!
!       lneo_opt             LOGICAL, T= run NEO to calculate effective ripple
!                            DEFAULT false
!       sigma_neo            array of sigmas for effective ripple for each surface
!
!       lorbit_opt           =T optimize fast ion confinement vis MC simulation
!       sigma_orbit          sigma for the fast ion loss
!
!       sigma_pseudo         Sigma for pseudo-symmetry 'water volume' of wells
!                            relative to toroidal well volume
!       sigma_pseudo2        Sigma for second 'water' measure: width of wells
!                            relative to field-line length
!
!       sigma_bmin           1D Array (nrad) of sigmas for BMIN
!       sigma_bmax           1D Array (nrad) of sigmas for BMAX
!       sigma_ripple         1D Array (nrad) of sigmas for ripple
!
!
!   Bootstrap current related
!---------------------------
!       lbootsj_opt          =T, add match to bootstrap current to chi-sq
!
!       jboot                Bootstrap model to use:
!                               =0  Tolliver 3D, asymptotic collisionless
!                               =1  Simple 2D, asymptotic collisionless
!       sigma_bootsj         Sigma for matching bootstrap current self-consistently
!       sigma_boot           Sigma for bootstrap current density, rel to MAX of
!                            bootstrap or equil. current, w/prog. spatial weight
!
!       fboot                array of s**j coefficients for multiplicative factor
!                            on bootstrap current to mimic nu-star effects
!       zeff_boot            Zeff value for use in bootstrap current calculation
!       at                   polynomial coefficients for specifying plasma temperature
!                            Te=Ti for bootstrap, (keV)
!       lseedcur             =T, add seed current to bootstrap current before matching
!       aseedcur             polynomial coefficients for seed current, similar to AC coefs.
!       lseed_opt            =F, the aseedcur coefficients are fixed
!                            =T, the non-zero aseedcur cofficients are optimized
!                            ***not yet implemented***
!
!
!    Flux surface quality related
!--------------------------------
!
!       n_jac                   array of n-number for resonantBoozer jacobian components
!       m_jac                   array of m-numbers for resonant Boozer jacobian components
!       sigma_jac               sigmas for supressing resonant Boozer jacobian components
!
!       n_vac_island            array of n-number for vacuum island targetting
!       m_vac_island            array of m-numbers for vacuum island targetting
!       sigma_vac_island        sigmas for supressing vacuum islands
!
!       ldsubr_opt              LOGICAL, T= run JMC to calculate Dr
!       sigma_dsubr             array of sigmas for resistive interchange (Dr) stability
!
!
!    Vacuum vessel and limiter matching related
!--------------------------------
!       target_vv               Target for closest approach distance
!       target_vv_rms           Target for RMS average distance
!       sigma_vv                sigma for closest distance to specified
!                                  3D boundary shape (if <~1.e10, else ignore)
!       sigma_vv_rms            sigma for RMS average distance from plasma to
!                                  3D boundary (if <~ 1.e10, else ignore)
!       mpol_vv                 maximum m for 3D boundary
!       ntor_vv                 maximum n for 3D boundary
!       rbc_vv                  R-COS(theta) array for 3D boundary
!       zbs_vv                  Z-SIN(theta) array for 3D boundary
!
!       shapeweight             If F boundary weighting is turned off and
!                               the VV routine is used. Default = .T.
!       planes_bw(3) (radians)  planes in which shape is evaluated
!                               Do not normalize to Nfp
!                               0, Pi/2 and Pi are recommended
!       amplw_bw(3) (radians)   relative weighting of planes(3)
!       theta0_bw(3) (radians)  poloidal angle at which weight FUNCTION maiimizes
!                               Gaussian Weigts: 15 terms
!                               theta0 centers are at theta0, theta0+Pi/2,
!                               theta0+Pi, theta0+3Pi/2 &  theta0+2 Pi
!       phi0_bw (radians)       toroidal angles at which weights maximize
!                               centers at phi0/Nfp, (phi0+Pi/2/)Nfp, & (phi0+ Pi)/Nfp
!       wtheta_bw (radians)     poloidal half-width of shaping weights; e.g.,Pi/8
!       wphi_bw (radians)       toroidal half-width of shaping weights
!
!    equilibrium reconstruction  from magnetic sensors
!------------------------------------------------------------------------------
!	lphiedge 	if true phiedge will be varied to minimize 
!			chisq_diagnostics
!	phiedge_min	if set and lphiedg=T, the search on phiedge is bounded
!			below
!	phiedge_max     if set and lphiedg=T, the search on phiedge is bounded
!			above
!	ledge_current	if true total plasma current will be varied to minimize
!                       chisq_diagnostics 
!     lreconp         if true am(0:kpp-1)  will be varied to minimize
!                       chisq_diagnostics with am(kjj:10)=0 and a constraint
!       lp1zero         if(lp1zero) am(kpp)=-sum(am(0:kpp-1)
!	kpp		number of terms in plasma pressure polynomial
!	lreconj        if true ac(0:kjj-1)  will be varied to minimize
!                       chisq_diagnostics with ac(kjj:10)=0 and a constraint 
!	lj1zero		if(lj1zero) ac(kjj)=-sum(ac(0:kjj-1)
!	kjj             number of terms in plasma current polynomial
!
!
!
!-----------------------------------------------------------------------------------------------
!     Load information from wout file
!-----------------------------------------------------------------------------------------------
      IF (nopt .le. 0) THEN
         i = 0
         CALL read_wout_opt(i, extension, jcount, iflag)                !loads boozer parameters 1st time thru
         IF (iflag .ne. 0) THEN
            iflag = -1
         ELSE IF (jcount .ne. 0) THEN
            iflag = -8
         END IF
         IF (iflag .ne. 0) RETURN

      ELSE
         dtheta = one
         IF (nu2_b .gt. 1) dtheta = one/REAL(nu2_b - 1,rprec)
         ALLOCATE (index_array(nopt),
     1             wegt(nopt), chisq_match(nopt), chisq_target(nopt),
     2             stat=n)
         IF (n .ne. 0) THEN
            iflag = -9
            RETURN
         END IF

         index_array = 0
         chisq_target = 0
         chisq_match = 0
         wegt = bigno

      ENDIF

      ALLOCATE (bmin_b(nu2_b), bmax_b(nu2_b), bmod_b(nv_boz,nu_boz),
     1         stat=n)

      IF (nrad .gt. 1) hs = one/REAL((nrad - 1),rprec)
      nmodB = nu2_b

!-----------------------------------------------------------------------------------------------
!     COMPUTE BOOZER TRANSFORMATION HERE ON RADIAL SURFACES
!     DETERMINED BY LSURF_MASK LOGICAL ARRAY. IF NO RADIAL DERIVATIVES
!     OF Bmn REQUIRED, ONLY NEED ONE SURFACE AT A TIME (OTHERWISE, COMPUTE
!     AND STORE ADJACENT SURFACES).
!-----------------------------------------------------------------------------------------------

      mskip     = nu2_b/mboz
      mboz_surf = 0
      DO n = 1, nu2_b, mskip
         mboz_surf = mboz_surf + 1
      END DO

      sqrt_nsurf = SQRT(REAL(mboz_surf,rprec))
      WRITE (comm1,'(1x,l1)') lscreen            !!Needed for command line arg
      version = TRIM(home_dir) // 'xbooz_xform' // TRIM(exe_suffix)
      lscreen_master = lscreen .and. (myid.eq.master)

      IF (nopt.gt.0) THEN
!-----------------------------------------------------------------------------------------------
!     OPEN INPUT FILES NEEDED FOR COMMUNICATING WITH EXECUTABLES
!-----------------------------------------------------------------------------------------------
         CALL open_comm_files(mboz, nboz, ns_booz, ns_booz_max,
     1      ns_surf, ns_surf_max, lbootsj_opt, extension, iflag)

!-----------------------------------------------------------------------------------------------
!        Write and read in bmn from boozmn.extension file 
!-----------------------------------------------------------------------------------------------
         IF (ns_booz_max .ge. 1) THEN ! SAL - So boozer doesn't get called if not needed
            CALL load_physics_codes (version, 'in_booz', comm1,
     1                               'boozmn',extension, iunit, iflag)
            IF (iflag .eq. 0) THEN
               IF (mnboz.ne.mnboz_b .or.  mboz.ne.mboz_b .or.
     1             nboz.ne.nboz_b) THEN
                   iflag = -10
               END IF
            END IF
         END IF

!-----------------------------------------------------------------------------------------------
!        Open DKES_OPT file to accumulate OPT_DKES data for all surfaces
!-----------------------------------------------------------------------------------------------
         IF (ldkes_opt) THEN
            iunit_dkes = 19
            CALL safe_open(iunit_dkes, k, 'dkes_opt.'//TRIM(extension),
     1                  'replace', 'formatted')
            IF (k .ne. 0) THEN
               STOP 'Error opening dkes_opt file in load_target'
            ELSE
               WRITE(iunit_dkes,'(2a,/)')
     1           'SUMMARY OF DKES TRANSPORT COEFFICIENTS FOR CASE: ',
     2           TRIM(extension)
               WRITE(iunit_dkes,"(78('='),/,3x,3(9x,a,9x),/,78('='))")
     1           'L11(+/-)', 'L13(+/-)', 'L33(+/-)'
            END IF
         END IF
!-----------------------------------------------------------------------------------------------
!     initialize AJAX for inverse coordinate mapping, if needed
!-----------------------------------------------------------------------------------------------
         IF( lajax) THEN
            CALL init_ajax
         END IF
         
      ELSE

         INQUIRE(file=TRIM(version), exist=ex, iostat=istat)
         IF ((istat.ne.0 .or. .not.ex) .and. myid.eq.master) THEN
            PRINT *, 'xbooz_xform file not found in ' // TRIM(home_dir)
            iflag = -16
            RETURN
         END IF
         
         lajax = .false.                                 ! initialize flag on AJAX inverse mapping

      ENDIF

!-----------------------------------------------------------------------------------------------
!     BEGIN MATCHING TARGETS
!-----------------------------------------------------------------------------------------------
      no_elements = 0
      iflag = 0

!-----------------------------------------------------------------------------------------------
!     COMPUTE DIAGNOSTIC RESPONSES TARGETS
!-----------------------------------------------------------------------------------------------
      IF (lv3post) CALL chisq_diagnostics (no_elements, nopt, iflag, 
     1                           extension, lscreen_master)
      IF (iflag .ne. 0) RETURN
!-----------------------------------------------------------------------------------------------
!
!     MATCH ASPECT RATIO
!
      CALL chisq_aspect (aspect_opt, target_aspectratio,
     1     sigma_aspect, no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH MAXIMUM TOROIDAL CURRENT
!    (DO NOT USE ABSOLUTE VALUE: OPTIMIZER UNSTABLE THAT WAY)
!
      CALL chisq_maxcurrent (target_maxcurrent, sigma_maxcurrent,
     1                       no_elements, nrad, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH DESIRED POLOIDAL FLUX               MCZ  Sep. 98
!
      IF (sigma_fluxp .lt. bigno) CALL chisq_polflux (target_fluxp,
     1    sigma_fluxp, hs, nrad, no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH VOLUME-AVERAGED BETA
!
!      CALL chisq_beta (Target_Beta, sigma_beta, no_elements, nopt)
       CALL chisq_beta (Target_Beta, sigma_beta,
     1                 Target_eplasma, sigma_eplasma,
     2                 no_elements, nopt)

!------------------------------------------------------------------------------
!
!     MATCH Total current
!
      call chisq_curtor (Target_curtor, sigma_curtor,
     1                   no_elements, nrad, nopt)
     
!------------------------------------------------------------------------------
!
!     MATCH edge current density (PPPL)
!
      call chisq_jedge (Target_jedge, sigma_jedge,
     1                   no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH R-Btor AT S=1 (USED IN FREE-BDY MODE TO SET PHIEDGE)             MCZ  July 00
!
      CALL chisq_rbtor (rbtor_opt, target_rbtor,
     1     sigma_rbtor, no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MATCH IOTA PROFILE (IF NCURR=1)
!
      CALL chisq_iota(iota_opt, sigma_iota, hs, no_elements, nrad, nopt)

! -----------------------------------------------------------------------
!     iota-prime                             M.Zarnstorff
!
      CALL chisq_iota_p(iota_opt, sigma_iota_pmax, sigma_iota_pmin,
     1                 hs, no_elements, nrad, nopt)

!-----------------------------------------------------------------------------------------------
!
!     MERCIER CRITERION (DETERMINES STABLE MAGNETIC WELL = vp(s)-vp(0), normed to vp(0))
!     Dmerc < 0 is UNSTABLE
!
      CALL chisq_mercier (Dmerc_opt, sigma_mercier, no_elements,
     1     nrad, nopt)

!-----------------------------------------------------------------------------------------------
!
!     EXPLICITLY DETERMINE MAGNETIC WELL = vp(s)-vp(0), normed to vp(0)
!
      IF (ANY(ABS(sigma_vp(2:nrad)) .lt. bigno)) CALL chisq_magwell
     1   (hs, sigma_vp, no_elements, nrad, nopt)

!-----------------------------------------------------------------------------------------------
!     CODE added by R.SANCHEZ (01/19/99). Modified (02/01/99) to include
!     multiple initial toroidal and poloidal positions.
!-----------------------------------------------------------------------------------------------

      IF (lballoon_opt) THEN
         CALL chisq_ballooning (sigma_balloon, target_balloon,
     1       no_elements,  ns_ball, ns_ball_max,                !!VMEC COBRA (RS)
     2       nopt, iflag, extension, comm1)
         IF (iflag .ne. 0) RETURN
      END IF

      IF(     (np_prof > 0 .and. ANY(sigma_p_prof<bigno))
     1   .or. (nne_prof > 0 .and. ANY(sigma_ne_prof < bigno))
     2   .or. (nte_prof > 0 .and. ANY(sigma_te_prof < bigno))
     3   .or. (nti_prof > 0 .and. ANY(sigma_ti_prof < bigno)) ) THEN
         CALL chisq_p_prof(pres_opt, ivar_p_prof, no_elements, nrad,
     1                        nopt, iflag, extension)
         IF (iflag .ne. 0) RETURN
         IF (isote) THEN
            CALL chisq_isote(ivar_isote, no_elements, nrad,
     1                        nopt, iflag, extension)
            IF (iflag .ne. 0) RETURN
         END IF
      END IF
      IF (lpres_prof_opt
     1   .or. (nne_prof > 0 .and. ANY(sigma_ne_prof < bigno))
     2   .or. (nte_prof > 0 .and. ANY(sigma_te_prof < bigno))
     3   .or. (nti_prof > 0 .and. ANY(sigma_ti_prof < bigno)) ) THEN
         pres_opt = pres_opt
         CALL chisq_bpres (pres_opt, sigma_pedge, ivar_pedge,
     1        no_elements, nrad, nopt, 0)
         CALL chisq_bpres (pres_opt, sigma_pgrad, ivar_pgrad,
     1        no_elements, nrad, nopt, 1)
      END IF

!-----------------------------------------------------------------------------------------------
!
!     CODE added by R. SANCHEZ, L. BERRY, and S. HIRSHMAN (01/19/99) to MATCH
!     BOOTSTRAP CURRENT <J dot B> PROFILE COMPUTED FROM VMEC
!-----------------------------------------------------------------------------------------------

      IF(lbootsj_opt) THEN
         CALL chisq_bootsj (sigma_bootsj, no_elements,
     1      ns_surf, ns_surf_max, nrad, nopt, iflag, extension,
     2      lscreen_master)
         IF (iflag .ne. 0) RETURN
      END IF

!-----------------------------------------------------------------------------------------
!     Motional Stark Effect Diagnostic (MSE) (PPPL)
!-----------------------------------------------------------------------------------------

      IF(ANY(sigma_mse_pol < bigno)) THEN
          call chisq_mse (ivar_mse, ivar_iota_extrap, no_elements,nopt, 
     1                    iflag, extension)

          IF (iflag .ne. 0) RETURN
      END IF


!-----------------------------------------------------------------------------------------
!     DIAGNO targetting magnetic diagnostics (PPPL)
!-----------------------------------------------------------------------------------------

      if( ldiagno_opt .and. ndiagno_seg+ndiagno_flx+ndiagno_bp > 0) then
          call chisq_diagno (no_elements, nopt, iflag,
     1                            extension, lscreen_master)

          if (iflag .ne. 0) return
      endif

!-----------------------------------------------------------------------------------------------
!     EXTERNAL KINK MODE CRITERION (LPK, MCZ, PPPL - Feb 1999, SPH Feb 2000)
!-----------------------------------------------------------------------------------------------
      IF (lkink_opt) CALL chisq_kink (no_elements,
     1      nopt, iflag, extension, lscreen_master)
      IF (iflag .ne. 0) RETURN

!-----------------------------------------------------------------------------------------------
!     COIL TARGETS FROM XCOILGEOM, AS CONTROLED IN THE COILSIN NAMELIST
!-----------------------------------------------------------------------------------------------
!      IF (lcoil_geom) CALL chisq_coilgeom (no_elements, nopt, iflag,
!     1    extension, lscreen_master)
!      IF (iflag .ne. 0) RETURN
      IF (lcoil_geom .or. sigma_berr_avg < bigno .or.
     1    sigma_berr_max < bigno  )
     2       CALL chisq_coilgeom (no_elements, nopt, iflag,
     3                            extension, lscreen_master)
      IF (iflag .ne. 0) RETURN

!-----------------------------------------------------------------------------------------------
!     COMPUTE TARGET CRITERIA THAT DEPEND ON BOOZER-TRANSFORMED QUANTITIES
!     SUCH AS: J*, Jinvariant, Bmin, Bmax, Ripple, Bmn
!-----------------------------------------------------------------------------------------------
      IF (lprint_bmin) THEN
         iunit = 18
         CALL safe_open (iunit, k, 'bmin_bmax.' // TRIM(extension),
     1                  'replace', 'formatted')
         WRITE (iunit, *) ns_surf_max, nu2_b
      END IF

      icontrol = 0
      SURF: DO jcount = 1, ns_booz_max
         nsval = ns_booz(jcount)
!-----------------------------------------------------------------------------------------------
!        COMPUTE BOOZER-SPACE MOD-B vs theta-booz, zeta-booz IF NECESSARY
!-----------------------------------------------------------------------------------------------
         IF (nopt .gt. 0) THEN
            IF (jcount .eq. ns_booz_max) icontrol = -1    !!release memory here
            IF (lneed_modB(nsval) .or. (icontrol .eq. -1)) THEN
               CALL modbooz(bmod_b, nfp_opt, icontrol, nsval)
               CALL bextrema (bmod_b, bmin_b, bmax_b, nv_boz, nu2_b)
            END IF

            IF (lprint_bmin) THEN
               WRITE (iunit, *) nsval
               DO k = 1,nu2_b
                    WRITE (iunit, *) k, bmin_b(k), bmax_b(k)
               END DO
               IF (jcount .eq. ns_booz_max) CLOSE (iunit)
            END IF
         END IF                !!End nopt > 0

!-----------------------------------------------------------------------------------------------
!     added by R. FOWLER (03/03/00)
!-----------------------------------------------------------------------------------------------
         IF (ldkes_opt .and. ldkes_mask(nsval)) THEN
            CALL chisq_dkes(sigma_dkes(nsval), no_elements,
     1          nopt, iflag, nsval, extension, comm1)
            IF (iflag .ne. 0) RETURN
         END IF

         IF (.not.lsurf_mask(nsval)) CYCLE
   
!-----------------------------------------------------------------------------------------------
!     MATCH BMN-BOOZER CRITERION. DUE TO THE VARIETY OF POSSIBLE
!     WEIGHTS (SIGMAS) AND TARGET COMBINATIONS (BMN_B=0, BMN_B(m,n) = BMN_B(m1,n1), etc),
!     THE USER MAY DECIDE TO HAND-CODE THE WEIGHTS AND TARGETS IN THE FOLLOWING LOOP
!
!     modified to allow targetting specifict Bmn values
!     M. Zarnstorff              July 2006
!-----------------------------------------------------------------------------------------------
         if (lbmn) then

         IF (lbmn) CALL chisq_bmn(sigma_bmn(nsval), hs, no_elements, 
     1                            mnboz, nsval, nopt)

            if( any(sigma_bmn_tgt(nsval,:) .lt. bigno) ) then
               CALL chisq_bmn_tgt(target_bmn_tgt(nsval,:),
     1                            sigma_bmn_tgt(nsval,:),
     2                            n_bmn_tgt,
     3                            m_bmn_tgt,size(sigma_bmn_tgt,2),
     4                            hs, no_elements, mnboz, nsval, nopt,
     5                            lscreen)
            endif

         endif

!-----------------------------------------------------------------------------------------------
!   Pseudo-symmetry criteria; try and discourage local wells. Modified by MCZ, Feb. 99, SPH Feb. 00
!-----------------------------------------------------------------------------------------------
         if (sigma_pseudo(nsval).lt.bigno .or.
     1       sigma_pseudo2(nsval).lt.bigno)
     2      CALL chisq_pseudo (sigma_pseudo, sigma_pseudo2, hs,
     3                         no_elements, nopt, nsval,
     4                         lpseudo_sin, lscreen_master)

!-----------------------------------------------------------------------------------------------
!        MINIMIZE THE VARIATION OF THE SECOND ADIABATIC INVARIANT, J-INVARIANT,
!        AROUND A FLUX SURFACE FOR SELECTED RADII AND VALUES
!        OF EPSILON/mu (=Bturning). THIS CALLS xj_invariant WHICH COMPUTES |B|
!        ITSELF FROM THE BOOZMN FILE, SO IT DOES NOT REQUIRE A CALL HERE TO MODBOOZ
!        Added by D. A. Spong, 1/00
!-----------------------------------------------------------------------------------------------
         IF (lj_invariant) CALL chisq_jinvar
     1     (sigma_jinvariant(nsval,1:NumJinvariant),
     2      nsval, no_elements, nopt, nu2_b/mskip, sqrt_nsurf, iflag,
     3      extension, comm1, lscreen_master)
         IF (iflag .ne. 0) RETURN

!-----------------------------------------------------------------------------------------------
!
!        I M P O R T A N T    D E V E L O P E R    N O T E
!
!        <<ALL>> ROUTINES HEREAFTER (IN SURF LOOP) REQUIRE BMOD_B, Bmin_B, AND/OR Bmax_B.
!        IF ANY ARE USED (FINITE SIGMA), THEN THE CALL TO MODBOOZ ABOVE MUST HAVE BEEN
!        MADE. THIS IS DETERMINED BY THE LOGICAL LNEED_MODB, WHICH IS COMPUTE THE FIRST TIME
!        THROUGH (FOR NOPT <= 0).
!
!-----------------------------------------------------------------------------------------------
         IF (nopt .le. 0) THEN
            n = no_elements
         ELSE IF (.not.lneed_modB(nsval)) THEN
            CYCLE
         END IF

!-----------------------------------------------------------------------------------------------
!        Minimize (BMIN - <Bmin>)
!-----------------------------------------------------------------------------------------------
         CALL chisq_bmin (sigma_bmin(nsval), ivar_bmin, bmin_b,
     1          no_elements, nopt, mskip, sqrt_nsurf, dtheta)

!-----------------------------------------------------------------------------------------------
!        Minimize (BMAX - <Bmax>)
!-----------------------------------------------------------------------------------------------
         CALL chisq_bmin (sigma_bmax(nsval), ivar_bmax, bmax_b,
     1          no_elements, nopt, mskip, sqrt_nsurf, dtheta)

!-----------------------------------------------------------------------------------------------
!        Minimize magnetic ripple on each surface
!        7/98: added SIN(u/2) weight to emphaSIZE outboard ripple reduction
!-----------------------------------------------------------------------------------------------
         CALL chisq_bripple (sigma_ripple(nsval), bmax_b,
     1          bmin_b, no_elements, nopt, mskip, sqrt_nsurf, dtheta)

!-----------------------------------------------------------------------------------------------
!        MINIMIZE THE THE VARIATION OF THE TRAPPED PARTICLE J*
!        AROUND A FLUX SURFACE FOR SELECTED RADII AND VALUES
!        OF EPSILON/mu (=Bturning).
!-----------------------------------------------------------------------------------------------
         IF (lj_star) CALL chisq_jstar (sigma_jstar(nsval,1:NumJstar),
     1      bmax_b, bmin_b, bmod_b, no_elements, nopt,
     2      mskip, sqrt_nsurf)
!-----------------------------------------------------------------------------------------------

         IF (nopt.le.0 .and. no_elements.ne.n) lneed_modB(nsval)=.true.

      END DO SURF

!
!     Clean up: close any open files, etc.
!
      IF (ldkes_opt .and. nopt.gt.0) CLOSE(unit=iunit_dkes)

!-----------------------------------------------------------------------------------------------
!!
!!    NEO  LPK 01-19-01
!!
         IF (lneo_opt) CALL chisq_neo(sigma_neo, ivar_neo, no_elements,
     1          ns_neo, ns_neo_max, nopt, iflag, extension,
     2          lscreen_master)
         IF (iflag .ne. 0) RETURN

         IF (ldsubr_opt) CALL chisq_dsubr(sigma_dsubr, ivar_dsubr,
     1          no_elements, nopt, iflag, extension, lscreen_master)
         IF (iflag .ne. 0) RETURN

         IF (lorbit_opt) CALL chisq_orbit(sigma_orbit, ivar_orbit,
     1          no_elements, nopt, iflag, extension, lscreen_master)
         IF (iflag .ne. 0) RETURN

!------------------------------------------------------------------------------
!
!    JConf    M.Zarnstorff  Aug. 2002
!
!------------------------------------------------------------------------------

         call chisq_jconf(sigma_jconf, nrad, no_elements, nopt,
     1      nu2_b/mskip, iflag, extension, comm1, lscreen_master)
         if (iflag .ne. 0) return

!-----------------------------------------------------------------------------------------------
!
!     Minimize peak curvature of outer flux surfaces
!     at four different toroidal planes
!     (phi = 0, 90, 180, and 270 degrees)
!     added by d.a. spong
!     PPPL - Enclosed in IF statement
      IF ( sigma_curv < bigno ) then
         CALL chisq_curvature(sigma_curv, no_elements, nopt)
      END IF

!-----------------------------------------------------------------------------------------------
!
!     Maintain boundary between Target_rmax, Target_rmin, and below Target_zmax
!
      CALL chisq_rzbdy(Target_Rmax, rmax_opt, sigma_rmax,
     1         ivar_center, no_elements, nopt, .true.)
      CALL chisq_rzbdy(Target_Rmin, rmin_opt, sigma_rmin,
     1         ivar_center, no_elements, nopt, .false.)
      CALL chisq_rzbdy(Target_Zmax, zmax_opt, sigma_zmax,
     1         ivar_zmax, no_elements, nopt, .true.)


! -----------------------------------------------------------------------
!
!     3D boundary limiting                      R.Hatcher & M.Zarnstorff
!
!
!      IF ((sigma_vv < bigno) .or. (sigma_vv_rms < bigno)) THEN
!         CALL chisq_3dbd(sigma_vv, ivar_vv, sigma_vv_rms,
!     1        target_vv, target_vv_rms,
!     2        ivar_vv_rms, no_elements, nopt, iflag, lscreen_master)
!      END IF
! -----------------------------------------------------------------------
!
!     3D boundary limiting  (PPPL)              R.Hatcher & M.Zarnstorff
!

      if (sigma_vv < bigno .or. sigma_vv_rms < bigno .or.
     1    sigma_vv_max < bigno ) then
         call chisq_3dbd(sigma_vv, ivar_vv, sigma_vv_rms,
     1        sigma_vv_max, target_vv, target_vv_rms, lvv_tgt_min,
     2        ivar_vv_rms, ivar_vv_max, no_elements, nopt, iflag,
     3        lscreen_master)
      end if

      if (sigma_bd < bigno .or. sigma_bd_rms < bigno .or.
     1    sigma_bd_max < bigno ) then
         call chisq_3dbd_bd(sigma_bd, ivar_vv, sigma_bd_rms,
     1        sigma_bd_max, target_bd, target_bd_rms,
     2        ivar_vv_rms, ivar_vv_max, no_elements, nopt, iflag,
     3        lscreen_master)
      end if

! -----------------------------------------------------------------------
!     Regularize Free-Bdry Coil Currents            M.Zarnstorff
!
      CALL chisq_extcur(sigma_extcur, target_extcur, ivar_extcur,
     1                   no_elements, nopt)

! -----------------------------------------------------------------------
!     Linear regularization constraint on Coil Currents     M.Zarnstorff
!
      CALL chisq_OH(sigma_oh, ivar_oh, no_elements, nopt)

! -----------------------------------------------------------------------
!     <kappa> : surface-averaged elongation         M.Zarnstorff
!
      CALL chisq_kappa(target_kappa, sigma_kappa, ivar_kappa,
     1                 no_elements, nopt)

!------------------------------------------------------------------------------
!
!     Suppress targeted resonant Boozer Jacobian components on their
!     respective resonant surfaces.
!     The resonances are indicated by each pair of entries (n_jac, m_jac)
!     with the sigmas being the respective sigma_jac entry
!
!     MCZ  Jan. 01

      CALL chisq_jac (iota_opt, hs, ivar_jac, no_elements, ns_booz, 
     1                ns_booz_max, mnboz, nopt, lscreen)

! -----------------------------------------------------------------------
!
!     Limit ellipticity of the phi=0 cross section   A. Ware, D. Spong, S. Hirshman
!
      CALL chisq_ellipticity (Target_ellipticity, sigma_ellipticity,
     1        no_elements, nopt)

!-----------------------------------------------------------------------------------------------
!
!     Minimize coil-complexity measure (>1, the minimum value)
!     and maximum current density
!
      IF (lnescoil_opt) CALL chisq_nescoil
     1    (no_elements, nopt, iflag, extension, lscreen_master)

!-----------------------------------------------------------------------------------------------
!
!     Island due to 3D Shaping                       L-P. Ku
!
      CALL chisq_vac_island(ivar_vac_island, no_elements, nopt, iflag,
     1     extension, lscreen_master)
!-----------------------------------------------------------------------------------------------
!
!     Vacuum Islands                                 S.P.Hirshman/D.J.Strickler
!
      IF (lvac_opt) CALL chisq_vacopt (no_elements, nopt, iflag, 
     1              extension, lscreen_master)
!-----------------------------------------------------------------------------------------------
!
!     Free memory associated with boozer transform arrays
!
      CALL read_boozer_deallocate

!-----------------------------------------------------------------------------------------------
      IF (nopt .le. 0) THEN
         IF (nopt .eq. -1) 
     1      ALLOCATE (chisq_descript(no_elements), stat=n)
         nopt = no_elements
         DEALLOCATE (bmin_b, bmax_b, bmod_b)
         RETURN
      ELSE IF (nopt .ne. no_elements) THEN
         IF (myid .eq. master)                                           !MPI
     1       PRINT *,' NOPT = ',nopt,'!= NO_ELEMENTS = ',no_elements
         iflag = -11
         RETURN
      ELSE IF (.not.lfreeb) THEN
!        Check that boundary data are correctly read in
         CALL chk_rzmnb(xc_opt, rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy,
     1        iflag)
         IF (iflag .ne. 0) RETURN
      ENDIF

!
!     COMPUTE TOTAL CHISQ. NOTE THE "chisq" VARIABLES HERE HAVE
!     NO BEEN SQUARED YET (AS REQUIRED BY LEVENBERG...)
!
!      epsm = EPSILON(wegt(1))                                           !ORNL
      epsm = EPSILON(wegt(1))**2                                         !PPPL
      wegt(:nopt) = MAX (ABS(wegt(:nopt)), epsm)                         !PPPL
      WHERE (ABS(wegt(:nopt)) .gt. epsm)
        fvec(:nopt) =
     1      (chisq_target(:nopt) - chisq_match(:nopt))/wegt(:nopt)
      ELSEWHERE
        fvec(:nopt) = 0
      END WHERE


      chisq_tot = SUM(fvec(:nopt)**2)
      
!-----------------------------------------------------------------------------------------------
!     COMPUTE Error information
!-----------------------------------------------------------------------------------------------
      CALL safe_open(iunit_err,k,'errvals.'//TRIM(extension),
     1               'replace','formatted')
      WRITE(iunit_err,'(A)') 'chisq_target(1:nopt)-chisq_match(1:nopt)'
      WRITE(iunit_err,'(I5)') nopt
      WRITE(iunit_err,'(E20.10)') 
     1                 (chisq_target(i)-chisq_match(i),i=1,nopt)
      WRITE(iunit_err,'(A)') 'wegt(1:nopt)'
      WRITE(iunit_err,'(E20.10)') wegt(1:nopt)
      CLOSE(iunit_err)

!-----------------------------------------------------------------------------------------------
!     COMPUTE INDIVIDUAL CHISQs
!-----------------------------------------------------------------------------------------------

      DO i = 1, ntargets
         nweight(i) = COUNT(index_array.eq.i .and. wegt.gt.epsm)
         IF (nweight(i) .eq. 0) THEN
            chisq(i) = 0
            mean_weight(i) = 0
         ELSE
!            chisq(i) = SUM (fvec(:nopt)**2, mask = index_array.eq.i)    !ORNL
!            mean_weight(i) = SUM (ABS(one/wegt(:nopt)**2), mask =
!     1           (index_array.eq.i .and. wegt.gt.epsm))
            chisq(i) = sum (fvec(:nopt)**2,
     1           mask = (index_array.eq.i .and. abs(wegt).gt.epsm))      !PPPL
            mean_weight(i) = sum (abs(one/wegt(:nopt)**2),
     1           mask = (index_array.eq.i .and. abs(wegt).gt.epsm))      !PPPL
         END IF
      END DO

      bmin_0 = 0
      bmin_pi2 = 0
      bmin_pi = 0
      bmax_0 = 0
      bmax_pi2 = 0
      bmax_pi = 0

      IF (ns_booz_max .ge. 1) THEN
         nsval = ns_booz(ns_booz_max)

         if( any(lneed_modB(1:ns_booz_max))) then
            bmin_0 = bmin_b(1)
            bmin_pi2 = bmin_b(nu2_b/2)
            bmin_pi = bmin_b(nu2_b)
            bmax_0 = bmax_b(1)
            bmax_pi2 = bmax_b(nu2_b/2)
            bmax_pi = bmax_b(nu2_b)
         endif
      ELSE
         nsval = 0
      END IF

!      bmin_0 = bmin_b(1)
!      bmin_pi2 = bmin_b(nu2_b/2)
!      bmin_pi = bmin_b(nu2_b)
!      bmax_0 = bmax_b(1)
!      bmax_pi2 = bmax_b(nu2_b/2)
!      bmax_pi = bmax_b(nu2_b)
      dmax_j = twopi * MAXVAL(ABS(buco_opt(2:nrad))) / mu0

!
!     Begin writing to output record file
!
      IF (lscreen_master) THEN
         WRITE (*, 544) aspect_opt, chisq_tot, ncnt, nsval,
     1                  nrad, bmin_0, bmin_pi2, bmin_pi, bmax_0,
     2                  bmax_pi2, bmax_pi,  dmax_j, wp_opt/wb_opt
         IF (sigma_vv < bigno .or. sigma_vv_rms < bigno)
     1      WRITE (*, 545) vv_dist, vv_dist_rms, vv_dist_max
         if (sigma_bd < bigno .or. sigma_bd_rms < bigno .or.
     1      sigma_bd_max < bigno )
     1      write (*, 547) bd_dist, bd_dist_rms, bd_dist_max
      ENDIF
      WRITE (iunit_opt_local, 544) aspect_opt, chisq_tot, ncnt, nsval,
     1    nrad, bmin_0, bmin_pi2, bmin_pi, bmax_0, bmax_pi2, bmax_pi,
     2    dmax_j, wp_opt/wb_opt
      IF (sigma_vv < bigno .or. sigma_vv_rms < bigno)
     1   WRITE (iunit_opt_local, 545) vv_dist, vv_dist_rms, vv_dist_max
      if (sigma_bd < bigno .or. sigma_bd_rms < bigno)
     1   write (iunit_opt_local, 547) bd_dist,bd_dist_rms,bd_dist_max

!
!     output either coil current or boundary coefficients for intermediate restart
!
      IF (ldiag_opt) THEN
        WRITE(iunit_opt_local,*)
        IF( lfreeb ) THEN
            WRITE(iunit_opt_local,'(a)') " EXTCUR ="
            WRITE(iunit_opt_local,*) extcur
            write(iunit_opt_local,'(a,es14.7)') " PHIEDGE = ",phiedge
        ELSE
          write(iunit_opt_local,'(a,es14.7)') " PHIEDGE = ",phiedge
          CALL write_rbzb(iunit_opt_local, i)
        ENDIF
      END IF

      DO i = 1, ntargets
         IF (nweight(i) .gt. 0)
     1      WRITE (iunit_opt_local, 646)
     2      descript(i), chisq(i), SQRT(chisq(i)/mean_weight(i)),
     3      SQRT(mean_weight(i)/nweight(i)), nweight(i)
      END DO
 646  FORMAT(' Chi-Sq(',a,')=',1p,e12.5, '  RMS Err=', e12.5,
     1     '  <Wgt>=', e12.5, ' <Constraints> = ', i5)

!     DO i = 1, nopt
!        WRITE (iunit_opt_local, 546) chisq_descript(i), i, fvec(i),
!    1      chisq_target(i), chisq_match(i), wegt(i)
!     END DO

 544  FORMAT(/' Aspect Ratio = ',f12.5,' Chi-Sq = ',1p,e12.5,
     1   ' Function Evaluations = ',i5,/,' At js = ',i3,
     2   ' out of ',i3,' radial points:',/,' Bmin_0 = ',e12.5,
     3   ' Bmin_pi/2 = ',e12.5,' Bmin_pi = ',e12.5,/,' Bmax_0 = ',
     4     e12.5,' Bmax_pi/2 = ',e12.5,' Bmax_pi = ',e12.5,/,
     5   ' Maximum Toroidal Current = ',e12.5,' Average Beta = ',
     6     e12.5)

 545  format(' Min Plasma-VV Sep =', 1pe12.5,
     1          '  RMS dist =', 1pe12.5,
     2          '  Max dist =', 1pe12.5)


 546  FORMAT(a,3x,' fvec(',i3,') = ',d12.5,' Target = ',d12.5,
     1   ' Match = ',d12.5,' Sigma(with units) = ',d12.5)

 547  format(' Min Plasma-BD Sep =', 1pe12.5,
     1          '  RMS dist =', 1pe12.5,
     2          '  Max dist =', 1pe12.5)


      DEALLOCATE (bmin_b, bmax_b, index_array, wegt,
     1            chisq_match, chisq_target, bmod_b, stat=n)


      END SUBROUTINE load_target
