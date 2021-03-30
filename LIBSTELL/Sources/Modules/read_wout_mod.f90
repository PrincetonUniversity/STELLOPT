      MODULE read_wout_mod
!
!     USE READ_WOUT_MOD to include variables dynamically allocated
!     in this module
!     Call DEALLOCATE_READ_WOUT to free this memory when it is no longer needed
! 
!     Reads in output from VMEC equilibrium code(s), contained in wout file
!
!     Contained subroutines:
!
!     read_wout_file      wrapper alias called to read/open wout file      
!     read_wout_text      called by read_wout_file to read text file wout
!     read_wout_nc        called by read_wout_file to read netcdf file wout
!
!     Post-processing routines
!
!     mse_pitch           user-callable function to compute mse pitch angle
!                         for the computed equilibrium
!

      USE v3_utilities
      USE vmec_input, ONLY: lrfp, lmove_axis, nbfld
      USE mgrid_mod

      IMPLICIT NONE
#if defined(NETCDF)
!------------------------------------------------
!   L O C A L   P A R A M E T E R S
!------------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER(LEN=*), PARAMETER :: vn_version = 'version_',           &
        vn_extension = 'input_extension', vn_mgrid = 'mgrid_file',      &
        vn_magen = 'wb', vn_therm = 'wp', vn_gam = 'gamma',             &
        vn_maxr = 'rmax_surf', vn_minr = 'rmin_surf',                   &
        vn_maxz = 'zmax_surf', vn_fp = 'nfp',                           &
        vn_radnod = 'ns', vn_polmod = 'mpol', vn_tormod = 'ntor',       &
        vn_maxmod = 'mnmax', vn_maxit = 'niter', vn_actit = 'itfsq',    &
        vn_maxpot ='mnmaxpot', vn_potsin = 'potsin', vn_potcos = 'potcos', &
        vn_xmpot = 'xmpot', vn_xnpot='xnpot',                           &
        vn_asym = 'lasym', vn_recon = 'lrecon',                         &

        vn_free = 'lfreeb', vn_moveaxis = 'lmove_axis',                 &
        vn_error = 'ier_flag', vn_aspect = 'aspect', vn_rfp = 'lrfp',   &
        vn_maxmod_nyq = 'mnmax_nyq',                                    &
        vn_beta = 'betatotal', vn_pbeta = 'betapol',                    &
        vn_tbeta = 'betator', vn_abeta = 'betaxis',                     &
        vn_b0 = 'b0', vn_rbt0 = 'rbtor0', vn_rbt1 = 'rbtor',            &
        vn_sgs = 'signgs', vn_lar = 'IonLarmor', vn_modB = 'volavgB',   &
        vn_ctor = 'ctor', vn_amin = 'Aminor_p', vn_Rmaj = 'Rmajor_p',   &
        vn_vol = 'volume_p', vn_am = 'am', vn_ai = 'ai', vn_ac = 'ac',  &
        vn_ah = 'hot particle fraction', vn_atuname = 'T-perp/T-par',   &
        vn_pmass_type = 'pmass_type', vn_piota_type = 'piota_type',     &
        vn_pcurr_type = 'pcurr_type',                                   &
        vn_am_aux_s = 'am_aux_s', vn_am_aux_f = 'am_aux_f',             &
        vn_ai_aux_s = 'ai_aux_s', vn_ai_aux_f = 'ai_aux_f',             &
        vn_ac_aux_s = 'ac_aux_s', vn_ac_aux_f = 'ac_aux_f',             &
        vn_mse = 'imse', vn_thom = 'itse',                              &
        vn_pmod = 'xm', vn_tmod = 'xn', vn_pmod_nyq = 'xm_nyq',         &
        vn_tmod_nyq = 'xn_nyq',                                         &
        vn_racc = 'raxis_cc', vn_zacs = 'zaxis_cs',                     &
        vn_racs = 'raxis_cs', vn_zacc = 'zaxis_cc', vn_iotaf = 'iotaf', &

        vn_qfact='q-factor', vn_chi='chi', vn_chipf='chipf',            &
        vn_presf = 'presf', vn_phi = 'phi', vn_phipf = 'phipf',         &
        vn_jcuru = 'jcuru', vn_jcurv = 'jcurv', vn_iotah = 'iotas',     &
        vn_mass = 'mass', vn_presh = 'pres', vn_betah = 'beta_vol',     &
        vn_buco = 'buco', vn_bvco = 'bvco', vn_vp = 'vp',               &
        vn_specw = 'specw', vn_phip = 'phips', vn_jdotb = 'jdotb',      &
        vn_bdotb = 'bdotb', vn_overr = 'over_r',                        &
        vn_bgrv = 'bdotgradv', vn_merc = 'DMerc', vn_mshear = 'DShear', &
        vn_mwell = 'DWell', vn_mcurr = 'DCurr', vn_mgeo = 'DGeod',      &
        vn_equif = 'equif', vn_fsq = 'fsqt', vn_wdot = 'wdot',          &
        vn_ftolv = 'ftolv', vn_fsql= 'fsql', vn_fsqr = 'fsqr',          &
        vn_fsqz = 'fsqz',                                               &
        vn_extcur = 'extcur', vn_curlab = 'curlabel', vn_rmnc = 'rmnc', &
        vn_zmns = 'zmns', vn_lmns = 'lmns', vn_gmnc = 'gmnc',           &
        vn_bmnc = 'bmnc', vn_bsubumnc = 'bsubumnc',                     &
        vn_bsubvmnc = 'bsubvmnc', vn_bsubsmns = 'bsubsmns',             &
        vn_bsupumnc = 'bsupumnc', vn_bsupvmnc = 'bsupvmnc',             &
        vn_rmns = 'rmns', vn_zmnc = 'zmnc',                             &
        vn_lmnc = 'lmnc', vn_gmns = 'gmns', vn_bmns = 'bmns',           &
        vn_bsubumns = 'bsubumns', vn_bsubvmns = 'bsubvmns',             &
        vn_bsubsmnc = 'bsubsmnc', vn_bsupumns = 'bsupumns',             &
        vn_currumnc = 'currumnc', vn_currumns = 'currumns',             &
        vn_currvmnc = 'currvmnc', vn_currvmns = 'currvmns',             &
        vn_bsupvmns = 'bsupvmns',                                       &
        vn_bsubumnc_sur = 'bsubumnc_sur',                               &
        vn_bsubvmnc_sur = 'bsubvmnc_sur',                               &
        vn_bsupumnc_sur = 'bsupumnc_sur',                               &
        vn_bsupvmnc_sur = 'bsupvmnc_sur',                               &
        vn_bsubumns_sur = 'bsubumns_sur',                               &
        vn_bsubvmns_sur = 'bsubvmns_sur',                               &
        vn_bsupumns_sur = 'bsupumns_sur',                               &
        vn_bsupvmns_sur = 'bsupvmns_sur',                               &
        vn_rbc = 'rbc', vn_zbs = 'zbs', vn_rbs = 'rbs', vn_zbc = 'zbc', &

        vn_wpar = 'wpar', vn_pparmnc = 'pparmnc', vn_ppermnc ='ppermnc',&
        vn_hotdmnc = 'hotdmnc', vn_pbprmnc = 'pbprmnc',                 &
        vn_ppprmnc = 'ppprmnc', vn_sigmnc  = 'sigmnc',                  &
        vn_taumnc  = 'taumnc',                                          &
        vn_pparmns = 'pparmns', vn_ppermns = 'ppermns',                 &
        vn_hotdmns = 'hotdmns', vn_pbprmns = 'pbprmns',                 &
        vn_ppprmns = 'ppprmns', vn_sigmns  = 'sigmns',                  &
        vn_taumns  = 'taumns',                                          &

        vn_machsq = 'machsq',                                           &
        vn_protmnc = 'protmnc', vn_protrsqmnc = 'protrsqmnc',           &
        vn_prprmnc = 'prprmnc',                                         &
        vn_protmns = 'protmns', vn_protrsqmns = 'protrsqmns',           &
        vn_prprmns = 'prprmns',                                         &
        vn_pmap = 'pmap', vn_omega = 'omega', vn_tpotb = 'tpotb'        

! Long names (ln_...)
      CHARACTER(LEN=*), PARAMETER :: ln_version = 'VMEC Version',       &
        ln_extension = 'Input file extension', ln_mgrid = 'MGRID file', &
        ln_magen = 'Magnetic Energy', ln_therm = 'Thermal Energy',      &
        ln_gam = 'Gamma', ln_maxr = 'Maximum R', ln_minr = 'Minimum R', &
        ln_maxz = 'Maximum Z', ln_fp = 'Field Periods',                 &
        ln_radnod = 'Radial nodes', ln_polmod = 'Poloidal modes',       &
        ln_tormod = 'Toroidal modes', ln_maxmod = 'Fourier modes',      &
        ln_maxmod_nyq = 'Max # Fourier modes (Nyquist)',                &
        ln_maxpot = 'Max # Fourier modes (vacuum potential)',           &
        ln_potsin = 'Vacuum potential sin modes',                       &
        ln_potcos = 'Vacuum potential cos modes',                       &
        ln_xmpot = 'Vacuum potential poloidal modes',                   &
        ln_xnpot = 'Vacuum potential toroidal modes',                   &
        ln_maxit = 'Max iterations', ln_actit = 'Actual iterations',    &
        ln_asym = 'Asymmetry', ln_recon = 'Reconstruction',             &
        ln_free = 'Free boundary',                                      &
        ln_error = 'Error flag', ln_aspect = 'Aspect ratio',            &
        ln_beta = 'Total beta', ln_pbeta = 'Poloidal beta',             &
        ln_tbeta = 'Toroidal beta', ln_abeta = 'Beta axis',             &
        ln_b0 = 'RB-t over R axis', ln_rbt0 = 'RB-t axis',              &
        ln_rbt1 = 'RB-t edge', ln_sgs = 'Sign jacobian',                &
        ln_lar = 'Ion Larmor radius', ln_modB = 'avg mod B',            &
        ln_ctor = 'Toroidal current', ln_amin = 'minor radius',         &
        ln_Rmaj = 'major radius', ln_vol = 'Plasma volume',             &
        ln_mse = 'Number of MSE points',                                &
        ln_thom = 'Number of Thompson scattering points',               &
        ln_am = 'Specification parameters for mass(s)',                 &
        ln_ac = 'Specification parameters for <J>(s)',                  &
        ln_ai = 'Specification parameters for iota(s)',                 &
        ln_pmass_type = 'Profile type specifier for mass(s)',           &
        ln_pcurr_type = 'Profile type specifier for <J>(s)',            &
        ln_piota_type = 'Profile type specifier for iota(s)',           &
        ln_am_aux_s = 'Auxiliary-s parameters for mass(s)',             &
        ln_am_aux_f = 'Auxiliary-f parameters for mass(s)',             &
        ln_ac_aux_s = 'Auxiliary-s parameters for <J>(s)',              &
        ln_ac_aux_f = 'Auxiliary-f parameters for <J>(s)',              &
        ln_ai_aux_s = 'Auxiliary-s parameters for iota(s)',             &
        ln_ai_aux_f = 'Auxiliary-f parameters for iota(s)',             &
        ln_pmod = 'Poloidal mode numbers',                              &
        ln_tmod = 'Toroidal mode numbers',                              &
        ln_pmod_nyq = 'Poloidal mode numbers (Nyquist)',                &
        ln_tmod_nyq = 'Toroidal mode numbers (Nyquist)',                &
        ln_racc = 'raxis (cosnv)', ln_racs = 'raxis (sinnv)',           &
        ln_zacs = 'zaxis (sinnv)', ln_zacc = 'zaxis (cosnv)',           &
        ln_iotaf = 'iota on full mesh',                                 &
        ln_qfact = 'q-factor on full mesh',                             &

        ln_presf = 'pressure on full mesh',                             &
        ln_phi = 'Toroidal flux on full mesh',                          &
        ln_phipf = 'd(phi)/ds: Toroidal flux deriv on full mesh',       &
        ln_chi = 'Poloidal flux on full mesh',                          &

        ln_chipf = 'd(chi)/ds: Poroidal flux deriv on full mesh',       &

        ln_jcuru = 'j dot gradu full',                                  &
        ln_jcurv = 'j dot gradv full', ln_iotah = 'iota half',          &
        ln_mass = 'mass half', ln_presh = 'pressure half',              &
        ln_betah = 'beta half', ln_buco = 'bsubu half',                 &
        ln_bvco = 'bsubv half', ln_vp = 'volume deriv half',            &
        ln_specw = 'Spectral width half',                               &
        ln_phip = 'tor flux deriv over 2pi half',                       &
        ln_jdotb = 'J dot B', ln_bgrv = 'B dot grad v',                 &
        ln_bdotb = 'B dot B',                                           &
        ln_merc = 'Mercier criterion', ln_mshear = 'Shear Mercier',     &
        ln_mwell = 'Well Mercier', ln_mcurr = 'Current Mercier',        &
        ln_mgeo = 'Geodesic Mercier', ln_equif='Average force balance', &
        ln_fsq = 'Residual decay',                                      &
        ln_wdot = 'Wdot decay', ln_extcur = 'External coil currents',   &
        ln_fsqr = 'Residual decay - radial',                            &
        ln_fsqz = 'Residual decay - vertical',                          &
        ln_fsql = 'Residual decay - hoop',                              &
        ln_ftolv = 'Residual decay - requested',                        &
        ln_curlab = 'External current names',                           &

        ln_rmnc = 'cosmn component of cylindrical R, full mesh',        &
        ln_zmns = 'sinmn component of cylindrical Z, full mesh',        &
        ln_lmns = 'sinmn component of lambda, half mesh',               &
        ln_gmnc = 'cosmn component of jacobian, half mesh',             &
        ln_bmnc = 'cosmn component of mod-B, half mesh',                &
        ln_bsubumnc = 'cosmn covariant u-component of B, half mesh',    &
        ln_bsubvmnc = 'cosmn covariant v-component of B, half mesh',    &
        ln_bsubsmns = 'sinmn covariant s-component of B, half mesh',    &

        ln_bsubumnc_sur = 'cosmn bsubu of B, surface',                  &
        ln_bsubvmnc_sur = 'cosmn bsubv of B, surface',                  &
        ln_bsupumnc_sur = 'cosmn bsupu of B, surface',                  &
        ln_bsupvmnc_sur = 'cosmn bsupv of B, surface',                  &

        ln_bsupumnc = 'BSUPUmnc half', ln_bsupvmnc = 'BSUPVmnc half',   & 

        ln_rmns = 'sinmn component of cylindrical R, full mesh',        &
        ln_zmnc = 'cosmn component of cylindrical Z, full mesh',        &
        ln_lmnc = 'cosmn component of lambda, half mesh',               &
        ln_gmns = 'sinmn component of jacobian, half mesh',             &
        ln_bmns = 'sinmn component of mod-B, half mesh',                &
        ln_bsubumns = 'sinmn covariant u-component of B, half mesh',    &
        ln_bsubvmns = 'sinmn covariant v-component of B, half mesh',    &
        ln_bsubsmnc = 'cosmn covariant s-component of B, half mesh',    &

        ln_currumnc = 'cosmn covariant u-component of J, full mesh',    &
        ln_currumns = 'sinmn covariant u-component of J, full mesh',    &
        ln_currvmnc = 'cosmn covariant v-component of J, full mesh',    &
        ln_currvmns = 'sinmn covariant v-component of J, full mesh',    &

        ln_bsubumns_sur = 'sinmn bsubu of B, surface',                  &
        ln_bsubvmns_sur = 'sinmn bsubv of B, surface',                  &
        ln_bsupumns_sur = 'sinmn bsupu of B, surface',                  &
        ln_bsupvmns_sur = 'sinmn bsupv of B, surface',                  &

        ln_bsupumns = 'BSUPUmns half', ln_bsupvmns = 'BSUPVmns half',   &
        ln_rbc = 'Initial boundary R cos(mu-nv) coefficients',          &
        ln_zbs = 'Initial boundary Z sin(mu-nv) coefficients',          &
        ln_rbs = 'Initial boundary R sin(mu-nv) coefficients',          &
        ln_zbc = 'Initial boundary Z cos(mu-nv) coefficients',          &

        ln_wpar = 'Energy',                                             &
        ln_pparmnc = 'cosmn compoents of hot part. para. pressure',     &
        ln_ppermnc = 'cosmn compoents of hot part. perp. pressure',     &
        ln_hotdmnc = 'cosmn compoents of hot part. density',            &
        ln_pbprmnc = 'cosmn compoents of hot part. para. pres. grad.',  &
        ln_ppprmnc = 'cosmn compoents of hot part. perp. pres. grad.',  &
        ln_sigmnc  = 'cosmn firehose stability variable',               &
        ln_taumnc  = 'cosmn mirror stability variable',                 &
        ln_pparmns = 'sinmn compoents of hot part. para. pressure',     &
        ln_ppermns = 'sinmn compoents of hot part. perp. pressure',     &
        ln_hotdmns = 'sinmn compoents of hot part. density',            &
        ln_pbprmns = 'sinmn compoents of hot part. para. pres. grad.',  &
        ln_ppprmns = 'sinmn compoents of hot part. perp. pres. grad.',  &
        ln_sigmns  = 'sinmn firehose stability variable',               &
        ln_taumns  = 'sinmn mirror stability variable',                 &

        ln_machsq = 'Mach # on axis (squared)',                         &
        ln_protmnc = 'cosmn components of pressure',                    &
        ln_protrsqmnc = 'cosmn component of rotational energy',         &
        ln_prprmnc = 'cosmn components of radial pressure gradient',    &
        ln_protmns = 'sinmn components of pressure',                    &
        ln_protrsqmns = 'sinmn component of rotational energy',         &
        ln_prprmns = 'sinmn components of radial pressure gradient',    &
        ln_pmap = '<p(s,R)>', ln_omega = 'Toroidal Angular Freq.',      &
        ln_tpotb = 'T_perp/T_parallel or T(flow)'                   
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nfp, ns, mpol, ntor, mnmax, mnmax_nyq, itfsq, niter,   &
          iasym, ireconstruct, ierr_vmec, imse, itse, nstore_seq,       &
          isnodes, ipnodes, imatch_phiedge, isigng, mnyq, nnyq, ntmax,  &
          mnmaxpot, vmec_type
      REAL(rprec) :: wb, wp, gamma, pfac, rmax_surf, rmin_surf,         &
          zmax_surf, aspect, betatot, betapol, betator, betaxis, b0,    &
          tswgt, msewgt, flmwgt, bcwgt, phidiam, version_,              &
          delphid, IonLarmor, VolAvgB,                                  &
          fsql, fsqr, fsqz, ftolv,                                      &
          Aminor, Rmajor, Volume, RBtor, RBtor0, Itor, machsq
      REAL(rprec), ALLOCATABLE :: rzl_local(:,:,:,:)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
          rmnc, zmns, lmns, rmns, zmnc, lmnc, bmnc, gmnc, bsubumnc,     &
          bsubvmnc, bsubsmns, bsupumnc, bsupvmnc, currvmnc,             &
          currumnc, bbc, raxis, zaxis 
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
          bmns, gmns, bsubumns, bsubvmns, bsubsmnc,                     &
          bsupumns, bsupvmns, currumns, currvmns
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
          pparmnc, ppermnc, hotdmnc, pbprmnc, ppprmnc, sigmnc, taumnc,  & ! SAL - ANIMEC
          pparmns, ppermns, hotdmns, pbprmns, ppprmns, sigmns, taumns,  & ! SAL - ANIMEC
          protmnc, protrsqmnc, prprmnc,                                 & ! SAL - FLOW
          protmns, protrsqmns, prprmns                                    ! SAL - FLOW
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::                         &
          bsubumnc_sur, bsubumns_sur, bsubvmnc_sur, bsubvmns_sur,       &
          bsupumnc_sur, bsupumns_sur, bsupvmnc_sur, bsupvmns_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::                         &
         iotas, iotaf, presf, phipf, mass, pres, beta_vol, xm, xn,      &
         potsin, potcos, xmpot, xnpot, qfact, chipf, phi, chi,          &
         xm_nyq, xn_nyq, phip, buco, bvco, vp, overr, jcuru, jcurv,     & 
         specw, jdotb, bdotb, bdotgradv, fsqt, wdot, am, ac, ai,        & 
         am_aux_s, am_aux_f, ac_aux_s, ac_aux_f, ai_aux_s, ai_aux_f,    &
         Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, extcur,             &
         sknots, ystark, y2stark, pknots, ythom, y2thom,                &
         anglemse, rmid, qmid, shear, presmid, alfa, curmid, rstark,    &
         qmeas, datastark, rthom, datathom, dsiobt      
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: pmap, omega, tpotb        ! SAL -FLOW
      LOGICAL :: lasym, lthreed, lwout_opened=.false.
      CHARACTER :: mgrid_file*200, input_extension*100
      CHARACTER :: pmass_type*20, pcurr_type*20, piota_type*20

      INTEGER, PARAMETER :: norm_term_flag=0,                           &
         bad_jacobian_flag=1, more_iter_flag=2, jac75_flag=4

!     OVERLOAD SUBROUTINE READ_WOUT_FILE TO ACCEPT BOTH UNIT NO. (OPENED EXTERNALLY)
!     OR FILENAME (HANDLE OPEN/CLOSE HERE)
      INTERFACE read_wout_file
          MODULE PROCEDURE readw_and_open, readw_only
      END INTERFACE

#if defined(NETCDF)
      PRIVATE :: read_wout_text, read_wout_nc, write_wout_text, write_wout_nc
#else
      PRIVATE :: read_wout_text, write_wout_text
#endif
      PRIVATE :: norm_term_flag, bad_jacobian_flag,                     &
                 more_iter_flag, jac75_flag

      CONTAINS

      SUBROUTINE readw_and_open(file_or_extension, ierr, iopen)
      USE safe_open_mod
      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      INTEGER, INTENT(out) :: ierr
      INTEGER, OPTIONAL :: iopen
      CHARACTER(LEN=*), INTENT(in) :: file_or_extension
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER, PARAMETER :: iunit_init = 10
      INTEGER :: iunit
      LOGICAL :: isnc
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
!------------------------------------------------
!
!     THIS SUBROUTINE READS THE WOUT FILE CREATED BY THE VMEC CODE
!     AND STORES THE DATA IN THE READ_WOUT MODULE
!
!     FIRST, CHECK IF THIS IS A FULLY-QUALIFIED PATH NAME
!     MAKE SURE wout IS NOT EMBEDDED IN THE NAME (PERVERSE USER...)
!
      filename = 'wout'
      CALL parse_extension(filename, file_or_extension, isnc)
      CALL flush(6)
!SPH  IF (.not.isnc) STOP 'ISNC ERR IN READ_WOUT_MOD'
      IF (isnc) THEN
#if defined(NETCDF)
         CALL read_wout_nc(filename, ierr)
#else
         PRINT *, "NETCDF wout file can not be opened on this platform"
         ierr = -100
#endif
      ELSE
         iunit = iunit_init
         CALL safe_open (iunit, ierr, filename, 'old', 'formatted')
         IF (ierr .eq. 0) CALL read_wout_text(iunit, ierr)
         CLOSE(unit=iunit)
      END IF
      
      IF (PRESENT(iopen)) iopen = ierr
      lwout_opened = (ierr .eq. 0)
      ! WHEN READING A NETCDF FILE, A BAD RUN MAY PREVENT XN FROM BEING
      ! READ, SUBSEQUENTLY WE MUST CHECK TO SEE IF XN HAS BEEN ALLOCATED
      ! BEFORE DOING ANYTHING WITH IT OTHERWISE WE DEFAULT LTHREED TO
      ! FALSE.  - SAL 09/07/11
      IF (ALLOCATED(XN)) THEN
         lthreed = ANY(NINT(xn) .ne. 0)
      ELSE
         lthreed = .FALSE.
      END IF

      END SUBROUTINE readw_and_open

      SUBROUTINE write_wout_file(file_or_extension, ierr)
      USE safe_open_mod
      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      INTEGER, INTENT(out) :: ierr
      CHARACTER(LEN=*), INTENT(in) :: file_or_extension
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER, PARAMETER :: iunit_init = 10
      INTEGER :: iunit
      LOGICAL :: isnc
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
!------------------------------------------------
!
!     THIS SUBROUTINE WRITES THE WOUT FILE FROM 
!     THE DATA IN THE READ_WOUT MODULE
!
!     FIRST, CHECK IF THIS IS A FULLY-QUALIFIED PATH NAME
!     MAKE SURE wout IS NOT EMBEDDED IN THE NAME (PERVERSE USER...)
!
      filename = 'wout'
      CALL parse_extension(filename, file_or_extension, isnc)
      IF (isnc) THEN
#if defined(NETCDF)
         CALL write_wout_nc(filename, ierr)
#else
         PRINT *, "NETCDF wout file can not be opened on this platform"
         ierr = -100
#endif
      ELSE
         IF (ierr .eq. 0) CALL write_wout_text(filename, ierr)
      END IF

      END SUBROUTINE write_wout_file


      SUBROUTINE readw_only(iunit, ierr, iopen)
      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(out):: ierr
      INTEGER, OPTIONAL :: iopen
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER :: istat
      CHARACTER(LEN=256) :: vmec_version
      LOGICAL :: exfile
!------------------------------------------------
!
!     User opened the file externally and has a unit number, iunit
!
      ierr = 0

      INQUIRE(unit=iunit, exist=exfile, name=vmec_version,iostat=istat)
      IF (istat.ne.0 .or. .not.exfile) THEN
        PRINT *,' In READ_WOUT_FILE, Unit = ',iunit,                    &
                ' File = ',TRIM(vmec_version),' DOES NOT EXIST'
        IF (PRESENT(iopen)) iopen = -1
        ierr = -1
        RETURN
      ELSE
        IF (PRESENT(iopen)) iopen = 0
      END IF

      CALL read_wout_text(iunit, ierr)
      lwout_opened = (ierr .eq. 0)
      lthreed = ANY(NINT(xn) .ne. 0)

      END SUBROUTINE readw_only


      SUBROUTINE read_wout_text(iunit, ierr)
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      INTEGER :: iunit, ierr
!------------------------------------------------
!   L o c a l   P a r a m e t e r s
!------------------------------------------------
      REAL(rprec), PARAMETER :: eps_w = 1.e-4_dp
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER :: istat(15), i, j, k, js, m, n, n1, mn, nparts_in, &
                 i_animec, i_flow
      CHARACTER(LEN=256) :: vmec_version
      LOGICAL            :: lcurr
!------------------------------------------------
!
!     THIS SUBROUTINE READS THE TEXT FILE WOUT CREATED BY THE VMEC CODE
!     AND STORES THE INFORMATION IN THE read_WOUT MODULE
!
!     CALL read_wout_file - GENERIC INTERFACE - CAN BE CALLED WITH EITHER UNIT NO. OR FILENAME
!
!     RMNC, ZMNS: FULL-GRID
!     LMNS      : HALF-GRID
!
      istat = 0
      ierr = 0
      nextcur = 0

      READ (iunit, '(a)', iostat=istat(2), err=1000) vmec_version

      i = INDEX(vmec_version,'=')
      !!!! ADDED BY SAL
      i_animec = INDEX(vmec_version,'_ANIMEC')
      i_flow = INDEX(vmec_version,'_FLOW')
      vmec_type = 0
      IF (i_animec > 0) THEN
         vmec_type = 1 ! ANIMEC
         vmec_version = vmec_version(1:i_animec-1)
      END IF
      IF (i_flow > 0) THEN
         vmec_type = 2 ! FLOW
         vmec_version = vmec_version(1:i_flow-1)
      END IF
      !!!! END SAL Addition


      IF (i .ge. 0) THEN
         READ(vmec_version(i+1:LEN_TRIM(vmec_version)),*) version_
      ELSE
         version_ = -1.0
      END IF

      ierr_vmec = norm_term_flag

      IF (version_ .le. (5.10 + eps_w)) THEN
         READ (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma,      &
            pfac, nfp, ns, mpol, ntor, mnmax, itfsq, niter, iasym,      &
            ireconstruct
      ELSE
         IF (version_ .lt. 6.54) THEN
            READ (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma,   &
              pfac, rmax_surf, rmin_surf
         ELSE
            READ (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma,   &
               pfac, rmax_surf, rmin_surf, zmax_surf
         END IF
         IF (vmec_type == 2) THEN !SAL
            READ(iunit, *,iostat=istat(2), err=1000) machsq
         END IF
         IF (version_ .le. (8.0+eps_w)) THEN
            READ (iunit, *, iostat=istat(2), err=1000) nfp, ns, mpol,   &
            ntor, mnmax, itfsq, niter, iasym, ireconstruct, ierr_vmec
            mnmax_nyq = mnmax
         ELSE
            READ (iunit, *, iostat=istat(2), err=1000) nfp, ns, mpol,   &
            ntor, mnmax, mnmax_nyq, itfsq, niter, iasym, ireconstruct,  &
            ierr_vmec
         END IF
      END IF

      lasym = (iasym .gt. 0)

      IF (version_ .gt. (6.20+eps_w)) THEN
         READ (iunit, *, iostat=istat(1), err=1000)imse, itse, nbsets,  &
            nobd, nextcur, nstore_seq
      ELSE
         READ (iunit, *, iostat=istat(1), err=1000)imse, itse, nbsets,  &
               nobd, nextcur
         nstore_seq = 100
      END IF
     
      IF ((ierr_vmec .ne. norm_term_flag) .and.                         &
          (ierr_vmec .ne. more_iter_flag) .and.                         &
          (ierr_vmec .ne. jac75_flag)) THEN
        ierr = -2
        GOTO 1000
      ENDIF

      IF (nextcur .gt. nigroup) istat(15) = -1

      IF (ALLOCATED(xm)) CALL read_wout_deallocate

      ALLOCATE (xm(mnmax), xn(mnmax), xm_nyq(mnmax_nyq),                &
        xn_nyq(mnmax_nyq),rmnc(mnmax,ns), zmns(mnmax,ns),               &
        lmns(mnmax,ns), bmnc(mnmax_nyq,ns), gmnc(mnmax_nyq,ns),         &
        bsubumnc(mnmax_nyq,ns), bsubvmnc(mnmax_nyq,ns),                 &
        bsubsmns(mnmax_nyq,ns), bsupumnc(mnmax_nyq,ns),                 &
        bsupvmnc(mnmax_nyq,ns), currvmnc(mnmax_nyq,ns),                 &
        currumnc(mnmax_nyq,ns),                                         &
        iotas(ns), mass(ns), pres(ns), beta_vol(ns), phip(ns),          &
        buco(ns), bvco(ns), phi(ns), iotaf(ns), presf(ns), phipf(ns),   &
        chipf(ns),                                                      &
        vp(ns), overr(ns), jcuru(ns), jcurv(ns), specw(ns), Dmerc(ns),  &
        Dshear(ns), Dwell(ns), Dcurr(ns), Dgeod(ns), equif(ns),         &
        raxis(0:ntor,2), zaxis(0:ntor,2), jdotb(ns), bdotb(ns),         &
        bdotgradv(ns), am(0:20), ac(0:20),  ai(0:20),                   &
        fsqt(nstore_seq), wdot(nstore_seq), stat = istat(6))

      IF (nextcur .GT. 0) ALLOCATE(extcur(nextcur), curlabel(nextcur),  &
        stat = istat(6))

      IF (lasym)                                                        &
         ALLOCATE (rmns(mnmax,ns), zmnc(mnmax,ns), lmnc(mnmax,ns),      &
                   bmns(mnmax_nyq,ns), gmns(mnmax_nyq,ns),              &
                   bsubumns(mnmax_nyq,ns),                              &
                   bsubvmns(mnmax_nyq,ns), bsubsmnc(mnmax_nyq,ns),      &
                   bsupumns(mnmax_nyq,ns), bsupvmns(mnmax_nyq,ns),      &
                   currumns(mnmax_nyq,ns), currvmns(mnmax_nyq,ns),      &
                   stat=istat(6))
     
      IF (vmec_type == 1) THEN   ! SAL
         ALLOCATE (pparmnc(mnmax_nyq,ns),ppermnc(mnmax_nyq,ns),         &
                   hotdmnc(mnmax_nyq,ns),pbprmnc(mnmax_nyq,ns),         &
                   ppprmnc(mnmax_nyq,ns),sigmnc(mnmax_nyq,ns),          &
                   taumnc(mnmax_nyq,ns),stat=istat(6))
         IF (lasym)                                                     &
            ALLOCATE (pparmns(mnmax_nyq,ns),ppermns(mnmax_nyq,ns),      &
                      hotdmns(mnmax_nyq,ns),pbprmns(mnmax_nyq,ns),      &
                      ppprmns(mnmax_nyq,ns),sigmns(mnmax_nyq,ns),       &
                      taumns(mnmax_nyq,ns),stat=istat(6))
      ELSE IF (vmec_type == 2) THEN
         ALLOCATE (pmap(ns), omega(ns), tpotb(ns),stat=istat(6))
         ALLOCATE (protmnc(mnmax_nyq,ns),protrsqmnc(mnmax_nyq,ns),      &
                   prprmnc(mnmax_nyq,ns),stat=istat(6))
         IF (lasym)                                                     &
            ALLOCATE (protmns(mnmax_nyq,ns),protrsqmns(mnmax_nyq,ns),   &
                      prprmns(mnmax_nyq,ns),stat=istat(6))
      END IF

      fsqt = 0; wdot = 0; raxis = 0; zaxis = 0

      IF (nbsets .gt. 0) READ (iunit, *, iostat=istat(4), err=1000)     &
         (nbfld(i),i=1,nbsets)
      READ (iunit, '(a)', iostat=istat(5), err=1000) mgrid_file

      DO js = 1, ns
         DO mn = 1, mnmax
            IF(js .eq. 1) THEN
               READ (iunit, *, iostat=istat(7), err=1000) m, n
               xm(mn) = REAL(m,rprec)
               xn(mn) = REAL(n,rprec)
            END IF
            IF (version_ .le. (6.20+eps_w)) THEN
              READ (iunit, 730, iostat=istat(8), err=1000)              &
              rmnc(mn,js), zmns(mn,js), lmns(mn,js), bmnc(mn,js),       &
              gmnc(mn,js), bsubumnc(mn,js), bsubvmnc(mn,js),            &
              bsubsmns(mn,js), bsupumnc(mn,js), bsupvmnc(mn,js),        &
              currvmnc(mn,js)
            ELSE IF (version_ .le. (8.0+eps_w)) THEN
              READ (iunit, *, iostat=istat(8), err=1000)                &
              rmnc(mn,js), zmns(mn,js), lmns(mn,js), bmnc(mn,js),       &
              gmnc(mn,js), bsubumnc(mn,js), bsubvmnc(mn,js),            &
              bsubsmns(mn,js), bsupumnc(mn,js), bsupvmnc(mn,js),        &
              currvmnc(mn,js)
            ELSE
              READ (iunit, *, iostat=istat(8), err=1000)                &
              rmnc(mn,js), zmns(mn,js), lmns(mn,js)
            END IF
            IF (lasym) THEN
               IF (version_ .le. (8.0+eps_w)) THEN
                  READ (iunit, *, iostat=istat(8), err=1000)            &
                  rmns(mn,js), zmnc(mn,js), lmnc(mn,js),                &
                  bmns(mn,js), gmns(mn,js), bsubumns(mn,js),            &
                  bsubvmns(mn,js), bsubsmnc(mn,js),                     &
                  bsupumns(mn,js), bsupvmns(mn,js)
               ELSE
                  READ (iunit, *, iostat=istat(8), err=1000)            &
                  rmns(mn,js), zmnc(mn,js), lmnc(mn,js)
               END IF
            END IF
            IF (js.eq.1 .and. m.eq.0) THEN
               n1 = ABS(n/nfp)
               IF (n1 .le. ntor) THEN
                  raxis(n1,1) = rmnc(mn,1)
                  zaxis(n1,1) = zmns(mn,1)
                  IF (lasym) THEN
                     raxis(n1,2) = rmns(mn,1)
                     zaxis(n1,2) = zmnc(mn,1)
                  END IF
               END IF
            END IF
         END DO

         IF (version_ .le. (8.0+eps_w)) CYCLE
         DO mn = 1, mnmax_nyq
            IF(js .eq. 1) THEN
               READ (iunit, *, iostat=istat(7), err=1000) m, n
               xm_nyq(mn) = REAL(m,rprec)
               xn_nyq(mn) = REAL(n,rprec)
            END IF
            IF (vmec_type == 1) THEN  !SAL (ELSE statement below is orriginal)
               READ (iunit, *, iostat=istat(8), err=1000)               &
               bmnc(mn,js), gmnc(mn,js), bsubumnc(mn,js),               &
               bsubvmnc(mn,js), bsubsmns(mn,js),                        &
               bsupumnc(mn,js), bsupvmnc(mn,js),                        &
               pparmnc (mn,js), ppermnc (mn,js), hotdmnc (mn,js),       &
               pbprmnc (mn,js), ppprmnc (mn,js), sigmnc  (mn,js),       &
               taumnc  (mn,js)
               IF (lasym) THEN
                  READ (iunit, *, iostat=istat(8), err=1000)            &
                  bmns(mn,js), gmns(mn,js), bsubumns(mn,js),            &
                  bsubvmns(mn,js), bsubsmnc(mn,js),                     &
                  bsupumns(mn,js), bsupvmns(mn,js),                     &
                  pparmns (mn,js), ppermns (mn,js), hotdmns (mn,js),    &
                  pbprmns (mn,js), ppprmns (mn,js), sigmns  (mn,js),    &
                  taumns  (mn,js)
               END IF
            ELSE IF (vmec_type ==2) THEN
               READ (iunit, *, iostat=istat(8), err=1000)               &
               bmnc(mn,js), gmnc(mn,js), bsubumnc(mn,js),               &
               bsubvmnc(mn,js), bsubsmns(mn,js),                        &
               bsupumnc(mn,js), bsupvmnc(mn,js),                        &
               protmnc (mn,js), protrsqmnc(mn,js), prprmnc(mn,js)
               IF (lasym) THEN
                  READ (iunit, *, iostat=istat(8), err=1000)            &
                  bmns(mn,js), gmns(mn,js), bsubumns(mn,js),            &
                  bsubvmns(mn,js), bsubsmnc(mn,js),                     &
                  bsupumns(mn,js), bsupvmns(mn,js),                     &
                  protmns (mn,js), protrsqmns(mn,js), prprmns(mn,js)
               END IF
            ELSE
               READ (iunit, *, iostat=istat(8), err=1000)               &
                 bmnc(mn,js), gmnc(mn,js), bsubumnc(mn,js),             &
                 bsubvmnc(mn,js), bsubsmns(mn,js),                      &
                 bsupumnc(mn,js), bsupvmnc(mn,js)
               IF (lasym) THEN
                  READ (iunit, *, iostat=istat(8), err=1000)            &
                  bmns(mn,js), gmns(mn,js), bsubumns(mn,js),            &
                  bsubvmns(mn,js), bsubsmnc(mn,js),                     &
                  bsupumns(mn,js), bsupvmns(mn,js)
            END IF
            END IF
         END DO

      END DO

!     Compute current coefficients on full mesh
      IF (version_ .gt. (8.0+eps_w)) THEN
         CALL Compute_Currents(bsubsmnc, bsubsmns,                             &
     &                         bsubumnc, bsubumns,                             &
     &                         bsubvmnc, bsubvmns,                             &
     &                         xm_nyq, xn_nyq, mnmax_nyq, lasym, ns,           &
     &                         currumnc, currvmnc,                             &
     &                         currumns, currvmns)
      END IF

      mnyq = INT(MAXVAL(xm_nyq));  nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp

!
!     Read FULL AND HALF-MESH QUANTITIES 
!
!     NOTE: In version_ <= 6.00, mass, press were written out in INTERNAL (VMEC) units
!     and are therefore multiplied here by 1/mu0 to transform to pascals. Same is true
!     for ALL the currents (jcuru, jcurv, jdotb). Also, in version_ = 6.10 and
!     above, PHI is the true (physical) toroidal flux (has the sign of jacobian correctly
!     built into it)
!
      iotas(1) = 0; mass(1) = 0; pres(1) = 0; phip(1) = 0; 
      buco(1) = 0; bvco(1) = 0; vp(1) = 0; overr(1) = 0;  specw(1) = 1
      beta_vol(1) = 0

      IF (version_ .le. (6.05+eps_w)) THEN
         READ (iunit, 730, iostat=istat(9), err=1000)                   &
           (iotas(js), mass(js), pres(js),                              &
            phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),   &
            jcuru(js), jcurv(js), specw(js),js=2,ns)
         READ (iunit, 730, iostat=istat(10), err=1000)                  &
            aspect, betatot, betapol, betator, betaxis, b0
      ELSE IF (version_ .le. (6.20+eps_w)) THEN
         READ (iunit, 730, iostat=istat(9), err=1000)                   &
           (iotas(js), mass(js), pres(js), beta_vol(js),                &
            phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),   &
            jcuru(js), jcurv(js), specw(js),js=2,ns)
         READ (iunit, 730, iostat=istat(10), err=1000)                  &
            aspect, betatot, betapol, betator, betaxis, b0
      ELSE IF (version_ .le. (6.95+eps_w)) THEN
         READ (iunit, *, iostat=istat(9), err=1000)                     &
           (iotas(js), mass(js), pres(js), beta_vol(js),                &
            phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),   &
            jcuru(js), jcurv(js), specw(js),js=2,ns)
         READ (iunit, *, iostat=istat(10), err=1000)                    &
            aspect, betatot, betapol, betator, betaxis, b0
      ELSE
         READ (iunit, *, iostat=istat(9), err=1000)                     &
         (iotaf(js), presf(js), phipf(js), phi(js),                     &
         jcuru(js), jcurv(js), js=1,ns)
         READ (iunit, *, iostat=istat(9), err=1000)                     &
         (iotas(js), mass(js), pres(js),                                &
         beta_vol(js), phip(js), buco(js), bvco(js), vp(js),            &
         overr(js), specw(js),js=2,ns)
         READ (iunit, *, iostat=istat(10), err=1000)                    &
            aspect, betatot, betapol, betator, betaxis, b0
      END IF


      IF (version_ .gt. (6.10+eps_w)) THEN
         READ (iunit, *, iostat=istat(10), err=1000) isigng
         READ (iunit, *, iostat=istat(10), err=1000) input_extension
         READ (iunit, *, iostat=istat(10), err=1000) IonLarmor,         &
           VolAvgB, RBtor0, RBtor, Itor, Aminor, Rmajor, Volume
      END IF

!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      IF (version_.gt.(5.10+eps_w) .and. version_.lt.(6.20-eps_w)) THEN
         READ (iunit, 730, iostat=istat(11), err=1000)                  &
            (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),               &
             Dgeod(js), equif(js), js=2,ns-1)
      ELSE IF (version_ .ge. (6.20-eps_w)) THEN
         READ (iunit, *, iostat=istat(11), err=1000)                    &
            (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),               &
             Dgeod(js), equif(js), js=2,ns-1)
      END IF

      IF (nextcur .gt. 0) THEN
         IF (version_ .le. (6.20+eps_w)) THEN
            READ (iunit, 730, iostat=istat(12), err=1000)               &
            (extcur(i),i=1,nextcur)
         ELSE
            READ (iunit, *, iostat=istat(12), err=1000)                 &
            (extcur(i),i=1,nextcur)
         END IF
         READ (iunit, *, iostat=istat(13)) lcurr
         IF (lcurr) READ (iunit, *, iostat=istat(13), err=1000)         &
            (curlabel(i),i=1,nextcur)
      END IF

      IF (version_ .le. (6.20+eps_w)) THEN
         READ (iunit, 730, iostat=istat(14))                            &
            (fsqt(i), wdot(i), i=1,nstore_seq)
      ELSE
         READ (iunit, *, iostat=istat(14))                              &
           (fsqt(i), wdot(i), i=1,nstore_seq)
      END IF

      IF ((version_.ge.6.20-eps_w) .and. (version_ .lt. (6.50-eps_w))   &
         .and. (istat(14).eq.0)) THEN  
         READ (iunit, 730, iostat=istat(14), err=1000)                  &
           (jdotb(js), bdotgradv(js), bdotb(js), js=1,ns)
      ELSE IF (version_ .ge. (6.50-eps_w)) THEN
         READ (iunit, *, iostat=istat(14), err=1000)                    &
           (jdotb(js), bdotgradv(js), bdotb(js), js=1,ns)
      ELSE
         istat(14) = 0
      END IF

      chipf = iotaf*phipf
!
!     CONVERT FROM INTERNAL UNITS TO PHYSICAL UNITS IF NEEDED
!
      IF (version_ .le. (6.05+eps_w)) THEN
         mass = mass/mu0
         pres = pres/mu0
         jcuru = jcuru/mu0
         jcurv = jcurv/mu0
         jdotb = jdotb/mu0
         phi   = -phi
      END IF

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      IF (ireconstruct .gt. 0) THEN

        n1 = MAXVAL(nbfld(:nbsets))
        ALLOCATE (sknots(isnodes), ystark(isnodes), y2stark(isnodes),   &
           pknots(ipnodes), ythom(ipnodes), y2thom(ipnodes),            &
           anglemse(2*ns), rmid(2*ns), qmid(2*ns), shear(2*ns),         &
           presmid(2*ns), alfa(2*ns), curmid(2*ns), rstark(imse),       &
           datastark(imse), rthom(itse), datathom(itse),                &
           dsiext(nobd), plflux(nobd), dsiobt(nobd), bcoil(n1,nbsets),  &
           plbfld(n1,nbsets), bbc(n1,nbsets))
         IF (imse.ge.2 .or. itse.gt.0) THEN
            READ (iunit, *) tswgt, msewgt
            READ (iunit, *) isnodes, (sknots(i),ystark(i),y2stark(i),   &
               i=1,isnodes)
            READ (iunit, *) ipnodes, (pknots(i), ythom(i),              &
               y2thom(i),i=1,ipnodes)
            READ(iunit, *)(anglemse(i),rmid(i),qmid(i),shear(i),        &
            presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
            READ(iunit, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
            READ(iunit, *)(rthom(i),datathom(i),i=1,itse)
         END IF

         IF (nobd .gt. 0) THEN
            READ (iunit, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
            READ (iunit, *) flmwgt
         END IF

         nbfldn = SUM(nbfld(:nbsets))
         IF (nbfldn .gt. 0) THEN
            DO n = 1, nbsets
               READ (iunit, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),        &
                  i=1,nbfld(n))
            END DO
            READ (iunit, *) bcwgt
         END IF

         READ (iunit, *) phidiam, delphid
!
!     READ Limiter & Prout plotting specs
!
         READ (iunit, *) nsets, nparts_in, nlim

         ALLOCATE (nsetsn(nsets))
         READ (iunit, *) (nsetsn(i),i=1,nsets)

         n1 = MAXVAL(nsetsn(:nsets))
         ALLOCATE (pfcspec(nparts_in,n1,nsets), limitr(nlim))

         READ (iunit, *) (((pfcspec(i,j,k),i=1,nparts_in),              &
            j=1,nsetsn(k)),k=1,nsets)

         READ (iunit, *) (limitr(i), i=1,nlim)

         m  = MAXVAL(limitr(:nlim))
         ALLOCATE (rlim(m,nlim), zlim(m,nlim))

         READ (iunit, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),          &
            j=1,nlim)
         READ (iunit, *) nrgrid, nzgrid
         READ (iunit, *) tokid
         READ (iunit, *) rx1, rx2, zy1, zy2, condif
         READ (iunit, *) imatch_phiedge

      END IF

 1000 CONTINUE

      READ (iunit, iostat=ierr) mgrid_mode
      IF (ierr .ne. 0) THEN
         ierr = 0; mgrid_mode = 'N'
      END IF

      IF (istat(2) .ne. 0) ierr_vmec = 1

      DO m = 1,15
        IF (istat(m) .gt. 0) THEN
           PRINT *,' Error No. ',m,' in READ_WOUT, iostat = ',istat(m)
           ierr = m
           EXIT
        END IF
      END DO


  720 FORMAT(8i10)
  730 FORMAT(5e20.13)
  740 FORMAT(a)
  790 FORMAT(i5,/,(1p,3e12.4))

      END SUBROUTINE read_wout_text


#if defined(NETCDF)
      SUBROUTINE read_wout_nc(filename, ierr)
      USE ezcdf
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      INTEGER, INTENT(out) :: ierr
      CHARACTER(LEN=*), INTENT(in) :: filename
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER :: nwout, ierror, i_animec, i_flow
      INTEGER, DIMENSION(3)   :: dimlens
!      REAL(rprec) :: ohs
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_cc, raxis_cs,     &
                                                zaxis_cs, zaxis_cc
!------------------------------------------------
! Open cdf File
      CALL cdf_open(nwout,filename,'r', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening wout .nc file'
         RETURN
      END IF

! Be sure all arrays are deallocated
      CALL read_wout_deallocate
      
      ! ANIMEC/FLOW -SAL
      i_animec = 0
      i_flow   = 0
      vmec_type = 0
      
      CALL cdf_inquire(nwout, vn_pparmnc, dimlens, ier=ierror)
      IF (ierror.eq.0) vmec_type = 1
      CALL cdf_inquire(nwout, vn_omega, dimlens, ier=ierror)
      IF (ierror.eq.0) vmec_type = 2

! Read in scalar variables
      CALL cdf_read(nwout, vn_error, ierr_vmec)
            
      IF ((ierr_vmec .ne. norm_term_flag) .and.                         &
          (ierr_vmec .ne. more_iter_flag) .and.                         &
          (ierr_vmec .ne. jac75_flag)) THEN
        ierr = -2
        GOTO 1000
      ENDIF

      CALL cdf_read(nwout, vn_version, version_)
      CALL cdf_read(nwout, vn_extension, input_extension)
      CALL cdf_read(nwout, vn_mgrid, mgrid_file)
      CALL cdf_read(nwout, vn_magen, wb)
      CALL cdf_read(nwout, vn_therm, wp)
      CALL cdf_read(nwout, vn_gam, gamma)
      CALL cdf_read(nwout, vn_maxr, rmax_surf)
      CALL cdf_read(nwout, vn_minr, rmin_surf)
      CALL cdf_read(nwout, vn_maxz, zmax_surf)
      CALL cdf_read(nwout, vn_fp, nfp)
      CALL cdf_read(nwout, vn_radnod, ns)
      CALL cdf_read(nwout, vn_polmod, mpol)
      CALL cdf_read(nwout, vn_tormod, ntor)
      CALL cdf_read(nwout, vn_maxmod, mnmax)
      mnmaxpot=0
      CALL cdf_read(nwout, vn_maxpot, mnmaxpot)
      mnmax_nyq = -1
      CALL cdf_read(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_read(nwout, vn_maxit, niter)
      CALL cdf_read(nwout, vn_actit, itfsq)
      CALL cdf_read(nwout, vn_asym, lasym)
      IF (lasym) iasym = 1
      CALL cdf_read(nwout, vn_recon, lrecon)
      IF (lrecon) ireconstruct = 1
      CALL cdf_read(nwout, vn_free, lfreeb)
      CALL cdf_read(nwout, vn_moveaxis, lmove_axis)

      CALL cdf_read(nwout, vn_rfp, lrfp)
      CALL cdf_read(nwout, vn_aspect, aspect)
      CALL cdf_read(nwout, vn_beta, betatot)
      CALL cdf_read(nwout, vn_pbeta, betapol)
      CALL cdf_read(nwout, vn_tbeta, betator)
      CALL cdf_read(nwout, vn_abeta, betaxis)
      CALL cdf_read(nwout, vn_b0, b0)
      CALL cdf_read(nwout, vn_rbt0, rbtor0)
      CALL cdf_read(nwout, vn_rbt1, rbtor)
      CALL cdf_read(nwout, vn_sgs, isigng)
      CALL cdf_read(nwout, vn_lar, IonLarmor)
      CALL cdf_read(nwout, vn_modB, volAvgB)
      CALL cdf_read(nwout, vn_ctor, Itor)
      CALL cdf_read(nwout, vn_amin, Aminor)
      CALL cdf_read(nwout, vn_rmaj, Rmajor)
      CALL cdf_read(nwout, vn_vol, volume)
      CALL cdf_read(nwout, vn_ftolv, ftolv)
      CALL cdf_read(nwout, vn_fsqr, fsqr)
      CALL cdf_read(nwout, vn_fsqz, fsqz)
      CALL cdf_read(nwout, vn_fsql, fsql)
      CALL cdf_read(nwout, vn_pcurr_type, pcurr_type)
      CALL cdf_read(nwout, vn_piota_type, piota_type)
      CALL cdf_read(nwout, vn_pmass_type, pmass_type)
      imse = -1
      IF (lrecon) THEN
         CALL cdf_read(nwout, vn_mse, imse)
         CALL cdf_read(nwout, vn_thom, itse)
      END IF
      CALL cdf_read(nwout, vn_nextcur, nextcur)

      mgrid_mode = 'N'
      CALL cdf_inquire(nwout, vn_mgmode, dimlens, ier=ierror)
      IF (ierror.eq.0) CALL cdf_read(nwout, vn_mgmode, mgrid_mode)
      IF (lfreeb) THEN
         CALL cdf_read(nwout, vn_flp, nobser)
         CALL cdf_read(nwout, vn_nobd, nobd)
         CALL cdf_read(nwout, vn_nbset, nbsets)
      END IF
      
      ! ANIMEC/FLOW -SAL
      CALL cdf_inquire(nwout, vn_machsq, dimlens, ier=ierror)
      IF (ierror.eq.0) CALL cdf_read(nwout, vn_machsq, machsq)
      CALL cdf_inquire(nwout, vn_wpar, dimlens, ier=ierror)
      IF (ierror.eq.0) CALL cdf_read(nwout, vn_wpar, wp)  ! This overwrites wp with wpar

! Inquire existence, dimensions of arrays for allocation
! 1D Arrays
      IF (lfreeb .and. nbsets.gt.0) THEN
         CALL cdf_read(nwout, vn_nbfld, nbfld)
      END IF

      CALL cdf_inquire(nwout, vn_pmod, dimlens)
      ALLOCATE (xm(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_tmod, dimlens)
      ALLOCATE (xn(dimlens(1)), stat = ierror)
      IF (mnmax_nyq .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_pmod_nyq, dimlens)
         ALLOCATE (xm_nyq(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_tmod_nyq, dimlens)
         ALLOCATE (xn_nyq(dimlens(1)), stat = ierror)
      END IF

      IF (mnmaxpot .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_potsin, dimlens)
         ALLOCATE (potsin(dimlens(1)), xmpot(dimlens(1)),               &
                   xnpot(dimlens(1)), stat = ierror)
         IF (lasym) ALLOCATE (potcos(dimlens(1)), stat = ierror)
      END IF

      CALL cdf_inquire(nwout, vn_racc, dimlens)
      ALLOCATE (raxis_cc(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_zacs, dimlens)
      ALLOCATE (zaxis_cs(0:dimlens(1)-1), stat = ierror)
      IF (lasym) THEN
         CALL cdf_inquire(nwout, vn_racs, dimlens) 
         ALLOCATE (raxis_cs(0:dimlens(1)-1), stat = ierror)
         CALL cdf_inquire(nwout, vn_zacc, dimlens)
         ALLOCATE (zaxis_cc(0:dimlens(1)-1), stat = ierror)
      END IF

!  Profile coefficients, dimensioned from 0
      CALL cdf_inquire(nwout, vn_am, dimlens)
      ALLOCATE (am(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_ac, dimlens)
      ALLOCATE (ac(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai, dimlens)
      ALLOCATE (ai(0:dimlens(1)-1), stat = ierror)
      
      CALL cdf_inquire(nwout, vn_ac_aux_s, dimlens)
      ALLOCATE (ac_aux_s(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ac_aux_f, dimlens)
      ALLOCATE (ac_aux_f(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai_aux_s, dimlens)
      ALLOCATE (ai_aux_s(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai_aux_f, dimlens)
      ALLOCATE (ai_aux_f(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_am_aux_s, dimlens)
      ALLOCATE (am_aux_s(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_am_aux_f, dimlens)
      ALLOCATE (am_aux_f(dimlens(1)), stat = ierror)
      
      CALL cdf_inquire(nwout, vn_iotaf, dimlens)
      ALLOCATE (iotaf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_qfact, dimlens)

      ALLOCATE (qfact(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_presf, dimlens) 
      ALLOCATE (presf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phi, dimlens)
      ALLOCATE (phi(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_chi, dimlens)

      ALLOCATE (chi(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_phipf, dimlens)
      ALLOCATE (phipf(dimlens(1)), stat = ierror)
!OLD VERSION MAY NOT HAVE THIS!      CALL cdf_inquire(nwout, vn_chipf, dimlens)
      ALLOCATE (chipf(dimlens(1)), stat = ierror)
      chipf = 1.E30_dp

      CALL cdf_inquire(nwout, vn_jcuru, dimlens)
      ALLOCATE (jcuru(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_jcurv, dimlens)
      ALLOCATE (jcurv(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_iotah, dimlens)
      ALLOCATE (iotas(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mass, dimlens) 
      ALLOCATE (mass(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_presh, dimlens)
      ALLOCATE (pres(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_betah, dimlens)
      ALLOCATE (beta_vol(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_buco, dimlens) 
      ALLOCATE (buco(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bvco, dimlens)
      ALLOCATE (bvco(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_vp, dimlens)
      ALLOCATE (vp(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_specw, dimlens) 
      ALLOCATE (specw(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phip, dimlens)
      ALLOCATE (phip(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_overr, dimlens)
      ALLOCATE (overr(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_jdotb, dimlens)
      ALLOCATE (jdotb(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bdotb, dimlens)
      ALLOCATE (bdotb(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bgrv, dimlens)
      ALLOCATE (bdotgradv(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_merc, dimlens)
      ALLOCATE (Dmerc(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mshear, dimlens)
      ALLOCATE (Dshear(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mwell, dimlens)
      ALLOCATE (Dwell(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mcurr, dimlens)
      ALLOCATE (Dcurr(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mgeo, dimlens)
      ALLOCATE (Dgeod(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_equif, dimlens)
      ALLOCATE (equif(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_fsq, dimlens)
      ALLOCATE (fsqt(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_wdot, dimlens)
      ALLOCATE (wdot(dimlens(1)), stat = ierror)

      IF (nextcur .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_extcur, dimlens)
         ALLOCATE (extcur(dimlens(1)), stat = ierror)
!NOTE: curlabel is an array of CHARACTER(30) strings - defined in mgrid_mod 
!      so dimlens(1) == 30 (check this) and dimlens(2) is the number of strings in the array
         CALL cdf_inquire(nwout, vn_curlab, dimlens)
         ALLOCATE (curlabel(dimlens(2)), stat = ierror)
      ENDIF
      
      ! ANIMEC/FLOW -SAL
      CALL cdf_inquire(nwout, vn_pmap, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (pmap(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_omega, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (omega(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_tpotb, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (tpotb(dimlens(1)), stat = ierror)

! 2D Arrays
      CALL cdf_inquire(nwout, vn_rmnc, dimlens)
      ALLOCATE (rmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_zmns, dimlens)
      ALLOCATE (zmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_lmns, dimlens)
      ALLOCATE (lmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_gmnc, dimlens)
      ALLOCATE (gmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bmnc, dimlens)
      ALLOCATE (bmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubumnc, dimlens)
      ALLOCATE (bsubumnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubvmnc, dimlens)
      ALLOCATE (bsubvmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubsmns, dimlens)
      ALLOCATE (bsubsmns(dimlens(1),dimlens(2)), stat = ierror)

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
      CALL cdf_inquire(nwout, vn_bsupumnc, dimlens)
      ALLOCATE (bsupumnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsupvmnc, dimlens)
      ALLOCATE (bsupvmnc(dimlens(1),dimlens(2)), stat = ierror)

!     ANIMEC/FLOW -SAL
      CALL cdf_inquire(nwout, vn_pparmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (pparmnc(dimlens(1),dimlens(2)),        &
                                stat = ierror)
      CALL cdf_inquire(nwout, vn_ppermnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (ppermnc(dimlens(1),dimlens(2)),        &
                                stat = ierror)
      CALL cdf_inquire(nwout, vn_hotdmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (hotdmnc(dimlens(1),dimlens(2)),        &
                                stat = ierror)
      CALL cdf_inquire(nwout, vn_pbprmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (pbprmnc(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_ppprmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (ppprmnc(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_sigmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (sigmnc(dimlens(1),dimlens(2)),         &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_taumnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (taumnc(dimlens(1),dimlens(2)),         &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_protmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (protmnc(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_protrsqmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (protrsqmnc(dimlens(1),dimlens(2)),     &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_prprmnc, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (prprmnc(dimlens(1),dimlens(2)),        &
                                 stat = ierror)

!  The curr*mn* arrays have the same dimensions as the bsu**mn* arrays. No need
!  to inquire about the dimension sizes.
      ALLOCATE (currumnc(dimlens(1),dimlens(2)), stat = ierror)
      ALLOCATE (currvmnc(dimlens(1),dimlens(2)), stat = ierror)

      IF (lfreeb) THEN
         CALL cdf_inquire(nwout, vn_bsubumnc_sur, dimlens)
         ALLOCATE (bsubumnc_sur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_bsubvmnc_sur, dimlens)
         ALLOCATE (bsubvmnc_sur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_bsupumnc_sur, dimlens)
         ALLOCATE (bsupumnc_sur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_bsupvmnc_sur, dimlens)
         ALLOCATE (bsupvmnc_sur(dimlens(1)), stat = ierror)
      END IF

      IF (.NOT. lasym) GO TO 800

      CALL cdf_inquire(nwout, vn_rmns, dimlens)
      ALLOCATE (rmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_zmnc, dimlens)
      ALLOCATE (zmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_lmnc, dimlens)
      ALLOCATE (lmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_gmns, dimlens)
      ALLOCATE (gmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bmns, dimlens)
      ALLOCATE (bmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubumns, dimlens)
      ALLOCATE (bsubumns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubvmns, dimlens)
      ALLOCATE (bsubvmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubsmnc, dimlens)
      ALLOCATE (bsubsmnc(dimlens(1),dimlens(2)), stat = ierror)

!     ELIMINATE THESE EVENTUALLY: DO NOT NEED THEM
      CALL cdf_inquire(nwout, vn_bsupumns, dimlens)
      ALLOCATE (bsupumns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsupvmns, dimlens)
      ALLOCATE (bsupvmns(dimlens(1),dimlens(2)), stat = ierror)
      
      ! ANIMEC/FLOW -SAL
      CALL cdf_inquire(nwout, vn_pparmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (pparmns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_ppermns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (ppermns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_hotdmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (hotdmns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_pbprmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (pbprmns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_ppprmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (ppprmns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_sigmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (sigmns(dimlens(1),dimlens(2)),         &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_taumns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (taumns(dimlens(1),dimlens(2)),         &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_protmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (protmns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_protrsqmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (protrsqmns(dimlens(1),dimlens(2)),     &
                                 stat = ierror)
      CALL cdf_inquire(nwout, vn_prprmns, dimlens, ier=ierror)
      IF (ierror == 0) ALLOCATE (prprmns(dimlens(1),dimlens(2)),        &
                                 stat = ierror)

!  The curr*mn* arrays have the same dimensions as the bsu**mn* arrays. No need
!  to inquire about the dimension sizes.
      ALLOCATE (currumns(dimlens(1),dimlens(2)), stat = ierror)
      ALLOCATE (currvmns(dimlens(1),dimlens(2)), stat = ierror)

      IF (lfreeb) THEN
         CALL cdf_inquire(nwout, vn_bsubumns_sur, dimlens)
         ALLOCATE (bsubumns_sur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_bsubvmns_sur, dimlens)
         ALLOCATE (bsubvmns_sur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_bsupumns_sur, dimlens)
         ALLOCATE (bsupumns_sur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_bsupvmns_sur, dimlens)
         ALLOCATE (bsupvmns_sur(dimlens(1)), stat = ierror)
      END IF

 800  CONTINUE

! Read Arrays
      CALL cdf_read(nwout, vn_pmod, xm)
      CALL cdf_read(nwout, vn_tmod, xn)
      IF (mnmax_nyq .le. 0) THEN
         mnmax_nyq = mnmax
         ALLOCATE (xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq), stat=ierror)
         xm_nyq = xm;  xn_nyq = xn
      ELSE
         CALL cdf_read(nwout, vn_pmod_nyq, xm_nyq)
         CALL cdf_read(nwout, vn_tmod_nyq, xn_nyq)
      END IF

      IF (mnmaxpot .GT. 0) THEN
         CALL cdf_read(nwout, vn_potsin, potsin)
         CALL cdf_read(nwout, vn_xmpot, xmpot)
         CALL cdf_read(nwout, vn_xnpot, xnpot)           !includes nfp factor
         IF (lasym) CALL cdf_read(nwout, vn_potcos, potcos)
      END IF

      mnyq = INT(MAXVAL(xm_nyq));  nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp

      CALL cdf_read(nwout, vn_racc, raxis_cc)
      CALL cdf_read(nwout, vn_zacs, zaxis_cs) 

      IF (SIZE(raxis_cc) .ne. ntor+1)                                   &
         STOP 'WRONG SIZE(raxis_cc) in READ_WOUT_NC'
      ALLOCATE (raxis(0:ntor,2), zaxis(0:ntor,2), stat=ierror)
      raxis(:,1) = raxis_cc(0:ntor);   zaxis(:,1) = zaxis_cs(0:ntor)
      raxis(:,2) = 0;                  zaxis(:,2) = 0
      DEALLOCATE (raxis_cc, zaxis_cs, stat=ierror)

      CALL cdf_read(nwout, vn_rmnc, rmnc)
      CALL cdf_read(nwout, vn_zmns, zmns)
      CALL cdf_read(nwout, vn_lmns, lmns)
      CALL cdf_read(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_read(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_read(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_read(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_read(nwout, vn_bsubsmns, bsubsmns)      !Full mesh
!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
      CALL cdf_read(nwout, vn_bsupumnc, bsupumnc)
      CALL cdf_read(nwout, vn_bsupvmnc, bsupvmnc)

      IF (version_ .ge. 9.0) THEN
         CALL cdf_read(nwout, vn_currumnc, currumnc)
         CALL cdf_read(nwout, vn_currvmnc, currvmnc)
      END IF

      IF (lfreeb) THEN
         CALL cdf_read(nwout, vn_bsubumnc_sur, bsubumnc_sur)
         CALL cdf_read(nwout, vn_bsubvmnc_sur, bsubvmnc_sur)
         CALL cdf_read(nwout, vn_bsupumnc_sur, bsupumnc_sur)
         CALL cdf_read(nwout, vn_bsupvmnc_sur, bsupvmnc_sur)
      END IF

      IF (lasym) THEN
         CALL cdf_read(nwout, vn_racs, raxis_cs)
         CALL cdf_read(nwout, vn_zacc, zaxis_cc) 
         raxis(:,2) = raxis_cs;   zaxis(:,2) = zaxis_cc
         DEALLOCATE (raxis_cs, zaxis_cc, stat=ierror)
         CALL cdf_read(nwout, vn_rmns, rmns)
         CALL cdf_read(nwout, vn_zmnc, zmnc)
         CALL cdf_read(nwout, vn_lmnc, lmnc)
         CALL cdf_read(nwout, vn_gmns, gmns)
         CALL cdf_read(nwout, vn_bmns, bmns) 
         CALL cdf_read(nwout, vn_bsubumns, bsubumns)
         CALL cdf_read(nwout, vn_bsubvmns, bsubvmns)
         CALL cdf_read(nwout, vn_bsubsmnc, bsubsmnc)
!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_read(nwout, vn_bsupumns, bsupumns)
         CALL cdf_read(nwout, vn_bsupvmns, bsupvmns)

         IF (version_ .ge. 9.0) THEN
            CALL cdf_read(nwout, vn_currumns, currumns)
            CALL cdf_read(nwout, vn_currvmns, currvmns)
         END IF

         IF (lfreeb) THEN
            CALL cdf_read(nwout, vn_bsubumns_sur, bsubumns_sur)
            CALL cdf_read(nwout, vn_bsubvmns_sur, bsubvmns_sur)
            CALL cdf_read(nwout, vn_bsupumns_sur, bsupumns_sur)
            CALL cdf_read(nwout, vn_bsupvmns_sur, bsupvmns_sur)
         END IF
      END IF
      
      ! ANIMEC/FLOW -SAL
      IF (vmec_type == 1) THEN
         CALL cdf_read(nwout, vn_pparmnc, pparmnc)
         CALL cdf_read(nwout, vn_ppermnc, ppermnc)
         CALL cdf_read(nwout, vn_hotdmnc, hotdmnc)
         CALL cdf_read(nwout, vn_pbprmnc, pbprmnc)
         CALL cdf_read(nwout, vn_ppprmnc, ppprmnc)
         CALL cdf_read(nwout,  vn_sigmnc,  sigmnc)
         CALL cdf_read(nwout,  vn_taumnc,  taumnc)
         IF (lasym) THEN
            CALL cdf_read(nwout, vn_pparmns, pparmns)
            CALL cdf_read(nwout, vn_ppermns, ppermns)
            CALL cdf_read(nwout, vn_hotdmns, hotdmns)
            CALL cdf_read(nwout, vn_pbprmns, pbprmns)
            CALL cdf_read(nwout, vn_ppprmns, ppprmns)
            CALL cdf_read(nwout,  vn_sigmns,  sigmns)
            CALL cdf_read(nwout,  vn_taumns,  taumns)
         END IF
      ELSE IF (vmec_type == 2) THEN
         CALL cdf_read(nwout,    vn_protmnc,    protmnc)
         CALL cdf_read(nwout,    vn_prprmnc,    prprmnc)
         CALL cdf_read(nwout, vn_protrsqmnc, protrsqmnc)
         IF (lasym) THEN
            CALL cdf_read(nwout,    vn_protmns,    protmns)
            CALL cdf_read(nwout,    vn_prprmns,    prprmns)
            CALL cdf_read(nwout, vn_protrsqmns, protrsqmns)
         END IF
      END IF

      CALL cdf_read(nwout, vn_am, am)
      CALL cdf_read(nwout, vn_ac, ac)
      CALL cdf_read(nwout, vn_ai, ai)
      
      CALL cdf_read(nwout, vn_am_aux_s, am_aux_s)
      CALL cdf_read(nwout, vn_am_aux_f, am_aux_f)
      CALL cdf_read(nwout, vn_ac_aux_s, ac_aux_s)
      CALL cdf_read(nwout, vn_ac_aux_f, ac_aux_f)
      CALL cdf_read(nwout, vn_ai_aux_s, ai_aux_s)
      CALL cdf_read(nwout, vn_ai_aux_f, ai_aux_f)

      CALL cdf_read(nwout, vn_iotaf, iotaf) 
      CALL cdf_read(nwout, vn_qfact, qfact) 

      CALL cdf_read(nwout, vn_presf, presf) 
      CALL cdf_read(nwout, vn_phi, phi) 
      CALL cdf_read(nwout, vn_phipf, phipf)
      CALL cdf_read(nwout, vn_chi, chi) 

      CALL cdf_read(nwout, vn_chipf, chipf)
      IF (ALL(chipf .EQ. 1.E30_dp)) THEN
         chipf = iotaf*phipf
      END IF

      CALL cdf_read(nwout, vn_jcuru, jcuru)
      CALL cdf_read(nwout, vn_jcurv, jcurv)
 
!     HALF-MESH quantities
!     NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
      CALL cdf_read(nwout, vn_iotah, iotas)
      CALL cdf_read(nwout, vn_mass, mass) 
      CALL cdf_read(nwout, vn_presh, pres)
      CALL cdf_read(nwout, vn_betah, beta_vol)
      CALL cdf_read(nwout, vn_buco, buco)
      CALL cdf_read(nwout, vn_bvco, bvco) 
      CALL cdf_read(nwout, vn_vp, vp)
      CALL cdf_read(nwout, vn_specw, specw)
      CALL cdf_read(nwout, vn_phip, phip)
      CALL cdf_read(nwout, vn_jdotb, jdotb)
      CALL cdf_read(nwout, vn_bdotb, bdotb)
      CALL cdf_read(nwout, vn_bgrv, bdotgradv)

!     MERCIER_CRITERION
      CALL cdf_read(nwout, vn_merc, Dmerc)
      CALL cdf_read(nwout, vn_mshear, Dshear)
      CALL cdf_read(nwout, vn_mwell, Dwell)
      CALL cdf_read(nwout, vn_mcurr, Dcurr)
      CALL cdf_read(nwout, vn_mgeo, Dgeod)
      CALL cdf_read(nwout, vn_equif, equif)

      CALL cdf_read(nwout, vn_fsq, fsqt)
      CALL cdf_read(nwout, vn_wdot, wdot)
     
      IF (nextcur .gt. 0) THEN
         CALL cdf_read(nwout, vn_extcur, extcur)
         CALL cdf_read(nwout, vn_curlab, curlabel)
      ENDIF

 1000 CONTINUE

      CALL cdf_close(nwout, ierr)

      IF (.not.ALLOCATED(bsubumnc)) RETURN                              !Moved this here because ns may not be set. SAL -09/07/11
!
!     COMPUTE CONTRAVARIANT CURRENT COMPONENTS IN AMPS 
!     ON THE FULL RADIAL MESH, WHERE JACOBIAN = SQRT(G)
!
!     CURRU = SQRT(G) * J dot grad(u)
!     CURRV = SQRT(G) * J dot grad(v)
!
      IF (version_ .lt. 9.0) THEN
         IF (ierror .eq. 0) THEN
            CALL Compute_Currents(bsubsmnc, bsubsmns,                          &
     &                            bsubumnc, bsubumns,                          &
     &                            bsubvmnc, bsubvmns,                          &
     &                            xm_nyq, xn_nyq, mnmax_nyq, lasym, ns,        &
     &                            currumnc, currvmnc,                          &
     &                            currumns, currvmns)
         END IF
      END IF

      IF (ierr .ne. 0)   PRINT *,'in read_wout_nc ierr=',ierr
      IF (ierror .ne. 0) PRINT *,'in read_wout_nc ierror=',ierror

      END SUBROUTINE read_wout_nc
#endif

      SUBROUTINE write_wout_text(filename, ierr)
      USE vsvd0, ONLY: nparts
      USE safe_open_mod
      USE stel_constants, ONLY: mu0

      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      CHARACTER (len=*)    :: filename
      INTEGER, INTENT(out) :: ierr
!------------------------------------------------
!   L o c a l   P a r a m e t e r s
!------------------------------------------------
      REAL(rprec), PARAMETER :: eps_w = 1.e-4_dp
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER              :: iounit, js, mn, i, j, k, m, n, iasymm
      LOGICAL              :: lcurr
!------------------------------------------------
!
!     THIS SUBROUTINE WRITES A TEXT FILE WOUT CREATED BY STORED THE INFORMATION 
!     IN THE read_WOUT MODULE. This routine can only be called if the wout has 
!     already been read in.

      iounit = 0
      ierr = 0
      CALL safe_open(iounit, ierr,                                             &
     &               TRIM(filename),                      &
     &               'replace', 'formatted')

      CALL assert_eq(0, ierr, 'Error opening text wout file in ' //            &
     &               'write_wout_text of read_wout_mod.')


!  Write version info
      WRITE (iounit, '(a15,f4.2)') 'VMEC VERSION = ', version_

!  Check version numbers since values change.
      IF (lasym) THEN
         iasymm = 1
      ELSE
         iasym = 0
      END IF

      IF (version_ .le. (5.10 + eps_w)) THEN
         WRITE (iounit, *) wb, wp, gamma, pfac, nfp, ns, mpol, ntor,           &
     &      mnmax, itfsq, niter, iasymm, ireconstruct
      ELSE
         IF (version_ .lt. 6.54) THEN
            WRITE (iounit, *) wb, wp, gamma, pfac, rmax_surf, rmin_surf
         ELSE
            WRITE (iounit, *) wb, wp, gamma, pfac, rmax_surf, rmin_surf,       &
     &                        zmax_surf
         END IF
         IF (version_ .le. (8.0 + eps_w)) THEN
            WRITE (iounit, *) nfp, ns, mpol, ntor, mnmax, itfsq, niter,        &
     &                        iasym, ireconstruct, ierr_vmec
         ELSE
            WRITE (iounit, *) nfp, ns, mpol, ntor, mnmax, mnmax_nyq,           &
     &                        itfsq, niter, iasym, ireconstruct,               &
     &                        ierr_vmec
         END IF
      END IF

      IF (version_ .gt. (6.20 + eps_w)) THEN
         WRITE (iounit, *) imse, itse, nbsets, nobd, nextcur, nstore_seq
      ELSE
         WRITE (iounit, *) imse, itse, nbsets, nobd, nextcur
      END IF
      
      IF ((ierr_vmec .ne. norm_term_flag) .and.                                &
          (ierr_vmec .ne. more_iter_flag) .and.                                &
          (ierr_vmec .ne. jac75_flag)) THEN
        ierr = -2
        GOTO 1000
      ENDIF

      IF (nbsets .gt. 0) THEN
         WRITE (iounit, *) nbfld(1:nbsets)
      END IF
      WRITE (iounit, *) TRIM(mgrid_file)

      DO js = 1, ns
         DO mn = 1, mnmax
            IF (js .eq. 1) THEN
               WRITE (iounit, *) NINT(xm(mn)), NINT(xn(mn)/nfp)
            END IF
            IF (version_ .le. (6.20 + eps_w)) THEN
               WRITE (iounit, 730) rmnc(mn,js), zmns(mn,js),                   &
     &                             lmns(mn,js), bmnc(mn,js),                   &
     &                             gmnc(mn,js), bsubumnc(mn,js),               &
     &                             bsubvmnc(mn,js), bsubsmns(mn,js),           &
     &                             bsupumnc(mn,js), bsupvmnc(mn,js),           &
     &                             currvmnc(mn,js)
            ELSE IF (version_ .le. (8.0 + eps_w)) THEN
               WRITE (iounit, *) rmnc(mn,js), zmns(mn,js), lmns(mn,js),        &
     &                           bmnc(mn,js), gmnc(mn,js),                     &
     &                           bsubumnc(mn,js), bsubvmnc(mn,js),             &
     &                           bsubsmns(mn,js), bsupumnc(mn,js),             &
     &                           bsupvmnc(mn,js), currvmnc(mn,js)
            ELSE
               WRITE (iounit, *) rmnc(mn,js), zmns(mn,js), lmns(mn,js)
            END IF

!  Write asymmetric components.
            IF (lasym) THEN
               IF (version_ .le. (8.0 + eps_w)) THEN
                  WRITE (iounit, *) rmns(mn,js), zmnc(mn,js),                  &
     &                              lmnc(mn,js), bmns(mn,js),                  &
     &                              gmns(mn,js), bsubumns(mn,js),              &
     &                              bsubvmns(mn,js), bsubsmnc(mn,js),          &
     &                              bsupumns(mn,js), bsubvmns(mn,js)
               ELSE
                  WRITE (iounit, *) rmns(mn,js), zmnc(mn,js),                  &
     &                              lmnc(mn,js)
               END IF
            END IF
         END DO

         IF (version_ .le. (8.0 + eps_w)) THEN
            CYCLE
         END IF

         DO mn = 1, mnmax_nyq
            IF (js .eq. 1) THEN
               WRITE (iounit, *) NINT(xm_nyq(mn)),                             &
     &                           NINT(xn_nyq(mn)/nfp)
            END IF
            WRITE (iounit, *) bmnc(mn,js), gmnc(mn,js),                        &
     &                        bsubumnc(mn,js), bsubvmnc(mn,js),                &
     &                        bsubsmns(mn,js), bsupumnc(mn,js),                &
     &                        bsupvmnc(mn,js)
            IF (lasym) THEN
               WRITE (iounit, *) bmns(mn,js), gmns(mn,js),                     &
     &                           bsubumns(mn,js), bsubvmns(mn,js),             &
     &                           bsubsmnc(mn,js), bsupumns(mn,js),             &
     &                           bsupvmns(mn,js)
            END IF
         END DO
      END DO

!
!     Write FULL AND HALF-MESH QUANTITIES
!
!     NOTE: In version_ <= 6.00, mass, press were written out in INTERNAL (VMEC) units
!     and are therefore multiplied here by 1/mu0 to transform to pascals. Same is true
!     for ALL the currents (jcuru, jcurv, jdotb). Also, in version_ = 6.10 and
!     above, PHI is the true (physical) toroidal flux (has the sign of jacobian correctly
!     built into it)
!

      IF (version_ .le. (6.05 + eps_w)) THEN
         WRITE (iounit, 730) (iotas(js), mass(js)*mu0, pres(js)*mu0,           &
     &                        phip(js), buco(js), bvco(js), -phi(js),          &
     &                        vp(js), overr(js), jcuru(js)*mu0,                &
     &                        jcurv(js)*mu0, specw(js), js=2, ns)
         WRITE (iounit, 730) aspect, betatot, betapol, betaxis, b0
      ELSE IF (version_ .le. (6.20 + eps_w)) THEN
         WRITE (iounit, 730) (iotas(js), mass(js), pres(js),                   &
     &                        beta_vol(js), phip(js), buco(js),                &
     &                        bvco(js), phi(js), vp(js), overr(js),            &
     &                        jcuru(js), jcurv(js), specw(js),                 &
     &                        js=2, ns)
         WRITE (iounit, 730) aspect, betatot, betapol, betaxis, b0
      ELSE IF (version_ .le. (6.95 + eps_w)) THEN
         WRITE (iounit, *) (iotas(js), mass(js), pres(js),                     &
     &                      beta_vol(js), phip(js), buco(js),                  &
     &                      bvco(js), phi(js), vp(js), overr(js),              &
     &                      jcuru(js), jcurv(js), specw(js),                   &
     &                      js=2, ns)
         WRITE (iounit, *) aspect, betatot, betapol, betaxis, b0
      ELSE
         WRITE (iounit, *) (iotaf(js), presf(js), phipf(js), phi(js),          &
     &                       jcuru(js), jcurv(js), js=1, ns)
         WRITE (iounit, *) (iotas(js), mass(js), pres(js),                     &
     &                      beta_vol(js), phip(js), buco(js),                  &
     &                      bvco(js), vp(js), overr(js), specw(js),            &
     &                      js = 2, ns)
         WRITE (iounit, *) aspect, betatot, betapol, betaxis, b0
      END IF

      IF (version_ .gt. (6.10 + eps_w)) THEN
         WRITE (iounit, *) isigng
         WRITE (iounit, *) TRIM(input_extension)
         WRITE (iounit, *) IonLarmor, VolAvgB, RBtor0, RBtor, Itor,            &
     &                     Aminor, Rmajor, Volume
      END IF

!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      IF (version_ .gt. (5.10 + eps_w) .and.                                   &
     &    version_ .lt. (6.20 - eps_w)) THEN
         WRITE (iounit, 730) (Dmerc(js), Dshear(js), Dwell(js),                &
     &                        Dcurr(js), Dgeod(js), equif(js),                 &
     &                        js=2, ns - 1)
      ELSE IF (version_ .ge. (6.20 - eps_w)) THEN
         WRITE (iounit, *) (Dmerc(js), Dshear(js), Dwell(js),                  &
     &                      Dcurr(js), Dgeod(js), equif(js),                   &
     &                      js=2, ns - 1)
      END IF

      IF (nextcur .gt. 0) THEN
         IF (version_ .le. (6.20 + eps_w)) THEN
            WRITE (iounit, 730) (extcur(js), js=1, nextcur)
         ELSE
            WRITE (iounit, *) (extcur(js), js=1, nextcur)
         END IF

         lcurr = LEN_TRIM(curlabel(1)) .gt. 0
         WRITE (iounit, *) lcurr
         IF (lcurr) THEN
            WRITE (iounit, *) (TRIM(curlabel(js)), js=1, nextcur)
         END IF
      END IF

      IF (version_ .le. (6.20 + eps_w)) THEN
         WRITE (iounit, 730) (fsqt(js), wdot(js), js = 1, nstore_seq)
      ELSE
         WRITE (iounit, *) (fsqt(js), wdot(js), js = 1, nstore_seq)
      END IF

      IF (version_ .ge. (6.20 - eps_w) .and.                                   &
     &    version_ .lt. (6.50 - eps_w)) THEN
         WRITE (iounit, 730) (jdotb(js), bdotgradv(js), js=1, ns)
      ELSE IF (version_ .ge. (6.50 - eps_w)) THEN
         WRITE (iounit, *) (jdotb(js), bdotgradv(js), js=1, ns)
      END IF

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      IF (ireconstruct .gt. 0) THEN
         IF (imse .ge. 2 .or. itse .gt. 0) THEN
            WRITE (iounit, *) tswgt, msewgt
            WRITE (iounit, *) isnodes, (sknots(js), ystark(js),                &
     &                                  y2stark(js), js=1, isnodes)
            WRITE (iounit, *) ipnodes, (pknots(js), ythom(js),                 &
     &                                  y2thom(js), js=1, ipnodes)
            WRITE (iounit, *) (anglemse(js), rmid(js), qmid(js),               &
     &                         shear(js), presmid(js), alfa(js),               &
     &                         curmid(js), js=1, 2*ns - 1)
            WRITE (iounit, *) (rstark(js), datastark(js), qmeas(js),           &
     &                         js=1, imse)
            WRITE (iounit, *) (rthom(js), datathom(i), js=1, itse)
         END IF

         IF (nobd .gt. 0) THEN
            WRITE (iounit, *) (dsiext(js), plflux(js), dsiobt(js),             &
     &                         js=1, nobd)
            WRITE (iounit, *) flmwgt
         END IF

         IF (nbfldn .gt. 0) THEN
            DO n = 1, nbsets
               READ (iounit, *) (bcoil(i,n), plbfld(i,n), bbc(i,n),            &
     &                           i=1,nbfld(n))
            END DO
            WRITE (iounit, *) bcwgt
         END IF

         WRITE (iounit, *) phidiam, delphid
!
!     READ Limiter & Prout plotting specs
!
         WRITE (iounit, *) nsets, nparts, nlim

         WRITE (iounit, *) (nsetsn(js), js=1, nsets)

         WRITE (iounit, *) (((pfcspec(i,j,k), i=1, nparts),                    &
     &                       j=1, nsetsn(k)), k=1, nsets)

         WRITE (iounit, *) (limitr(i), i=1, nlim)

         WRITE (iounit, *) ((rlim(i,j), zlim(i,j), i=1, limitr(j)),            &
     &                      j=1, nlim)
         WRITE (iounit, *) nrgrid, nzgrid
         WRITE (iounit, *) tokid
         WRITE (iounit, *) rx1, rx2, zy1, zy2, condif
         WRITE (iounit, *) imatch_phiedge
      END IF

1000  CONTINUE

      WRITE (iounit, *) mgrid_mode

  730 FORMAT(5e20.13)

      CLOSE (iounit, iostat = ierr)
      CALL assert_eq(0, ierr, 'Error closing text wout file in ' //            &
     &               'write_wout_text of read_wout_mod.')

      END SUBROUTINE

#if defined(NETCDF)
      SUBROUTINE write_wout_nc(filename, ierr)
      USE vsvd0, ONLY: nparts
      USE stel_constants, ONLY: mu0, zero, one
      USE ezcdf

      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      CHARACTER (len=*)    :: filename
      INTEGER, INTENT(out) :: ierr
!------------------------------------------------
!   L o c a l   P a r a m e t e r s
!------------------------------------------------
      REAL(rprec), PARAMETER :: eps_w = 1.e-4_dp
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER              :: nwout, js, mn, i, j, k, m, n, iasymm, iwout0
      LOGICAL              :: lcurr
#ifdef NETCDF
      CHARACTER(LEN=*), PARAMETER, DIMENSION(1) ::  &
                   r1dim = (/'radius'/), mn1dim = (/'mn_mode'/), &
                   mn2dim = (/'mn_mode_nyq'/), &
                   mnpotdim = (/'mn_mode_pot'/), &
                   currg = (/'ext_current'/), &
                   currl = (/'current_label'/)
      CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: &
                   r2dim = (/'mn_mode','radius '/), &
                   r3dim = (/'mn_mode_nyq','radius     '/)
#endif
!------------------------------------------------
!
!     THIS SUBROUTINE WRITES A netCDF FILE WOUT CREATED BY STORED THE INFORMATION 
!     IN THE read_WOUT MODULE. This routine can only be called if the wout has 
!     already been read in.
      ierr   = 0
      iwout0 = 0     
      CALL cdf_open(nwout,TRIM(filename),'w',iwout0)
      IF (iwout0 .ne. 0) STOP 'Error opening wout.nc file VMEC WROUT'


      nstore_seq = 100
!      IF (.not. lrecon) THEN
!         itse = 0
!         imse2 = 0
!      END IF
      
!================================
! Define Variables
!================================
!  Scalars 
      CALL cdf_define(nwout, vn_version, version_)
      CALL cdf_define(nwout, vn_extension, input_extension)
      CALL cdf_define(nwout, vn_mgrid, mgrid_file)
      CALL cdf_define(nwout, vn_pcurr_type, pcurr_type)
      CALL cdf_define(nwout, vn_pmass_type, pmass_type)
      CALL cdf_define(nwout, vn_piota_type, piota_type)
      CALL cdf_define(nwout, vn_magen, wb)
      CALL cdf_define(nwout, vn_therm, wp)
      CALL cdf_define(nwout, vn_gam, gamma)
      CALL cdf_define(nwout, vn_maxr, rmax_surf)
      CALL cdf_define(nwout, vn_minr, rmin_surf)
      CALL cdf_define(nwout, vn_maxz, zmax_surf)
      CALL cdf_define(nwout, vn_fp, nfp)
      CALL cdf_define(nwout, vn_radnod, ns)
      CALL cdf_define(nwout, vn_polmod, mpol)
      CALL cdf_define(nwout, vn_tormod, ntor)
      CALL cdf_define(nwout, vn_maxmod, mnmax)
      CALL cdf_define(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_define(nwout, vn_maxit, niter)
      CALL cdf_define(nwout, vn_actit, itfsq)
      CALL cdf_define(nwout, vn_asym, lasym)
      CALL cdf_define(nwout, vn_recon, lrecon)
      CALL cdf_define(nwout, vn_free, lfreeb)
      CALL cdf_define(nwout, vn_rfp, lrfp)
      CALL cdf_define(nwout, vn_error, ierr_vmec)
      CALL cdf_define(nwout, vn_aspect, aspect)
      CALL cdf_define(nwout, vn_beta, betatot)
      CALL cdf_define(nwout, vn_pbeta, betapol)
      CALL cdf_define(nwout, vn_tbeta, betator)
      CALL cdf_define(nwout, vn_abeta, betaxis)
      CALL cdf_define(nwout, vn_b0, b0)
      CALL cdf_define(nwout, vn_rbt0, rbtor0)
      CALL cdf_define(nwout, vn_rbt1, rbtor)
      CALL cdf_define(nwout, vn_sgs, isigng)
      CALL cdf_define(nwout, vn_lar, IonLarmor)
      CALL cdf_define(nwout, vn_modB, volAvgB)
      CALL cdf_define(nwout, vn_ctor, Itor)
      CALL cdf_define(nwout, vn_amin, Aminor)
      CALL cdf_define(nwout, vn_Rmaj, Rmajor)
      CALL cdf_define(nwout, vn_vol, Volume)
      CALL cdf_define(nwout, vn_ftolv, ftolv)
      CALL cdf_define(nwout, vn_fsql, fsql)
      CALL cdf_define(nwout, vn_fsqr, fsqr)
      CALL cdf_define(nwout, vn_fsqz, fsqz)

!      IF (lrecon) THEN
!         CALL cdf_define(nwout, vn_mse, imse2)
!         CALL cdf_define(nwout, vn_thom, itse)
!      END IF

      CALL cdf_define(nwout, vn_nextcur, nextcur)
      CALL cdf_define(nwout, vn_extcur, extcur(1:nextcur),dimname=currg)
      CALL cdf_define(nwout, vn_mgmode, mgrid_mode)
!      IF (lfreeb) THEN
!         CALL cdf_define(nwout, vn_maxpot, mnpd)
!         CALL cdf_define(nwout, vn_flp, nobser)
!         CALL cdf_define(nwout, vn_nobd, nobd)
!         CALL cdf_define(nwout, vn_nbset, nbsets)
!         IF (nbsets .gt. 0) CALL cdf_define(nwout,vn_nbfld,nbfld(1:nbsets))
!      END IF
      CALL cdf_define(nwout, vn_pmod, xm, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_pmod, ln_pmod)
      CALL cdf_define(nwout, vn_tmod, xn, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_tmod, ln_tmod)
      CALL cdf_define(nwout, vn_pmod_nyq, xm_nyq, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_pmod_nyq, ln_pmod_nyq)
      CALL cdf_define(nwout, vn_tmod_nyq, xn_nyq, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_tmod_nyq, ln_tmod_nyq)

      CALL cdf_define(nwout, vn_racc, raxis(0:ntor,1),dimname=(/'n_tor'/))
      CALL cdf_setatt(nwout, vn_racc, ln_racc)
      CALL cdf_define(nwout, vn_zacs, zaxis(0:ntor,1),dimname=(/'n_tor'/))
      CALL cdf_setatt(nwout, vn_zacs, ln_zacs)
      IF (lasym) THEN
         CALL cdf_define(nwout, vn_racs, raxis(0:ntor,2),dimname=(/'n_tor'/))
         CALL cdf_setatt(nwout, vn_racs, ln_racs)
         CALL cdf_define(nwout, vn_zacc, zaxis(0:ntor,2),dimname=(/'n_tor'/))
         CALL cdf_setatt(nwout, vn_zacc, ln_zacc)
      END IF

      j = SIZE(am)-1
      CALL cdf_define(nwout, vn_am, am(0:j),dimname=(/'preset'/))
      j = SIZE(ac)-1
      CALL cdf_define(nwout, vn_ac, ac(0:j),dimname=(/'preset'/))
      j = SIZE(ai)-1
      CALL cdf_define(nwout, vn_ai, ai(0:j),dimname=(/'preset'/))
     
      j = SIZE(am_aux_s)
      CALL cdf_define(nwout, vn_am_aux_s, am_aux_s(1:j),dimname=(/'ndfmax'/))
      j = SIZE(am_aux_f)
      CALL cdf_define(nwout, vn_am_aux_f, am_aux_f(1:j),dimname=(/'ndfmax'/))
      j = SIZE(ai_aux_s)
      CALL cdf_define(nwout, vn_ai_aux_s, ai_aux_s(1:j),dimname=(/'ndfmax'/))
      j = SIZE(ai_aux_f)
      CALL cdf_define(nwout, vn_ai_aux_f, ai_aux_f(1:j),dimname=(/'ndfmax'/))
      j = SIZE(ac_aux_s)
      CALL cdf_define(nwout, vn_ac_aux_s, ac_aux_s(1:j),dimname=(/'ndfmax'/))
      j = SIZE(ac_aux_f)
      CALL cdf_define(nwout, vn_ac_aux_f, ac_aux_f(1:j),dimname=(/'ndfmax'/))


      CALL cdf_define(nwout, vn_iotaf, iotaf(1:ns),dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_iotaf)

      qfact=HUGE(qfact)
      WHERE (iotaf(1:ns) .NE. zero) qfact=one/iotaf(1:ns)

      CALL cdf_define(nwout, vn_qfact, qfact(1:ns),dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_qfact)
      CALL cdf_define(nwout, vn_presf, presf,dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presf, ln_presf, units='Pa')
      CALL cdf_define(nwout, vn_phi, phi,dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phi, ln_phi, units='wb')
      CALL cdf_define(nwout, vn_phipf,phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phipf, ln_phipf)
      CALL cdf_define(nwout, vn_chi, chi,dimname=r1dim)
      CALL cdf_setatt(nwout, vn_chi, ln_chi, units='wb')
      CALL cdf_define(nwout, vn_chipf,phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_chipf, ln_chipf)
      CALL cdf_define(nwout, vn_jcuru,jcuru, dimname=r1dim)
      CALL cdf_define(nwout, vn_jcurv,jcurv, dimname=r1dim)
 
      CALL cdf_define(nwout, vn_iotah, iotas(1:ns),dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotah, ln_iotah)
      CALL cdf_define(nwout, vn_mass, mass,dimname=r1dim)
      CALL cdf_setatt(nwout, vn_mass, ln_mass)
      CALL cdf_define(nwout, vn_presh, pres(1:ns),dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presh, ln_presh, units='Pa')
      CALL cdf_define(nwout, vn_betah, beta_vol,dimname=r1dim)
      CALL cdf_define(nwout, vn_buco, buco,dimname=r1dim)
      CALL cdf_define(nwout, vn_bvco, bvco,dimname=r1dim)
      CALL cdf_define(nwout, vn_vp, vp(1:ns),dimname=r1dim)
      CALL cdf_define(nwout, vn_specw, specw,dimname=r1dim)
      CALL cdf_define(nwout, vn_phip,phip(1:ns), dimname=r1dim)
      CALL cdf_define(nwout, vn_overr,overr(1:ns), dimname=r1dim)

      CALL cdf_define(nwout, vn_jdotb, jdotb,dimname=r1dim)
      CALL cdf_define(nwout, vn_bdotb, bdotb,dimname=r1dim)
      CALL cdf_define(nwout, vn_bgrv, bdotgradv,dimname=r1dim)

      CALL cdf_define(nwout, vn_merc, Dmerc,dimname=r1dim)
      CALL cdf_define(nwout, vn_mshear, Dshear,dimname=r1dim)
      CALL cdf_define(nwout, vn_mwell, Dwell,dimname=r1dim)
      CALL cdf_define(nwout, vn_mcurr, Dcurr,dimname=r1dim)
      CALL cdf_define(nwout, vn_mgeo,Dgeod, dimname=r1dim)
      CALL cdf_define(nwout, vn_equif,equif, dimname=r1dim)

      CALL cdf_define(nwout, vn_fsq, fsqt(1:nstore_seq),dimname=(/'time'/))
      CALL cdf_define(nwout, vn_wdot, wdot(1:nstore_seq),dimname=(/'time'/))

!      IF (lfreeb) THEN
!         CALL cdf_define(nwout, vn_potsin, potvac(1:mnpd),dimname=mnpotdim)
!         CALL cdf_setatt(nwout, vn_potsin, ln_potsin)
!         CALL cdf_define(nwout, vn_xmpot, xmpot(1:mnpd),dimname=mnpotdim)
!         CALL cdf_define(nwout, vn_xnpot, xnpot(1:mnpd),dimname=mnpotdim)
!         IF (lasym) THEN 
!            CALL cdf_define(nwout, vn_potcos,potvac(1+mnpd:2*mnpd), dimname=mnpotdim)
!            CALL cdf_setatt(nwout, vn_potcos, ln_potcos)
!         END IF
!
!         IF (nextcur.gt.0 .and. ALLOCATED(curlabel)) CALL cdf_define(nwout, vn_curlab,curlabel(1:nextcur), dimname=currl)
!      ENDIF
! 2D Arrays
      CALL cdf_define(nwout, vn_rmnc, rmnc, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_rmnc, ln_rmnc, units='m')
      CALL cdf_define(nwout, vn_zmns, zmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_zmns, ln_zmns, units='m')
      CALL cdf_define(nwout, vn_lmns, lmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_lmns, ln_lmns)
      CALL cdf_define(nwout, vn_gmnc, gmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_gmnc, ln_gmnc)
      CALL cdf_define(nwout, vn_bmnc, bmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bmnc, ln_bmnc)
      CALL cdf_define(nwout, vn_bsubumnc, bsubumnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubumnc, ln_bsubumnc)
      CALL cdf_define(nwout, vn_bsubvmnc, bsubvmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubvmnc, ln_bsubvmnc)
      CALL cdf_define(nwout, vn_bsubsmns, bsubsmns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubsmns, ln_bsubsmns)

      CALL cdf_define(nwout, vn_currumnc, currumnc, dimname=r3dim)    !MRC 8-12-16
      CALL cdf_setatt(nwout, vn_currumnc, ln_currumnc)
      CALL cdf_define(nwout, vn_currvmnc, currvmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_currvmnc, ln_currvmnc)

      IF (lfreeb) THEN
         CALL cdf_define(nwout, vn_bsubumnc_sur, bsubumnc_sur,dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsubumnc_sur, ln_bsubumnc_sur)
         CALL cdf_define(nwout, vn_bsubvmnc_sur, bsubvmnc_sur,dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsubvmnc_sur, ln_bsubvmnc_sur)
         CALL cdf_define(nwout, vn_bsupumnc_sur, bsupumnc_sur,dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsupumnc_sur, ln_bsupumnc_sur)
         CALL cdf_define(nwout, vn_bsupvmnc_sur, bsupvmnc_sur,dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsupvmnc_sur, ln_bsupvmnc_sur)
      END IF

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM - CAN COMPUTE FROM GSQRT
      CALL cdf_define(nwout, vn_bsupumnc, bsupumnc, dimname=r3dim)
      CALL cdf_define(nwout, vn_bsupvmnc, bsupvmnc, dimname=r3dim)

      IF (lasym) THEN
         CALL cdf_define(nwout, vn_rmns, rmns, dimname=r2dim)
         CALL cdf_setatt(nwout, vn_rmns, ln_rmns, units='m')
         CALL cdf_define(nwout, vn_zmnc, zmnc, dimname=r2dim)
         CALL cdf_setatt(nwout, vn_zmnc, ln_zmnc, units='m')
         CALL cdf_define(nwout, vn_lmnc, lmnc, dimname=r2dim)
         CALL cdf_setatt(nwout, vn_lmnc, ln_lmnc)
         CALL cdf_define(nwout, vn_gmns, gmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_gmns, ln_gmns)
         CALL cdf_define(nwout, vn_bmns, bmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bmns, ln_bmns)
         CALL cdf_define(nwout, vn_bsubumns, bsubumns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bsubumns, ln_bsubumns)
         CALL cdf_define(nwout, vn_bsubvmns, bsubvmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bsubvmns, ln_bsubvmns)
         CALL cdf_define(nwout, vn_bsubsmnc, bsubsmnc, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bsubsmnc, ln_bsubsmnc)

         CALL cdf_define(nwout, vn_currumns, currumns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_currumns, ln_currumns)
         CALL cdf_define(nwout, vn_currvmns, currvmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_currvmns, ln_currvmns)

         IF (lfreeb) THEN
            CALL cdf_define(nwout, vn_bsubumns_sur, bsubumns_sur,dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsubumns_sur, ln_bsubumns_sur)
            CALL cdf_define(nwout, vn_bsubvmns_sur, bsubvmns_sur,dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsubvmns_sur, ln_bsubvmns_sur)
            CALL cdf_define(nwout, vn_bsupumns_sur, bsupumns_sur,dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsupumns_sur, ln_bsupumns_sur)
            CALL cdf_define(nwout, vn_bsupvmns_sur, bsupvmns_sur,dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsupvmns_sur, ln_bsupvmns_sur)
         END IF

         CALL cdf_define(nwout, vn_bsupumns, bsupumns, dimname=r3dim)
         CALL cdf_define(nwout, vn_bsupvmns, bsupvmns, dimname=r3dim)
      END IF

! Write Vars
! Scalars
      CALL cdf_write(nwout, vn_version, version_)
      CALL cdf_write(nwout, vn_extension, input_extension)
      CALL cdf_write(nwout, vn_mgrid, mgrid_file)
      CALL cdf_write(nwout, vn_pcurr_type, pcurr_type)
      CALL cdf_write(nwout, vn_piota_type, piota_type)
      CALL cdf_write(nwout, vn_pmass_type, pmass_type)
      CALL cdf_write(nwout, vn_magen, wb)
      CALL cdf_write(nwout, vn_therm, wp)
      CALL cdf_write(nwout, vn_gam, gamma)
      CALL cdf_write(nwout, vn_maxr, rmax_surf)
      CALL cdf_write(nwout, vn_minr, rmin_surf)
      CALL cdf_write(nwout, vn_maxz, zmax_surf)
      CALL cdf_write(nwout, vn_fp, nfp)
      CALL cdf_write(nwout, vn_radnod, ns)
      CALL cdf_write(nwout, vn_polmod, mpol)
      CALL cdf_write(nwout, vn_tormod, ntor)
      CALL cdf_write(nwout, vn_maxmod, mnmax)
      CALL cdf_write(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_write(nwout, vn_maxit, niter)
      CALL cdf_write(nwout, vn_actit, itfsq)
      CALL cdf_write(nwout, vn_asym, lasym)
      CALL cdf_write(nwout, vn_recon, lrecon)
      CALL cdf_write(nwout, vn_free, lfreeb)
      CALL cdf_write(nwout, vn_rfp, lrfp)
      CALL cdf_write(nwout, vn_error, ierr_vmec)
!
      CALL cdf_write(nwout, vn_aspect, aspect)
      CALL cdf_write(nwout, vn_beta, betatot)
      CALL cdf_write(nwout, vn_pbeta, betapol)
      CALL cdf_write(nwout, vn_tbeta, betator)
      CALL cdf_write(nwout, vn_abeta, betaxis)
      CALL cdf_write(nwout, vn_b0, b0)
      CALL cdf_write(nwout, vn_rbt0, rbtor0)
      CALL cdf_write(nwout, vn_rbt1, rbtor)
      CALL cdf_write(nwout, vn_sgs, isigng)
      CALL cdf_write(nwout, vn_lar, IonLarmor)
      CALL cdf_write(nwout, vn_modB, volAvgB)
      CALL cdf_write(nwout, vn_ctor, Itor)
      CALL cdf_write(nwout, vn_amin, Aminor)
      CALL cdf_write(nwout, vn_rmaj, Rmajor)
      CALL cdf_write(nwout, vn_vol, Volume)
      CALL cdf_write(nwout, vn_ftolv, ftolv)
      CALL cdf_write(nwout, vn_fsql, fsql)
      CALL cdf_write(nwout, vn_fsqr, fsqr)
      CALL cdf_write(nwout, vn_fsqz, fsqz)

!      IF (lrecon) THEN
!         CALL cdf_write(nwout, vn_mse, imse2-1)
!         CALL cdf_write(nwout, vn_thom, itse)
!      END IF

      CALL cdf_write(nwout, vn_nextcur, nextcur)
      IF (nextcur .gt. 0) THEN
         CALL cdf_write(nwout, vn_extcur, extcur(1:nextcur))
         CALL cdf_write(nwout, vn_mgmode, mgrid_mode)
      ENDIF
!      IF (lfreeb) THEN
!         CALL cdf_write(nwout, vn_flp, nobser)
!         CALL cdf_write(nwout, vn_maxpot, mnpd)
!         CALL cdf_write(nwout, vn_nobd, nobd)
!         CALL cdf_write(nwout, vn_nbset, nbsets)
!         IF (nextcur.gt.0 .and. ALLOCATED(curlabel)) CALL cdf_write(nwout, vn_curlab, curlabel(1:nextcur))
!      END IF
! 1D Arrays
      IF (nbsets .gt. 0) CALL cdf_write(nwout,vn_nbfld,nbfld(1:nbsets))

      CALL cdf_write(nwout, vn_pmod, xm)
      CALL cdf_write(nwout, vn_tmod, xn)
      CALL cdf_write(nwout, vn_pmod_nyq, xm_nyq)
      CALL cdf_write(nwout, vn_tmod_nyq, xn_nyq)

!      IF (lfreeb) THEN
!         CALL cdf_write(nwout, vn_potsin, potvac(1:mnpd))
!         IF (lasym) CALL cdf_write(nwout, vn_potcos, potvac(1+mnpd:2*mnpd))
!         CALL cdf_write(nwout, vn_xmpot, xmpot)
!         CALL cdf_write(nwout, vn_xnpot, xnpot)
!      END IF

      CALL cdf_write(nwout, vn_racc, raxis(0:ntor,1))
      CALL cdf_write(nwout, vn_zacs, zaxis(0:ntor,1)) 
      CALL cdf_write(nwout, vn_rmnc, rmnc)
      CALL cdf_write(nwout, vn_zmns, zmns)
      CALL cdf_write(nwout, vn_lmns, lmns)
      CALL cdf_write(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubsmns, bsubsmns)      !Full mesh

      CALL cdf_write(nwout, vn_currumnc, currumnc)      !MRK 8-12-16
      CALL cdf_write(nwout, vn_currvmnc, currvmnc)

!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
      CALL cdf_write(nwout, vn_bsupumnc, bsupumnc)
      CALL cdf_write(nwout, vn_bsupvmnc, bsupvmnc)

      IF (lfreeb) THEN        !MRC    10-15-15
         CALL cdf_write(nwout, vn_bsubumnc_sur, bsubumnc_sur)
         CALL cdf_write(nwout, vn_bsubvmnc_sur, bsubvmnc_sur)
         CALL cdf_write(nwout, vn_bsupumnc_sur, bsupumnc_sur)
         CALL cdf_write(nwout, vn_bsupvmnc_sur, bsupvmnc_sur)
      END IF

!     FULL-MESH quantities
!     NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...

      j = SIZE(am)-1
      CALL cdf_write(nwout, vn_am, am(0:j))
      j = SIZE(ac)-1
      CALL cdf_write(nwout, vn_ac, ac(0:j))
      j = SIZE(ai)-1
      CALL cdf_write(nwout, vn_ai, ai(0:j))

      j = SIZE(am_aux_s)
      CALL cdf_write(nwout, vn_am_aux_s, am_aux_s(1:j))
      j = SIZE(am_aux_f)
      CALL cdf_write(nwout, vn_am_aux_f, am_aux_f(1:j))
      j = SIZE(ac_aux_s)
      CALL cdf_write(nwout, vn_ac_aux_s, ac_aux_s(1:j))
      j = SIZE(ac_aux_f)
      CALL cdf_write(nwout, vn_ac_aux_f, ac_aux_f(1:j))
      j = SIZE(ai_aux_s)
      CALL cdf_write(nwout, vn_ai_aux_s, ai_aux_s(1:j))
      j = SIZE(ai_aux_f)
      CALL cdf_write(nwout, vn_ai_aux_f, ai_aux_f(1:j))

      CALL cdf_write(nwout, vn_iotaf, iotaf(1:ns))
      CALL cdf_write(nwout, vn_qfact, qfact(1:ns))
      CALL cdf_write(nwout, vn_presf, presf) 
      CALL cdf_write(nwout, vn_phi, phi) 
      CALL cdf_write(nwout, vn_phipf, phipf)
      CALL cdf_write(nwout, vn_chi, chi) 
      CALL cdf_write(nwout, vn_chipf, chipf)
      CALL cdf_write(nwout, vn_jcuru, jcuru)
      CALL cdf_write(nwout, vn_jcurv, jcurv)
      CALL cdf_write(nwout, vn_jdotb, jdotb)
      CALL cdf_write(nwout, vn_bdotb, bdotb)
      CALL cdf_write(nwout, vn_bgrv, bdotgradv)
 
!     HALF-MESH quantities
      iotas(1) = 0; mass(1) = 0; pres(1) = 0; phip(1) = 0; 
      buco(1) = 0; bvco(1) = 0; vp(1) = 0; overr(1) = 0;  specw(1) = 1
      beta_vol(1) = 0
      CALL cdf_write(nwout, vn_iotah, iotas(1:ns))
      CALL cdf_write(nwout, vn_mass, mass/mu0) 
      CALL cdf_write(nwout, vn_presh, pres(1:ns)/mu0)
      CALL cdf_write(nwout, vn_betah, beta_vol)
      CALL cdf_write(nwout, vn_buco, buco)
      CALL cdf_write(nwout, vn_bvco, bvco) 
      CALL cdf_write(nwout, vn_vp, vp(1:ns))
      CALL cdf_write(nwout, vn_specw, specw)
      CALL cdf_write(nwout, vn_phip, phip(1:ns))
      CALL cdf_write(nwout, vn_overr, overr(1:ns))

!     MERCIER_CRITERION
      CALL cdf_write(nwout, vn_merc, Dmerc)
      CALL cdf_write(nwout, vn_mshear, Dshear)
      CALL cdf_write(nwout, vn_mwell, Dwell)
      CALL cdf_write(nwout, vn_mcurr, Dcurr)
      CALL cdf_write(nwout, vn_mgeo, Dgeod)
      CALL cdf_write(nwout, vn_equif, equif)

      CALL cdf_write(nwout, vn_fsq, fsqt(1:nstore_seq))
      CALL cdf_write(nwout, vn_wdot, wdot(1:nstore_seq))  
      IF (lasym) THEN
         CALL cdf_write(nwout, vn_racs, raxis(0:ntor,2))
         CALL cdf_write(nwout, vn_zacc, zaxis(0:ntor,2)) 
         CALL cdf_write(nwout, vn_rmns, rmns)
         CALL cdf_write(nwout, vn_zmnc, zmnc)
         CALL cdf_write(nwout, vn_lmnc, lmnc)
         CALL cdf_write(nwout, vn_gmns, gmns)
         CALL cdf_write(nwout, vn_bmns, bmns) 
         CALL cdf_write(nwout, vn_bsubumns, bsubumns)
         CALL cdf_write(nwout, vn_bsubvmns, bsubvmns)
         CALL cdf_write(nwout, vn_bsubsmnc, bsubsmnc)

         CALL cdf_write(nwout, vn_currumns, currumns)     !MRC  8-12-16
         CALL cdf_write(nwout, vn_currvmns, currvmns)

!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_write(nwout, vn_bsupumns, bsupumns)
         CALL cdf_write(nwout, vn_bsupvmns, bsupvmns)

         IF (lfreeb) THEN     !MRC    10-15-15
            CALL cdf_write(nwout, vn_bsubumns_sur, bsubumns_sur)
            CALL cdf_write(nwout, vn_bsubvmns_sur, bsubvmns_sur)
            CALL cdf_write(nwout, vn_bsupumns_sur, bsupumns_sur)
            CALL cdf_write(nwout, vn_bsupvmns_sur, bsupvmns_sur)
         END IF
      END IF
      CALL cdf_close(nwout)

      END SUBROUTINE
#endif


      SUBROUTINE Compute_Currents(bsubsmnc_, bsubsmns_,                          &
     &                            bsubumnc_, bsubumns_,                          &
     &                            bsubvmnc_, bsubvmns_,                          &
     &                            xm_nyq_, xn_nyq_, mnmax_nyq_,                  &
     &                            lasym_, ns_,     &
     &                            currumnc_, currvmnc_,                          &
     &                            currumns_, currvmns_)
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE

      REAL(rprec), DIMENSION(:,:), INTENT(in)  :: bsubsmnc_
      REAL(rprec), DIMENSION(:,:), INTENT(in)  :: bsubsmns_
      REAL(rprec), DIMENSION(:,:), INTENT(in)  :: bsubumnc_
      REAL(rprec), DIMENSION(:,:), INTENT(in)  :: bsubumns_
      REAL(rprec), DIMENSION(:,:), INTENT(in)  :: bsubvmnc_
      REAL(rprec), DIMENSION(:,:), INTENT(in)  :: bsubvmns_

      REAL(rprec), DIMENSION(:), INTENT(in)    :: xm_nyq_
      REAL(rprec), DIMENSION(:), INTENT(in)    :: xn_nyq_

      INTEGER, INTENT(in)                      :: mnmax_nyq_
      LOGICAL, INTENT(in)                      :: lasym_
      INTEGER, INTENT(in)                      :: ns_

      REAL(rprec), DIMENSION(:,:), INTENT(out) :: currumnc_
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: currvmnc_
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: currumns_
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: currvmns_

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js
      REAL(rprec) :: ohs, hs, shalf(ns_), sfull(ns_)
      REAL(rprec), DIMENSION(mnmax_nyq_) :: bu1, bu0, bv1, bv0, t1, t2, t3
!-----------------------------------------------
!
!     Computes current harmonics for currXmn == sqrt(g)*JsupX, X = u,v
!     [Corrected above "JsubX" to "JsupX", JDH 2010-08-16]

!     NOTE: bsub(s,u,v)mn are on HALF radial grid
!          (in earlier versions, bsubsmn was on FULL radial grid)

!     NOTE: near the axis, b_s is dominated by the m=1 component of gsv ~ cos(u)/sqrt(s)
!           we average it with a weight-factor of sqrt(s)
!

      ohs = (ns_-1)
      hs  = 1._dp/ohs

      DO js = 2, ns_
         shalf(js) = SQRT(hs*(js-1.5_dp))
         sfull(js) = SQRT(hs*(js-1))
      END DO

      DO js = 2, ns_-1
         WHERE (MOD(INT(xm_nyq_),2) .EQ. 1)
            t1 = 0.5_dp*(shalf(js+1)*bsubsmns_(:,js+1)                         &
               +         shalf(js)  *bsubsmns_(:,js)) /sfull(js)
            bu0 = bsubumnc_(:,js  )/shalf(js)
            bu1 = bsubumnc_(:,js+1)/shalf(js+1)
            t2 = ohs*(bu1-bu0)*sfull(js)+0.25_dp*(bu0+bu1)/sfull(js)
            bv0 = bsubvmnc_(:,js  )/shalf(js)
            bv1 = bsubvmnc_(:,js+1)/shalf(js+1)
            t3 = ohs*(bv1-bv0)*sfull(js)+0.25_dp*(bv0+bv1)/sfull(js)
         ELSEWHERE
            t1 = 0.5_dp*(bsubsmns_(:,js+1)+bsubsmns_(:,js))
            t2 = ohs*(bsubumnc_(:,js+1)-bsubumnc_(:,js))
            t3 = ohs*(bsubvmnc_(:,js+1)-bsubvmnc_(:,js))
         ENDWHERE
         currumnc_(:,js) = -xn_nyq_(:)*t1 - t3
         currvmnc_(:,js) = -xm_nyq_(:)*t1 + t2
      END DO         
   
      WHERE (xm_nyq_ .LE. 1)
         currvmnc_(:,1) =  2*currvmnc_(:,2) - currvmnc_(:,3)
         currumnc_(:,1) =  2*currumnc_(:,2) - currumnc_(:,3)
      ELSEWHERE
         currvmnc_(:,1) = 0
         currumnc_(:,1) = 0
      ENDWHERE

      currumnc_(:,ns_) = 2*currumnc_(:,ns_-1) - currumnc_(:,ns_-2)
      currvmnc_(:,ns_) = 2*currvmnc_(:,ns_-1) - currvmnc_(:,ns_-2)
      currumnc_ = currumnc_ /mu0;   currvmnc_ = currvmnc_/mu0

      IF (.NOT.lasym_) RETURN

      DO js = 2, ns_-1
         WHERE (MOD(INT(xm_nyq_),2) .EQ. 1)
            t1 = 0.5_dp*(shalf(js+1)*bsubsmnc_(:,js+1)                         &
               +         shalf(js)  *bsubsmnc_(:,js)) / sfull(js)
            bu0 = bsubumns_(:,js  )/shalf(js+1)
            bu1 = bsubumns_(:,js+1)/shalf(js+1)
            t2 = ohs*(bu1-bu0)*sfull(js) + 0.25_dp*(bu0+bu1)/sfull(js)
            bv0 = bsubvmns_(:,js  )/shalf(js)
            bv1 = bsubvmns_(:,js+1)/shalf(js+1)
            t3 = ohs*(bv1-bv0)*sfull(js)+0.25_dp*(bv0+bv1)/sfull(js)
         ELSEWHERE
            t1 = 0.5_dp*(bsubsmnc_(:,js+1) + bsubsmnc_(:,js))
            t2 = ohs*(bsubumns_(:,js+1)-bsubumns_(:,js))
            t3 = ohs*(bsubvmns_(:,js+1)-bsubvmns_(:,js))
         END WHERE
         currumns_(:,js) =  xn_nyq_(:)*t1 - t3
         currvmns_(:,js) =  xm_nyq_(:)*t1 + t2
      END DO         

      WHERE (xm_nyq_ .LE. 1)
         currvmns_(:,1) =  2*currvmns_(:,2) - currvmns_(:,3)
         currumns_(:,1) =  2*currumns_(:,2) - currumns_(:,3)
      ELSEWHERE
         currvmns_(:,1) = 0
         currumns_(:,1) = 0
      END WHERE
      currumns_(:,ns_) = 2*currumns_(:,ns_-1) - currumns_(:,ns_-2)
      currvmns_(:,ns_) = 2*currvmns_(:,ns_-1) - currvmns_(:,ns_-2)
      currumns_ = currumns_/mu0;   currvmns_ = currvmns_/mu0

      END SUBROUTINE Compute_Currents

      SUBROUTINE read_wout_deallocate
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat(10)
!-----------------------------------------------
      istat=0
      lwout_opened=.false.

      IF (ALLOCATED(extcur)) DEALLOCATE (extcur, curlabel,              &
              stat = istat(1))
      IF (ALLOCATED(overr)) DEALLOCATE (overr, stat = istat(2))

      IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, xm_nyq, xn_nyq,            &
        rmnc, zmns, lmns, bmnc, gmnc, bsubumnc, iotaf, presf, phipf,    &
        chipf,                                                          &
        bsubvmnc, bsubsmns, bsupumnc, bsupvmnc, currvmnc, iotas, mass,  &
        pres, beta_vol, phip, buco, bvco, phi, vp, jcuru, am, ac, ai,   &
        jcurv, specw, Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, jdotb, &
        bdotb, bdotgradv, raxis, zaxis, fsqt, wdot, stat=istat(3))
 
      IF (ALLOCATED(potsin)) DEALLOCATE (potsin)
      IF (ALLOCATED(potcos)) DEALLOCATE (potcos)

      IF (ALLOCATED(chipf)) DEALLOCATE (chipf, chi)

      IF (ALLOCATED(am_aux_s)) DEALLOCATE (am_aux_s, am_aux_f,          &
          ac_aux_s, ac_aux_f, ai_aux_s, ai_aux_f, stat=istat(6)) 

      IF (ireconstruct.gt.0 .and. ALLOCATED(sknots)) DEALLOCATE (       &
          ystark, y2stark, pknots, anglemse, rmid, qmid, shear,         &
          presmid, alfa, curmid, rstark, datastark, rthom, datathom,    &
          ythom, y2thom, plflux, dsiobt, bcoil, plbfld, bbc, sknots,    &
          pfcspec, limitr, rlim, zlim, nsetsn, stat = istat(4))

      IF (ALLOCATED(rmns)) DEALLOCATE (rmns, zmnc, lmnc,                &
          bmns, gmns, bsubumns, bsubvmns, bsubsmnc,                     &
          bsupumns, bsupvmns, stat=istat(5))

      IF (ALLOCATED(bsubumnc_sur)) THEN
         DEALLOCATE(bsubumnc_sur, bsubvmnc_sur,                         &
                    bsupumnc_sur, bsupvmnc_sur)
      END IF

      IF (ALLOCATED(bsubumns_sur)) THEN
         DEALLOCATE(bsubumns_sur, bsubvmns_sur,                         &
                    bsupumns_sur, bsupvmns_sur)
      END IF

!  Note currvmnc deallocated above.
      IF (ALLOCATED(currumnc)) DEALLOCATE (currumnc)
      IF (ALLOCATED(currumns)) DEALLOCATE (currumns, currvmns)
      IF (ALLOCATED(rzl_local)) DEALLOCATE (rzl_local)

      IF (ANY(istat .ne. 0))                                            &
            STOP 'Deallocation error in read_wout_deallocate'

      END SUBROUTINE read_wout_deallocate

      SUBROUTINE tosuvspace (s_in, u_in, v_in, gsqrt,                   &
                             bsupu, bsupv, jsupu, jsupv, lam)
      USE stel_constants, ONLY: zero, one
      IMPLICIT NONE
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!------------------------------------------------
      REAL(rprec), INTENT(in) :: s_in, u_in, v_in
      REAL(rprec), INTENT(out), OPTIONAL :: gsqrt, bsupu, bsupv,        &
                                            jsupu, jsupv, lam
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      INTEGER :: m, n, n1, mn, ipresent, jslo, jshi
      REAL(rprec) :: hs1, wlo, whi, wlo_odd, whi_odd
      REAL(rprec), DIMENSION(mnmax_nyq) :: gmnc1, gmns1, bsupumnc1,     &
         bsupumns1, bsupvmnc1, bsupvmns1, jsupumnc1, jsupumns1,         &
         jsupvmnc1, jsupvmns1, wmins, wplus,  lammns1, lammnc1
      REAL(rprec) :: cosu, sinu, cosv, sinv, tcosmn, tsinmn, sgn
      REAL(rprec) :: cosmu(0:mnyq), sinmu(0:mnyq),                      &
                     cosnv(0:nnyq), sinnv(0:nnyq)
      LOGICAL :: lgsqrt, lbsupu, lbsupv, ljsupu, ljsupv, llam
!------------------------------------------------
!
!     COMPUTE VARIOUS HALF/FULL-RADIAL GRID QUANTITIES AT THE INPUT POINT
!     (S, U, V) , WHERE 
!        S = normalized toroidal flux (0 - 1),
!        U = poloidal angle 
!        V = N*phi = toroidal angle * no. field periods
!
!     HALF-RADIAL GRID QUANTITIES
!     gsqrt, bsupu, bsupv
!   
!     FULL-RADIAL GRID QUANTITIES
!     dbsubuds, dbsubvds, dbsubsdu, dbsubsdv
!
!------------------------------------------------
      IF (s_in.lt.zero .or. s_in.gt.one) THEN
         WRITE(6, *) ' In tosuvspace, s(flux) must be between 0 and 1'
         RETURN
      END IF

      IF (.not.lwout_opened) THEN
         WRITE(6, *)' tosuvspace can only be called AFTER opening wout file!'
         RETURN
      END IF

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(u_in);   sinu = SIN(u_in)
      cosv = COS(v_in);   sinv = SIN(v_in)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mnyq
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      DO n = 1, nnyq
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO


!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE HALF-RADIAL GRID...
!     s-half(j) = (j-1.5)*hs, for j = 2,...ns
!
      hs1 = one/(ns-1)
      jslo = INT(c1p5 + s_in/hs1)
      jshi = jslo+1
      wlo = (hs1*(jshi-c1p5) - s_in)/hs1
      whi = 1 - wlo
      IF (jslo .eq. ns) THEN
!        USE Xhalf(ns+1) = 2*Xhalf(ns) - Xhalf(ns-1) FOR "GHOST" POINT VALUE 1/2hs OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns) + whi*Xhalf(ns+1) == Xhalf(ns) + whi*(Xhalf(ns) - Xhalf(ns-1)) 
         jshi = jslo-1
         wlo = 1+whi; whi = -whi
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
! 
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-c1p5)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-c1p5)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-c1p5)))
      END IF

      WHERE (MOD(NINT(xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      lgsqrt = PRESENT(gsqrt)
      IF (lgsqrt) THEN
         gsqrt = 0 ;  ipresent = ipresent+1
         gmnc1 = wmins*gmnc(:,jslo) + wplus*gmnc(:,jshi)
         IF (lasym) gmns1 = wmins*gmns(:,jslo) + wplus*gmns(:,jshi)
      END IF
      lbsupu = PRESENT(bsupu)
      IF (lbsupu) THEN
         bsupu = 0 ;  ipresent = ipresent+1
         bsupumnc1 = wmins*bsupumnc(:,jslo) + wplus*bsupumnc(:,jshi)
         IF (lasym) bsupumns1 = wmins*bsupumns(:,jslo) + wplus*bsupumns(:,jshi)
      END IF
      lbsupv = PRESENT(bsupv)
      IF (lbsupv) THEN
         bsupv = 0 ;  ipresent = ipresent+1
         bsupvmnc1 = wmins*bsupvmnc(:,jslo) + wplus*bsupvmnc(:,jshi)
         IF (lasym) bsupvmns1 = wmins*bsupvmns(:,jslo) + wplus*bsupvmns(:,jshi)
      END IF
      llam = PRESENT(lam)
      IF (llam) THEN
         lam = zero ;  ipresent = ipresent+1
         lammns1 = wmins*lmns(:,jslo) + wplus*lmns(:,jshi)
         IF (lasym) lammnc1 = wmins*lmnc(:,jslo) + wplus*lmnc(:,jshi)
      END IF

      IF (ipresent .eq. 0) GOTO 1000

!
!     COMPUTE LAM, ... IN REAL SPACE 
!     tcosmn = cos(mu - nv);  tsinmn = sin(mu - nv)
!     Half Grid but only mnmax sized
      IF (llam) THEN
         DO mn = 1, mnmax
            m = NINT(xm(mn));  n = NINT(xn(mn))/nfp
            n1 = ABS(n);   sgn = SIGN(1,n)   
            tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
            lam   = lam   + lammns1(mn)*tsinmn
         END DO
         IF (lasym) THEN
            DO mn = 1, mnmax
               m = NINT(xm(mn));  n = NINT(xn(mn))/nfp
               n1 = ABS(n);   sgn = SIGN(1,n)
               tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1) 
               lam   = lam   + lammnc1(mn)*tcosmn
            END DO
         END IF
      END IF
!
!     COMPUTE GSQRT, ... IN REAL SPACE
!     tcosmn = cos(mu - nv);  tsinmn = sin(mu - nv)
!

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)   
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmnc1(mn)*tcosmn
         IF (lbsupu) bsupu = bsupu + bsupumnc1(mn)*tcosmn
         IF (lbsupv) bsupv = bsupv + bsupvmnc1(mn)*tcosmn
         !IF (llam)   lam   = lam   + lammns1(mn)*tsinmn
      END DO

      IF (.not.lasym) GOTO 1000

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)   
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmns1(mn)*tsinmn
         IF (lbsupu) bsupu = bsupu + bsupumns1(mn)*tsinmn
         IF (lbsupv) bsupv = bsupv + bsupvmns1(mn)*tsinmn
         !IF (llam)   lam   = lam   + lammnc1(mn)*tcosmn
      END DO

 1000 CONTINUE

!     FULL-MESH QUANTITIES
!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE FULL-RADIAL GRID...
!     s-full(j) = (j-1)*hs, for j = 1,...ns
!
      hs1 = one/(ns-1)
      jslo = 1+INT(s_in/hs1)
      jshi = jslo+1
      IF (jslo .eq. ns) jshi = ns
      wlo = (hs1*(jshi-1) - s_in)/hs1
      whi = 1 - wlo
!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
! 
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-1)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-1)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-1)))
      END IF

      WHERE (MOD(NINT(xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      ljsupu = PRESENT(jsupu)
      IF (ljsupu) THEN
         IF (.not.lgsqrt) STOP 'MUST compute gsqrt for jsupu'
         jsupu = 0 ;  ipresent = ipresent+1
         jsupumnc1 = wmins*currumnc(:,jslo) + wplus*currumnc(:,jshi)
         IF (lasym) jsupumns1 = wmins*currumns(:,jslo) + wplus*currumns(:,jshi)
      END IF

      ljsupv = PRESENT(jsupv)
      IF (ljsupv) THEN
         IF (.not.lgsqrt) STOP 'MUST compute gsqrt for jsupv'
         jsupv = 0 ;  ipresent = ipresent+1
         jsupvmnc1 = wmins*currvmnc(:,jslo) + wplus*currvmnc(:,jshi)
         IF (lasym) jsupvmns1 = wmins*currvmns(:,jslo) + wplus*currvmns(:,jshi)
      END IF

      IF (ipresent .eq. 0) RETURN

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)   
         IF (ljsupu) jsupu = jsupu + jsupumnc1(mn)*tcosmn
         IF (ljsupv) jsupv = jsupv + jsupvmnc1(mn)*tcosmn
      END DO

      IF (.not.lasym) GOTO 2000

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (ljsupu) jsupu = jsupu + jsupumns1(mn)*tsinmn
         IF (ljsupv) jsupv = jsupv + jsupvmns1(mn)*tsinmn
      END DO

 2000 CONTINUE

      IF (ljsupu) jsupu = jsupu/gsqrt
      IF (ljsupv) jsupv = jsupv/gsqrt

      END SUBROUTINE tosuvspace

      SUBROUTINE LoadRZL
      IMPLICIT NONE
!------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------
      INTEGER     :: rcc, rss, zsc, zcs, rsc, rcs, zcc, zss
      INTEGER     :: mpol1, mn, m, n, n1
      REAL(rprec) :: sgn
!------------------------------------------------
!
!     Arrays must be stacked (and ns,ntor,mpol ordering imposed)
!     as coefficients of cos(mu)*cos(nv), etc
!     Only need R, Z components(not lambda, for now anyhow)
!      
      IF (ALLOCATED(rzl_local)) RETURN

      mpol1 = mpol-1
      rcc = 1;  zsc = 1
      IF (.not.lasym) THEN
         IF (lthreed) THEN
            ntmax = 2
            rss = 2;  zcs = 2
         ELSE
            ntmax = 1
         END IF
      ELSE
         IF (lthreed) THEN
            ntmax = 4
            rss = 2;  rsc = 3;  rcs = 4
            zcs = 2;  zcc = 3;  zss = 4
         ELSE
            ntmax = 2
            rsc = 2;  zcc = 2
         END IF
      END IF

!     really only need to ALLOCATE 2*ntmax (don't need lambdas)
!     for consistency, we'll allocate 3*ntmax and set lambdas = 0
      zsc = 1+ntmax; zcs = zcs+ntmax; zcc = zcc+ntmax; zss = zss+ntmax
      ALLOCATE(rzl_local(ns,0:ntor,0:mpol1,3*ntmax), stat=n)
      IF (n .ne. 0) STOP 'Allocation error in LoadRZL'
      rzl_local = 0

      DO mn = 1, mnmax
         m = NINT(xm(mn));  n = NINT(xn(mn))/nfp; n1 = ABS(n)
         sgn = SIGN(1, n)
         rzl_local(:,n1,m,rcc) = rzl_local(:,n1,m,rcc) + rmnc(mn,:)
         rzl_local(:,n1,m,zsc) = rzl_local(:,n1,m,zsc) + zmns(mn,:)
         IF (lthreed) THEN
            rzl_local(:,n1,m,rss) = rzl_local(:,n1,m,rss) + sgn*rmnc(mn,:)
            rzl_local(:,n1,m,zcs) = rzl_local(:,n1,m,zcs) - sgn*zmns(mn,:)
         END IF
         IF (lasym) THEN
            rzl_local(:,n1,m,rsc) = rzl_local(:,n1,m,rsc) + rmns(mn,:)
            rzl_local(:,n1,m,zcc) = rzl_local(:,n1,m,zcc) + zmnc(mn,:)
            IF (lthreed) THEN
                rzl_local(:,n1,m,rcs) = rzl_local(:,n1,m,rcs)           &
                                      - sgn*rmns(mn,:)
                rzl_local(:,n1,m,zss) = rzl_local(:,n1,m,zss)           &
                                      + sgn*zmnc(mn,:)
            END IF
         END IF
      END DO

      END SUBROUTINE LoadRZL

      END MODULE read_wout_mod

