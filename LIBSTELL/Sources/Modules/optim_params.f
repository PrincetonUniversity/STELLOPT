      MODULE optim_params
      USE vparams, ONLY: rprec, dp, nsd, ntord, ntor1d, mpol1d
      USE vsvd0, ONLY: nigroup
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NumJstard = 10
      INTEGER, PARAMETER :: ini_max=100                                  !! COBRA
      REAL(rprec) :: bigno = 1.e10_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: niter_opt, NumJstar, NumJinvariant,
     1   mboz_opt, nboz_opt,  nbmn, nproc, jboot,                        !LPK
     2   NS_JConf_Src, NS_JConf_Tgt, NPitch_Jconf                        !PPPL
      INTEGER :: num_processors, num_levmar_params
      INTEGER, DIMENSION(20) :: n_jac, m_jac
      INTEGER, DIMENSION(20) :: n_vac_island, m_vac_island
      INTEGER, DIMENSION(20) :: n_bmn_tgt, m_bmn_tgt                     !PPPL
      REAL(rprec) :: r00_scale, b00_scale, epsfcn, rgrid_min, rgrid_max,
     1               zgrid_min, zgrid_max, r00_opt
      REAL(rprec) :: sigma_jstar(nsd, NumJstard),
     1   sigma_jinvariant(nsd, NumJstard), sigma_jconf
      REAL(rprec), DIMENSION(nsd) ::
     1   sigma_iota, sigma_mercier, sigma_vp
      REAL(rprec), DIMENSION(nsd) :: sigma_bmin,
     1   sigma_bmax, sigma_ripple, sigma_bmn, nsurf_mask
      REAL(rprec) :: sigma_aspect, sigma_ellipticity, sigma_maxcurrent,
     1   sigma_coil_complex, sigma_curv, sigma_beta, target_aspectratio,
     2   target_maxcurrent, target_beta, target_rmax, target_rmin,
     3   target_ellipticity, sigma_rmax, sigma_rmin, sigma_iota_max,
     4   sigma_iota_min, target_iota_max, target_iota_min,
     5   sigma_iota_max_min, target_iota_max_min,                        !PPPL
     6   sigma_berr_avg, sigma_berr_max, target_eplasma, sigma_eplasma,  !PPPL
     7   target_curtor, sigma_curtor,                                    !PPPL
     8   sigma_centering                                                 !!Obsolete, replaced by sigma_rmax, sigma_rmin
      REAL(rprec) :: coil_separation, phiedge_min, phiedge_max
      REAL(rprec), DIMENSION(0:10) :: target_iota, target_well,
     1                                aseedcur, fboot                    !!LPK (note: at removed due to conflict with at in VMEC namelist)
!     1   at, aseedcur, fboot     !!LPK
      REAL(rprec) :: sigma_bal, sigma_boot, zeff_boot,                   !!LPK
     1    sigma_fluxp, target_fluxp, sigma_zmax, target_zmax,
     2    sigma_rbtor, target_rbtor                                      !PPPL
!     2    sigma_pseudo, sigma_pseudo2, sigma_rbtor, target_rbtor        !ORNL
      real(rprec), dimension(nsd) :: sigma_pseudo, sigma_pseudo2         !PPPL
      REAL(rprec), DIMENSION(9) :: sigma_kink, target_kink
      REAL(rprec), DIMENSION(20) :: sigma_jac, sigma_vac_island
      real(rprec), dimension(128) :: phi_lim                              !PPPL
      real(rprec), dimension(128,128) :: r_lim, z_lim                      !PPPL
      real(rprec), dimension(nsd, 20) :: sigma_bmn_tgt, target_bmn_tgt   !PPPL

      COMPLEX :: helicity
      LOGICAL :: lfix_ntor(-ntord:ntord), 
     1           lfix_rhob(-ntord:ntord,0:mpol1d+1)
      LOGICAL, DIMENSION(nsd) :: lsurf_mask, ldkes_mask
      LOGICAL :: lextcur(nigroup), lcoilp_sep
      LOGICAL :: lreset_opt, lbmn, lcoil_complex, lboundary,             !PPPL
     1   lcur_prof_opt, liota_prof_opt, lbootsj_opt, lj_star,
     2   lj_invariant, lcoil_opt, lcoil_geom, laspect_max, lbeta_min,
     3   lbal_opt, lbootstrap, lseedcur, lpress_opt, lkink_opt,          !!LPK
     4   lnescoil_opt, lvv_tgt_min, lballoon_flip, lpseudo_sin,
     5   lcurprof_opt, lprof_opt,                                         !!Obsolete, replaced by lcur_prof_opt, liota_prof_opt
     6   ledge_current, lvac_opt, lphiedge,phiedge_diode
!      LOGICAL :: lreset_opt, lbmn, lcoil_complex, ledge_current,         !ORNL
!     1   lcur_prof_opt, liota_prof_opt, lbootsj_opt, lj_star,
!     2   lj_invariant, lcoil_opt, lcoil_geom, laspect_max, lbeta_min,
!     3   lvac_opt, lphiedge,
!     4   lbal_opt, lbootstrap, lseedcur, lkink_opt, lnescoil_opt         !!LPK
!     5  ,lcurprof_opt, lprof_opt, lpress_opt, lboundary                  !!Obsolete, replaced by lcur_prof_opt, liota_prof_opt, lpres_prof_opt
      LOGICAL :: lreconp, lreconj, lp1zero, lj1zero
      INTEGER :: kpp, kjj
      CHARACTER(len=100) :: seq_ext, opt_ext
      CHARACTER(len=200) :: v3rfun_dir, v3post_in
      LOGICAL :: lv3post
      CHARACTER(len=4) :: sym_type

      REAL(rprec) :: target_coil_complex, target_coil_jmax               !!LPK
      REAL(rprec) :: sigma_coil_jmax, sigma_berr_ave                     !!LPK, SPH

      REAL(rprec), DIMENSION(nsd) :: sigma_neo, nneo_mask                !! NEO
      LOGICAL, DIMENSION(nsd) :: lneo_mask                               !! NEO
      LOGICAL :: lneo_opt                                                !! NEO
      REAL(rprec), DIMENSION(nsd) :: sigma_dsubr
      LOGICAL :: ldsubr_opt
      REAL(rprec) :: sigma_orbit
      LOGICAL :: lorbit_opt
      INTEGER :: nopt_alg, nopt_boundary
      LOGICAL :: ldiag_opt, lkeep_mins                                   !! Diagnostic output

      REAL(rprec) :: sigma_kappa, target_kappa, sigma_oh
!      REAL(rprec), DIMENSION(nigroup) :: sigma_extcur, oh_coefs         !ORNL
      REAL(rprec), DIMENSION(nigroup) :: sigma_extcur, target_extcur,    !PPPL
     1                                   oh_coefs                        !PPPL
      REAL(rprec), DIMENSION(0:10) :: target_iota_p
      REAL(rprec), DIMENSION(nsd) :: sigma_iota_pmax, sigma_iota_pmin
      REAL(rprec) :: sigma_jedge, target_jedge                           !PPPL

      INTEGER :: nini_theta, nini_zeta, nini_tot                         !! COBRA
      REAL(rprec), DIMENSION(nsd) :: sigma_balloon, sigma_pgrad,
     1                               target_balloon                      !! COBRA
      REAL(rprec), DIMENSION(ini_max):: bal_theta0, bal_zeta0            !! COBRA

      REAL(rprec) :: sigma_pedge(1)                                      !! COBRA
      REAL(rprec) :: sigma_bootsj(nsd)
      LOGICAL :: lballoon_mask(nsd)                                      !! COBRA
!      LOGICAL :: lballoon_opt, lpres_prof_opt                            !! COBRA (ORNL)
      LOGICAL :: lballoon_opt, lpres_prof_opt, lpres_opt_edge0,          !! COBRA (PPPL)
     1           lpres_opt_edgegr0, lcur_opt_edge0
      INTEGER :: ac_mask(0:10)                                           !! SAL AC masking
      INTEGER :: pres_opt_nmax
      LOGICAL :: lpres_prof_fit                                          !! SAL AM fitting

      REAL(rprec) :: nballoon_mask(nsd)                                  !! VMECCOBRA (RS)
      LOGICAL :: l_legendre                                              !! LEGENDRE (RS)

      REAL(rprec), DIMENSION(nsd) :: ndkes_mask, dkes_nu, dkes_efield,
     1     sigma_dkes                                                    !! RHF
      LOGICAL :: ldkes_opt                                               !! RHF

!      REAL(rprec) :: sigma_vv, sigma_vv_rms, vv_dist, vv_dist_rms,
!     1               target_vv, target_vv_rms                            !! RH & MZ (ORNL)
      REAL(rprec) :: sigma_vv, sigma_vv_rms, vv_dist, vv_dist_rms,
     1           target_vv, target_vv_rms, sigma_vv_max, vv_dist_max     !! RH & MZ (PPPL)
      INTEGER :: mpol_vv, ntor_vv, nu_vv, nv_vv
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) ::
     1   rbc_vv, zbs_vv
     
      REAL(rprec) :: sigma_bd, sigma_bd_rms, bd_dist, bd_dist_rms,
     1           target_bd, target_bd_rms, sigma_bd_max, bd_dist_max     !! MZ (PPPL)
      INTEGER :: mpol_bd, ntor_bd, nu_bd, nv_bd                          !PPPL
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: rbc_bd, zbs_bd    !PPPL
     
      LOGICAL :: shapeweight
CEAL    !       deviation weighting defaults
      REAL(rprec) :: theta0_bw(3), phi0_bw,
     1   wtheta_bw, wphi_bw, amplw_bw(3),  planes_bw(3)

      REAL(rprec) :: sigma_diagno(1000), data_diagno(1000)              !ORNL
      CHARACTER(len=30) :: name_diagno(1000)                            !ORNL
      LOGICAL :: ldiagno_opt                                             !! DIAGNO mcz
      CHARACTER(256) :: diagno_control
      REAL(rprec), DIMENSION(512) :: sigma_diagno_seg,                   !! diagno mcz
     1         target_diagno_seg, sigma_diagno_flx, target_diagno_flx    !! diagno mcz
      REAL(rprec), DIMENSION(512) :: target_diagno_bp, sigma_diagno_bp   !! diagno SAL
      INTEGER :: ndiagno_seg, ndiagno_flx, ndiagno_bp
      CHARACTER(256) :: diagno_coil                                     !SAL DIAGNO v3 coil file (for total response)

      
      INTEGER :: np_prof                        ! MCZ  pressure profile matching
      INTEGER :: nne_prof
      INTEGER :: nte_prof
      INTEGER :: nti_prof
      REAL(rprec), DIMENSION(512) :: r_te_prof, phi_te_prof, z_te_prof,
     1                               te_prof, sigma_te_prof
      REAL(rprec), DIMENSION(512) :: r_ne_prof, phi_ne_prof, z_ne_prof,
     1                               ne_prof, sigma_ne_prof
      REAL(rprec), DIMENSION(512) :: r_ti_prof, phi_ti_prof, z_ti_prof,
     1                               ti_prof, sigma_ti_prof
      REAL(rprec), DIMENSION (200)   :: ti_aux_s,ti_aux_f,
     1                                  te_aux_s,te_aux_f,
     2                                  ne_aux_s,ne_aux_f
      REAL(rprec), DIMENSION(512) :: sigma_p_prof, r_p_prof, z_p_prof,
     1           phi_p_prof, p_prof
      REAL(rprec) :: sigma_p_damp, factor_p_prof
      LOGICAL :: lp_prof_incl_edge
      LOGICAL :: isote                                                  ! SAL - ISO-Te contraint
      
      INTEGER :: n_xrc, m_xrc, nint_xrc
      REAL(rprec), DIMENSION(64,32) :: r_xrc, phi_xrc, z_xrc
      REAL(rprec), DIMENSION(64) :: signal_xrc

      INTEGER :: n_emis                                                  !! SXR chords
      REAL(rprec), DIMENSION(0:10) :: aemis                              !! mcz
      CHARACTER(256) :: emis_file
      REAL(rprec) :: sigma_emis_damp
      REAL(rprec), DIMENSION(512) :: sigma_emis, emis_chord
      
      LOGICAL :: lanimec                                                 !! PPPL - Use ANIMEC over VMEC
      LOGICAL :: lani_bcrit, lani_tperp, lani_phot
      INTEGER, DIMENSION(0:10) :: at_mask, ah_mask
      
      ! MSE Related (SAL 08/26/11)
      LOGICAL :: lmse_er
      INTEGER :: nmse_cams, nmse_chords
      REAL(rprec), DIMENSION (200)   :: er_aux_s,er_aux_f,
     1                                  ez_aux_s,ez_aux_f
      REAL(rprec), DIMENSION (8,128) :: mse_r, mse_phi, mse_z,
     1                                   mse_alpha,mse_beta,mse_theta,
     2                                   mse_pol,sigma_mse_pol,mse_vac,
     3                                   mse_a1_coef,mse_a2_coef,
     4                                   mse_a3_coef,mse_a4_coef,
     5                                   mse_a5_coef,mse_a6_coef,
     6                                   mse_a7_coef,mse_er,mse_ez
      

!     PPPL+ORNL Alphabetized Namelist
      NAMELIST /optimum/ ac_mask, aemis, amplw_bw, aseedcur,             &
     &   at_mask, ah_mask,                                               &
     &   b00_scale, bal_theta0, bal_zeta0, coil_separation,              &
     &   diagno_control, data_diagno, dkes_efield, dkes_nu, emis_chord,  &
     &   emis_file, epsfcn, factor_p_prof, fboot, helicity, jboot, kjj,  &
     &   kpp, l_legendre, laspect_max, lanimec, lani_bcrit, lani_tperp,  &
     &   lani_phot, lbal_opt, lballoon_flip,                             &
     &   lballoon_mask, lballoon_opt, lbeta_min, lbmn, lbootsj_opt,      &
     &   lbootstrap, lboundary, lcoil_complex, lcoil_geom, lcoil_opt,    &
     &   lcur_opt_edge0, lcur_prof_opt, lcurprof_opt, ldiag_opt,         &
     &   ldiagno_opt, ldkes_mask, ldkes_opt, ldsubr_opt, ledge_current,  &
     &   lextcur, lfix_ntor, lfix_rhob, liota_prof_opt, lj1zero,         &
     &   lj_invariant, lj_star, lkeep_mins, lkink_opt, lneo_mask,        &
     &   lneo_opt, lnescoil_opt, lorbit_opt, lp_prof_incl_edge, lp1zero, &
     &   lphiedge, lpress_opt, lpres_opt_edge0, lpres_opt_edgegr0,       &
     &   lpres_prof_opt, lpres_prof_fit, lprof_opt, lpseudo_sin,         &
     &   lreconj, lreconp, lreset_opt, lseedcur, lsurf_mask, lv3post,    &
     &   lvac_opt, lvv_tgt_min, m_bmn_tgt, m_jac, m_vac_island,          &
     &   mboz_opt, mpol_bd, mpol_vv, n_bmn_tgt, n_emis, n_jac,           &
     &   n_vac_island, name_diagno, nballoon_mask, nbmn, nboz_opt,       &
     &   ndkes_mask, niter_opt, nneo_mask, nopt_alg, nopt_boundary,      &
     &   np_prof, NPitch_Jconf, nproc, NS_JConf_Src, NS_JConf_Tgt,       &
     &   nsurf_mask, ntor_bd, ntor_vv, nu_vv, num_levmar_params,         &
     &   num_processors, NumJInvariant, NumJstar, nv_bd, nv_vv,          &
     &   oh_coefs, p_prof, phi_lim, phi_p_prof, phi0_bw, phiedge_max,    &
     &   phiedge_min, planes_bw, pres_opt_nmax, r00_opt, r00_scale,      &
     &   r_lim, r_p_prof, rbc_bd, rbc_vv, rgrid_max, rgrid_min,          &
     &   shapeweight, sigma_aspect, sigma_bal, sigma_balloon, sigma_bd,  &
     &   sigma_bd_max, sigma_bd_rms, sigma_berr_avg, sigma_berr_ave,     &
     &   sigma_berr_max, sigma_beta, sigma_bmax, sigma_bmin, sigma_bmn,  &
     &   sigma_bmn_tgt, sigma_boot, sigma_bootsj, sigma_centering,       &
     &   sigma_coil_complex, sigma_coil_jmax, sigma_curtor, sigma_curv,  &
     &   sigma_diagno, sigma_diagno_flx, sigma_diagno_seg,               &
     &   sigma_diagno_bp, sigma_dkes,                                    &
     &   sigma_dsubr, sigma_emis, sigma_emis_damp, sigma_ellipticity,    &
     &   sigma_eplasma, sigma_extcur, sigma_fluxp, sigma_iota,           &
     &   sigma_iota_max, sigma_iota_max_min, sigma_iota_min,             &
     &   sigma_iota_pmax, sigma_iota_pmin, sigma_jac, sigma_jconf,       &
     &   sigma_jedge, sigma_jinvariant, sigma_jstar, sigma_kappa,        &
     &   sigma_kink, sigma_maxcurrent, sigma_mercier, sigma_neo,         &
     &   sigma_oh, sigma_orbit, sigma_p_damp, sigma_p_prof, sigma_pedge, &
     &   sigma_pgrad, sigma_pseudo, sigma_pseudo2, sigma_rbtor,          &
     &   sigma_ripple, sigma_rmax, sigma_rmin, sigma_vac_island,         &
     &   sigma_vp, sigma_vv, sigma_vv_max, sigma_vv_rms, sigma_zmax,     &
     &   sym_type, target_aspectratio, target_balloon, target_bd,        &
     &   target_bd_rms, target_beta, target_bmn_tgt,                     &
     &   target_coil_complex, target_coil_jmax, target_curtor,           &
     &   target_diagno_flx, target_diagno_seg,                           &
     &   target_diagno_bp, target_ellipticity,                           &
     &   target_eplasma,target_fluxp, target_iota, target_iota_max,      &
     &   target_iota_max_min, target_iota_min, target_iota_p,            &
     &   target_jedge, target_kappa, target_kink, target_maxcurrent,     &
     &   target_rbtor, target_rmax, target_rmin, target_vv,              &
     &   target_vv_rms, target_well, target_zmax, theta0_bw, v3post_in,  &
     &   v3rfun_dir, wphi_bw, wtheta_bw, z_lim, z_p_prof, zbs_bd,        &
     &   zbs_vv, zeff_boot, zgrid_max, zgrid_min,                        &
     &   nmse_cams,nmse_chords,mse_r,mse_phi,mse_z,mse_alpha,mse_beta,   &
     &   mse_theta,mse_pol,mse_vac,sigma_mse_pol,mse_a1_coef,            &
     &   mse_a2_coef, mse_a3_coef,mse_a4_coef,mse_a5_coef,mse_a6_coef,   &
     &   mse_a7_coef, mse_er, mse_ez,                                    &
     &   phiedge_diode, isote,diagno_coil,target_extcur,                 &
     &   lmse_er, er_aux_s, er_aux_f, ez_aux_s, ez_aux_f,                &
     &   nne_prof, nte_prof, nti_prof,                                   &
     &   r_te_prof, phi_te_prof, z_te_prof, te_prof, sigma_te_prof,      &
     &   r_ne_prof, phi_ne_prof, z_ne_prof, ne_prof, sigma_ne_prof,      &
     &   r_ti_prof, phi_ti_prof, z_ti_prof, ti_prof, sigma_ti_prof,      &
     &   ti_aux_s, ti_aux_f, te_aux_s, te_aux_f, ne_aux_s, ne_aux_f
      

!     BELOW is the ORNL Input Namelist
!      NAMELIST /optimum/ epsfcn, niter_opt, num_processors,             & 
!     &   num_levmar_params, nsurf_mask, nopt_alg, nopt_boundary,        &
!     &   lreset_opt, ldiag_opt, lkeep_mins, lbmn, lj_star, laspect_max, &
!     &   lbeta_min, lj_invariant, liota_prof_opt, lcur_prof_opt,        &
!     &   ledge_current, lphiedge, lbootsj_opt, lkink_opt, lballoon_opt, &
!     &   l_legendre, ldkes_opt, lneo_opt, ldsubr_opt, lorbit_opt,       &
!     &   lpres_prof_opt, lnescoil_opt, lcoil_geom, lv3post, lvac_opt,   &
!     &   lfix_ntor, lfix_rhob, lextcur, sigma_extcur, oh_coefs,sigma_oh,&
!     &   r00_opt, r00_scale, b00_scale, rgrid_min, rgrid_max,           &
!     &   zgrid_min, zgrid_max, mboz_opt, nboz_opt, phiedge_min,         &
!     &   phiedge_max, coil_separation, target_aspectratio, sigma_aspect,&
!     &   target_beta, sigma_beta, target_kink, sigma_kink,              &
!     &   target_maxcurrent, sigma_maxcurrent, target_rmax, sigma_rmax,  &
!     &   target_rmin, sigma_rmin, target_zmax, sigma_zmax, target_kappa,&
!     &   sigma_kappa,target_ellipticity,sigma_ellipticity, target_fluxp,&
!     &   sigma_fluxp, target_rbtor, sigma_rbtor, target_coil_complex,   &
!     &   sigma_coil_complex, target_coil_jmax, sigma_coil_jmax,         &
!     &   target_iota, sigma_iota, target_iota_p, sigma_iota_pmax,       &
!     &   sigma_iota_pmin, target_iota_min, sigma_iota_min,              &
!     &   target_iota_max, sigma_iota_max, sigma_curv, sigma_berr_ave,   &
!     &   sigma_pseudo, sigma_pseudo2, sigma_mercier, sigma_jac, n_jac,  &
!     &   m_jac, n_vac_island, m_vac_island, sigma_vac_island,  helicity,&
!     &   sigma_bmin, sigma_bmax, sigma_bmn, sigma_ripple, NumJstar,     &
!     &   NumJInvariant, sigma_jstar, sigma_jinvariant, nballoon_mask,   &
!     &   target_balloon, sigma_balloon, sigma_pgrad, sigma_pedge,       &
!     &   bal_theta0, bal_zeta0, fboot, aseedcur,sigma_bootsj,ndkes_mask,&
!     &   dkes_nu, dkes_efield, sigma_dkes, nneo_mask, sigma_neo,        &
!     &   sigma_dsubr, sigma_orbit, v3rfun_dir, v3post_in, name_diagno,  &
!     &   data_diagno, sigma_diagno, target_vv, sigma_vv, target_vv_rms, &
!     &   lreconp, lp1zero, kpp, lreconj, lj1zero, kjj,                  &
!     &   sigma_vv_rms, mpol_vv, ntor_vv, nu_vv, nv_vv, rbc_vv, zbs_vv,  &
!     &   shapeweight, planes_bw, amplw_bw, theta0_bw,phi0_bw, wtheta_bw,&
!     &   wphi_bw,  
!     &   lcoil_complex, lballoon_mask, ldkes_mask, lseedcur, lneo_mask  &
!     &   ,lsurf_mask, lcurprof_opt, lprof_opt, target_well, lcoil_opt   & !!Obsolete, retain for consistency
!     &   ,sigma_centering, lpress_opt, lbootstrap, lbal_opt, lboundary  &
!     &   ,sigma_vp, sym_type, sigma_boot, sigma_bal, nproc, nbmn        &
!     &   ,zeff_boot, jboot, at
!
!        VARIABLE DESCRIPTIONS
!------------------------------------
!        SCALARS/VECTORS
!------------------------------------
!        nopt_alg             determines the algorithm for optimization
!                             0 = Levenberg-Marquardt (default)
!                             1 = Genetic (GA)
!                             2 = Differential Evolution (DE)
!        nopt_boundary        determines internal representation used for optimization 
!                             in fixed-boundary runs (LFREEB=.f)
!                             0 = Hirshman-Breslau (default)
!                             1 = Advanced HB
!                             2 = Garabedian delta_mn
!        niter_opt            Maximum number of optimization iterations
!        num_processors       Maximum number of processors to try and use simultaneously
!                             (obsolescent: nproc)
!        num_levmar_params
!        nproc
!        epsfcn               'Annealing' parameter, 1.E-2 or less. If
!                             lreset_opt=F, epsfcn <= 1.E-3; if lreset_opt=T,
!                             epsfcn <= 1.E-4 is recommended
!        r00_scale 
!        b00_scale
!        r00_opt
!        mboz_opt             User-input number of poloidal boozer modes
!        nboz_opt             User-input number of toroidal boozer modes
!        numJstar
!        numJInvariant
!        sym_type 
!        m_vac_island, n_vac_island 
!        m_jac, n_jac
!        jboot
!        fboot
!        at
!        aseedcur
!        dkes_nu
!        dkes_efield
!        mpol_vv, ntor_vv
!        coil_separation,
!        bal_theta0, bal_zeta0
!        rgrid_min, rgrid_max
!        zgrid_min, zgrid_max
!        nbmn
!        zeff_boot
!        oh_coefs
!        rbc_vv, zbs_vv
!        theta0_bw, phi0_bw, amplw_bw  
!        planes_bw
!------------------------------------
!        CONTROL-MASK ARRAYS
!------------------------------------
!        nsurf_mask           Real array (size=nrad) giving the fractional radial s-values
!                             at which optimizations are to be done for bmin, bmax, Jinvariant, 
!                             DKES, NEO calculations
!        lsurf_mask (obs)     Logical array (size=nrad); use nsurf_mask 
!        nballoon_mask        Real array (size=nrad) containing fractional radial s-values 
!                             where ballooning growth rates are to be computed
!        lballoon_mask (obs)  Logical array (size=nrad) designating surfaces where
!                             ballooning growth rates are to be computed (=T); use nballoon mask
!        lextcur              Logical array (size>=nextcur); if (i-th) element =T and lfreeb=True, 
!                             extcur(i) will be varied as an independent variable.
!        ndkes_mask
!        ldkes_mask
!        nneo_mask
!        lneo_mask
!------------------------------------
!        LOGICAL CONTROL VARIABLES
!------------------------------------
!        lreset_opt           =F, VMEC runs without resetting (after convergence)
!                             to a coarse radial grid
!                             =T (default), VMEC resets after each new iteration
!        lcur_prof_opt        =T, ncurr set to 1 and ac current series expansion coefficients
!                             are varied as independent variables
!                             =F (default), ac coefficients are fixed (if ncurr = 1)
!        liota_prof_opt       =T, ncurr set to 0 and ai iota series expansion coefficients 
!                             are varied as independent variables
!        ledge_current        =T, vary curtor (edge current), provided ncurr=1
!        lphiedge             =T, vary phiedge (edge toroidal flux)
!        lprof_opt (obs)      (see lcur_prof_opt and liota_prof_opt)
!        lcurprof_opt (obs)   (see lcur_prof_opt)
!        lpres_prof_opt       =T  optimize pressure profile shape (vary mass expansion coefficients
!                             as independent variables); used to achieve ballooning stability
!                             =F  keep pressure profile shape fixed
!        lpress_opt (obs)     (see lpres_prof_opt)
!        laspect_max          =T, then sigma_aspect = bigno if aspect ratio <= target_aspectratio.
!                             used to guarantee a no larger than the target.
!        lbeta_min            =T, then sigma_beta = bigno if beta >= Target_Beta. Used to
!                             guarantee a minimum beta
!        lkink_opt            =T, do global kink stability calculation
!        lballoon_opt         =T, do ballooning calculation
!                             =F, do not do ballooning calculation if nballoon_mask prescribed
!        lbal_opt (obs)       (use nballoon_mask)
!                             =F (or lfreeb=False), extcur(i) is fixed during optimization
!        lnescoil_opt         =T, evaluate NESCOIL coil optimization targets
!
!        lfix_ntor(n)         Logical array. =F, r0n's, z0n's are free to vary (i.e., they ARE NOT fixed)
!                             for the toroidal index 'n' of the arrays rbc(n,m=0), zbs(n,m=0)
!                             =T, r0n's, z0n's are fixed (not allowed to varied) during optimization
!                             Useful if certain n-s are to be kept fixed (such as the
!                             axi-symmetric boundary components: lfix_ntor(0) = T)
!                             NO EFFECT IN FREE-BDY OPTIMIZATION
!        lfix_rhob(n,m)       2D logical array. =F, rhobc(n,m) component is varied during optimization (default)
!                             =T, rhobc(n,m) compoents is fixed during optimization
!                             NO EFFECT IN FREE-BDY OPTIMIZATION
!        lbmn                 =T, add Boozer spectra target to chi-sq (for targetting symmetries: QA, QH, QP)
!        lbootsj_opt
!        lbootstrap (obs)     (see lbootsj_opt)
!        lseedcur
!        lj_star
!        lj_invariant
!        l_legendre
!        lkeep_mins
!        lcoil_complex
!        lcoil_opt            (obsolete)
!        lcoil_geom           =T, do coil geometry optimization
!                             (coil harmonics are independent variables in the optimization)
!                             Note that in this case, the initial coil shapes
!                             and the weights for the coil evaluations targets
!                             (e.g. coil curvatures and separations) are
!                             specified in the COILSIN NAMELIST (in coilsnamin.f in LIBSTELL).
!
!                             =F, coil geometry is FIXED and either the
!                             EXTERNAL currents are varied (lfreeb=.true., lextcur)
!                             or the bdy coefficients (rbc, zbs) are varied.
!        lcoilp_sep           a COMPUTED variable, which = F if lcoil_geom = T and
!                             MAY = T when lcoil_geom=F. Used to force computation of
!                             minimum coil-plasma separation even when the coil geometry is not
!                             evolving in the optimizer
!        lvac_opt             =T, use rotations and shifts of coils for independent optimization variables
!        ldkes_opt
!        lneo_opt
!        ldiag_opt
!        ldsubr_opt
!        lorbit_opt
!        phiedge_max,min      max/min allowable range for variation of phiedge, when lphiedge = TRUE
!------------------------------------
!        TARGET/WEIGHTS (SIGMA = 1/WEIGHT)
!------------------------------------
!        Equilibrium and Geometry
!------------------------------------
!        Target_AspectRatio   Aspect Ratio 
!        sigma_aspect         
!        Target_MaxCurrent    Maximum integrated toroidal current
!                             (bounding current, matched at each radial position)
!        sigma_maxcurrent
!        Target_Beta          Volume averaged beta to match
!        sigma beta
!        Target_Iota          Coefficients of power series (in flux s)
!                             for iota profile
!        sigma_iota
!        Target_Iota_P        Coefficients of power series (in s) for iota-prime
!                             = d(iota)/ds
!        sigma_iota_pmax      Array (size nrad) of sigmas for target_iota_p as
!                                   upper bound on d-iota/ds as
!        sigma_iota_pmin      Array (size nrad) of sigmas for target_iota_p as
!                                   lower bound on d-iota/ds as
!        Target_iota_min      minimum iota value (as lower bound)
!        sigma_iota_min
!        Target_iota_max      maximum iota value (as upper bound)
!        sigma_iota_max
!        Target_rmin          minimum major radius over the boundary
!        sigma_rmin           
!        Target_rmax          maximum major radius over the boundary
!        sigma_rmax           
!        Target_zmax          maximum height over the boundary
!        sigma_zmax           
!        Target_ellipticity   desired elongation of phi=0 cross section
!        sigma_ellipticity
!        Target_kappa         desired <kappa> (n=0 component of elongation)
!        sigma_kappa          
!        Target_fluxp         Minimum poloidal flux (Wb)
!        sigma_fluxp          
!        sigma_curv           Forces curvature kurtosis to zero at phi=0,90,180,270
!
!        Stability
!------------------------------------
!        Target_Well          Coefficients Power series (in flux s)
!                             for the magnetic well
!                             [Vp(s) - Vp(0)]/Vp(0) [Replaced with Mercier criterion]
!        sigma_vp             Array (nrad) of sigmas for well
!        sigma_Mercier        Sigma for Mercier stability (not relative)
!        Target_balloon       Real array (size=nrad) of ballooning eigenvalue (=0 for marginal stability)
!        sigma_balloon        Real array (size=nrad) 
!        sigma_pgrad          Real array (size=nrad) forcing local pressure gradient to zero
!        sigma_pedge          Real array (size=1) for forcing pressure value at edge to zero
!        Target_kink          array of kink eigenvalues. Can be used to bias the chi-square
!                             to achieve stability
!        sigma_kink           Array of sigmas for full kink stability calc.,
!                             The first sigma value is used in a call to xtprp/ft5tpr
!                             the second value is used for xtprp2/ft5tpr2,and so-on for successive values.  
!                             This allows targetting of multiple mode families. The sigma values are  *not relative*
!
!        Coils Currents
!------------------------------------
!        sigma_extcur         Array (dim=nextcur) of sigmas for the coil-currents used for regularizing 
!                             (forcing to zero) the coil currents in the extcur array for which the mask array 
!                             lextcur = T
!        Target_rbtor         The effective poloidal coil currents (R*Bt) used to constrain sum of varied external 
!                             coil currents ([T] - [m]). Only needed if lcoil_geom=F (fixed coil geometry); otherwise,
!                             handled by xcoilgeom.
!        sigma_rbtor  
!        oh_coefs             Array (size=nextcur) for imposing a linear sum constraint on the coil currents, 
!                             e.g. for ensuring that a PF set does not generate net poloidal flux.
!                             Constraint is sum_i (oh_coefs(i)*extcur(i))
!        sigma_oh             sigma for the OH flux constraint (weighted linear sum
!                             coil currents, targetted to zero)
!
!        NESCOIL current-sheet (lnescoil_opt=T)
!------------------------------------
!        Target_coil_complex  desired maximum coil complexity measure (<M> for coils on sheet)
!        sigma_coil_complex   
!        Target_coil_jmax     desired maximum coil current density
!        sigma_coil_jmax
!        sigma_berr_ave       forces NESCOIL berr to zero
!
!        COILOPT
!------------------------------------
!
!        Transport (see lfix_ntor, lfix_rhob, lbmn described above)
!------------------------------------
!        helicity             If lbmn=T, used to select helicity of Boozer spectra that
!                             influence symmetry hence transport properties.
!                             For kh = real(helicity), lh = imag(helicity), then
!                             kh = 0         =>   quasi-poloidal symmetry
!                             lh = 0         =>   quasi-axisymmetry
!                             m*kh+n*lh = 0  =>  quasi-helical symmetry (lu + kv)
!        sigma_bmn            Real array (size=nrad) for forcing specific bmns to zero
!                             (relative to largest bmn). Used to obtain quasi-symmetric spectra
!                             in conjunction with helicity specification.
!    Vacuum vessel and limiter matching related
!--------------------------------
!    ability to target two boundary shapes is provided (_vv and _bd)
!    normally _VV is used to target an enclosing PFC shape, and _BD is used to
!    target a desired plasma shape
!
!       target_vv               target for closest approach distance
!       lvv_tgt_min             if T, then there is no penalty if closest
!                               approach distance is larger than target_vv
!
!       target_vv_rms           target for RMS average distance
!       sigma_vv                sigma for closest distance to specified -VV
!                                  3D boundary shape (if <~1.e10, else ignore)
!                                  and limiters
!       sigma_vv_rms            sigma for RMS average distance from plasma to
!                                _VV 3D boundary (if <~ 1.e10, else ignore)
!       sigma_vv_max            sigma for Maximum deviation distance from
!                                  plasma to _VV 3D boundary
!       mpol_vv                 maximum m for _VV 3D boundary
!       ntor_vv                 maximum n for _VV 3D boundary
!       rbc_vv                  R-cos(theta) array for _VV 3D boundary
!       zbs_vv                  Z-sin(theta) array for _VV 3D boundary
!
!
!       target_bd               target for closest approach distance
!       target_bd_rms           target for RMS average distance
!       sigma_bd                sigma for closest distance to specified _BD
!                                  3D boundary shape (if <~1.e10, else ignore)
!                                  and limiters
!       sigma_bd_rms            sigma for RMS average distance from plasma to
!                                 _BD 3D boundary (if <~ 1.e10, else ignore)
!       sigma_bd_max            sigma for Maximum deviation distance from
!                                  plasma to _BD 3D boundary
!       mpol_bd                 maximum m for _BD 3D boundary
!       ntor_bd                 maximum n for _BD 3D boundary
!       rbc_bd                  R-cos(theta) array for _BD 3D boundary
!       zbs_bd                  Z-sin(theta) array for _BD 3D boundary
!
!  Piece-wise linear limiter specification, used with _VV targets
!
!       phi_lim                 toroidal angles of limiters (1-D array)
!                               limiters are specified at each toroidal angle
!                               as an array of piecewise-linear segments,
!                               going counter-clockwise around the plasma.
!                               For each phi_lim(iphi), the i-th segment
!                               goes from (r_lim(i,iphi),z_lim(i,iphi)) to
!                               (r_lim(i+1,iphi),zlim(i+1,iphi))
!       r_lim                   2D-array of major radial locations of limiter
!                               segments  (i, iphi)
!       z_lim                   2D-array of vertical locations of limiter
!                               segments  (i, iphi)
!
!    Reconstruction related
!--------------------------------
!       ldiagno_opt             simulate & match magnetic diagnostics with diagno
!       diagno_control          name of diagno control file to use
!       sigma_diagno_seg        array of sigmas for segmented coils
!       sigma_diagno_flx        array of sigmas for flux loops
!       target_diagno_seg       array of target values for segmented coils
!       target_diagno_flx       array of target values for flux loops
!
!       np_prof                 number of points to match in pressure profile
!       p_prof                  array of pressure values (depricated)
!       ne_prof                 array of electron number density (m^-3)
!       te_prof                 array of electron temperature (eV)
!       r_p_prof                array of major radius values for pressure points (m)
!       z_p_prof                array of height values for pressure points (m)
!       phi_p_prof              array of toroidal angles for pressure points (degrees)
!       sigma_p_prof            array of sigmas for pressure points
!
!       target_curtor           target total current
!       sigma_curtor            sigma for total current
!
!       target_jedge            target for edge current density
!       sigma_jedge             sigma for edge current density        
!    
!       nmse_cams               Number of MSE Cameras
!       nmse_chords             Maximum number of chords per camera
!       mse_r                   Array of Radial datapoint locations size(nmse_cams,nmse_chords)
!       mse_phi                 Array of toroidal datapoint locations size(nmse_cams,nmse_chords)
!       mse_z                   Array of vertical datapoint locations size(nmse_cams,nmse_chords)
!       mse_alpha               Array of angles between toroidal direction and beam size(nmse_cams,nmse_chords)
!       mse_beta                Array of angles between chord and beam size(nmse_cams,nmse_chords)
!       mse_theta               Array of angles between chord and midplane size(nmse_cams,nmse_chords)
!       mse_pol                 Array of measured MSE polarizations size(nmse_cams,nmse_chords)
!       mse_vac                 Array of vacuum MSE polarizations (size(nmse_cams,nmse_chords)
!       sigma_mse_pol           Array of sigmas for measured MSE polarizations size(nmse_cams,nmse_chords)
! 
!        target_iota_max, sigma_iota_max
!        target_iota_min, sigma_iota_min
!        target_vv, sigma_vv
!        target_vv_rms, sigma_vv_rms
!        sigma_bootsj (sigma_boot)
!        sigma_pseudo
!        sigma_pseudo2
!        sigma_bal
!        sigma_vac_island,
!        sigma_jstar 
!        sigma_jinvariant
!        sigma_jac
!        sigma_oh         
!        sigma_bmin
!        sigma_bmax 
!        sigma_ripple
!        sigma_curv
!        sigma_pgrad 
!        sigma_pedge
!        sigma_neo
!        sigma_dsubr
!        sigma_orbit
!        sigma_dkes
!        shapeweight
!        wtheta_bw
!        wphi_bw
!        sigma_centering
!
!        Magnetic Diagnostics
!------------------------------------
!        lv3post              =T, match diagnostic signals, =F, ignore diagnostic signals
!                             in optimization
!        v3rfun_dir           fully-qualified path to directory containing pre-computed
!                             response tables for the diagnostic set used for matching
!        v3post_in            fully-qualified name of v3post input namelist file

      CONTAINS

      SUBROUTINE read_optimum_namelist (iunit, istat)
      INTEGER, INTENT(in)    :: iunit
      INTEGER, INTENT(inout) :: istat

      istat = 0

!
!     DEFAULT VALUES
!
      r00_opt = -1
      r00_scale = 1
      b00_scale = 1
      kpp = 10 ; kjj = 10
      lreconp = .false.  ;  lreconj = .false.
      lp1zero = .false.  ;  lj1zero = .false. 

      mboz_opt = 0; nboz_opt = 0; nbmn=0; nproc=0
      opt_ext = 'none'
      v3rfun_dir = ' '
      v3post_in = ' '
      lbeta_min = .false.
      laspect_max = .false.
      lcoil_geom = .false.
      lvac_opt = .false.
      lcoilp_sep = .false.
      lv3post = .false.
      phiedge_diode = .false.
      sigma_berr_avg = bigno
      sigma_berr_max = bigno
      coil_separation = 0
      NumJstar = 0
      NumJinvariant = 0
      NPitch_JConf = 4
      NS_JConf_Src = 0
      NS_JConf_Tgt = 0
      target_iota_max = 1
      target_iota_max_min = 1
      target_iota_min = 1
      sigma_iota_max = bigno
      sigma_iota_max_min = bigno
      sigma_iota_min = bigno
      sigma_iota = bigno
      sigma_jedge = bigno
      target_jedge = 0
      sigma_vp   = bigno
      sigma_mercier = bigno
      target_kink = 0
      sigma_kink = bigno
      sigma_curv = bigno
      sigma_aspect = bigno
      Target_AspectRatio = 3
      sigma_coil_complex = bigno
      sigma_MaxCurrent = bigno
      target_MaxCurrent = 0
      sigma_beta = bigno
      target_beta = 0
      sigma_centering = bigno
      sigma_rmin = bigno
      sigma_rmax = bigno
      sigma_zmax = bigno
      target_rmin = 0
      target_rmax = 0
      target_zmax = 0
      sigma_ellipticity = bigno
      target_ellipticity = 0
      sigma_bmin = bigno
      sigma_bmax = bigno
      sigma_bmn  = bigno
      sigma_bmn_tgt = bigno
      target_bmn_tgt = 0
      n_bmn_tgt = 0
      m_bmn_tgt = 0
      sigma_jstar = bigno
      sigma_jinvariant = bigno
      sigma_jconf = bigno
      sigma_ripple = bigno
      sigma_bootsj = bigno
      sigma_boot   = bigno
!      at = 0
      sigma_jac = bigno
      n_jac = 0
      m_jac =0
      sigma_vac_island = bigno
      n_vac_island = 0
      m_vac_island = 0
      sigma_diagno = bigno
      data_diagno = 0
      name_diagno = ""

!
!     INITIALIZATION OBSOLETE VARS, TOO
!
      lcurprof_opt = .false.; lcoil_opt = .false.; lboundary = .false.

      num_levmar_params = 1
      num_processors =    1
      
!
!  Experiment matching (pressure profile totoal current, DIAGNO, SXR,etc)
!
      sigma_eplasma = bigno           ! MCZ plasma kinetic stored energy
      target_eplasma = 0
      sigma_curtor = bigno            ! MCZ plasma current
      target_curtor = 0

      ! Pressure
      sigma_p_prof = bigno            ! MCZ experimental pressure profile
      sigma_p_damp = bigno
      p_prof = 0
      factor_p_prof = 0
      lp_prof_incl_edge = .true.
      lpres_prof_fit = .false.        ! Use AM fitting not optimization
      r_p_prof = 0
      z_p_prof = 0
      phi_p_prof = 0  
      np_prof = 0
      
      ! Ne
      nne_prof = 0
      sigma_ne_prof = bigno
      ne_prof = 0      
      r_ne_prof = 0
      z_ne_prof = 0
      phi_ne_prof = 0  
      ne_aux_s = 0.0
      ne_aux_f = 1.0                  ! Needs to be 1.0 so code doesn't puke
      
      ! Te
      nte_prof = 0
      sigma_te_prof = bigno
      te_prof = 0 
      r_te_prof = 0
      z_te_prof = 0
      phi_te_prof = 0 
      te_aux_s = 0.0
      te_aux_f = 0.0
      
      ! Ti
      nti_prof = 0
      sigma_ti_prof = bigno
      ti_prof = 0            
      r_ti_prof = 0
      z_ti_prof = 0
      phi_ti_prof = 0  
      ti_aux_s = 0.0
      ti_aux_f = 0.0
      
      ! ISOTE (doesn't work yet)
      isote = .false.                 ! Default to old method if isote=.false.

      ! SXR - (not implemented)
      sigma_emis = bigno              ! MCZ SXR emissivity chords
      sigma_emis_damp = bigno
      n_emis = 0
      emis_file = " "
      aemis = 0

      ! Magnetics DIAGNO
      ldiagno_opt = .false.           ! diagno simulations of mag. diags.
      diagno_control = 'diagno.control'
      sigma_diagno_seg = bigno
      target_diagno_seg = 0
      sigma_diagno_flx = bigno
      target_diagno_flx = 0
      sigma_diagno_bp = bigno
      target_diagno_bp = 0
      diagno_coil = ''                ! SAL diagno coils file for total response
      ac_mask = 0                     ! SAL AC array masking (defaults to using all)
      
      ! MSE
      nmse_cams=0
      nmse_chords=0
      mse_r=0
      mse_phi=0
      mse_z=0
      mse_alpha=0
      mse_beta=0
      mse_theta=0
      mse_pol=0
      mse_vac=0
      mse_a1_coef = 0             ! Do this so we can detect mse coefficients being set
      mse_a2_coef = 0
      mse_a3_coef = 0
      mse_a4_coef = 0
      mse_a5_coef = 0
      mse_a6_coef = 0
      mse_a7_coef = 0
      mse_er = 0.0
      mse_ez = 0.0
      lmse_er = .FALSE.
      er_aux_s = 0.0
      er_aux_f = 0.0
      ez_aux_s = 0.0
      ez_aux_f = 0.0
      sigma_mse_pol=bigno
      
      
!
!     ANIMEC Additions (SAL - PPPL)
!

      lanimec = .false.
      lani_bcrit = .false.
      lani_tperp = .false.
      lani_phot = .false.
      at_mask = 0
      ah_mask = 0
!
!     LPK ADDITIONS
!
      sigma_fluxp = bigno
      sigma_pseudo = bigno
      sigma_pseudo2 = bigno
      nproc = 1
      lkink_opt = .false.
      lbal_opt = .false.
      lbootstrap = .false.
      jboot = 0
      fboot = 0
      fboot(0) = 1
      zeff_boot = 1
      lseedcur = .false.
      lpress_opt = .false.
      lcoil_complex = .false.
      lnescoil_opt = .false.
      sigma_coil_jmax = bigno
      sigma_berr_ave = bigno
      target_coil_complex = 1
      target_coil_jmax = 0
      ldkes_opt = .false.
      ldkes_mask = .false.
      sigma_dkes = bigno
      ndkes_mask = 0
      dkes_efield = 0
      dkes_nu = 0.001_dp
      aseedcur = 0
      Target_RBtor = 1
      Sigma_RBtor = bigno
!!
!!    NEO
!!
      sigma_neo = bigno
      lneo_opt = .false.
      nneo_mask = 0
      lneo_mask = .false.
      ldiag_opt = .false.
      lkeep_mins = .false.
      ldsubr_opt = .false.
      sigma_dsubr = bigno
      lorbit_opt = .false.
      sigma_orbit = bigno
!!
!!    END LPK ADDITIONS
!!
      Target_AspectRatio = 3
      Target_MaxCurrent = 0
      Target_Iota = 0
      Target_Well = (/ 0.0_dp, -0.39_dp, 0.19_dp, (0.0_dp,istat=1,8) /)
      niter_opt = 1
      nopt_alg = 0
      nopt_boundary = 0
      mboz_opt = 0
      nboz_opt = 0
      lreset_opt = .true.
      lcur_prof_opt = .false.
      lcurprof_opt = .false.       ! obsolete !
      liota_prof_opt = .false.
      ledge_current = .false.
      lphiedge = .false.
      lbmn = .false.
      lj_star = .false.
      lj_invariant = .false.
      lprof_opt = .false.
      lextcur = .false.
      lbootsj_opt = .false.
      lsurf_mask = .false.
      lfix_ntor = .false.
      lfix_rhob = .false.
      phiedge_max = 0;  phiedge_min = 0
      nsurf_mask = 0
      epsfcn = -1
      helicity = CMPLX(0._dp, 0._dp)
      sym_type = 'NONE'

!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (01/19/99): ballooning related variables and sigmas. (PPPL)

      target_balloon = 0
      lballoon_flip = .false.
      sigma_balloon = bigno
      sigma_bal = bigno                                                  !Old style - LPK
      bal_zeta0 = 0
      bal_theta0 = 0
      sigma_pgrad = bigno
      sigma_pedge = bigno
      lballoon_opt = .false.
      lballoon_mask = .false.
      lpres_prof_opt = .false.
      lpres_opt_edge0 = .false.
      lpres_opt_edgegr0 = .false.
      pres_opt_nmax = 11
      
      nballoon_mask = 0                                                  !VMECCOBRA
      l_legendre = .false.                                               !LEGENDRE

!------------------------------------------------------------------------------
!  initialize 3D boundary limiter parameters,  RH & MZ  June 2000 (PPPL)

      lvv_tgt_min = .false.
      target_vv = 0
      target_vv_rms = 0
      sigma_vv = bigno
      sigma_vv_rms = bigno
      sigma_vv_max = bigno
      rbc_vv = 0
      zbs_vv = 0
      mpol_vv = 0
      ntor_vv = 0
      nu_vv = 5
      nv_vv = 2

      target_bd = 0
      target_bd_rms = 0
      sigma_bd = bigno
      sigma_bd_rms = bigno
      sigma_bd_max = bigno
      rbc_bd = 0
      zbs_bd = 0
      mpol_bd = 0
      ntor_bd = 0
      nu_bd = 5
      nv_bd = 2

      r00_opt = -1
      r00_scale = 1
      b00_scale = 1
!      rmax_opt = 0; rmin_opt = 0; zmax_opt = 0
      rgrid_max = 0; rgrid_min = 0; zgrid_max = 0; zgrid_min = 0

      shapeweight = .FALSE.
      theta0_bw(1) = 0 ; theta0_bw(2) = 0 ; theta0_bw(3) = 0
      phi0_bw = 0 ; wtheta_bw = 0.392699 ; wphi_bw = 0.392699
      planes_bw(1) = 0 ; planes_bw(2) = 1.5707 ; planes_bw(3) = 3.14159

      phi_lim = -400  ! piecewise linear limiters
      r_lim = 0
      z_lim = 0

!------------------------------------------------------------------------------

      sigma_kappa = bigno
      target_kappa = 0
      Target_Iota = 0
      target_iota_p = 0
      Target_Iota_min = 0
      sigma_iota_pmax = bigno
      sigma_iota_pmin = bigno
      sigma_extcur = bigno
      target_extcur = 0
      sigma_oh = bigno
      oh_coefs = 0
      bal_theta0 = 0; bal_zeta0 = 0

      READ (iunit, nml=optimum, iostat=istat)

      END SUBROUTINE read_optimum_namelist

      END MODULE optim_params
