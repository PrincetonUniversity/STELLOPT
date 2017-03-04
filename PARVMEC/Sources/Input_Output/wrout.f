      SUBROUTINE wrout(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu, 
     1                 rzl_array, gc_array, ier_flag, lwrite
#ifdef _ANIMEC
     2                ,tau_an, sigma_an, ppar, pperp, onembc, pbprim,
     3                 ppprim, densit
#endif
     4                 )
! ... from SPH 2009-10-05; changes for modB sine-harmonics included
      USE vmec_main, p5 => cp5, two => c2p0
      USE vmec_input, ONLY: ns_array, ftol_array, lwouttxt
      USE vmec_params
      USE vmercier
      USE vmec_persistent
      USE vsvd
      USE vspline
      USE xstuff
      USE vmec_io
      USE realspace, ONLY: phip, chip, gsqrta=>z1, z1=>z1
      USE totzsp_mod
      USE vforces, ONLY: bsupua=>brmn_e, bsupva=>czmn_o, bsqa=>bzmn_e, 
     1                   bsubsa=>armn_e, bsubua=>azmn_e, bsubva=>armn_o 
      USE vacmod, ONLY: potvac, mnpd,                                   !added for diagno, J.Geiger
     1                  bsubu_sur, bsubv_sur, bsupu_sur, bsupv_sur      !MRC  10-15-15
!     USE parallel_include_module 
#ifdef _HBANGLE
      USE angle_constraints, ONLY: getrz
#endif
!#undef NETCDF 
#ifdef NETCDF      
      USE ezcdf
      USE read_wout_mod, ONLY: vn_version, vn_extension, vn_mgrid,
     1  vn_magen, vn_therm, vn_gam, vn_maxr, vn_minr, vn_maxz, vn_fp,
     2  vn_radnod, vn_polmod, vn_tormod, vn_maxmod, vn_maxit, vn_actit,
     3  vn_asym, vn_recon, vn_free, vn_error, vn_aspect, vn_beta, 
     4  vn_pbeta, vn_tbeta, vn_abeta, vn_b0, vn_rbt0, vn_maxmod_nyq,
     5  vn_rbt1, vn_sgs, vn_lar, vn_modB, vn_ctor, vn_amin, vn_Rmaj, 
     6  vn_vol, vn_mse, vn_thom, vn_ac, vn_ai, vn_am, vn_rfp, 
     6  vn_pmass_type, vn_pcurr_type, vn_piota_type,
     6  vn_am_aux_s, vn_am_aux_f, vn_ac_aux_s, vn_ac_aux_f, 
     6  vn_ai_aux_s, vn_ai_aux_f, 
     6  vn_ftolv, vn_fsqr, vn_fsqz, vn_fsql,
     7  vn_pmod, vn_tmod, vn_pmod_nyq, vn_tmod_nyq,
     7  vn_racc, vn_zacs, vn_racs, vn_zacc, vn_iotaf, vn_qfact,
     8  vn_presf, vn_phi, vn_phipf, vn_jcuru, vn_jcurv, vn_iotah,
     8  vn_chi, vn_chipf, 
     9  vn_mass, vn_presh, vn_betah, vn_buco, vn_bvco, vn_vp, vn_specw, 
     A  vn_phip, vn_jdotb, vn_overr, vn_bgrv, vn_merc, vn_mshear,
     B  vn_mwell, vn_mcurr, vn_mgeo, vn_equif, vn_fsq, vn_wdot, 
     C  vn_extcur, vn_curlab, vn_rmnc, vn_zmns, vn_lmns, vn_gmnc, 
     D  vn_bmnc, vn_bsubumnc, vn_bsubvmnc, vn_bsubsmns, 
     E  vn_bsupumnc, vn_bsupvmnc, vn_rmns, vn_zmnc, vn_lmnc, vn_gmns,
     F  vn_bmns, vn_bsubumns, vn_bsubvmns, vn_bsubsmnc, vn_bsupumns, 
     G  vn_bsupvmns, vn_rbc, vn_zbs, vn_rbs, vn_zbc,
     H  ln_version, ln_extension, ln_mgrid,
     &  vn_bsubumnc_sur, vn_bsubvmnc_sur,     !MRK 10-15-15
     &  vn_bsupumnc_sur, vn_bsupvmnc_sur,
     &  vn_bsubumns_sur, vn_bsubvmns_sur,
     &  vn_bsupumns_sur, vn_bsupvmns_sur,
     1  ln_magen, ln_therm, ln_gam, ln_maxr, ln_minr, ln_maxz, ln_fp,
     2  ln_radnod, ln_polmod, ln_tormod, ln_maxmod, ln_maxit, ln_actit,
     3  ln_asym, ln_recon, ln_free, ln_error, ln_aspect, ln_beta, 
     4  ln_pbeta, ln_tbeta, ln_abeta, ln_b0, ln_rbt0, ln_maxmod_nyq,
     5  ln_rbt1, ln_sgs, ln_lar, ln_modB, ln_ctor, ln_amin, ln_Rmaj, 
     6  ln_mse, ln_thom, ln_flp, ln_nobd, ln_nbset, ln_next, ln_nbfld,
     7  ln_pmod, ln_tmod, ln_pmod_nyq, ln_tmod_nyq, ln_racc, ln_zacs, 
     7  ln_racs, ln_zacc, ln_iotaf, ln_qfact, ln_am, ln_ac, ln_ai,
     7  ln_pmass_type, ln_pcurr_type, ln_piota_type,
     7  ln_am_aux_s, ln_am_aux_f, ln_ac_aux_s, ln_ac_aux_f, 
     7  ln_ai_aux_s, ln_ai_aux_f, ln_chi, ln_chipf, 
     8  ln_presf, ln_phi, ln_phipf, ln_jcuru, ln_jcurv, ln_iotah,
     9  ln_mass, ln_presh, ln_betah, ln_buco, ln_bvco, ln_vp, ln_specw, 
     A  ln_vol, ln_phip, ln_jdotb, ln_bgrv, ln_merc, ln_mshear,
     B  ln_mwell, ln_mcurr, ln_mgeo, ln_equif, ln_fsq, ln_wdot, 
     C  ln_extcur, ln_curlab, ln_rmnc, ln_zmns, ln_lmns, ln_gmnc, 
     D  ln_bmnc, ln_bsubumnc, ln_bsubvmnc, ln_bsubsmns, 
     E  ln_bsupumnc, ln_bsupvmnc, ln_rmns, ln_zmnc, ln_lmnc, ln_gmns,
     F  ln_bmns, ln_bsubumns, ln_bsubvmns, ln_bsubsmnc, ln_bsupumns, 
     G  ln_bsupvmns, ln_rbc, ln_zbs, ln_rbs, ln_zbc,
     &  ln_bsubumnc_sur, ln_bsubvmnc_sur,          !MRC 10-15-15
     &  ln_bsupumnc_sur, ln_bsupvmnc_sur,
     &  ln_bsubumns_sur, ln_bsubvmns_sur,
     &  ln_bsupumns_sur, ln_bsupvmns_sur

!------------------DEC$ ELSE !to use safe_open_mod in any case (J.Geiger)
#endif
      USE safe_open_mod
      USE mgrid_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ier_flag
      REAL(rprec), DIMENSION(mnmax,ns,3*MAX(ntmax/2,1)),           !reverse ns, mnmax for backwards compatibility
     1   INTENT(inout), TARGET :: rzl_array, gc_array
      REAL(rprec), DIMENSION(ns,nznt), INTENT(inout) ::
     1   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu
#ifdef _ANIMEC
     2  ,tau_an, ppar, pperp, onembc, sigma_an
      REAL(rprec), DIMENSION(ns,nznt), INTENT(out) ::
     1   densit, pbprim, ppprim
#endif
      REAL(rprec) :: qfact(ns)
      LOGICAL :: lwrite
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      LOGICAL :: lnyquist = .FALSE.                               !=false, eliminate nyquist stuff
#ifdef NETCDF
      CHARACTER(LEN=*), PARAMETER, DIMENSION(1) ::
     1             r1dim = (/'radius'/), mn1dim = (/'mn_mode'/),
     2             mn2dim = (/'mn_mode_nyq'/),
     3             currg = (/'ext_current'/),
     4             currl = (/'current_label'/)
      CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: 
     1             r2dim = (/'mn_mode','radius '/),
     1             r3dim = (/'mn_mode_nyq','radius     '/)
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j, js, jlk, mn, lk, iasym, ireconstruct,
     1           m, n, k, iwout0, n1, nwout, istat, i, indx1(1),
     2           mnmax_nyq0, mnyq0, nnyq0, nwout2   ! nwout2 by J.Geiger
     3          ,isgn, js2, nfort      !for diagno 1.5
      REAL(rprec) :: dmult, tcosi, tsini, vversion, sgn, tmult,
     1               presfactor, ftolx1, d_bsupumn, d_bsupvmn   ! diagno 1.5
#ifdef _ANIMEC
     2              ,hotdam, omtbc, optbc, pdh, pmh, pde, pme, eps
#endif
      REAL(rprec), POINTER, DIMENSION(:,:) :: rmnc, rmns, zmns, 
     1   zmnc, lmns, lmnc
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: 
     1   gmnc, bmnc, gmns, bmns, 
     2   bsubumnc, bsubvmnc, bsubsmns, bsubumns, bsubvmns, bsubsmnc
#ifdef _ANIMEC
     3  ,sigmnc  , taumnc  , pparmnc , ppermnc , pbprmnc , ppprmnc ,
     4   hotdmnc , hotdmns ,
     5   sigmns  , taumns  , pparmns , ppermns , pbprmns , ppprmns
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sigma_ana, tau_ana,
     1                      ppara, pperpa, pbprima, ppprima, densita     
#endif
      REAL(rprec), DIMENSION(mnmax) :: rmnc1, zmns1, lmns1,
     1   rmns1, zmnc1, lmnc1, bmodmn, bmodmn1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gmn, bmn,
     1   bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
#ifdef _ANIMEC
     2  ,sigmn  , taumn  , pparmn , ppermn , pbprmn , ppprmn , hotdmn
#endif
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubumnc_sur  !MRK 10-15-15
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubvmnc_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsupumnc_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsupvmnc_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubumns_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubvmns_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsupumns_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsupvmns_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubua_sur, bsubva_sur
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsupua_sur, bsupva_sur

      CHARACTER(LEN=120) :: wout_file, wout2_file         ! wout2_file by J.Geiger
      CHARACTER(LEN=120) :: fort_file   ! fort_file for diagno 1.5
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xfinal
      REAL(rprec), DIMENSION(:), POINTER ::   xm_nyq0, xn_nyq0
!     ELIMINATE THESE EVENTUALLY
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: 
     1   bsupumnc, bsupumns, bsupvmnc, bsupvmns
      LOGICAL :: lcurr
      INTEGER :: nmin0     ! J Geiger:   Added for diagno-file

!-----------------------------------------------
!
!  Pointer assignments for storage arrays
!
      n1 = MAX(1,ntmax/2)
      rmnc => rzl_array(:,:,1)            !!store COS(mu-nv) components
      zmns => rzl_array(:,:,1+n1)         !!store SIN(mu-nv)
      lmns => rzl_array(:,:,1+2*n1)       !!store SIN(mu-nv)

      IF (lasym) THEN
         rmns => gc_array(:,:,1)            !!store SIN(mu-nv)
         zmnc => gc_array(:,:,1+n1)         !!store COS(mu-nv)
         lmnc => gc_array(:,:,1+2*n1)       !!store COS(mu-nv)
      END IF

!
!     THIS SUBROUTINE CREATES THE FILE WOUT.IT CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (full), LMN (half_mesh - CONVERTED FROM
!     INTERNAL full REPRESENTATION), AS WELL AS COEFFICIENTS (ON NYQ MESH) FOR COMPUTED
!     QUANTITIES:
!
!     BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF); BSUBS (FULL-CONVERTED IN JXBFORCE)
!
      IF (lnyquist) THEN
         mnmax_nyq0 = mnmax_nyq
         mnyq0 = mnyq
         nnyq0 = nnyq
         xm_nyq0 => xm_nyq; xn_nyq0 => xn_nyq
      ELSE
         mnmax_nyq0 = mnmax
         mnyq0 = mpol1
         nnyq0 = ntor
         xm_nyq0 => xm; xn_nyq0 => xn
      END IF

      ALLOCATE (gmn(mnmax_nyq0), bmn(mnmax_nyq0),
     1   bsubumn(mnmax_nyq0), bsubvmn(mnmax_nyq0), bsubsmn(mnmax_nyq0), 
     2   bsupumn(mnmax_nyq0), bsupvmn(mnmax_nyq0), 
#ifdef _ANIMEC
     3   sigmn(mnmax_nyq0)  ,
     4   taumn(mnmax_nyq0)  , pparmn(mnmax_nyq0) , ppermn(mnmax_nyq0) ,
     5   pbprmn(mnmax_nyq0) , ppprmn(mnmax_nyq0) , hotdmn(mnmax_nyq0) ,
#endif
     6   stat=istat)

      IF (lfreeb) THEN        !MRK 10-15-15
         ALLOCATE (bsubua_sur(nzeta*ntheta2), bsubva_sur(nzeta*ntheta2))
         ALLOCATE (bsupua_sur(nzeta*ntheta2), bsupva_sur(nzeta*ntheta2))

         ALLOCATE (bsubumnc_sur(mnmax_nyq0), bsubvmnc_sur(mnmax_nyq0))
         ALLOCATE (bsupumnc_sur(mnmax_nyq0), bsupvmnc_sur(mnmax_nyq0))
         IF (lasym) THEN
            ALLOCATE (bsubumns_sur(mnmax_nyq0),                                &
     &                bsubvmns_sur(mnmax_nyq0))
            ALLOCATE (bsupumns_sur(mnmax_nyq0),                                &
     &                bsupvmns_sur(mnmax_nyq0))
         END IF
      END IF

      ALLOCATE (gmnc(mnmax_nyq0,ns), bmnc(mnmax_nyq0,ns),
     1          bsubumnc(mnmax_nyq0,ns), bsubvmnc(mnmax_nyq0,ns),
     2          bsubsmns(mnmax_nyq0,ns), bsupumnc(mnmax_nyq0,ns),
     3          bsupvmnc(mnmax_nyq0,ns), 
#ifdef _ANIMEC
     4          sigmnc(mnmax_nyq0,ns)  ,
     5          taumnc(mnmax_nyq0,ns)  , pparmnc(mnmax_nyq0,ns) ,
     6          ppermnc(mnmax_nyq0,ns) , pbprmnc(mnmax_nyq0,ns) ,
     7          ppprmnc(mnmax_nyq0,ns) , hotdmnc(mnmax_nyq0,ns) ,
#endif
     8          stat=istat)
      IF (lasym) THEN
      ALLOCATE (gmns(mnmax_nyq0,ns), bmns(mnmax_nyq0,ns),
     1          bsubumns(mnmax_nyq0,ns), bsubvmns(mnmax_nyq0,ns),
     2          bsubsmnc(mnmax_nyq0,ns), bsupumns(mnmax_nyq0,ns),
     3          bsupvmns(mnmax_nyq0,ns), 
#ifdef _ANIMEC
     4          sigmns(mnmax_nyq0,ns)  ,
     5          taumns(mnmax_nyq0,ns)  , pparmns(mnmax_nyq0,ns) ,
     6          ppermns(mnmax_nyq0,ns) , pbprmns(mnmax_nyq0,ns) ,
     7          ppprmns(mnmax_nyq0,ns) , hotdmns(mnmax_nyq0,ns) ,
#endif
     8          stat=istat)
#ifdef _ANIMEC
      ALLOCATE (sigma_ana(ns,nznt) ,tau_ana(ns,nznt) ,densita(ns,nznt),
     1          ppara(ns,nznt)     ,pperpa(ns,nznt)  ,pbprima(ns,nznt),
     2          ppprima(ns,nznt), stat=istat)
#endif
      END IF
      IF (istat .ne. 0) STOP 'Error allocating arrays in VMEC WROUT'

!      IF (nextcur .eq. 0) THEN
!         DO j = SIZE(extcur), 1, -1
!           IF (extcur(j) .ne. zero) THEN
!               nextcur = j
!               EXIT
!            END IF
!         END DO
!      END IF

! ftol info evaluated here!
      indx1=MAXLOC(ns_array)
      ftolx1=ftol_array(indx1(1))

!     NYQUIST FREQUENCY REQUIRES FACTOR OF 1/2
      IF (lnyquist) THEN
         IF (mnyq .ne. 0) cosmui(:,mnyq) = p5*cosmui(:,mnyq)
         IF (nnyq .ne. 0) cosnv (:,nnyq) = p5*cosnv (:,nnyq)
      END IF
      wout_file = version_
      READ (wout_file, *) vversion

#ifdef NETCDF
      wout_file = 'wout_' // TRIM(input_extension) // '.nc'
      CALL cdf_open(nwout,wout_file,'w',iwout0)
      IF (iwout0 .ne. 0) STOP 'Error opening wout.nc file VMEC WROUT'

      IF (.not. lrecon) THEN
         itse = 0
         imse2 = 0
      END IF
      
!================================
! Define Variables
!================================
!  Scalars 
      CALL cdf_define(nwout, vn_version, vversion)
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
      CALL cdf_define(nwout, vn_maxmod_nyq, mnmax_nyq0)
      CALL cdf_define(nwout, vn_maxit, iter2)
      CALL cdf_define(nwout, vn_actit, itfsq)
      CALL cdf_define(nwout, vn_asym, lasym)
      CALL cdf_define(nwout, vn_recon, lrecon)
      CALL cdf_define(nwout, vn_free, lfreeb)
      CALL cdf_define(nwout, vn_rfp, lrfp)
      CALL cdf_define(nwout, vn_error, ier_flag)
      CALL cdf_define(nwout, vn_aspect, aspect)
      CALL cdf_define(nwout, vn_beta, betatot)
      CALL cdf_define(nwout, vn_pbeta, betapol)
      CALL cdf_define(nwout, vn_tbeta, betator)
      CALL cdf_define(nwout, vn_abeta, betaxis)
      CALL cdf_define(nwout, vn_b0, b0)
      CALL cdf_define(nwout, vn_rbt0, rbtor0)
      CALL cdf_define(nwout, vn_rbt1, rbtor)
      CALL cdf_define(nwout, vn_sgs, NINT(signgs))
      CALL cdf_define(nwout, vn_lar, IonLarmor)
      CALL cdf_define(nwout, vn_modB, volAvgB)
      CALL cdf_define(nwout, vn_ctor, ctor)
      CALL cdf_define(nwout, vn_amin, Aminor_p)
      CALL cdf_define(nwout, vn_Rmaj, Rmajor_p)
      CALL cdf_define(nwout, vn_vol, volume_p)
      CALL cdf_define(nwout, vn_ftolv, ftolx1)
      CALL cdf_define(nwout, vn_fsql, fsql)
      CALL cdf_define(nwout, vn_fsqr, fsqr)
      CALL cdf_define(nwout, vn_fsqz, fsqz)

      IF (lrecon) THEN
         CALL cdf_define(nwout, vn_mse, imse2)
         CALL cdf_define(nwout, vn_thom, itse)
      END IF

      CALL cdf_define(nwout, vn_nextcur, nextcur)
      CALL cdf_define(nwout, vn_extcur,
     1        extcur(1:nextcur), dimname=currg)
      CALL cdf_define(nwout, vn_mgmode, mgrid_mode)
      IF (lfreeb) THEN
         CALL cdf_define(nwout, vn_flp, nobser)
         CALL cdf_define(nwout, vn_nobd, nobd)
         CALL cdf_define(nwout, vn_nbset, nbsets)
         IF (nbsets .gt. 0) 
     1      CALL cdf_define(nwout,vn_nbfld,nbfld(1:nbsets))
      END IF

      IF (.not.lwrite) GO TO 800

! 1D Arrays

      CALL cdf_define(nwout, vn_pmod, xm, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_pmod, ln_pmod)
      CALL cdf_define(nwout, vn_tmod, xn, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_tmod, ln_tmod)
      CALL cdf_define(nwout, vn_pmod_nyq, xm_nyq0, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_pmod_nyq, ln_pmod_nyq)
      CALL cdf_define(nwout, vn_tmod_nyq, xn_nyq0, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_tmod_nyq, ln_tmod_nyq)

      CALL cdf_define(nwout, vn_racc, raxis_cc(0:ntor), 
     1                dimname=(/'n_tor'/))
      CALL cdf_setatt(nwout, vn_racc, ln_racc)
      CALL cdf_define(nwout, vn_zacs, zaxis_cs(0:ntor), 
     1                dimname=(/'n_tor'/))
      CALL cdf_setatt(nwout, vn_zacs, ln_zacs)
      IF (lasym) THEN
         CALL cdf_define(nwout, vn_racs, raxis_cs(0:ntor), 
     1                dimname=(/'n_tor'/))
         CALL cdf_setatt(nwout, vn_racs, ln_racs)
         CALL cdf_define(nwout, vn_zacc, zaxis_cc(0:ntor), 
     1                dimname=(/'n_tor'/))
         CALL cdf_setatt(nwout, vn_zacc, ln_zacc)
      END IF

      j = SIZE(am)-1
      CALL cdf_define(nwout, vn_am, am(0:j),
     1                dimname=(/'preset'/))
      j = SIZE(ac)-1
      CALL cdf_define(nwout, vn_ac, ac(0:j),
     1                dimname=(/'preset'/))
      j = SIZE(ai)-1
      CALL cdf_define(nwout, vn_ai, ai(0:j),
     1                dimname=(/'preset'/))
     
      j = SIZE(am_aux_s)
      CALL cdf_define(nwout, vn_am_aux_s, am_aux_s(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(am_aux_f)
      CALL cdf_define(nwout, vn_am_aux_f, am_aux_f(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ai_aux_s)
      CALL cdf_define(nwout, vn_ai_aux_s, ai_aux_s(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ai_aux_f)
      CALL cdf_define(nwout, vn_ai_aux_f, ai_aux_f(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ac_aux_s)
      CALL cdf_define(nwout, vn_ac_aux_s, ac_aux_s(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ac_aux_f)
      CALL cdf_define(nwout, vn_ac_aux_f, ac_aux_f(1:j),
     1                dimname=(/'ndfmax'/))


      CALL cdf_define(nwout, vn_iotaf, iotaf(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_iotaf)

      qfact=HUGE(qfact)
      WHERE (iotaf .NE. zero) qfact=one/iotaf

      CALL cdf_define(nwout, vn_qfact, qfact(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_qfact)
      CALL cdf_define(nwout, vn_presf, presf, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presf, ln_presf, units='Pa')
      CALL cdf_define(nwout, vn_phi, phi, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phi, ln_phi, units='wb')
      CALL cdf_define(nwout, vn_phipf, 
     1                phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phipf, ln_phipf)
      CALL cdf_define(nwout, vn_chi, chi, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_chi, ln_chi, units='wb')
      CALL cdf_define(nwout, vn_chipf, 
     1                phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_chipf, ln_chipf)
      CALL cdf_define(nwout, vn_jcuru, 
     1                jcuru, dimname=r1dim)
      CALL cdf_define(nwout, vn_jcurv, 
     1                jcurv, dimname=r1dim)
 
      CALL cdf_define(nwout, vn_iotah, iotas(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotah, ln_iotah)
      CALL cdf_define(nwout, vn_mass, mass, 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_mass, ln_mass)
      CALL cdf_define(nwout, vn_presh, pres(1:ns), 
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presh, ln_presh, units='Pa')
      CALL cdf_define(nwout, vn_betah, beta_vol, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_buco, buco, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bvco, bvco, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_vp, vp(1:ns), 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_specw, specw, 
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_phip, 
     1                phips(1:ns), dimname=r1dim)
      CALL cdf_define(nwout, vn_overr, 
     2                overr(1:ns), dimname=r1dim)

      CALL cdf_define(nwout, vn_jdotb, jdotb,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bgrv, bdotgradv,
     1                dimname=r1dim)

      CALL cdf_define(nwout, vn_merc, Dmerc,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mshear, Dshear,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mwell, Dwell,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mcurr, Dcurr,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mgeo,
     1                Dgeod, dimname=r1dim)
      CALL cdf_define(nwout, vn_equif,
     1                equif, dimname=r1dim)

      CALL cdf_define(nwout, vn_fsq, fsqt(1:nstore_seq),
     1                dimname=(/'time'/))
      CALL cdf_define(nwout, vn_wdot, wdot(1:nstore_seq),
     1                dimname=(/'time'/))

      IF (lfreeb .and. nextcur.gt.0 .and. ALLOCATED(curlabel)) THEN
         CALL cdf_define(nwout, vn_curlab,
     1        curlabel(1:nextcur), dimname=currl)
      ENDIF

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

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM - CAN COMPUTE FROM GSQRT
      CALL cdf_define(nwout, vn_bsupumnc, bsupumnc, dimname=r3dim)
      CALL cdf_define(nwout, vn_bsupvmnc, bsupvmnc, dimname=r3dim)
!     IF (lfreeb) THEN
!         CALL cdf_define(nwout, vn_rbc, rbc, 
!    1                dimname=(/'n_mode','m_mode'/))
!         CALL cdf_setatt(nwout, vn_rbc, ln_rbc, units='m')
!         CALL cdf_define(nwout, vn_zbs, zbs, 
!    1                dimname=(/'n_mode','m_mode'/))
!         CALL cdf_setatt(nwout, vn_zbs, ln_zbs, units='m')
!        IF (lasym) THEN
!           CALL cdf_define(nwout, vn_rbs, rbs, 
!    1                dimname=(/'n_mode','m_mode'/))
!           CALL cdf_define(nwout, vn_zbc, zbc, 
!    1                dimname=(/'n_mode','m_mode'/))
!        END IF
!     END IF

      IF (.NOT. lasym) GO TO 800

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

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
      CALL cdf_define(nwout, vn_bsupumns, bsupumns, dimname=r3dim)
      CALL cdf_define(nwout, vn_bsupvmns, bsupvmns, dimname=r3dim)

 800  CONTINUE

!================================
! Write Variables
!================================

! Scalars
      CALL cdf_write(nwout, vn_version, vversion)
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
      CALL cdf_write(nwout, vn_maxmod_nyq, mnmax_nyq0)
      CALL cdf_write(nwout, vn_maxit, iter2)
      CALL cdf_write(nwout, vn_actit, itfsq)
      CALL cdf_write(nwout, vn_asym, lasym)
      CALL cdf_write(nwout, vn_recon, lrecon)
      CALL cdf_write(nwout, vn_free, lfreeb)
      CALL cdf_write(nwout, vn_rfp, lrfp)
      CALL cdf_write(nwout, vn_error, ier_flag)
!
      CALL cdf_write(nwout, vn_aspect, aspect)
      CALL cdf_write(nwout, vn_beta, betatot)
      CALL cdf_write(nwout, vn_pbeta, betapol)
      CALL cdf_write(nwout, vn_tbeta, betator)
      CALL cdf_write(nwout, vn_abeta, betaxis)
      CALL cdf_write(nwout, vn_b0, b0)
      CALL cdf_write(nwout, vn_rbt0, rbtor0)
      CALL cdf_write(nwout, vn_rbt1, rbtor)
      CALL cdf_write(nwout, vn_sgs, NINT(signgs))
      CALL cdf_write(nwout, vn_lar, IonLarmor)
      CALL cdf_write(nwout, vn_modB, volAvgB)
      CALL cdf_write(nwout, vn_ctor, ctor/mu0)
      CALL cdf_write(nwout, vn_amin, Aminor_p)
      CALL cdf_write(nwout, vn_rmaj, Rmajor_p)
      CALL cdf_write(nwout, vn_vol, volume_p)
      CALL cdf_write(nwout, vn_ftolv, ftolx1)
      CALL cdf_write(nwout, vn_fsql, fsql)
      CALL cdf_write(nwout, vn_fsqr, fsqr)
      CALL cdf_write(nwout, vn_fsqz, fsqz)

      IF (lrecon) THEN
         CALL cdf_write(nwout, vn_mse, imse2-1)
         CALL cdf_write(nwout, vn_thom, itse)
      END IF

      CALL cdf_write(nwout, vn_nextcur, nextcur)
      IF (nextcur .gt. 0) THEN
         CALL cdf_write(nwout, vn_extcur, extcur(1:nextcur))
         CALL cdf_write(nwout, vn_mgmode, mgrid_mode)
      ENDIF
      IF (lfreeb) THEN
         CALL cdf_write(nwout, vn_flp, nobser)
         CALL cdf_write(nwout, vn_nobd, nobd)
         CALL cdf_write(nwout, vn_nbset, nbsets)
         IF (nextcur.gt.0 .and. ALLOCATED(curlabel))
     1   CALL cdf_write(nwout, vn_curlab, curlabel(1:nextcur))
      END IF

! 1D Arrays
      IF (nbsets .gt. 0) CALL cdf_write(nwout,vn_nbfld,nbfld(1:nbsets))

      IF (.not.lwrite) GO TO 940   !change 970 to 940, J.Geiger
                                   !skip the next cdf_write-calls
      CALL cdf_write(nwout, vn_pmod, xm)
      CALL cdf_write(nwout, vn_tmod, xn)
      CALL cdf_write(nwout, vn_pmod_nyq, xm_nyq0)
      CALL cdf_write(nwout, vn_tmod_nyq, xn_nyq0)

940   CONTINUE   ! before closing, write the initial part of the wouttxt-file
#endif
      IF(lwouttxt) THEN
!      IF(.TRUE.) THEN !SKS: Dec 5, 2014
         wout2_file = 'wout_'//TRIM(input_extension) // '.txt'
         nwout2 = nwout0
         CALL safe_open(nwout2, iwout0, wout2_file,
     1                 'replace', 'formatted')
         IF (iwout0 .ne. 0) STOP 'Error opening WOUT.txt file in WROUT'

         IF (lasym) THEN
            iasym = 1                                  ! asymmetric mode
         ELSE
            iasym = 0
         END IF
         IF (lrecon) THEN
            ireconstruct = 1
         ELSE
            itse = 0
            imse2 = 0
            ireconstruct = 0
         END IF

!
!     Insert version information into wout file. This will be parsed in
!     read_wout_file to return the real value version_ to check the version number.
!
#if defined(_ANIMEC)
         WRITE (nwout2, '(a15,a,a)') 'VMEC VERSION = ', version_,
     1                               '_ANIMEC'
#elif defined(_FLOW)
         WRITE (nwout2, '(a15,a,a)') 'VMEC VERSION = ', version_,'_FLOW'
#else
         WRITE (nwout2, '(a15,a)') 'VMEC VERSION = ', version_
#endif

#ifdef _ANIMEC
         WRITE (nwout2, *) wb, wpar, gamma, pfac,
#else
         WRITE (nwout2, *) wb, wp, gamma, pfac,
#endif
     1    rmax_surf, rmin_surf, zmax_surf

         WRITE (nwout2, *) nfp, ns, mpol, ntor, mnmax, mnmax_nyq0,
     1     itfsq, iter2, iasym, ireconstruct, ier_flag

         WRITE (nwout2, *) imse2 - 1, itse, nbsets, nobd, nextcur,
     1     nstore_seq
         IF (nbsets .gt. 0) WRITE (nwout2, *) (nbfld(i),i=1,nbsets)
         WRITE (nwout2, '(a)') mgrid_file

         IF (.not. lwrite) GOTO 950  ! J Geiger: At this point only closing of
                                     !           txt- and nc-file is needed
         pres(1) = pres(2)
      END IF
!---------------------DEC$ ENDIF

      IF (.not. lwrite) GOTO 970   ! J Geiger: in case lwouttxt is not true
                                   !           jump to close nc-file
      ALLOCATE (xfinal(neqs2), stat=js)
      IF (js .ne. 0) STOP 'Allocation error for xfinal in WROUT!'
      xfinal = xc
#ifdef _HBANGLE
      CALL getrz(xfinal)
#else
!
!     MUST CONVERT m=1 MODES... FROM INTERNAL TO PHYSICAL FORM
!     Extrapolation of m=0 Lambda (cs) modes, which are not evolved at j=1, done in CONVERT
!

      IF (lthreed) CALL convert_sym  (xfinal(1+mns*(rss-1)), 
     1                                xfinal(1+irzloff+mns*(zcs-1)))
      IF (lasym)   CALL convert_asym (xfinal(1+mns*(rsc-1)), 
     1                                xfinal(1+irzloff+mns*(zcc-1)))
#endif
!
!     CONVERT TO rmnc, zmns, lmns, etc EXTERNAL representation (without internal mscale, nscale)
!     IF B^v ~ phip + lamu, MUST DIVIDE BY phipf(js) below to maintain old-style format
!     THIS COULD BE A PROBLEM FOR RFP WHERE PHIPF->0 INSIDE THE PLASMA!
!
      RADIUS1: DO js = 1, ns

         CALL convert (rmnc1, zmns1, lmns1, rmns1, zmnc1, lmnc1, 
     1                 xfinal, js)

         rmnc(:,js) = rmnc1(:)
         zmns(:,js) = zmns1(:)
         lmns(:,js) = (lmns1(:)/phipf(js)) * lamscale
         IF (lasym) THEN
            rmns(:,js) = rmns1(:)
            zmnc(:,js) = zmnc1(:)
            lmnc(:,js) = (lmnc1(:)/phipf(js)) * lamscale
         END IF

      END DO RADIUS1

      DEALLOCATE (xfinal)

!
!     INTERPOLATE LAMBDA ONTO HALF-MESH FOR BACKWARDS CONSISTENCY WITH EARLIER VERSIONS OF VMEC
!     AND SMOOTHS POSSIBLE UNPHYSICAL "WIGGLE" ON RADIAL MESH
!

      WHERE (NINT(xm) .le. 1) lmns(:,1) = lmns(:,2)
      DO js = ns,2,-1
         WHERE (MOD(NINT(xm),2) .eq. 0) 
            lmns(:,js) = p5*(lmns(:,js) + lmns(:,js-1))
         ELSEWHERE
            lmns(:,js) = p5*(sm(js)*lmns(:,js) + sp(js-1)*lmns(:,js-1))
         END WHERE
      END DO

      lmns(:,1) = 0  
      raxis_cc(0:ntor) = rmnc(1:ntor+1,1)
      zaxis_cs(0:ntor) = zmns(1:ntor+1,1)

      IF (.not.lasym) GOTO 900

      WHERE (NINT(xm) .le. 1) lmnc(:,1) = lmnc(:,2)
      DO js = ns,2,-1
         WHERE (MOD(NINT(xm),2) .eq. 0) 
            lmnc(:,js) = p5*(lmnc(:,js) + lmnc(:,js-1))
         ELSEWHERE
            lmnc(:,js) = p5*(sm(js)*lmnc(:,js) + sp(js-1)*lmnc(:,js-1))
         END WHERE
      END DO

      lmnc(:,1) = 0;   
      raxis_cs(0:ntor) = rmns(1:ntor+1,1)
      zaxis_cc(0:ntor) = zmnc(1:ntor+1,1)

 900  CONTINUE

#ifdef _ANIMEC
!... CALCULATE RADIAL DERIVATIVES OF HOT PARTICLE PRESSURE TERMS
!... STORE IN ARRAYS pm AND pd PREVIOUSLY USED IN PRESSURE AND EQFOR
      eps = EPSILON(eps)
      DO js=2,ns-1
         pd(js) = ohs * (pres(js+1) * phot(js+1) - pres(js) * phot(js))
         pmap(js) = ohs * (tpotb(js+1) - tpotb(js))
      END DO
!... INTERPOLATE (EXTRAPOLATE) TO HALF INTEGER MESH
      pdh = c1p5 * pd(2) - p5 * pd(3)
      pmh = c1p5 * pmap(2) - p5 * pmap(3)
      pde = c1p5 * pd(ns-1) - p5 * pd(ns-2)
      pme = c1p5 * pmap(ns-1) - p5 * pmap(ns-2)
      DO js=ns-2,2,-1
         pd(js+1) = p5*(pd(js+1) + pd(js)) / (pres(js+1)*phot(js+1)+eps)
         pmap(js+1) = p5 * (pmap(js+1) + pmap(js))
      END DO
      pd(2)  = pdh / (pres(2)*phot(2)+eps)
      pd(ns) = pde / (pres(ns)*phot(ns)+eps)
      pmap(2)  = pmh
      pmap(ns) = pme
!ALTERNATE EXTRAPOLATION
      pd(2) = 2*pd(3) - pd(4)
      pd(ns) = 2*pd(ns-1) - pd(ns-2) 

!CALCULATE HOT PARTICLE PARALLEL AND PERPENDICULAR PRESSURE GRADIENT; DENSITY
      DO 20 js = 2, ns 
        hotdam = pres(js) * phot(js) / SQRT(tpotb(js)+eps)
        DO 10 lk = 1, nznt  
!  
           omtbc = one - tpotb(js) * onembc(js,lk)
           optbc = one + tpotb(js) * onembc(js,lk)
        IF (onembc(js,lk) <= zero) THEN
           densit(js,lk)= (ppar(js,lk) - pres(js))*hotdam / 
     &                    (pres(js)*phot(js)+eps)
           pbprim(js,lk) =  (ppar(js,lk) -pres(js)) *
     &             (pd(js) + onembc(js,lk) * pmap(js) / (omtbc+eps))
           ppprim(js,lk) =  (pperp(js,lk)-pres(js)) *
     &            (pd(js) + optbc * pmap(js) / (omtbc * tpotb(js)+eps))
        ELSE
          densit(js,lk) = hotdam * (one - onembc(js,lk)) *
     &  (optbc - 2*(tpotb(js)*onembc(js,lk))**c1p5) / (omtbc*optbc+eps)
          pbprim(js,lk) =  (ppar(js,lk) -pres(js)) * pd(js) +
     &    ( 2 * tpotb(js) * onembc(js,lk)**2 * (ppar(js,lk)-pres(js))
     &   + pres(js)*phot(js)*(one-onembc(js,lk))*onembc(js,lk)*(one -5
     & *(tpotb(js)*onembc(js,lk))**c1p5))* pmap(js) / (omtbc*optbc+eps)
          ppprim(js,lk) =  (pperp(js,lk)-pres(js)) * pd(js) +
     & ((pperp(js,lk)-pres(js))*(one+3*(tpotb(js)*onembc(js,lk))**2) /
     &  (tpotb(js)+eps)+ pres(js)*phot(js)*tpotb(js)
     &   *(one-onembc(js,lk))**2
     &   * onembc(js,lk)*(two*optbc-sqrt(tpotb(js)*onembc(js,lk))*(7.5
     &   - 3.5_dp*(tpotb(js)*onembc(js,lk))**2))/(omtbc*optbc+eps))
     &   * pmap(js)/ (omtbc * optbc + eps)
        END IF  
   10   END DO
   20  END DO
#endif
!SPH100209: COMPUTE |B| = SQRT(|B|**2) and store in bsq, bsqa
      DO js = 2, ns
         bsq(js,:nznt) = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
      END DO

      tmult = p5/r0scale**2
!SPH: FIXED THIS 03-05-07 TO CALL symmetrization routine
      IF (lasym) THEN
!Changed integration norm in fixaray, SPH012314
         tmult = 2*tmult
         bsubs(1,:) = 0
         CALL symoutput (bsq,   gsqrt,  bsubu,  bsubv,  bsupu,  
     1                   bsupv,  bsubs, 
#ifdef _ANIMEC
     2                   ppar   , pperp  , densit ,
     3                   sigma_an , tau_an , pbprim , ppprim ,
#endif
     4                   bsqa,  gsqrta, bsubua, bsubva, bsupua,
     5                   bsupva, bsubsa
#ifdef _ANIMEC
     6                  ,ppara  , pperpa , densita,
     7                   sigma_ana, tau_ana, pbprima, ppprima
#endif
     8                    )

         IF (lfreeb) THEN     !MRC  10-15-15
            CALL symoutput_sur(bsubu_sur, bsubv_sur,                           &
     &                         bsupu_sur, bsupv_sur,                           &
     &                         bsubua_sur, bsubva_sur,                         &
     &                         bsupua_sur, bsupva_sur)
         END IF
      END IF

!         DO js = 2, ns
!            WRITE (200, *) 'JS: ', js, 'BSUBU, BSUBV'
!            WRITE (200, '(1p,6e12.4)') bsubu(js,:), bsubv(js,:)
!         END DO

      RADIUS2: DO js = 2, ns
         gmn = 0
         bmn = 0
         bsubumn = 0
         bsubvmn = 0
         bsubsmn = 0
         bsupumn = 0
         bsupvmn = 0

         MN2: DO mn = 1, mnmax_nyq0
            n = NINT(xn_nyq0(mn))/nfp
            m = NINT(xm_nyq0(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1 
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))          !cos(mu - nv)
                  tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                       sgn*cosmui(j,m)*sinnv(k,n1))          !sin(mu - nv)
                  bmn(mn) = bmn(mn) + tcosi*bsq(js,lk)
                  gmn(mn) = gmn(mn) + tcosi*gsqrt(js,lk)
                  bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lk)
                  bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lk)
                  bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lk)
                  bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lk) 
                  bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lk) 
#ifdef _ANIMEC
                  pparmn(mn)  = pparmn(mn)  + tcosi*
     1                                       (ppar(js,lk)-pres(js))
                  ppermn(mn)  = ppermn(mn)  + tcosi*
     1                                       (pperp(js,lk)-pres(js))
                  sigmn(mn)   = sigmn(mn)   + tcosi*sigma_an(js,lk)
                  taumn(mn)   = taumn(mn)   + tcosi*tau_an(js,lk)
                  pbprmn(mn)  = pbprmn(mn)  + tcosi*pbprim(js,lk)
                  ppprmn(mn)  = ppprmn(mn)  + tcosi*ppprim(js,lk)
                  hotdmn(mn)  = hotdmn(mn)  + tcosi*densit(js,lk)
#endif
               END DO
            END DO
         END DO MN2

         IF (js .eq. ns/2) bmodmn = bmn(1:mnmax)
         IF (js .eq. ns) bmodmn1 = bmn(1:mnmax)
         gmnc(:,js) = gmn(:)
         bmnc(:,js) = bmn(:)
         bsubumnc(:,js) = bsubumn(:)
         bsubvmnc(:,js) = bsubvmn(:)
         bsubsmns(:,js) = bsubsmn(:)
         bsupumnc(:,js) = bsupumn(:)
         bsupvmnc(:,js) = bsupvmn(:)
#ifdef _ANIMEC
         pparmnc(:,js)  = pparmn(:)
         ppermnc(:,js)  = ppermn(:)
         sigmnc(:,js)   = sigmn(:)
         taumnc(:,js)   = taumn(:)
         pbprmnc(:,js)  = pbprmn(:)
         ppprmnc(:,js)  = ppprmn(:)
         hotdmnc(:,js)  = hotdmn(:)
#endif
      END DO RADIUS2

      IF (lfreeb) THEN    !MRC    10-15-15
         bsubumnc_sur = 0
         bsubvmnc_sur = 0
         bsupumnc_sur = 0
         bsupvmnc_sur = 0
         DO mn = 1, mnmax_nyq0
            n = NINT(xn_nyq0(mn))/nfp
            m = NINT(xm_nyq0(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))
                  bsubumnc_sur(mn) = bsubumnc_sur(mn)                          &
     &                             + tcosi*bsubu_sur(lk)
                  bsubvmnc_sur(mn) = bsubvmnc_sur(mn)                          &
     &                             + tcosi*bsubv_sur(lk)
                  bsupumnc_sur(mn) = bsupumnc_sur(mn)                          &
     &                             + tcosi*bsupu_sur(lk)
                  bsupvmnc_sur(mn) = bsupvmnc_sur(mn)                          &
     &                             + tcosi*bsupv_sur(lk)
               END DO
            END DO
         END DO
      END IF

      gmnc(:,1) = 0; bmnc(:,1) = 0; 
      bsubumnc(:,1) = 0
      bsubvmnc(:,1) = 0
      bsubsmns(:,1) = 2*bsubsmns(:,2) - bsubsmns(:,3)
      bsupumnc(:,1) = 0;  bsupvmnc(:,1) = 0

#ifdef _ANIMEC
      hotdmnc(:,1)  = 0;  pparmnc(:,1)  = 0;  ppermnc(:,1) = 0
      pbprmnc(:,1)  = 0;  ppprmnc(:,1)  = 0
      sigmnc(:,1)   = 0;  taumnc(:,1)   = 0   
#endif

      IF (.not.lasym) GO TO 200

      RADIUS3: DO js = 2, ns
         gmn = 0
         bmn = 0
         bsubumn = 0
         bsubvmn = 0
         bsubsmn = 0
         bsupumn = 0
         bsupvmn = 0
#ifdef _ANIMEC
         pparmn  = 0
         ppermn  = 0
         sigmn   = 0
         taumn   = 0
         pbprmn  = 0
         ppprmn  = 0
         hotdmn  = 0
#endif
         MN3: DO mn = 1, mnmax_nyq0
            n = NINT(xn_nyq0(mn))/nfp
            m = NINT(xm_nyq0(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            jlk = js
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))
                  tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                       sgn*cosmui(j,m)*sinnv(k,n1))
                  bmn(mn) = bmn(mn) + tsini*bsqa(jlk)
                  gmn(mn) = gmn(mn) + tsini*gsqrta(jlk,0)
                  bsubumn(mn) = bsubumn(mn) + tsini*bsubua(jlk)
                  bsubvmn(mn) = bsubvmn(mn) + tsini*bsubva(jlk)
                  bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubsa(jlk)
                  bsupumn(mn) = bsupumn(mn) + tsini*bsupua(jlk)
                  bsupvmn(mn) = bsupvmn(mn) + tsini*bsupva(jlk)
#ifdef _ANIMEC
                  pparmn(mn)  = pparmn(mn)  + tsini*ppara(js,lk)
                  ppermn(mn)  = ppermn(mn)  + tsini*pperpa(js,lk)
                  sigmn(mn)   = sigmn(mn)   + tsini*sigma_ana(js,lk)
                  taumn(mn)   = taumn(mn)   + tsini*tau_ana(js,lk)
                  pbprmn(mn)  = pbprmn(mn)  + tsini*pbprima(js,lk)
                  ppprmn(mn)  = ppprmn(mn)  + tsini*ppprima(js,lk)
                  hotdmn(mn)  = hotdmn(mn)  + tsini*densita(js,lk)
#endif
                  jlk = jlk+ns
               END DO
            END DO
         END DO MN3
   
         gmns(:,js) = gmn(:)
         bmns(:,js) = bmn(:)
         bsubumns(:,js) = bsubumn(:)
         bsubvmns(:,js) = bsubvmn(:)
         bsubsmnc(:,js) = bsubsmn(:)
         bsupumns(:,js) = bsupumn(:)
         bsupvmns(:,js) = bsupvmn(:)
#ifdef _ANIMEC
         pparmns(:,js)  = pparmn(:)
         ppermns(:,js)  = ppermn(:)
         sigmns(:,js)   = sigmn(:)
         taumns(:,js)   = taumn(:)
         pbprmns(:,js)  = pbprmn(:)
         ppprmns(:,js)  = ppprmn(:)
         hotdmns(:,js)  = hotdmn(:)
#endif
      END DO RADIUS3

      gmns(:,1) = 0; bmns(:,1) = 0
      bsubumns(:,1) = 0
      bsubvmns(:,1) = 0
      bsubsmnc(:,1) = 2*bsubsmnc(:,2) - bsubsmnc(:,3)
      bsupumns(:,1) = 0;  bsupvmns(:,1) = 0
#ifdef _ANIMEC
      hotdmns(:,1)  = 0;  pparmns(:,1)  = 0;  ppermns(:,1) = 0
      pbprmns(:,1)  = 0;  ppprmns(:,1)  = 0
      sigmns(:,1)   = 0;  taumns(:,1)   = 0   
#endif

      IF (lfreeb) THEN        !MRC  10-15-15
         bsubumns_sur = 0
         bsubvmns_sur = 0
         bsupumns_sur = 0
         bsupvmns_sur = 0

         DO mn = 1, mnmax_nyq0
            n = NINT(xn_nyq0(mn))/nfp
            m = NINT(xm_nyq0(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                       sgn*cosmui(j,m)*sinnv(k,n1))
                  bsubumns_sur(mn) = bsubumns_sur(mn)                          &
     &                             + tsini*bsubua_sur(lk)
                  bsubvmns_sur(mn) = bsubvmns_sur(mn)                          &
     &                             + tsini*bsubva_sur(lk)
                  bsupumns_sur(mn) = bsupumns_sur(mn)                          &
     &                             + tsini*bsupua_sur(lk)
                  bsupvmns_sur(mn) = bsupvmns_sur(mn)                          &
     &                             + tsini*bsupva_sur(lk)
               END DO
            END DO
         END DO
      END IF

 200  CONTINUE

#ifdef _DEBUG
      WRITE (333, *) '    JS     M*B_S     GRAD(B_U)    J^V'
      DO mn = 1, mnmax_nyq0
         WRITE (333,'(2(a,i4))') 'm=',  INT(xm_nyq0(mn)),
     1                           ' n=', INT(xn_nyq0(mn))/nfp
         DO js = 2,ns-1
            tmult=-xm_nyq0(mn)*bsubsmns(mn,js) +
     1                    ohs*(bsubumnc(mn,js+1)-bsubumnc(mn,js))
            WRITE (333,'(i6,1p,3e12.4)') js, 
     1                  bsubsmns(mn,js)*xm_nyq0(mn),
     2             ohs*(bsubumnc(mn,js+1)-bsubumnc(mn,js)),
     3             tmult
         END DO
      END DO
#endif
#ifdef _DEBUG
      WRITE(333,*)' VMEC2000 V8.51'
      WRITE(333,*)'    MN       GMNC         GMNS           BMNC' //
     1            '          BMNS'
      IF (.NOT.lasym) THEN
      DO mn = 1, mnmax_nyq0 
         WRITE (333,'(i6,1p,4e14.4)') mn, gmnc(mn,ns/2), 0., 
     1                                    bmnc(mn,ns/2), 0.
      END DO
      ELSE
      DO mn = 1, mnmax_nyq0 
         WRITE (333,'(i6,1p,4e14.4)') mn, gmnc(mn,ns/2), gmns(mn,ns/2),
     1                                bmnc(mn,ns/2), bmns(mn,ns/2)
      END DO
      END IF
#endif
!
!     WRITE OUT ARRAYS
!
#ifdef NETCDF
      CALL cdf_write(nwout, vn_racc, raxis_cc(0:ntor))
      CALL cdf_write(nwout, vn_zacs, zaxis_cs(0:ntor)) 
      CALL cdf_write(nwout, vn_rmnc, rmnc)
      CALL cdf_write(nwout, vn_zmns, zmns)
      CALL cdf_write(nwout, vn_lmns, lmns)
      CALL cdf_write(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubsmns, bsubsmns)      !Full mesh
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
      CALL cdf_write(nwout, vn_presf, presf/mu0) 
      CALL cdf_write(nwout, vn_phi, phi) 
      CALL cdf_write(nwout, vn_phipf, twopi*signgs*phipf)
      CALL cdf_write(nwout, vn_chi, chi) 
      CALL cdf_write(nwout, vn_chipf, twopi*signgs*chipf)
      CALL cdf_write(nwout, vn_jcuru, jcuru/mu0)
      CALL cdf_write(nwout, vn_jcurv, jcurv/mu0)
      CALL cdf_write(nwout, vn_jdotb, jdotb)
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
      CALL cdf_write(nwout, vn_phip, phips(1:ns))
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

!-----------------------------------------------
!     DATA AND MSE FITS : HAVE NOT CONVERTED TO NETCDF
!     SINCE THIS WILL BE REPLACED SOON
!-----------------------------------------------
!      IF (.not.lrecon) GOTO 925
! 925  CONTINUE

#endif
      IF(lwouttxt) THEN
         DO js = 1, ns
            WRITE(nwout2, *) "JS: ", js
            MN1_OUT: DO mn = 1, mnmax
               IF (js .eq. 1) THEN
                  WRITE (nwout2, *) NINT(xm(mn)), NINT(xn(mn))
               END IF

               WRITE (nwout2, *) rmnc(mn,js), zmns(mn,js), lmns(mn,js)
               IF (lasym) THEN
                  WRITE (nwout2, *)rmns(mn,js),zmnc(mn,js),lmnc(mn,js)
               ENDIF
            END DO MN1_OUT

            MN2_OUT: DO mn = 1, mnmax_nyq0
               IF (js .eq. 1) THEN
                  WRITE (nwout2, *) NINT(xm_nyq0(mn)), NINT(xn_nyq0(mn))
               END IF
               WRITE (nwout2, *) bmnc(mn,js), gmnc(mn,js),
     1               bsubumnc(mn,js), bsubvmnc(mn,js), bsubsmns(mn,js),
     2               bsupumnc(mn,js), bsupvmnc(mn,js)
#ifdef _ANIMEC
     3              ,pparmnc (mn,js), ppermnc (mn,js), hotdmnc (mn,js),
     4               pbprmnc (mn,js), ppprmnc (mn,js), sigmnc  (mn,js),
     5               taumnc  (mn,js)
#endif
               IF (lasym) THEN
                  WRITE (nwout2, *) bmns(mn,js), gmns(mn,js), 
     1               bsubumns(mn,js), bsubvmns(mn,js), bsubsmnc(mn,js), 
     2               bsupumns(mn,js), bsupvmns(mn,js)
#ifdef _ANIMEC
     3              ,pparmns (mn,js), ppermns (mn,js), hotdmns (mn,js),
     4               pbprmns (mn,js), ppprmns (mn,js), sigmns  (mn,js),
     5               taumns  (mn,js)
#endif
               ENDIF
            END DO MN2_OUT
         END DO

!
!     HALF-MESH QUANTITIES (except phi, jcuru, jcurv which are FULL MESH - computed in eqfor)
!     NOTE: jcuru, jcurv are local current densities, NOT integrated in s and normed to twopi
!     NOTE: In version <= 6.00, mass, press are written out in INTERNAL units
!     and should be multiplied by 1/mu0 to transform to pascals. In version > 6.00,
!     the pressure, mass are in correct (physical) units
!

!     NOTE: phipf has a twopi * signgs factor compared with phips...

         WRITE (nwout2, *) (iotaf(js), presf(js)/mu0,
     1       twopi*signgs*phipf(js),
     2       phi(js), jcuru(js)/mu0, jcurv(js)/mu0, js=1,ns)
         WRITE (nwout2, *) (iotas(js), mass(js)/mu0, pres(js)/mu0,
     1   beta_vol(js), phip(js), buco(js), bvco(js), vp(js),
     2   overr(js), specw(js),js=2,ns)
!-----------------------------------------------

         WRITE (nwout2, *) aspect, betatot, betapol, betator, betaxis,
     1       b0

!-----------------------------------------------
!     New output added to version 6.20
!-----------------------------------------------
         WRITE (nwout2, *) NINT(signgs)
         WRITE (nwout2, '(a)') input_extension
         WRITE (nwout2, *) IonLarmor, VolAvgB, rbtor0, rbtor, ctor/mu0,
     1       Aminor_p, Rmajor_p, volume_p
!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
         WRITE (nwout2, *) (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     1       Dgeod(js), equif(js), js=2,ns-1)

         IF (nextcur.gt.0) THEN
            WRITE (nwout2, *) (extcur(i),i=1,nextcur)
            lcurr = ALLOCATED(curlabel) .and. lfreeb
            WRITE (nwout2, *) lcurr
            IF (lcurr) WRITE (nwout2, *) (curlabel(i),i=1,nextcur)
         ENDIF

!-----------------------------------------------
!     NOTE: jdotb is in units of A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
!-----------------------------------------------
         WRITE (nwout2, *) (fsqt(i),wdot(i),i=1,nstore_seq)
         WRITE (nwout2, *) (jdotb(js),bdotgradv(js),js=1,ns)

!-----------------------------------------------
!   Modification to obtain old fort.8 file (depracated)
!   Write out only the stellarator symmetric parts
!   Only kept for old codes. (J. Geiger)
!-----------------------------------------------
      IF (loldout) THEN
         WRITE (nfort8, '(f10.3,7i6)')
     1      gamma, nfp, ns, mpol, ntor, mnmax, itfsq, iter2/100+1
         DO js = 1, ns
            mn = 0
            DO m = 0, mpol1
               nmin0 = -ntor
               IF (m .eq. 0) nmin0 = 0
               DO n = nmin0, ntor
                  mn = mn + 1
                  IF (js .eq. 1)
     1               WRITE (nfort8,'(2i10)') NINT(xm(mn)),NINT(xn(mn))
                  WRITE (nfort8,'(5e20.13)')
     1               rmnc(mn,js),zmns(mn,js),lmns(mn,js),
     2               bmnc(mn,js),gmnc(mn,js),
     3               bsubumnc(mn,js),bsubvmnc(mn,js),bsubsmns(mn,js),
     4               bsupumnc(mn,js),bsupvmnc(mn,js)
               END DO
            END DO
         END DO
         WRITE (nfort8,'(5e20.13)') (iotas(js),mass(js),pres(js),
     1          phips(js),buco(js),bvco(js),phi(js),vp(js),
     2          jcuru(js)/mu0,jcurv(js)/mu0,specw(js),js=2,ns)
         WRITE (nfort8,'(5e20.13)') (fsqt(i),wdot(i),i=1,100)
         CLOSE(nfort8)   !last write to nfort8
      END IF
!-----------------------------------------------
!   Write diagno file (J.Geiger)
!-----------------------------------------------
      IF ((.not.lasym).and. ldiagno .and.lfreeb) THEN
         open(unit=21,file='diagno_input.data',status='unknown',
     1            action='write')
         write(21,'(a8)') "vmec2000"
         write(21,*) "nfp  mpol  ntor"
         write(21,*) nfp, mpol, ntor
         write(21,*) "rmnc"
         write(21,'(1p,e24.16)') (rmnc(mn,ns),mn=1,mnmax)
         write(21,*) "zmns"
         write(21,'(1p,e24.16)') (zmns(mn,ns),mn=1,mnmax)
         write(21,*) "potsin"
         DO mn = 1, mnpd
               write(21,'(1p,e24.16)') potvac(mn)
         END DO
         write(21,*) "phiedge"
         write(21,*) phiedge
         write(21,*) "nextcur"
         write(21,*) nextcur
         write(21,*) "external currents"
         write(21,*) extcur(1:nextcur)
         write(21,*) "plasma current"
         write(21,*) ctor
         write(21,*) "plasma current filament fc R"
         write(21,*) rmnc(1:ntor+1,1)
         write(21,*) "plasma current filament fc z"
         write(21,*) zmns(1:ntor+1,1)
         close(unit=21)
      END IF
!-----------------------------------------------
!  for diagno version 1.5 written by Sam Lazerson (SAL) start
!-----------------------------------------------
      IF(ldiagno)THEN         
         IF(lfreeb .and. (.not.lasym))THEN
            nfort = 21
            fort_file = 'diagno1.5_in.'//input_extension
            call safe_open(nfort,iwout0,fort_file,'replace',
     1                     'formatted')
            if (iwout0 .ne. 0)
     1          stop 'Error writing diagno_in. file in VMEC WROUT'

            write(21,'(a)') "vmec2000_B"
            write(21,*) "nfp  mpol  ntor"
            write(21,*) nfp, mpol, ntor
            write(21,*) "rmnc"
            write(21,'(1p,e24.16)') (rmnc(mn,ns),mn=1,mnmax)
            write(21,*) "zmns"
            write(21,'(1p,e24.16)') (zmns(mn,ns),mn=1,mnmax)

            write(21,*) "potsin"
            DO i = 1, mnpd
               write(21,'(1p,e24.16)') potvac(i)
            END DO
            
!-----  Added by SAL 11/2010
            write(21,*) "bsupu"
            js=ns
            js2=ns-1
            do m = 0, mpol1
               nmin0 = -ntor
               if (m .eq. 0) nmin0 = 0
               do n = nmin0, ntor
                  dmult = two/(mscale(m)*nscale(abs(n)))
                  if (m .eq. 0 .and. n .eq. 0) dmult = p5*dmult
                  n1 = abs(n)
                  isgn = sign(1, n)
                  d_bsupumn = 0
                  do j = 1, ntheta2
                     do k = 1, nzeta
                        lk = k + nzeta*(j - 1)
                        tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                            isgn*sinmui(j,m)*sinnv(k,n1))
                        d_bsupumn = d_bsupumn + 1.5*tcosi*bsupu(js,lk) 
     1                                    - 0.5*tcosi*bsupu(js2,lk)
                     end do
                  end do
                  write (21,'(1p,e24.16)') d_bsupumn
               end do
            end do
!-----  Added by SAL 11/2010
            write(21,*) "bsupv"
            js=ns
            js2=ns-1
            do m = 0, mpol1
               nmin0 = -ntor
               if (m .eq. 0) nmin0 = 0
               do n = nmin0, ntor
                  dmult = two/(mscale(m)*nscale(abs(n)))
                  if (m .eq. 0 .and. n .eq. 0) dmult = p5*dmult
                  n1 = abs(n)
                  isgn = sign(1, n)
                  d_bsupvmn = 0
                  do j = 1, ntheta2
                     do k = 1, nzeta
                        lk = k + nzeta*(j - 1)
                        tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                            isgn*sinmui(j,m)*sinnv(k,n1))
                        d_bsupvmn = d_bsupvmn + 1.5*tcosi*bsupv(js,lk) 
     1                                    - 0.5*tcosi*bsupv(js2,lk)
                     end do
                  end do
                  write (21,'(1p,e24.16)') d_bsupvmn
               end do
            end do
            
            write(21,*) "phiedge"
            write(21,*) phiedge
            write(21,*) "nextcur"
            write(21,*) nextcur
            write(21,*) "external currents"
            write(21,*) extcur(1:nextcur)

            write(21,*) "plasma current"
            write(21,*) ctor
            write(21,*) "plasma current filament fc R"
            write(21,*) rmnc(1:ntor+1,1)
            write(21,*) "plasma current filament fc z"
            write(21,*) zmns(1:ntor+1,1)

            close(unit=21)
         ELSE
            write(6,*)"Diagno-file request not completed!"
            write(6,*)"VMEC2000 not running in free-boundary mode!"
            write(6,*)"-or- LASYM = .true. !"
            write(6,*)"LASYM  = ",lasym
            write(6,*)"LFREEB = ",lfreeb
            write(6,*)"Check mgrid-file and namelist!"
         ENDIF
      ENDIF                   !added for diagno version 1.5 end

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
         IF (.not.lrecon) GOTO 950

         IF (imse2 - 1.gt.2 .or. itse.gt.0) THEN
            WRITE (nwout2, *) tswgt, msewgt
            CALL smoothdata(nwout2)

!       These knot values are on SQRT(s) grid
            presfactor = mu0*pthommax             !!*pfac moved to getthom
            WRITE (nwout2, *) isnodes, (sknots(i),ystark(i),y2stark(i),
     1          i=1,isnodes)
            WRITE (nwout2, *) ipnodes, (pknots(i),presfactor*ythom(i),
     1          presfactor*y2thom(i),i=1,ipnodes)
            WRITE (nwout2, *)(datamse(i),rmid(i),qmid(i),shear(i),
     1          presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
            WRITE (nwout2, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
            WRITE (nwout2, *)(rthom(i),datathom(i),i=1,itse)
         ENDIF
         IF (nobd .gt. 0) THEN
            WRITE (nwout2, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
            WRITE (nwout2, *) flmwgt
         ENDIF
         IF (nbfldn .gt. 0) THEN
            DO n = 1, nbsets
               WRITE (nwout2, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),
     1             i=1,nbfld(n))
            END DO
            WRITE (nwout2, *) bcwgt
         ENDIF

         WRITE (nwout2, *) phidiam, delphid
!
!     Write Limiter & Prout plotting specs
!
         WRITE (nwout2, *) nsets, nparts, nlim
         WRITE (nwout2, *) (nsetsn(i),i=1,nsets)
         WRITE (nwout2, *) (((pfcspec(i,j,k),i=1,nparts),j=1,nsetsn(k)),
     1       k=1,nsets)
         WRITE (nwout2, *) (limitr(i), i=1,nlim)
         WRITE (nwout2, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),
     1       j=1,nlim)
         WRITE (nwout2, *) nrgrid, nzgrid
         WRITE (nwout2, *) tokid
         WRITE (nwout2, *) rx1, rx2, zy1, zy2, condif
         WRITE (nwout2, *) imatch_phiedge

      ENDIF

 950  CONTINUE

      IF (lwouttxt) CLOSE (unit=nwout2)   !J Geiger: Close only if open, i.e. lwouttxt==true
!--------------------DEC$ ENDIF
      IF (.not. lwrite) GOTO 970   ! J Geiger: in case lwouttxt is not true and netcdf-write is finished
#ifdef NETCDF
      IF (lasym) THEN
         CALL cdf_write(nwout, vn_racs, raxis_cs(0:ntor))
         CALL cdf_write(nwout, vn_zacc, zaxis_cc(0:ntor)) 
         CALL cdf_write(nwout, vn_rmns, rmns)
         CALL cdf_write(nwout, vn_zmnc, zmnc)
         CALL cdf_write(nwout, vn_lmnc, lmnc)
         CALL cdf_write(nwout, vn_gmns, gmns)
         CALL cdf_write(nwout, vn_bmns, bmns) 
         CALL cdf_write(nwout, vn_bsubumns, bsubumns)
         CALL cdf_write(nwout, vn_bsubvmns, bsubvmns)
         CALL cdf_write(nwout, vn_bsubsmnc, bsubsmnc)
!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_write(nwout, vn_bsupumns, bsupumns)
         CALL cdf_write(nwout, vn_bsupvmns, bsupvmns)
      END IF

         IF (lfreeb) THEN     !MRC    10-15-15
            CALL cdf_write(nwout, vn_bsubumns_sur, bsubumns_sur)
            CALL cdf_write(nwout, vn_bsubvmns_sur, bsubvmns_sur)
            CALL cdf_write(nwout, vn_bsupumns_sur, bsupumns_sur)
            CALL cdf_write(nwout, vn_bsupvmns_sur, bsupvmns_sur)
         END IF
#endif
 970  CONTINUE   ! J Geiger: need to keep label 970 out of NETCDF defines.

#ifdef NETCDF
      CALL cdf_close(nwout)
#endif
!
!     RESTORE nyq ENDPOINT VALUES
!
      IF (lnyquist) THEN
         IF (mnyq .ne. 0) cosmui(:,mnyq) = 2*cosmui(:,mnyq)
         IF (nnyq .ne. 0) cosnv (:,nnyq) = 2*cosnv (:,nnyq)
      END IF

!
! DEALLOCATIONS ! J Geiger: these have been moved downwards.
!
      IF (ALLOCATED(gmnc)) DEALLOCATE(gmnc, bmnc, bsubumnc, bsubvmnc,
     1                                bsubsmns, bsupumnc, bsupvmnc
#ifdef _ANIMEC
     2   ,sigmnc,taumnc,pparmnc,ppermnc,pbprmnc,ppprmnc,hotdmnc
#endif
     3                                )
      IF (ALLOCATED(gmns)) DEALLOCATE(gmns, bmns, bsubumns, bsubvmns,
     1                                bsubsmnc, bsupumns, bsupvmns
#ifdef _ANIMEC
     2   ,sigmns,taumns,pparmns,ppermns,pbprmns,ppprmns,hotdmns
#endif
     3                                )
#ifdef _ANIMEC
      IF (ALLOCATED(tau_ana)) DEALLOCATE(sigma_ana, tau_ana, ppara,
     1                        pperpa, pbprima, ppprima, densita)
#endif
! J Geiger: check also for allocation.
      IF (ALLOCATED(gmn)) DEALLOCATE (gmn, bmn, bsubumn, bsubvmn,
     1             bsubsmn, bsupumn, bsupvmn,
#ifdef _ANIMEC
     2             sigmn, taumn, pparmn, ppermn, pbprmn, ppprmn,
     3             hotdmn,
#endif
     4             stat=istat)

      IF (ALLOCATED(bsubumnc_sur)) THEN
         DEALLOCATE(bsubumnc_sur, bsubvmnc_sur)
         DEALLOCATE(bsupumnc_sur, bsupvmnc_sur)
      END IF
      IF (ALLOCATED(bsubumns_sur)) THEN
         DEALLOCATE(bsubumns_sur, bsubvmns_sur)
         DEALLOCATE(bsupumns_sur, bsupvmns_sur)
      END IF
      IF (ALLOCATED(bsubua_sur)) THEN
         DEALLOCATE(bsubua_sur, bsubva_sur)
         DEALLOCATE(bsupua_sur, bsupva_sur)
      END IF


!-----------------------------------------------
!     FREE BOUNDARY DATA
!-----------------------------------------------
      IF (lwrite )
     1   CALL freeb_data(rmnc1, zmns1, rmns1, zmnc1, bmodmn, bmodmn1)

 1000 CONTINUE

      rzl_array = 0

      END SUBROUTINE wrout
