!-----------------------------------------------------------------------
!     Module:        stellopt_vars
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/24/2012
!     Description:   This module contains the STELLOPT global variables
!                    which are associated with the optimizer input
!                    variables.  It also contains a utility subroutine
!                    which outputs a list of the variables to a unit
!                    number.
!-----------------------------------------------------------------------
      MODULE stellopt_vars
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE vparams, ONLY: ndatafmax, mpol1d, ntord, &
                         ntor_rcws, mpol_rcws
      USE vsvd0
!-----------------------------------------------------------------------
!     Module Variables
!            nfunc_max          Maximum number of function evaluations
!            lphiedge_opt       Logical to control PHIEDGE variation
!            lcurtor_opt        Logcaal to control CURTOR variation
!            lpscale_opt        Logical to control PRES_SCALE variation
!            lbcrit_opt         Logical to control BCRIT variation
!            lextcur_opt        Logical array to control EXTCUR varaition
!            laphi_opt          Logical array to control APHI variation
!            lam_opt            Logical array to control AM variation   
!            lac_opt            Logical array to control AC variation
!            lai_opt            Logical array to control AI variation
!            lah_opt            Logical array to control AH variation
!            lat_opt            Logical array to control AT variation
!            lte_opt            Logical array to control TE variation
!            lti_opt            Logical array to control TI variation
!            lne_opt            Logical array to control NE variation
!            lth_opt            Logical array to control TH variation
!            lam_s_opt          Logical array to control AM_AUX_S variation
!            lam_f_opt          Logical array to control AM_AUX_F variation
!            lac_s_opt          Logical array to control AC_AUX_S variation
!            lac_f_opt          Logical array to control AC_AUX_F variation
!            lai_s_opt          Logical array to control AI_AUX_S variation
!            lai_f_opt          Logical array to control AI_AUX_F variation
!            lne_f_opt          Logical array to control NE_AUX_F variation
!            lte_f_opt          Logical array to control TE_AUX_F variation
!            lti_f_opt          Logical array to control TI_AUX_F variation
!            lth_f_opt          Logical array to control TH_AUX_F variation
!            lphi_s_opt         Logical array to control PHI_AUX_S variation
!            lphi_f_opt         Logical array to control PHI_AUX_F variation
!            lbound_opt         Logical array to control Boudnary variation
!            equil_type         Name of Equilibrium Code
!            ne_aux_f           Spline Knots for NE Profile (normalized to 1E19)
!            te_aux_f           Spline Knots for TE Profile
!            ti_aux_f           Spline Knots for TI Profile
!            th_aux_f           Spline Knots for T-Hot Profile
!            beamj_aux_f        Spline Knots for Beam Current Profile
!            bootj_aux_f        Spline Knots for Bootstrap Current profile.
!            REGCOIL specific variables
!              regcoil_winding_surface_separation
!              regcoil_current_density
!              regcoil_nlambda
!              regcoil_num_field_periods
!              lregcoil_rcws_rbound_c_opt, lregcoil_rcws_rbound_s_opt,
!              lregcoil_rcws_zbound_c_opt, lregcoil_rcws_zbound_s_opt
!              dregcoil_rcws_rbound_c_opt, dregcoil_rcws_rbound_s_opt,
!              dregcoil_rcws_zbound_c_opt, dregcoil_rcws_zbound_s_opt
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL  ::  lphiedge_opt, lcurtor_opt, lpscale_opt, lbcrit_opt,&
                   lmix_ece_opt, lregcoil_winding_surface_separation_opt,&
                   lregcoil_current_density_opt
      LOGICAL, DIMENSION(nigroup)  ::  lextcur_opt
      LOGICAL, DIMENSION(1:20)  ::  laphi_opt
      LOGICAL, DIMENSION(0:20)  ::  lam_opt, lac_opt, lai_opt,&
                                    lah_opt, lat_opt, &
                                    lne_opt, lte_opt, lti_opt,&
                                    lth_opt, lzeff_opt
      LOGICAL, DIMENSION(ndatafmax)  ::  lam_s_opt, lam_f_opt,&
                                lac_s_opt, lac_f_opt,&
                                lai_s_opt, lai_f_opt,&
                                lne_f_opt, lte_f_opt,&
                                lti_f_opt, lth_f_opt,&
                                lphi_s_opt, lphi_f_opt,&
                                lzeff_f_opt, lemis_xics_f_opt, &
                                lbootj_f_opt, lbeamj_f_opt, &
                                lah_f_opt, lat_f_opt
      LOGICAL, DIMENSION(0:ntord)                ::  laxis_opt
      LOGICAL, DIMENSION(-ntord:ntord,0:mpol1d)  ::  lbound_opt, lrho_opt, lmode_opt
      LOGICAL, DIMENSION(-ntord:ntord,-mpol1d:mpol1d) :: ldeltamn_opt
      INTEGER, PARAMETER :: maxcoilknots=40
      LOGICAL, DIMENSION(nigroup,maxcoilknots)        ::  lcoil_spline
      INTEGER, DIMENSION(nigroup)                     ::  coil_nknots
      LOGICAL  ::  lwindsurf
      INTEGER  ::  nfunc_max
      REAL(rprec)     ::  dphiedge_opt, dcurtor_opt, dbcrit_opt, &
                          dpscale_opt, dmix_ece_opt, &
                          dregcoil_winding_surface_separation_opt, &
                          dregcoil_current_density_opt
      REAL(rprec)     ::  phiedge_min, curtor_min, bcrit_min, &
                          pscale_min, mix_ece_min, &
                          regcoil_winding_surface_separation_min, &
                          regcoil_current_density_min
      REAL(rprec)     ::  phiedge_max, curtor_max, bcrit_max, &
                          pscale_max, mix_ece_max, &
                          regcoil_winding_surface_separation_max, &
                          regcoil_current_density_max
      REAL(rprec), DIMENSION(nigroup)  ::  dextcur_opt,extcur_min,extcur_max
      REAL(rprec), DIMENSION(1:20)     ::  daphi_opt,aphi_min,aphi_max
      REAL(rprec), DIMENSION(0:20)     ::  dam_opt, dac_opt, dai_opt,&
                                           dah_opt, dat_opt,&
                                           dte_opt, dne_opt, dti_opt, dth_opt,&
                                           dzeff_opt, &
                                           am_min, ac_min, ai_min, &
                                           ah_min, at_min, &
                                           am_max, ac_max, ai_max, &
                                           ah_max, at_max,&
                                           te_min, ne_min, ti_min, th_min, &
                                           te_max, ne_max, ti_max, th_max, &
                                           zeff_max, zeff_min
      REAL(rprec)                       :: mix_ece
      REAL(rprec)                       :: regcoil_winding_surface_separation
      REAL(rprec)                       :: regcoil_current_density
      INTEGER :: regcoil_nlambda, regcoil_num_field_periods
      LOGICAL, DIMENSION(-mpol_rcws:mpol_rcws, &
                         -ntor_rcws:ntor_rcws) :: lregcoil_rcws_rbound_c_opt , &
                                                  lregcoil_rcws_rbound_s_opt, &
                                                  lregcoil_rcws_zbound_c_opt, &
                                                  lregcoil_rcws_zbound_s_opt
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, &
                             -ntor_rcws:ntor_rcws) :: dregcoil_rcws_rbound_c_opt , &
                                                      dregcoil_rcws_rbound_s_opt, &
                                                      dregcoil_rcws_zbound_c_opt, &
                                                      dregcoil_rcws_zbound_s_opt
      REAL(rprec), DIMENSION(0:20)      ::  te_opt, ti_opt, ne_opt, th_opt, zeff_opt                                     
      REAL(rprec), DIMENSION(ndatafmax) ::  ne_aux_s, te_aux_s, &
                                            ti_aux_s, th_aux_s, &
                                            zeff_aux_s, &
                                            phi_aux_s, beamj_aux_s,&
                                            bootj_aux_s, sfincs_s, emis_xics_s
      INTEGER                           ::  sfincs_min_procs
      CHARACTER(256)  ::  sfincs_Er_option
      REAL(rprec)                       ::  vboot_tolerance
      REAL(rprec), DIMENSION(ndatafmax) ::  ne_aux_f, te_aux_f, &
                                            ti_aux_f, th_aux_f,&
                                            zeff_aux_f, &
                                            phi_aux_f, beamj_aux_f, &
                                            bootj_aux_f, emis_xics_f
      REAL(rprec), DIMENSION(ndatafmax) ::  dam_s_opt, dam_f_opt, &
                                            dac_s_opt, dac_f_opt, &
                                            dai_s_opt, dai_f_opt, &
                                            dphi_s_opt, dphi_f_opt, &
                                            dne_f_opt, dte_f_opt, &
                                            dti_f_opt, dth_f_opt, &
                                            dzeff_f_opt, &
                                            dbeamj_f_opt, dbootj_f_opt,&
                                            dat_f_opt, dah_f_opt, &
                                            demis_xics_f_opt
      REAL(rprec), DIMENSION(ndatafmax) ::  am_f_min, ac_f_min, &
                                            ai_f_min, phi_f_min, &
                                            ne_f_min, te_f_min, &
                                            ti_f_min, th_f_min, &
                                            zeff_f_min, emis_xics_f_min, &
                                            beamj_f_min, bootj_f_min, &
                                            ah_f_min, at_f_min
      REAL(rprec), DIMENSION(ndatafmax) ::  am_f_max, ac_f_max, &
                                            ai_f_max, phi_f_max, &
                                            ne_f_max, te_f_max, &
                                            ti_f_max, th_f_max, &
                                            zeff_f_max,  emis_xics_f_max, &
                                            beamj_f_max, bootj_f_max,&
                                            ah_f_max, at_f_max
      REAL(rprec), DIMENSION(0:ntord)                     ::  daxis_opt
      REAL(rprec), DIMENSION(0:ntord)                     ::  raxis_min, raxis_max,&
                                                              zaxis_min, zaxis_max
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  rhobc, modemn
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  dbound_opt, drho_opt
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  rbc_min, rbc_max
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  zbs_min, zbs_max
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  bound_min, bound_max
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  zbc_min, zbc_max, rbs_min, rbs_max
      REAL(rprec), DIMENSION(-ntord:ntord,-mpol1d:mpol1d) ::  deltamn
      REAL(rprec), DIMENSION(-ntord:ntord,-mpol1d:mpol1d) ::  ddeltamn_opt
      REAL(rprec), DIMENSION(-ntord:ntord,-mpol1d:mpol1d) ::  delta_min, delta_max
      REAL(rprec), DIMENSION(nigroup,maxcoilknots) :: coil_splinesx,coil_splinesy,coil_splinesz,&
                                                      coil_splinefx,coil_splinefy,coil_splinefz
      REAL(rprec), DIMENSION(nigroup,maxcoilknots) :: dcoil_spline
      REAL(rprec), DIMENSION(nigroup,maxcoilknots) :: coil_splinefx_min,coil_splinefy_min,coil_splinefz_min,&
                                                      coil_splinefx_max,coil_splinefy_max,coil_splinefz_max

      ! Regcoil Winding Surface (rcws): Boundary+min/max
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, -ntor_rcws:ntor_rcws) :: regcoil_rcws_rbound_c, regcoil_rcws_rbound_s
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, -ntor_rcws:ntor_rcws) :: regcoil_rcws_rbound_c_min, regcoil_rcws_rbound_s_min
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, -ntor_rcws:ntor_rcws) :: regcoil_rcws_rbound_c_max, regcoil_rcws_rbound_s_max
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, -ntor_rcws:ntor_rcws) :: regcoil_rcws_zbound_c, regcoil_rcws_zbound_s
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, -ntor_rcws:ntor_rcws) :: regcoil_rcws_zbound_c_min, regcoil_rcws_zbound_s_min
      REAL(rprec), DIMENSION(-mpol_rcws:mpol_rcws, -ntor_rcws:ntor_rcws) :: regcoil_rcws_zbound_c_max, regcoil_rcws_zbound_s_max

      CHARACTER(256)  ::  equil_type, te_type, ne_type, ti_type, th_type, &
                          beamj_type, bootj_type, zeff_type, emis_xics_type, windsurfname, &
                          regcoil_nescin_filename, bootcalc_type
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sfincs_J_dot_B_flux_surface_average, sfincs_B_squared_flux_surface_average
      
      ! These are not really variable parameters as we don't vary them
      ! yet
      REAL(rprec), DIMENSION(ndatafmax) :: nustar_s, nustar_f
      
      CHARACTER, DIMENSION(nigroup) :: coil_type

      REAL(rprec) :: phiedge_old  ! For keeping track of phiedge
      
      INTEGER, PARAMETER ::  norm_dex   = -5
      
      INTEGER, PARAMETER ::  iphiedge   = 11
      INTEGER, PARAMETER ::  icurtor    = 12
      INTEGER, PARAMETER ::  ibcrit     = 13
      INTEGER, PARAMETER ::  ipscale    = 14
      INTEGER, PARAMETER ::  imixece    = 15
      INTEGER, PARAMETER ::  iextcur    = 21
      INTEGER, PARAMETER ::  iaphi      = 31
      INTEGER, PARAMETER ::  iam        = 32
      INTEGER, PARAMETER ::  iac        = 33
      INTEGER, PARAMETER ::  iai        = 34
      INTEGER, PARAMETER ::  iah        = 35
      INTEGER, PARAMETER ::  iat        = 36
      INTEGER, PARAMETER ::  iah_aux_f  = 361
      INTEGER, PARAMETER ::  iat_aux_f  = 362
      INTEGER, PARAMETER ::  ite        = 371
      INTEGER, PARAMETER ::  ine        = 372
      INTEGER, PARAMETER ::  iti        = 373
      INTEGER, PARAMETER ::  ith        = 374
      INTEGER, PARAMETER ::  izeff      = 375
      INTEGER, PARAMETER ::  iam_aux_s  = 41
      INTEGER, PARAMETER ::  iam_aux_f  = 51
      INTEGER, PARAMETER ::  iac_aux_s  = 42
      INTEGER, PARAMETER ::  iac_aux_f  = 52
      INTEGER, PARAMETER ::  ibeamj_aux_f = 521
      INTEGER, PARAMETER ::  ibootj_aux_f = 522
      INTEGER, PARAMETER ::  iemis_xics_f = 531
      INTEGER, PARAMETER ::  iai_aux_s  = 43
      INTEGER, PARAMETER ::  iai_aux_f  = 53
      INTEGER, PARAMETER ::  iphi_aux_s = 44
      INTEGER, PARAMETER ::  iphi_aux_f = 54
      INTEGER, PARAMETER ::  ine_aux_f  = 55
      INTEGER, PARAMETER ::  izeff_aux_f  = 551
      INTEGER, PARAMETER ::  ite_aux_f  = 56
      INTEGER, PARAMETER ::  iti_aux_f  = 57
      INTEGER, PARAMETER ::  ith_aux_f  = 58
      INTEGER, PARAMETER ::  icoil_splinefx  = 60
      INTEGER, PARAMETER ::  icoil_splinefy  = 61
      INTEGER, PARAMETER ::  icoil_splinefz  = 62
      INTEGER, PARAMETER ::  ibound_rbc = 91
      INTEGER, PARAMETER ::  ibound_rbs = 92
      INTEGER, PARAMETER ::  ibound_zbc = 93
      INTEGER, PARAMETER ::  ibound_zbs = 94
      INTEGER, PARAMETER ::  irhobc     = 95
      INTEGER, PARAMETER ::  ideltamn   = 96
      INTEGER, PARAMETER ::  imodemn    = 97
      INTEGER, PARAMETER ::  iraxis_cc  = 911
      INTEGER, PARAMETER ::  iraxis_cs  = 912
      INTEGER, PARAMETER ::  izaxis_cc  = 913
      INTEGER, PARAMETER ::  izaxis_cs  = 914
      INTEGER, PARAMETER ::  iregcoil_winding_surface_separation   = 5150
      INTEGER, PARAMETER ::  iregcoil_current_density   = 5151
      ! 5152 - 5159 reserved for future REGCOIIL options
      INTEGER, PARAMETER ::  iregcoil_rcws_rbound_c = 5160
      INTEGER, PARAMETER ::  iregcoil_rcws_rbound_s = 5161
      INTEGER, PARAMETER ::  iregcoil_rcws_zbound_c = 5162
      INTEGER, PARAMETER ::  iregcoil_rcws_zbound_s = 5163
      
      REAL(rprec), PARAMETER :: ne_norm = 1.0E18
      
      CONTAINS      
      
      SUBROUTINE write_vars(iunit,var_num,var_dex1,var_dex2)
      INTEGER, INTENT(in) :: iunit, var_num, var_dex1, var_dex2
      CHARACTER*(*), PARAMETER ::  out_format = '(5X,A)'
      CHARACTER*(*), PARAMETER ::  out_format_1D = '(5X,A,I3.3,A)'
      CHARACTER*(*), PARAMETER ::  out_format_2D = '(5X,A,I3.3,A,I3.3,A)'
      CHARACTER*(*), PARAMETER ::  out_format_2DB = '(5X,A,I4.3,A,I4.3,A)'
      SELECT CASE(var_num)

         CASE(iphiedge)
            WRITE(iunit,out_format) 'PHIEDGE:  Total Enclosed Toroidal Flux'
         CASE(imixece)
            WRITE(iunit,out_format) 'MIX_ECE:  O2/X2 ECE Mixing'
         CASE(icurtor)
            WRITE(iunit,out_format) 'CURTOR:   Total Toroidal Current'
         CASE(ibcrit)
            WRITE(iunit,out_format) 'BCRIT:    Aux Parameter (Critical Field/Mach)'
         CASE(ipscale)
            WRITE(iunit,out_format) 'PRES_SCALE:  Pressure Scaling Factor'
         CASE(iextcur)
            WRITE(iunit,out_format_1D) 'EXTCUR(',var_dex1,'):  Vacuum Field Current'
         CASE(iaphi)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'APHI_NORM:  Toroidal Flux Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'APHI(',var_dex1,'):  Toroidal Flux Coefficient'
            END IF
         CASE(iam)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AM:  Pressure Profile Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AM(',var_dex1,'):  Pressure Profile Coefficient'
            END IF
         CASE(iac)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AC:  Current Profile Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AC(',var_dex1,'):  Current Profile Coefficient'
            END IF
         CASE(iai)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AI:  Iota Profile Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AI(',var_dex1,'):  Iota Profile Coefficient'
            END IF
         CASE(iah)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AH:  Hot Particle Fraction Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AH(',var_dex1,'):  Hot Particle Fraction Coefficient'
            END IF
         CASE(iat)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AT:  Anisotropy Profile Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AT(',var_dex1,'):  Anisotropy Profile Coefficient'
            END IF
         CASE(ine)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'NE:  Electron Density Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'NE(',var_dex1,'):  Electron Density Coefficient'
            END IF
         CASE(izeff)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'ZEFF:  Z-Effective Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'ZEFF(',var_dex1,'):  ZEFF:  Z-Effective Coefficient'
            END IF
         CASE(ite)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'TE:  Electron Temp. Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'TE(',var_dex1,'):  Electron Temp Coefficient'
            END IF
         CASE(iti)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'TI:  Ion Temp. Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'TI(',var_dex1,'):  Ion Temp. Coefficient'
            END IF
         CASE(ith)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'TH:  Hot Particle Temp. Coefficient (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'TH(',var_dex1,'):  Hot Particle Temp. Coefficient'
            END IF
         CASE(iam_aux_s)
            WRITE(iunit,out_format_1D) 'AM_AUX_S(',var_dex1,'):  Pressure Profile Knot Location'
         CASE(iam_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AM_AUX_F:  Presure Profile Knot (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AM_AUX_F(',var_dex1,'):  Presure Profile Knot'
            END IF
         CASE(iac_aux_s)
            WRITE(iunit,out_format_1D) 'AC_AUX_S(',var_dex1,'):  Current Profile Knot Location'
         CASE(iac_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AC_AUX_F:  Current Profile Knot (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AC_AUX_F(',var_dex1,'):  Current Profile Knot'
            END IF
         CASE(ibeamj_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'BEAMJ_AUX_F:  Beam Driven Current Profile Knot (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'BEAMJ_AUX_F(',var_dex1,'):  Beam Driven Current Profile Knot'
            END IF
         CASE(ibootj_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'BOOTJ_AUX_F:  Bootstrap Current Profile Knot (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'BOOTJ_AUX_F(',var_dex1,'):  Bootstrap Current Profile Knot'
            END IF
         CASE(iemis_xics_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'EMIS_XICS_F:  XICS Emissivity Profile Coef (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'EMIS_XICS_F(',var_dex1,'):  XICS Emissivity Profile Coef'
            END IF
         CASE(iai_aux_s)
            WRITE(iunit,out_format_1D) 'AI_AUX_S(',var_dex1,'):  Iota Profile Knot Location'
         CASE(iai_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AI_AUX_F:  Iota Profile Knot (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AI_AUX_F(',var_dex1,'):  Iota Profile Knot'
            END IF
         CASE(iphi_aux_s)
            WRITE(iunit,out_format_1D) 'PHI_AUX_S(',var_dex1,'):  E-Static Potential Profile Knot Location'
         CASE(iphi_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'PHI_AUX_F:  E-Static Potential Profile Knot (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'PHI_AUX_F(',var_dex1,'):  E-Static Potential Profile Knot'
            END IF
         CASE(ine_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'NE_AUX_F:  Electron Density Profile (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'NE_AUX_F(',var_dex1,'):  Electron Density Profile Knot'
            END IF
         CASE(izeff_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'ZEFF_AUX_F:  Z-Effective Profile (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'ZEFF_AUX_F(',var_dex1,'):  Z-Effective Profile Knot'
            END IF
         CASE(ite_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'TE_AUX_F:  Electron Temperature Profile (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'TE_AUX_F(',var_dex1,'):  Electron Temperature Profile Knot'
            END IF
         CASE(iti_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'TI_AUX_F:  Ion Temperature Profile (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'TI_AUX_F(',var_dex1,'):  Ion Temperature Profile Knot'
            END IF
         CASE(ith_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'TH_AUX_F:  Hot Particle Temperature Profile (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'TH_AUX_F(',var_dex1,'):  Hot Particle Temperature Profile Knot'
            END IF
         CASE(iah_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AH_AUX_F:  Auxillary Array (hot temp/rotation) (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AH_AUX_F(',var_dex1,'):  Auxillary Array Knot (hot temp/rotation)'
            END IF
         CASE(iat_aux_f)
            IF (var_dex2 == norm_dex) THEN
               WRITE(iunit,out_format) 'AT_AUX_F:  Hot Particle Temperature Profile (normalization)'
            ELSE
               WRITE(iunit,out_format_1D) 'AT_AUX_F(',var_dex1,'):  Hot Particle Temperature Profile Knot'
            END IF
         CASE(iraxis_cc)
            WRITE(iunit,out_format_1D) 'RAXIS(',var_dex1,'):  Axis Specification (COS)'
         CASE(iraxis_cs)
            WRITE(iunit,out_format_1D) 'RAXIS(',var_dex1,'):  Axis Specification (SIN)'
         CASE(izaxis_cc)
            WRITE(iunit,out_format_1D) 'ZAXIS(',var_dex1,'):  Axis Specification (COS)'
         CASE(izaxis_cs)
            WRITE(iunit,out_format_1D) 'ZAXIS(',var_dex1,'):  Axis Specification (SIN)'
         CASE(ibound_rbc)
            WRITE(iunit,out_format_2DB) 'RBC(',var_dex1,',',var_dex2,'):  Radial Boundary Specification (COS)'
         CASE(ibound_rbs)
            WRITE(iunit,out_format_2DB) 'RBS(',var_dex1,',',var_dex2,'):  Radial Boundary Specification (SIN)'
         CASE(ibound_zbc)
            WRITE(iunit,out_format_2DB) 'ZBC(',var_dex1,',',var_dex2,'):  Vertical Boundary Specification (COS)'
         CASE(ibound_zbs)
            WRITE(iunit,out_format_2DB) 'ZBS(',var_dex1,',',var_dex2,'):  Vertical Boundary Specification (SIN)'
         CASE(irhobc)
            WRITE(iunit,out_format_2DB) 'RHO(',var_dex1,',',var_dex2,'):  Boundary Specifiction (Hirsh. -Bres.)'
         CASE(ideltamn)
            WRITE(iunit,out_format_2DB) 'DELTA(',var_dex1,',',var_dex2,'):  Boundary Specifiction (Garabedian)'
         CASE(imodemn)
            WRITE(iunit,out_format_2DB) 'MODE(',var_dex1,',',var_dex2,'):  Boundary Specifiction (Lazerson)'
         CASE(icoil_splinefx)
            WRITE(iunit,out_format_2DB) 'COIL_SPLINEX(',var_dex1,',',var_dex2,'):  Coil Spline Knots (X)'
         CASE(icoil_splinefy)
            WRITE(iunit,out_format_2DB) 'COIL_SPLINEY(',var_dex1,',',var_dex2,'):  Coil Spline Knots (Y)'
         CASE(icoil_splinefz)
            WRITE(iunit,out_format_2DB) 'COIL_SPLINEZ(',var_dex1,',',var_dex2,'):  Coil Spline Knots (Z)'

         ! REGCOIL cases
         CASE(iregcoil_winding_surface_separation)
            WRITE(iunit,out_format) 'REGCOIL_SEPARATION: Coil winding surface separation'
         CASE(iregcoil_current_density)
            WRITE(iunit,out_format) 'REGCOIL_SEPARATION: Winding surface current density'
         CASE(iregcoil_rcws_rbound_c)
            WRITE(iunit,out_format_2DB) 'REGCOIL_RCWS_rbound_c(',var_dex1,',',var_dex2,'):  REGCOIL Winding Surface Boundary Radial Specification (COS MN)'
         CASE(iregcoil_rcws_rbound_s)
            WRITE(iunit,out_format_2DB) 'REGCOIL_RCWS_rbound_s(',var_dex1,',',var_dex2,'):  REGCOIL Winding Surface Boundary Radial Specification (SIN MN)'
         CASE(iregcoil_rcws_zbound_c)
            WRITE(iunit,out_format_2DB) 'REGCOIL_RCWS_zbound_c(',var_dex1,',',var_dex2,'):  REGCOIL Winding Surface Boundary Vertical Specification (COS MN)'
         CASE(iregcoil_rcws_zbound_s)
            WRITE(iunit,out_format_2DB) 'REGCOIL_RCWS_zbound_s(',var_dex1,',',var_dex2,'):  REGCOIL Winding Surface Boundary Vertical Specification (SIN MN)'
         ! END of REGCOIL cases
      END SELECT
      END SUBROUTINE write_vars

      SUBROUTINE bcast_vars(loc_master,loc_comm,ierr)
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h' 
!DEC$ ENDIF        
      INTEGER :: loc_master, loc_comm, ierr
      INTEGER :: n_temp
      ierr = 0
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(phiedge_min,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(phiedge_max,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(mix_ece_min,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(mix_ece_max,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(curtor_min,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(curtor_max,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(pscale_min,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(pscale_max,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(bcrit_min,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(bcrit_max,1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(extcur_min,nigroup,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(extcur_max,nigroup,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(aphi_min,20,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(aphi_max,20,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(am_min,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(am_max,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ac_min,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ac_max,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ai_min,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ai_max,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ah_min,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ah_max,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(at_min,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(at_max,21,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(am_f_min,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(am_f_max,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ac_f_min,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ac_f_max,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ai_f_min,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ai_f_max,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ah_f_min,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(ah_f_max,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(at_f_min,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(at_f_max,ndatafmax,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(raxis_min,ntord+1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(raxis_max,ntord+1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(zaxis_min,ntord+1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(zaxis_max,ntord+1,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      n_temp = (2*ntord+1)*(mpol1d+1)
      CALL MPI_BCAST(rbc_min,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(rbc_max,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(zbs_min,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(zbs_max,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      n_temp = (2*ntord+1)*(mpol1d+1)
      CALL MPI_BCAST(bound_min,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(bound_max,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      n_temp = (2*ntord+1)*(2*mpol1d+1)
      CALL MPI_BCAST(delta_min,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN
      CALL MPI_BCAST(delta_max,n_temp,MPI_REAL8, loc_master, loc_comm,ierr)
      IF (ierr /= MPI_SUCCESS) RETURN


!DEC$ ENDIF
      RETURN

      END SUBROUTINE bcast_vars


      END MODULE stellopt_vars
