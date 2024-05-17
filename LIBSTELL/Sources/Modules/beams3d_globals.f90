!-----------------------------------------------------------------------
!     Module:        beams3d_globals
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          05/17/2024
!     Description:   This module contains the BEAMS3D global variables
!                    needed by the input namelist.
!-----------------------------------------------------------------------
      MODULE beams3d_globals
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! Moved from beams3d_runtime
      INTEGER, PARAMETER :: MAXPARTICLES = 2**18
      INTEGER, PARAMETER :: MAXBEAMS = 32
      INTEGER, PARAMETER :: MAXPROFLEN = 512
      INTEGER, PARAMETER :: NION = 4
      LOGICAL :: lverb, lcollision, lrestart_particles, ldebug, &
                 lfusion, lfusion_alpha, lfusion_He3, lfusion_proton, &
                 lfusion_tritium, lkick, lgcsim, lbeam, lbbnbi
      INTEGER :: npoinc, nbeams, nparticles_start, duplicate_factor
      INTEGER, DIMENSION(MAXBEAMS) :: Dex_beams
      REAL(rprec) :: follow_tol, pi2, ne_scale, te_scale, ti_scale, &
                     zeff_scale, fusion_scale, lendt_m, te_col_min
      REAL(rprec), DIMENSION(MAXBEAMS) :: Adist_beams, Asize_beams, Div_beams, E_beams, mass_beams, &
                                        charge_beams, Zatom_beams, P_beams
      REAL(rprec), DIMENSION(MAXBEAMS, 2) :: r_beams, z_beams, phi_beams
      REAL(rprec), DIMENSION(MAXPROFLEN) :: TE_AUX_S, TE_AUX_F, NE_AUX_S, NE_AUX_F, TI_AUX_S, TI_AUX_F,&
                                            POT_AUX_S, POT_AUX_F, ZEFF_AUX_S, ZEFF_AUX_F                                  
      INTEGER, DIMENSION(NION) :: NI_AUX_Z   
      REAL(rprec), DIMENSION(MAXPROFLEN) :: NI_AUX_S
      REAL(rprec), DIMENSION(NION,MAXPROFLEN) :: NI_AUX_F
      REAL(rprec), DIMENSION(NION) :: NI_AUX_M
      REAL(rprec), DIMENSION(MAXPARTICLES) :: r_start_in, phi_start_in, z_start_in, vll_start_in, &
                                            & mu_start_in, charge_in, Zatom_in, mass_in, t_end_in, &
                                            vr_start_in, vphi_start_in, vz_start_in, weight_in
      CHARACTER(256) :: id_string, int_type

      ! moved from beams3d_lines
      INTEGER  ::  ns_prof1, ns_prof2, ns_prof3, ns_prof4, ns_prof5, nsh_prof4, nparticles
      REAL(rprec) :: partvmax, plasma_Zmean, plasma_mass

      ! moved from beams3d_grid
      INTEGER  ::    nr, nphi, nz, nte, nne, nti, nzeff, npot, nr_fida, nphi_fida, nz_fida, nenergy_fida, npitch_fida, dexionT, dexionD, dexionHe3
      REAL :: rmin_fida, rmax_fida, zmin_fida, zmax_fida, phimin_fida, phimax_fida, emin_fida, pimin_fida
      REAL(rprec) :: rmin, rmax, zmin, zmax, phimin, phimax, &
                     vc_adapt_tol, therm_factor, B_kick_min, &
                     B_kick_max, E_kick, freq_kick, rho_fullorbit, &
                     t_fida, s_max,s_max_te, s_max_ne,s_max_zeff,s_max_ti, s_max_pot



      END MODULE beams3d_globals
