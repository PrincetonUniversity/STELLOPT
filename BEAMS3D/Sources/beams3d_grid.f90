!-----------------------------------------------------------------------
!     Module:        beams3d_grid
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This module contains the BEAMS3D grid variables.
!-----------------------------------------------------------------------
      MODULE beams3d_grid
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Module Variables
!         nr             Number of radial gridpoints
!         nphi           Number of toroidal gridpoints
!         nz             Number of vertical gridpoints
!         rmin           Minimum radial extent of grid [m]
!         rmax           Maximum radial extent of grid [m]
!         phimin         Minimum toroidal extent of grid [radians]
!         phimax         Maximum toroidal extent of grid [radians]
!         zmin           Minimum vertical extent of grid [m]
!         zmax           Maximum vertical extent of grid [m]
!         delta_t        Stepsize [seconds]
!         vc_adapt_tol   Adaptive Integration tolerance for virtual casing
!         raxis          Radial Grid
!         phiaxis        Toroidal Grid
!         zaxis          Vertical Grid
!         B_R            Radial Magnetic Field [T]
!         B_PHI          Toroidal Magnetic Field [T]
!         B_Z            Vertical Magnetic Field [T]
!         MODB           Magnitude of Magnetic Field (|B|) [T]
!         BR_spl         EZSpline Object for B_R
!         BPHI_spl       EZSpline Object for B_PHI
!         BZ_spl         EZSpline Object for B_Z
!         MODB_spl       EZSpline Object for MODB
!         plasma_Zavg    <Z> = sum(n*Z*Z)/sum(n*Z) sum over j ion species
!         plasma_Zmean   [Z] = sum(n*Z*Z*mi/mj)/sum(n*Z) sum over j ion species (mi: plasma mass)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  ::    nr, nphi, nz, nte, nne, nti, nzeff, npot
      INTEGER  ::    win_raxis, win_phiaxis, win_zaxis, win_B_R, win_B_PHI, win_B_Z,&
                     win_MODB, win_TE, win_NE, win_TI, win_ZEFF_ARR,&
                     win_S_ARR, win_U_ARR,win_X_ARR,win_Y_ARR, win_POT_ARR, win_BR4D, win_BPHI4D, &
                     win_BZ4D, win_MODB4D, win_TE4D, win_NE4D, win_TI4D, win_ZEFF4D, &
                     win_S4D, win_U4D,  win_X4D, win_Y4D, win_POT4D, win_req_axis, win_zeq_axis, &
                     win_wall_load, win_wall_shine, win_hr, win_hp, win_hz, &
                     win_hri, win_hpi, win_hzi, win_NI5D, win_NI, &
                     win_raxis_fida, win_phiaxis_fida, win_zaxis_fida, win_energy_fida, win_pitch_fida, &
                     nr_fida, nphi_fida, nz_fida, nenergy_fida, npitch_fida
      REAL(rprec) :: rmin, rmax, zmin, zmax, phimin, phimax, tmin, tmax, delta_t, &
                     vc_adapt_tol, psiedge_eq, phiedge_eq, plasma_Zmean, plasma_mass, &
                     reff_eq, therm_factor, B_kick_min, B_kick_max, &
                     E_kick, freq_kick, t_fida
      REAL(rprec), POINTER :: raxis(:), zaxis(:), phiaxis(:)
      REAL(rprec), POINTER :: req_axis(:), zeq_axis(:)
      REAL :: rmin_fida, rmax_fida, zmin_fida, zmax_fida, phimin_fida, phimax_fida, emin_fida, pimin_fida
      REAL(rprec), POINTER :: raxis_fida(:), zaxis_fida(:), phiaxis_fida(:), energy_fida(:), pitch_fida(:)
      REAL(rprec), POINTER :: wall_load(:,:), wall_shine(:,:)
      REAL(rprec), POINTER :: B_R(:,:,:),B_PHI(:,:,:), B_Z(:,:,:), MODB(:,:,:),&
                                  TE(:,:,:), NE(:,:,:), TI(:,:,:), ZEFF_ARR(:,:,:), &
                                  S_ARR(:,:,:), U_ARR(:,:,:), X_ARR(:,:,:), Y_ARR(:,:,:), POT_ARR(:,:,:)
      REAL(rprec), POINTER :: NI(:,:,:,:)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                                                  NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET
      REAL(rprec), DIMENSION(:,:,:,:), POINTER :: BR4D, BPHI4D, BZ4D, MODB4D, &
                                  TE4D, NE4D, TI4D, ZEFF4D, &
                                  S4D, U4D, X4D, Y4D, POT4D
      REAL(rprec), DIMENSION(:,:,:,:,:), POINTER :: NI5D
      REAL*8 ::      eps1, eps2, eps3
      REAL*8, parameter :: small = 1.e-10_ezspline_r8
      REAL*8, POINTER :: hr(:), hp(:), hz(:)
      REAL*8, POINTER :: hri(:), hpi(:), hzi(:)
      TYPE(EZspline3_r8) :: BR_spl, BPHI_spl, BZ_spl, MODB_spl, TE_spl, NE_spl, &
                            TI_spl, ZEFF_spl, S_spl, U_spl, X_spl, Y_spl,POT_spl
      TYPE(EZspline1_r8) :: TE_spl_s, NE_spl_s, TI_spl_S, ZEFF_spl_s, Vp_spl_s, POT_spl_s
      TYPE(EZspline1_r8), DIMENSION(4) :: NI_spl_s

      END MODULE beams3d_grid
