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
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  ::    nr, nphi, nz, nte, nne, nti, npot
      REAL(rprec) :: rmin, rmax, zmin, zmax, phimin, phimax, tmin, tmax, delta_t, &
                     vc_adapt_tol
      REAL(rprec), ALLOCATABLE :: raxis(:), zaxis(:), phiaxis(:)
      REAL(rprec), ALLOCATABLE :: B_R(:,:,:), B_PHI(:,:,:), B_Z(:,:,:), MODB(:,:,:),&
                                  TE(:,:,:), NE(:,:,:), TI(:,:,:), &
                                  S_ARR(:,:,:), U_ARR(:,:,:), POT_ARR(:,:,:)
      TYPE(EZspline3_r8) :: BR_spl, BPHI_spl, BZ_spl, MODB_spl, TE_spl, NE_spl, &
                            TI_spl, S_spl, U_spl, POT_spl
      TYPE(EZspline1_r8) :: TE_spl_s, NE_spl_s, TI_spl_S, Vp_spl_s, POT_spl_s

      END MODULE beams3d_grid
