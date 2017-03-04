!-----------------------------------------------------------------------
!     Module:        fieldlines_grid
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the FIELDLINES grid variables.
!-----------------------------------------------------------------------
      MODULE fieldlines_grid
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
!         delta_phi      Stepsize in toroidal angle [radians]
!         vc_adapt_tol   Adaptive Integration tolerance for virtual casing
!         raxis          Radial Grid
!         phiaxis        Toroidal Grid
!         zaxis          Vertical Grid
!         B_R            Radial Magnetic Field [T]
!         B_PHI          Toroidal Magnetic Field [T]
!         B_Z            Vertical Magnetic Field [T]
!         BR_spl         EZSpline Object for B_R/B_PHI
!         BZ_spl         EZSpline Object for B_Z/B_PHI
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  ::    nr, nphi, nz
      REAL(rprec) :: rmin, rmax, zmin, zmax, phimin, phimax, delta_phi,&
                     vc_adapt_tol
      REAL(rprec), ALLOCATABLE :: raxis(:), zaxis(:), phiaxis(:)
      REAL(rprec), ALLOCATABLE :: B_R(:,:,:), B_Z(:,:,:), B_PHI(:,:,:)
      REAL(rprec), ALLOCATABLE :: MU3D(:,:,:)
      TYPE(EZspline3_r8) :: BR_spl, BZ_spl, MU_spl, MODB_spl

      ! For the boundary scaling routine
      INTEGER :: mnmax_m, nfp_m, isgn_m
      INTEGER, ALLOCATABLE :: ixm_m(:), ixn_m(:)
      REAL(rprec) ::  vb_ws, u0, v0, dr_m, rb_ws, zb_ws
      REAL(rprec), ALLOCATABLE :: rmnc_m(:), zmns_m(:), rmns_m(:), zmnc_m(:)
      
      ! The following are for using the TOMS661 interpolation algorithm
      INTEGER :: n_qshep, nq_qshep, nw_qshep, nr_qshep
      INTEGER, ALLOCATABLE :: lnext_qshep(:)
      INTEGER, ALLOCATABLE :: lcell_qshep(:,:,:)
      REAL(rprec) :: rmin_qshep(3), rdel_qshep(3), rmax_qshep
      REAL(rprec), ALLOCATABLE :: rsq1_qshep(:), rsq2_qshep(:)
      REAL(rprec), ALLOCATABLE :: x_qshep(:), y_qshep(:), z_qshep(:)
      REAL(rprec), ALLOCATABLE :: v_qshep(:,:)
      REAL(rprec), ALLOCATABLE :: b1_qshep(:), b2_qshep(:)
      REAL(rprec), ALLOCATABLE :: a1_qshep(:,:), a2_qshep(:,:)
      

      END MODULE fieldlines_grid
