!-----------------------------------------------------------------------
!     Module:        torlines_background
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This module contains the background coordinate
!                    variables used by PIES.
!-----------------------------------------------------------------------
      MODULE torlines_background
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          nfp        Number of Field Periods
!          m          Number of Poloidal Modes
!          n          Number of Toroidal Modes
!          k          Number of Radial Surfaces
!          nu         Number of real space poloidal points
!          nv         Number of real space toroidal points
!          mnmax      Total number of Fourier Modes
!          xm         Poloidal Mode Number Array
!          xn         Toroidal Mode Number Array
!          iotamn     Minimum value of iota
!          iotamx     Maximum value of iota
!          rho        Radial Points Array [0,1]
!          xu         Poloidal Points Number Array [0,1]
!          xv         Toroidal Points Number Array [0,1]
!          rmnc       Even R Fourier Modes (cos)
!          rmns       Odd R Fourier Modes (sin)
!          zmnc       Even Z Fourier Modes (cos)
!          zmns       Odd Z Fourier Modes (sin)
!          bsmnc      Even B^s Fourier Modes
!          bsmns      Odd B^s Fourier Modes
!          bumnc      Even B^u Fourier Modes
!          bumns      Odd B^u Fourier Modes
!          bvmnc      Even B^v Fourier Modes
!          bvmns      Odd B^v Fourier Modes
!          jsmnc      Even J^s Fourier Modes
!          jsmns      Odd J^s Fourier Modes
!          jumnc      Even J^u Fourier Modes
!          jumns      Odd J^u Fourier Modes
!          jvmnc      Even J^v Fourier Modes
!          jvmns      Odd J^v Fourier Modes         
!----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8, parameter :: small = 1.e-10_ezspline_r8
      
      INTEGER ::  mnmax_m, nfp
      REAL(rprec) :: vc_adapt_tol, bound_separation, isgn
      REAL(rprec) :: u0, v0, dr_m, vb_ws, rb_ws, zb_ws, phimin, phimax, delta_phi
      INTEGER, ALLOCATABLE :: ixm_m(:), ixn_m(:)
      REAL(rprec), ALLOCATABLE :: rmnc_m(:),zmns_m(:)
      REAL(rprec), ALLOCATABLE :: rmns_m(:),zmnc_m(:)
      INTEGER, ALLOCATABLE :: xm_sav(:), xn_sav(:)
      REAL(rprec), ALLOCATABLE :: rmnc_sav(:,:),zmns_sav(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_sav(:,:),zmnc_sav(:,:)
      REAL(rprec), ALLOCATABLE :: bmnc_sav(:,:),bmns_sav(:,:)

      INTEGER                  :: win_R4D, win_Z4D, win_BXSI4D, &
                                  win_BETA4D, win_B4D, win_rho
      REAL(rprec), POINTER :: rho(:)
      REAL(rprec), DIMENSION(:,:,:,:), POINTER :: R4D, Z4D, BXSI4D, &
                                                  BETA4D, B4D
      
      REAL(rprec)              :: thmx, phmx
      TYPE(EZspline3_r8)       :: beta_spl
      TYPE(EZspline3_r8)       :: bxsi_spl
      TYPE(EZspline3_r8)       :: R_spl
      TYPE(EZspline3_r8)       :: Z_spl
      TYPE(EZspline3_r8)       :: B_spl

      ! Aliased VARS (did this so input namelist is clear
      INTEGER :: NZONE_EMC3, NRADI_EMC3, NPOLO_EMC3, NTORO_EMC3
      REAL(rprec) :: S0_EMC3, S1_EMC3
      ! EMC3 vars
      INTEGER :: NZONET, SRF_RADI, SRF_POLO, SRF_TORO
      REAL(rprec) :: s_inner, s_outer, s_inbg
      
      
      
      END MODULE torlines_background
