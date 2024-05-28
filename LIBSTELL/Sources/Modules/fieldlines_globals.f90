!-----------------------------------------------------------------------
!     Module:        fieldlines_globals
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          05/28/2024
!     Description:   This module contains the FIELDLINES global 
!                    variables needed by the input namelist.
!-----------------------------------------------------------------------
      MODULE fieldlines_globals
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE

      ! From fiellines_runtime
      INTEGER, PARAMETER ::  MAXLINES   = 2**19
      LOGICAL :: lerror_field
      INTEGER :: npoinc, num_hcp
      REAL(rprec) :: dphi, follow_tol, delta_hc, mu
      REAL(rprec), DIMENSION(20)           :: errorfield_amp, errorfield_phase
      REAL(rprec), DIMENSION(MAXLINES)     :: r_start, phi_start, &
                                              z_start, phi_end, &
                                              r_hc, z_hc, phi_hc
      CHARACTER(256) :: int_type


      ! From fieldlines_lines
      INTEGER :: nlines
      ! From fieldlines_grid
      INTEGER :: nr, nphi, nz
      REAL(rprec) :: rmin, rmax, zmin, zmax, phimin, phimax, &
                     vc_adapt_tol


      CONTAINS

      ! These expose the global variables through ctypes
      INTEGER FUNCTION getmaxlines()
      IMPLICIT NONE
      getmaxlines = MAXLINES
      END FUNCTION getmaxlines

      END MODULE fieldlines_globals
