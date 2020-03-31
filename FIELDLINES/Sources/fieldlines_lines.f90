!-----------------------------------------------------------------------
!     Module:        fieldlines_lines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the FIELDLINES field line
!                    variables.
!-----------------------------------------------------------------------
      MODULE fieldlines_lines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!          myline    Dummy index
!          nlines    Number of Fieldlines
!          nsteps    Number of integration steps along fieldline
!          R_lines   Radial locations along fieldline [m] (npoinc per field period)
!          Z_lines   Vertical locations along field line [m]
!          PHI_lines Toroidal locations along field line [radians]
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  :: nlines, nsteps, myline, myldex, myend

      REAL(rprec) :: xlast,ylast,zlast ! for storing position
      REAL(rprec), ALLOCATABLE :: R_lines(:,:),Z_lines(:,:),PHI_lines(:,:),&
                                  B_lines(:,:),L_lines(:)
      
      REAL(rprec), ALLOCATABLE :: Rhc_lines(:,:),Zhc_lines(:,:)

      END MODULE fieldlines_lines
