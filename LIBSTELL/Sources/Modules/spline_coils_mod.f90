!-----------------------------------------------------------------------
!     Module:        spline_coils_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/10/2025
!     Description:   This module contains routines for representing
!                    coils as splines in 3D space.
!-----------------------------------------------------------------------
      MODULE spline_coils_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Module SUBROUTINES/FUNCTIONS
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE init_splines(ns,k,s)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ns ! number of points to evaluate splines
      INTEGER, INTENT(IN) :: k  ! Order of splines
      DOUBLE PRECISION, DIMENSION(ns) :: s ! Where splines are evaluated
      END SUBROUTINE init_splines



      SUBROUTINE bspline_to_xyz(t,xr,xu,xv,s,x,y,z)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION()
      END SUBROUTINE bspline_to_ruv

      SUBROUTINE initialize_coils
      USE biotsavart
      IMPLICIT NONE
      ! We need to allocate coil_group from biotsavart
      ! 
      END SUBROUTINE initialize_coils

      
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE spline_coils_mod