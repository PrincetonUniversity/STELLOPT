!-----------------------------------------------------------------------
!     Subroutine:    current_density
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/14/2011
!     Description:   This subroutine calculates the total current
!                    density from the pressure and toroidal current.
!-----------------------------------------------------------------------
      SUBROUTINE current_density
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_runtime, ONLY: iter
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      ! Calculate Enclosed Toroidal flux
      CALL calc_metrics
      ! Calculate perpendicular current
      CALL jperp
      ! Calculate parallel current
      !CALL jpara
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE current_density
