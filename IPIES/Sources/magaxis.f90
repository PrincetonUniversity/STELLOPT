!-----------------------------------------------------------------------
!     Subroutine:    magaxis
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/8/2011
!     Description:   This subroutine calculates the location of the
!                    magnetic axis by following field lines and looking
!                    for the one which ends on itself.
!-----------------------------------------------------------------------
      SUBROUTINE magaxis
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          ierr        Error flag
!          xlb         X Lower Bound
!          xub         X upper Bound
!          x0          Initial guess for xaxis
!          z0          Initial guess for zaxis
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,u,v
      REAL    :: pi2, xlb, xub, x0, z0
!-----------------------------------------------------------------------
!     External Functions
!-----------------------------------------------------------------------
      EXTERNAL faxis
      EXTERNAL C05ADF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      x0 = 1e-6
      ! If symmetry use bisection method
      IF (lasym) THEN
      ELSE
      	xlb = -0.6
      	xub =  0.6
      	ier = 1
         CALL C05ADF(xlb, xub, magaxis_abs, magaxis_eta, faxis, x0, ier)
         IF (ier > 1) CALL handle_err(D05ADF_ERR,'magaxis',ier)
         IF (ier == 0) THEN
            r_axis = x0
            z_axis = 0.0
         END IF
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE magaxis
