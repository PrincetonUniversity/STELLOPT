!-----------------------------------------------------------------------
!     Function:      faxis_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/30/2012
!     Description:   This subroutine calculates the distance between
!                    the starting and ending point of a field line.
!-----------------------------------------------------------------------
      SUBROUTINE faxis_nag(iflag, m, n, xc, fvecc, iw, liw, w, lw)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid, ONLY: phimin, phimax
!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iflag, m, n, iw(liw), liw, lw
      DOUBLE PRECISION    :: xc(n), fvecc(m), w(lw)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      REAL(rprec) :: x0, phi0, z0, x1, phi1, z1
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      iflag = 1
      x0   = xc(1)
      z0   = xc(2)
      phi0 = phimin
      phi1 = phimax
      x1   = x0
      z1   = z0
      PRINT *,'     faxis_nag(0):  x0,phi0,z0  ', x0,phi0,z0
      PRINT *,'     faxis_nag(0):  x1,phi0,z1,phi1  ', x1,phi0,z1,phi1
      CALL follow_single(x1,phi0,z1,phi1)
      PRINT *,'     faxis_nag(1):  x1,phi0,z1,phi1  ', x1,phi0,z1,phi1
      fvecc(1) = x1 - x0
      fvecc(2) = z1 - z0
      PRINT *,'               fvecc ',fvecc

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE faxis_nag
