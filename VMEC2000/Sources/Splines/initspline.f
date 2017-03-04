      SUBROUTINE initspline(amat, splnot, h, weight, nots)
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nots
      REAL(rprec), DIMENSION(nots,nots) :: amat
      REAL(rprec), DIMENSION(*) :: splnot, h, weight
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER i
      REAL(rprec) :: eps
C-----------------------------------------------

      IF (nots .lt. 3) STOP 'nots<3'
      eps = 1.0/(splnot(nots)-splnot(1))

      amat = 0.
      DO i = 1, nots
         amat(i,i) = weight(i)
      END DO

      DO i = 1, nots - 1
         h(i) = splnot(i+1) - splnot(i)
         IF (eps*h(i) .le. 1.e-8_dp) STOP 'h(i)<1.e-8'
      END DO

      END SUBROUTINE initspline
