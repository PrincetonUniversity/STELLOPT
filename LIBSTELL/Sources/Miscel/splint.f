      SUBROUTINE splint(xa, ya, y2a, n, x, y)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER     :: n
      REAL(rprec) :: xa(n), ya(n), y2a(n), x, y
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: klo, khi, i, k
      REAL(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi
C-----------------------------------------------
!
!       SPLINE INTERPOLATION ROUTINE (Numerical Recipes, pg. 89)
!       XA: ordered array of length N of ordinates at which function YA=F(XA)
!           is tabulated
!       YA: array of length N , = F(XA)
!       Y2A: array of second derivatives at XA points
!       computed from call to SPLINE
!       X : value at which Y = F(X) is to be computed from splines
!       YP = dY/dX at X
!       NDIM: dimension of X, Y, YP arrays


      klo = 1
      khi = n
      DO WHILE(khi - klo .gt. 1)
         k = (khi + klo)/2
         IF (xa(k) .gt. x) THEN
            khi = k
         ELSE
            klo = k
         ENDIF
      END DO

      h = xa(khi) - xa(klo)
      a = xa(khi) - x
      b = x - xa(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2a(klo)
      y26hi = c1o6*y2a(khi)
      y     = (a*(ya(klo)+(a2-h2)*y26lo)+b*(ya(khi)+(b2-h2)*y26hi))/h

      END SUBROUTINE splint
