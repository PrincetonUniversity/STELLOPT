      FUNCTION splints (x)
      USE vspline
      USE vmec_input, ONLY: isnodes
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) x
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: klo, khi, k
      REAL(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi, yx, splints
C-----------------------------------------------

      klo = 1
      khi = isnodes

    1 CONTINUE
      IF (khi - klo .gt. 1) THEN
         k = (khi + klo)/2
         IF (sknots(k) .gt. x) THEN
            khi = k
         ELSE
            klo = k
         ENDIF
         GOTO 1
      ENDIF

      h = sknots(khi) - sknots(klo)
      a = sknots(khi) - x
      b = x - sknots(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2stark(klo)
      y26hi = c1o6*y2stark(khi)
      yx = (a*(ystark(klo)+(a2-h2)*y26lo)+b*(ystark(khi)+(b2-h2)*y26hi))
     1   /h
      splints = yx

      END FUNCTION splints
