      FUNCTION splintx(x)
      USE vparams
      USE vsvd0
      USE csplinx
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
      REAL(rprec) :: h, a, b, h2, a2, b2, y26lo, y26hi, qmidx0,
     1   splintx
C-----------------------------------------------

      CALL setspline (rmidx, wmidx, qmidx, hmidx, ymidx, y2midx,
     1    tenmidx, tenmidx(1), nptsx, natur)

      klo = 1
      khi = nptsx

    1 CONTINUE
      IF (khi - klo .gt. 1) THEN
         k = (khi + klo)/2
         IF (rmidx(k) .gt. x) THEN
            khi = k
         ELSE
            klo = k
         ENDIF
         GOTO 1
      ENDIF

      h = rmidx(khi) - rmidx(klo)
      IF( h.eq.zero )then
        splintx = zero
        RETURN
      END IF
      a = rmidx(khi) - x
      b = x - rmidx(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2midx(klo)
      y26hi = c1o6*y2midx(khi)
      qmidx0 = (a*(ymidx(klo)+(a2-h2)*y26lo)+b*(ymidx(khi)+
     1   (b2-h2)*y26hi))/h
      splintx = qmidx0

      END FUNCTION splintx
