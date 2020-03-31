      FUNCTION pythag (a, b)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) a, b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: absa, absb, pythag
C-----------------------------------------------
c computes SQRT(a^2+b^2) without destructive underflow or overflow
      absa = ABS(a)
      absb = ABS(b)
      IF (absa .gt. absb) THEN
         pythag = absa*SQRT(one + (absb/absa)**2)
      ELSE
         IF (absb .eq. zero) THEN
            pythag = zero
         ELSE
            pythag = absb*SQRT(one + (absa/absb)**2)
         END IF
      END IF

      END FUNCTION pythag
