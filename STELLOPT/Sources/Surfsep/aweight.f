      FUNCTION aweight(x, y, a, b, v, w)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: x, y, a, b, w, v
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: pi, aweight, n12, n4
      REAL(rprec), PARAMETER :: n3 = -1.e10_dp
C-----------------------------------------------
      Pi=4*ATAN(1._dp)
      n12=(-((-a + x)**2/v**2) - (-b + y)**2/w**2)
      IF(n12 > n3)then
         n4=EXP(n12)
      ELSE
         aweight=0
         RETURN
      ENDIF
      aweight = n4

      END FUNCTION aweight
