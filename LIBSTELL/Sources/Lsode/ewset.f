      SUBROUTINE ewset(n, itol, rtol, atol, ycur, ewt)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n, itol
      REAL(rprec), DIMENSION(*) :: rtol, atol
      REAL(rprec), DIMENSION(n) :: ycur, ewt
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c this subroutine sets the error weight vector ewt according to
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
c with the subscript on rtol and/or atol possibly replaced by 1 above,
c depending on the value of itol.
c-----------------------------------------------------------------------
c
      SELECT CASE (itol)
      CASE DEFAULT
         ewt = rtol(1)*ABS(ycur) + atol(1)
         RETURN
      CASE (2)
         ewt = rtol(1)*ABS(ycur) + atol(:n)
         RETURN
      CASE (3)
         ewt = rtol(:n)*ABS(ycur) + atol(1)
         RETURN
      CASE (4)
         ewt = rtol(:n)*ABS(ycur) + atol(:n)
         RETURN
      END SELECT

      END SUBROUTINE ewset
