      SUBROUTINE multxy(x,y,n)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: n
      REAL(rprec), DIMENSION(n) :: x, y
c-----------------------------------------------

      x = x * y

      END SUBROUTINE multxy
