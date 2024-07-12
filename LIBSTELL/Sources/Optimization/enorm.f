      FUNCTION enorm (n, x)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      REAL(rprec), DIMENSION(n) :: x
      REAL(rprec) :: enorm
!-----------------------------------------------
!
!     FUNCTION enorm
!
!     given an n-vector x, this FUNCTION calculates the
!     euclidean norm of x.
!
!
!     the FUNCTION statement is
!
!       FUNCTION enorm(n,x)
!
!     WHERE
!
!       n is a positive INTEGER input variable.
!
!       x is an input array of length n.
!
      enorm = SQRT (SUM(x(:n)*x(:n)))

      END FUNCTION enorm
