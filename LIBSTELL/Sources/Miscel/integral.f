      FUNCTION integral (n, x, y1, y2)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL(rprec), INTENT(IN), DIMENSION(n) :: y1, y2 , x
      REAL(rprec), PARAMETER :: one = 1
      REAL(rprec) :: h, integral
!------------------------------------------------------------------------------
!      WARNING: The integral formula assumes values on the half MESH,
!      that is, the independent variable is given from x(3/2) to x(M-1/2),
!      but the limits of the integral are x(1)=-1 and x(M)=1. N = M-1 is the
!      number of points on the half mesh.
!      In the current version, fulfilment of this condition is ONLY checked
!      for the case in which the interval of integration is the interval
!      of orthogonality of Legendre polynomias, i.e., [-1,1].
!      Integration formula is a 2-nd order HALF-MESH formula from
!      "Numerical Recipes", W.Press et al, Chapter 4, page 110.
!------------------------------------------------------------------------------
      IF (n < 10) STOP 'Too few points to carry out integration.!'
      IF (x(n) < x(1)) STOP ' B < A in INTEGRAL!'
      IF (x(n) == one .or. x(1) == -one) STOP 'HALF MESH INTEGRAL!'

      h = (x(n)-x(1))/(n-1)
      integral = h* SUM(y1*y2)

      END FUNCTION integral
