
SUBROUTINE neo_zeros2d(x, y, eps, iter_ma, iter, error)

  USE neo_precision

  IMPLICIT NONE

  INTEGER,        INTENT(out)   :: error, iter
  INTEGER,        INTENT(in)    :: iter_ma
  REAL (kind=dp), INTENT(in)    :: eps
  REAL (kind=dp), INTENT(inout) :: x, y

  REAL (kind=dp) :: x_n, y_n
  REAL (kind=dp) :: f, dfdx, dfdy, g, dgdx, dgdy
  REAL (kind=dp) :: f_n,g_n
  REAL (kind=dp) :: det
  REAL (kind=dp) :: x_err, y_err

  error = 0

! compute f(x,y), g(x,y) and all first derivatives
  CALL neo_bderiv(x, y, f, g, dfdx, dfdy, dgdx, dgdy)
  DO iter = 1, iter_ma

     det = dfdx * dgdy - dfdy * dgdx
!!$     PRINT *, 'f,    g    ', f, g
!!$     PRINT *, 'dfdx, dgdx ', dfdx, dgdx
!!$     PRINT *, 'dfdy, dgdy ', dfdy, dgdy
!!$     PRINT *, 'det        ', det
     x_n = x + ( dfdy *  g   -  f   * dgdy ) / det
     y_n = y + (  f   * dgdx - dfdx *  g   ) / det
!!$     PRINT *, 'x,    y    ', x, y
!!$     PRINT *, 'x_n,  y_n  ', x_n, y_n

!    compute f(x,y), g(x,y) and all first derivatives
!    at the new positions
!     PRINT *, x_n, y_n
     CALL neo_bderiv(x_n, y_n, f_n, g_n, dfdx, dfdy, dgdx, dgdy)

!    remaining relatve errors
!     IF (x_n .NE. ZERO) THEN
     IF (ABS(x_n) .GT. eps) THEN
       x_err = ABS ( (x_n - x) / x_n )
     ELSE
       x_err = ABS ( x_n - x )
     END IF
!     IF (y_n .NE. ZERO) THEN
     IF (ABS(y_n) .GT. eps) THEN
       y_err = ABS ( (y_n - y) / y_n )
     ELSE
       y_err = ABS ( y_n - y )
     END IF

!    new values
     f = f_n
     g = g_n
     x = x_n
     y = y_n

!    exit if error is small enough
     IF ( MAX ( x_err, y_err ) < eps ) RETURN

  END DO

  error = 1

  RETURN
END SUBROUTINE neo_zeros2d
