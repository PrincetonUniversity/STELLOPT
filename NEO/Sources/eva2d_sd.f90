!=====================================================
SUBROUTINE eva2d_sd(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates the second derivatives of 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval(3)            evaluated function values
!                             spval(1) = d^2f/dx^2
!                             spval(2) = d^2f/(dxdy)
!                             spval(3) = d^2f/dy^2

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp), DIMENSION(3),         INTENT(out) :: spval

  INTEGER                                          :: i,j
  REAL(kind=dp)                                    :: muli, mulj

  spval = 0

! d^2f/dx^2
  DO i=3,4
     IF (i == 3) THEN
        muli = 1
     ELSE
        muli = dx**(i-3)
     END IF
     muli = muli * (i-1) * (i-2)
     DO j=1,4
        IF (j == 1) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-1)
        END IF
        spval(1) = spval(1) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! d^2f/(dxdy)
  DO i=2,4
     IF (i == 2) THEN
        muli = 1
     ELSE
        muli = dx**(i-2)
     END IF
     muli = muli * (i-1)
     DO j=2,4
        IF (j == 2) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-2)
        END IF
        mulj = mulj * (j-1)
        spval(2) = spval(2) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! d^2f/dy^2
  DO i=1,4
     IF (i == 1) THEN
        muli = 1
     ELSE
        muli = dx**(i-1)
     END IF
     DO j=3,4
        IF (j == 3) THEN
           mulj = 1
        ELSE
           mulj = dy**(j-3)
        END IF
        mulj = mulj * (j-1) * (j-2)
        spval(3) = spval(3) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

  RETURN
END SUBROUTINE eva2d_sd
