!=====================================================
SUBROUTINE eva2d(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates a 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval               evaluated function value

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp),                       INTENT(out) :: spval

  REAL(kind=dp), DIMENSION(4)                      :: a
  INTEGER                                          :: l

  DO l=1,4
     a(l) = spl(1,l,ix,iy) + dx*(spl(2,l,ix,iy) +              &
          dx*(spl(3,l,ix,iy) + dx* spl(4,l,ix,iy)))
  END DO
  spval = a(1)+dy*(a(2)+dy*(a(3)+dy*a(4)))

  RETURN
END SUBROUTINE eva2d
