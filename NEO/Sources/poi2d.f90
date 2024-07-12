!=====================================================
SUBROUTINE poi2d(hx,hy,mx,my,                         &
     xmin,xmax,ymin,ymax,                             &
     x,y,ix,iy,dx,dy,ierr)
! Creates Pointers for eva2d
!
! Input:  hx, hy              increment in x and y
!         mx, my              standard (0) or periodic (1) spline
!         xmin, xmax          Minimum and maximum x
!         ymin, ymax          Minimum and maximum y
!         x, y                x and y values for spline avaluation
! Output: spval               evaluated function value
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         ierr                error (> 0)

  USE neo_precision

  IMPLICIT NONE

  REAL(kind=dp),                       INTENT(in)  :: hx, hy
  INTEGER,                             INTENT(in)  :: mx, my
  REAL(kind=dp),                       INTENT(in)  :: xmin, xmax, ymin, ymax
  REAL(kind=dp),                       INTENT(in)  :: x, y

  INTEGER,                             INTENT(out) :: ix, iy
  REAL(kind=dp),                       INTENT(out) :: dx, dy
  INTEGER,                             INTENT(out) :: ierr

  REAL(kind=dp)                                    :: dxx, x1, dyy, y1
  REAL(kind=dp)                                    :: dxmax, dymax

  ierr = 0

  dxx = x-xmin
  IF (mx .EQ. 0) THEN
     IF (dxx .LT. ZERO) THEN
        ierr = 1
        RETURN
     END IF
     IF (x .GT. xmax) THEN
        ierr = 2
        RETURN
     END IF
  ELSE
     dxmax = xmax - xmin
     IF(dxx .LT. ZERO) THEN
        dxx = dxx+(1+INT(abs(dxx/dxmax)))*dxmax
     ELSE IF(dxx .GT. dxmax) THEN
        dxx = dxx-(INT(abs(dxx/dxmax)))*dxmax
     END IF
  END IF
  x1 = dxx/hx
  ix = INT(x1)
  dx = hx*(x1-ix)
  ix = ix+1

  dyy = y-ymin
  IF (my .EQ. 0) THEN
     IF (dyy .LT. ZERO) THEN
        ierr = 3
        RETURN
     END IF
     IF (y .GT. ymax) THEN
        ierr = 4
        RETURN
     END IF
  ELSE
     dymax = ymax - ymin
     IF(dyy .LT. ZERO) THEN
        dyy = dyy+(1+INT(abs(dyy/dymax)))*dymax
     ELSE IF(dyy .GT. dymax) THEN
        dyy = dyy-(INT(abs(dyy/dymax)))*dymax
     END IF
  END IF
  y1 = dyy/hy
  iy = INT(y1)
  dy = hy*(y1-iy)
  iy = iy+1

  RETURN
END SUBROUTINE poi2d
