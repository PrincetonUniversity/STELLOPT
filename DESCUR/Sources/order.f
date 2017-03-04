      SUBROUTINE order(rval, zval, xaxis, yaxis, inside)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inside
      REAL(rprec) xaxis, yaxis
      REAL(rprec), DIMENSION(*) :: rval, zval
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, ip1, i1, j, next
      REAL(rprec), DIMENSION(nu) :: tempr, tempz
      REAL(rprec) ::
     1   newdist,olddist,shortest,saver,savez,residue,x,y,dx,dy
C-----------------------------------------------
c**********************************************************************
c       Program ORDER : orders points on a magnetic surface at a
c       fixed toroidal plane and assigns right-handed circulation
c       around flux surface. XAXIS, YAXIS:  Polar-TYPE axis (must lie
c       inside curve to check sign of rotation)
c
c Modified 8/13/97, J. Breslau. No longer quits IF the axis is not
c  contained by the points. instead, returns inside=0 so inguess can
c  try again.
c**********************************************************************
      olddist = 1.E20_dp
      DO i = 1, ntheta - 1
         ip1 = i + 1
         i1 = i
         shortest = 1.E20_dp
   15    CONTINUE
         DO j = ip1, ntheta
            IF (i1 > 1) olddist = (rval(i1-1)-rval(j))**2 + (zval(i1-1)-
     1         zval(j))**2
            newdist = (rval(i1)-rval(j))**2 + (zval(i1)-zval(j))**2
            IF (newdist.le.olddist .and. newdist<shortest) THEN
               next = j
               shortest = newdist
            ENDIF
         END DO
c**********************************************************************
c       Swap nearest point (next) with current point (ip1)
c**********************************************************************
         IF (shortest .ge. 1.E10_dp) THEN
            SAVEr = rval(i-1)
            rval(i-1) = rval(i)
            rval(i) = SAVEr
            SAVEz = zval(i-1)
            zval(i-1) = zval(i)
            zval(i) = SAVEz
            i1 = i1 - 1
            ip1 = ip1 - 1
            GOTO 15
         ENDIF
         SAVEr = rval(ip1)
         rval(ip1) = rval(next)
         rval(next) = SAVEr
         SAVEz = zval(ip1)
         zval(ip1) = zval(next)
         zval(next) = SAVEz
      END DO
c**********************************************************************
c       Check that xaxis,yaxis is inside surface and
c       ascertain that the angle rotates counterclockwise
c       using Cauchy theorem in "complex"-plane
c**********************************************************************
      inside = 1
      residue = zero
      DO i = 1, ntheta - 1
         x = 0.5_dp*(rval(i)+rval(i+1)) - xaxis
         y = 0.5_dp*(zval(i)+zval(i+1)) - yaxis
         dx = rval(i+1) - rval(i)
         dy = zval(i+1) - zval(i)
         residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.E-10_dp)
      END DO
      x = 0.5_dp*(rval(1)+rval(ntheta)) - xaxis
      y = 0.5_dp*(zval(1)+zval(ntheta)) - yaxis
      dx = rval(1) - rval(ntheta)
      dy = zval(1) - zval(ntheta)
      residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.E-10_dp)

      IF (residue < (-0.9_dp*twopi)) THEN
         DO i = 2, ntheta
            j = ntheta - i + 2
            tempr(i) = rval(j)
            tempz(i) = zval(j)
         END DO
         rval(2:ntheta) = tempr(2:ntheta)
         zval(2:ntheta) = tempz(2:ntheta)
      ELSE IF (ABS(residue) < 0.9_dp*twopi) THEN
         PRINT *, ' mag. axis not enclosed by bndry; trying again'
         WRITE (3, *) ' mag. axis not enclosed by bndry; trying again'
         inside = 0
c        STOP
      ENDIF

      END SUBROUTINE order
