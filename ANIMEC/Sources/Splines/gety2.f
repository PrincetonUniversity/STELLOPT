      SUBROUTINE gety2(y, y2, h, nots, nb)
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nots, nb
      REAL(rprec), DIMENSION(*) :: y, y2, h
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER jmax
      REAL(rprec), DIMENSION(nots) :: aspline, bspline, dspline
C-----------------------------------------------

      jspmin(1) = 2
      IF (nb .eq. ideriv) jspmin(1) = 1
      jmax = nots - 1
      aspline(1) = h(1)
      dspline(1) = 2.0*h(1)
      y2(1) = 0.
      y2(nots) = 0.
      IF (nb .eq. ideriv) y2(1) = 6.0*(y(2)-y(1))/h(1)
      aspline(2:jmax) = h(2:jmax)
      bspline(2:jmax) = h(:jmax-1)
      dspline(2:jmax) = 2.0*(h(2:jmax)+h(:jmax-1))
      y2(2:jmax) = 6.0*((y(3:jmax+1)-y(2:jmax))/h(2:jmax)
     1   -(y(2:jmax)-y(:jmax-1))/h(:jmax-1))

      CALL tridslv(aspline,dspline,bspline,y2,jspmin,jmax,0,nots,1)

      END SUBROUTINE gety2
