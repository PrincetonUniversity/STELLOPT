      SUBROUTINE setspline(x,weight,y,h,yfit,y2,wten,tens,nots,nb)
      USE vspline
      USE vparams, ONLY: zero
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nots, nb
      REAL(rprec) tens
      REAL(rprec), DIMENSION(*) :: x, weight
      REAL(rprec), DIMENSION(nots) :: y
      REAL(rprec), DIMENSION(*) :: h
      REAL(rprec), DIMENSION(nots) :: yfit
      REAL(rprec), DIMENSION(*) :: y2, wten
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ioff, info
      REAL(rprec) :: alsq(nots,nots)
C-----------------------------------------------
!
!       x:        independent coordinate array
!       y:        dependent y(x) array
!       yfit:     fitted values (under tension) to y array
!       h:        x(i+1) - x(i) array
!       y2:       y'' array used for splines
!       wten:        weight array for tension (changed on EXIT)
!       alsq:        matrix elements for least squares fit (from s-integrations)
!       nots:        number of independent coordinates (knots)
!       nb:        = NATUR, USE natural boundary condition at left knot
!                  = IDERIV, USE derivative (dy/dx =0) boundary condition at left knot

!
!       IT IS ASSUMED THAT X,Y,WTEN ARE ALL SORTED (ON X(I) < X(I+1))
!

!
!       INITIALIZE ALSQ TO ZERO, COMPUTE H ELEMENTS
!
      CALL initspline (alsq, x, h, weight, nots)

!
!       SET UP SPLINE MATRIX ASPLINE AND NON-DIMENSIONLIZE TENSION
!
      ioff = 0
      CALL add_tension (alsq, wten, h, tens, zero, zero, nots, nb,
     1   ioff, nots)

!
!       SOLVE FOR COEFFICIENTS
!
      yfit(:nots) = y(:nots)
      CALL solver (alsq, yfit, nots, 1, info)
!
!       OBTAIN Y'' COEFFICIENTS AND STORE IN Y2
!
      CALL gety2 (yfit, y2, h, nots, nb)

      END SUBROUTINE setspline
