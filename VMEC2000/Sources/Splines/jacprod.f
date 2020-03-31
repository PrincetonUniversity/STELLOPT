      SUBROUTINE jacprod(c, h, nots, nb)
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nots, nb, jmax
      REAL(rprec), DIMENSION(*) :: c, h
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nots) ::
     1   aspline, bspline, dspline, dum1
C-----------------------------------------------

!
!       THIS ROUTINE COMPUTES THE INNER PRODUCT COUT(I) = CIN(J)*JACOBIAN(J,I)
!       WHERE JACOBIAN(J,I) = D[G(J)]/D[F(I)]
!       HERE, G(J) ARE SECOND-DERIVATIVE KNOTS, F(I) FUNCTION KNOTS
!
!       COMPUTE COEFFICIENT ARRAY ELEMENTS A*X(I+1) + D*X(I) + B*X(I-1)
!       (TO BE SAFE, RECOMPUTE EACH TIME, SINCE IOTA, P SPLINES MAY
!        DIFFER FROM CALL TO CALL)
!
      aspline(1) = h(1)
      dspline(1) = 2.0*h(1)
      aspline(2:nots-1) = h(2:nots-1)
      bspline(2:nots-1) = h(:nots-2)
      dspline(2:nots-1) = 2.0*(h(2:nots-1)+h(:nots-2))

      jspmin(1) = 2
      IF (nb .eq. ideriv) jspmin(1) = 1
      jmax = nots - 1
      CALL tridslv(aspline,dspline,bspline,c,jspmin,jmax,0,nots,1)
      dum1(1) = 6.0*(c(2)-c(1))/h(1)
      dum1(2:nots) = 6.0*(c(:nots-1)-c(2:nots))/h(:nots-1)
      c(2:nots-1) = dum1(2:nots-1) - dum1(3:nots)
      c(1) = dum1(1)
      c(nots) = dum1(nots)

      END SUBROUTINE jacprod
