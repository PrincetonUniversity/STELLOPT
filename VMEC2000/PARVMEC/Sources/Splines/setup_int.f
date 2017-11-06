      SUBROUTINE setup_int(xknots,smesh,hx,w,w1,u,u1,nk,nots,nmesh)
      USE stel_kinds
      USE vparams, ONLY: epstan
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nots, nmesh
      INTEGER, DIMENSION(nots) :: nk
      REAL(rprec), DIMENSION(nots) :: xknots
      REAL(rprec), DIMENSION(nmesh) :: smesh
      REAL(rprec), DIMENSION(nots) :: hx
      REAL(rprec), DIMENSION(nmesh) :: w, w1, u, u1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ksp1, k, i, k1
      REAL(rprec) :: smesh1, hk6
C-----------------------------------------------
!
!       FOR THE SPLINED FUNCTIONS CORRESPONDING TO XKNOTS (PRESSURE,
!       IOTA), THIS ROUTINE COMPUTES THE A AND B MATRIX ELEMENTS
!       (STORED IN W,U,W1,U1) WHICH ARE NEEDED TO EVALUATE THE FUNCTIONS
!       IN REAL-SPACE IN TERMS OF THEIR (VARIABLE) SPLINE KNOT VALUES.
!       THIS 'UNDOES' WHAT SPLINT ROUTINE DOES. LET Y(I) DENOTE THE
!       FUNCTION AT THE POINT SMESH(I) SUCH THAT
!
!                   XKNOTS(K) < SMESH(I) <= XKNOTS(K)
!
!       THEN,  Y(I) = W(I)*YK  + U(I)*GK  + W1(I)*YK1   + U1(I)*GK1
!
!       WHERE YK, GK ARE THE SPLINE AND 2ND DERIVATIVES AT KNOT K
!             YK1,GK1 ARE THE SAME AT KNOT K+1
!
      ksp1 = nots - 1
      smesh1 = smesh(1)
      IF (smesh1 .le. xknots(1)) smesh(1) = xknots(1) + epstan

      nk = 0

      k = 1
      DO i = 1, nmesh
  140    CONTINUE
         k1 = k + 1
!
!       XKNOTS = SQRT(HS*(JS-1)) DEFINED IN STARK,PRESSURE ROUTINE
!       (THIS CORRESPONDS TO APPROXIMATELY EQUAL SPACING ALONG MIDPLANE)
!
         IF (smesh(i).gt.xknots(k) .and. smesh(i).le.xknots(k1)) THEN
            nk(k) = nk(k) + 1
            hk6 = hx(k)*hx(k)/6.0
            w1(i) = (smesh(i)-xknots(k))/hx(k)
            IF (w1(i)<(-epstan) .or. w1(i)>1.0+epstan) STOP 'w1(i)'
            w(i) = 1.0 - w1(i)
            u(i) = hk6*w(i)*(w(i)*w(i)-1.0)
            u1(i) = hk6*w1(i)*(w1(i)*w1(i)-1.0)
         ELSE
            k = k + 1
            IF (k .gt. ksp1) STOP 'K>KSP1'
            GOTO 140
         ENDIF
      END DO

      smesh(1) = smesh1

      END SUBROUTINE setup_int
