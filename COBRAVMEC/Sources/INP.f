      SUBROUTINE INP(A, X, Y, Z, U, V, W, L, R)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) :: A, X, Y, Z, U, V, W, L, R
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: B, C, P, Q, S, T
!-----------------------------------------------
      S = Z - X
      T = (Y - X)/S
      A = (U - V)/T + (W - V)/(1. - T)
      IF (A .NE. ZERO) THEN
         B = .5_dp*(W - U)/A - .5_dp
         C = U/A
         T = SQRT(ABS(C))
         IF (ABS(B) < SIGN(T,C)) GO TO 60
         T = MAX(T,ABS(B))
         IF (T .EQ. ZERO) GO TO 50
         Q = ONE/T
         P = SQRT((Q*B)**2 - Q*C*Q)
         P = T*P
         IF (ABS(P + B) .LE. ABS(P - B)) THEN
            Q = P - B
         ELSE
            Q = -(B + P)
         ENDIF
         P = C/Q
         Q = X + S*Q
         P = X + S*P
         IF (Q .GE. L) THEN
            IF (Q .LE. R) THEN
               A = Q
               RETURN
            ENDIF
         ENDIF
         A = P
         RETURN
      ENDIF
      IF (U .NE. W) THEN
         A = X + S*U/(U - W)
         RETURN
      ENDIF
   50 CONTINUE
      A = L
      RETURN
   60 CONTINUE
      A = X - S*B
      RETURN

      END SUBROUTINE INP
