      SUBROUTINE CP(T, B, V, L, D, P, N)
!------------------------------------------------
!*** EVALUATE CHARACTERISTIC POLYNOMIAL AND ***
!              SIGN ALTERNATIONS
!------------------------------------------------
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: L, N
      REAL(rprec) ::  T, B
      REAL(rprec) :: V
      REAL(rprec), DIMENSION(N) :: D, P
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, M
      REAL(rprec) :: A, C, E, F, U, Z
!-----------------------------------------------
      L = 0
      Z = 0
      V = Z
      IF (N .LE. 1) THEN
         B = D(1) - T
         IF (B .LE. Z) L = 1
         RETURN
      ENDIF
      F = 65536
      F = F**4
      E = ONE/F
      U = 16
      M = 0
      I = 1
   20 CONTINUE
      C = ONE
      B = D(I) - T
   30 CONTINUE
      IF (B .LE. Z) GO TO 70
      IF (I .GE. N) GO TO 130
   40 CONTINUE
      J = I
      I = I + 1
      A = (D(I)-T)*B - P(J)*C
      C = B
      B = A
   50 CONTINUE
      A = ABS(B)
      IF (A > F) GO TO 60
      IF (A > E) GO TO 30
      IF (A .EQ. Z) GO TO 70
      C = C*F
      B = B*F
      V = V - U
      GO TO 50
   60 CONTINUE
      C = C*E
      B = B*E
      V = V + U
      GO TO 50
   70 CONTINUE
      L = L + 1
      IF (I .GE. N) GO TO 130
      IF (B < Z) GO TO 90
      IF (P(I) > ZERO) GO TO 90
      I = I + 1
      M = 1
      V = Z
      GO TO 20
   80 CONTINUE
      IF (B .GE. Z) GO TO 120
      IF (I .GE. N) GO TO 130
   90 CONTINUE
      J = I
      I = I + 1
      A = (D(I)-T)*B - P(J)*C
      C = B
      B = A
  100 CONTINUE
      A = ABS(B)
      IF (A > F) GO TO 110
      IF (A > E) GO TO 80
      IF (A .EQ. Z) GO TO 120
      C = C*F
      B = B*F
      V = V - U
      GO TO 100
  110 CONTINUE
      C = C*E
      B = B*E
      V = V + U
      GO TO 100
  120 CONTINUE
      L = L + 1
      IF (I .GE. N) GO TO 130
      IF (B > Z) GO TO 40
      IF (P(I) > ZERO) GO TO 40
      I = I + 1
      M = 1
      V = Z
      GO TO 20
  130 CONTINUE
      IF (M .NE. 1) THEN
         IF (B .NE. ZERO) THEN
            A = ONE/U
            IF (ABS(B) .GE. A) THEN
  140          CONTINUE
               IF (ABS(B) < ONE) RETURN
               B = B*A
               V = V + ONE
               GO TO 140
            ENDIF
            B = B*U
            V = V - ONE
 1003       CONTINUE
            IF (ABS(B) .GE. A) GO TO 1002
            B = B*U
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002
            B = B*U
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002
            B = B*U
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002
            B = B*U
            V = V - ONE
            GO TO 1003
 1002       CONTINUE
            RETURN
         ENDIF
      ENDIF
      V = Z
      B = Z
      RETURN

      END SUBROUTINE CP
