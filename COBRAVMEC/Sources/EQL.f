      SUBROUTINE EQL(X, F, G, H, N)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER N
      REAL(rprec), DIMENSION(N) :: X, F, H
      REAL(rprec), DIMENSION(N) :: G
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I1, J2, I, J, K, L, M, NM1
      REAL(rprec) :: R, S, Z, U
!-----------------------------------------------
      Z = 0
      I = 1
      R = X(N)
      S = F(N)
      U = G(N)
      IF (S .EQ. Z) GO TO 30
      I1 = I
      J2 = MAX(N - 1,I1)
      DO I = I1, J2
         IF (F(I) .EQ. Z) CYCLE
         IF (U < G(I)) GO TO 30
         IF (U > G(I)) CYCLE
         IF (ABS(F(I)) .GE. ABS(S)) GO TO 30
      END DO
      GO TO 50
   30 CONTINUE
      M = N + I
      NM1 = N - 1
      DO J = I, NM1
         K = M - J
         L = K - 1
         X(K) = X(L)
         F(K) = F(L)
         G(K) = G(L)
      END DO
      X(I) = R
      F(I) = S
      G(I) = U
   50 CONTINUE
      U = G(N)
      WHERE ((G(:N)-U).LE.(-99._dp) .OR. F(:N).EQ.Z) H(:N) = Z
      WHERE (F(:N).NE.Z .AND. (G(:N)-U)>(-99._dp))
     1    H(:N) = F(:N)*16._dp**(G(:N)-U)

      RETURN
      END SUBROUTINE EQL
