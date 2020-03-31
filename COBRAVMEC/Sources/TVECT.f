      SUBROUTINE TVECT(E, X, L, D, U, N, W)
!      ________________________________________________________
!     |                                                        |
!     |     COMPUTE EIGENVECTOR CORRESPONDING TO GIVEN REAL    |
!     |        EIGENVALUE FOR A REAL TRIDIAGONAL MATRIX        |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |         L     --SUBDIAGONAL                            |
!     |                                                        |
!     |         D     --DIAGONAL                               |
!     |                                                        |
!     |         U     --SUPERDIAGONAL                          |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |         W     --WORK ARRAY (LENGTH AT LEAST 4N)        |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --IMPROVED ESTIMATE FOR EIGENVALUE       |
!     |                                                        |
!     |         X     --EIGENVECTOR                            |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,SQRT                         |
!     |________________________________________________________|
!
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N
      REAL(rprec):: E
      REAL(rprec), DIMENSION(N) :: X, L, D, U
      REAL(rprec), DIMENSION(4*N) :: W
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: G1, J1, J2, J4, J3, F, G, H, I, J, K, M
      REAL(rprec) :: O, P, Q, R, S, T, V, Y, Z
!-----------------------------------------------
      IF (N .LE. 1) THEN
         E = D(1)
         X(1) = ONE
         RETURN
      ENDIF
      M = N - 1
      J = 2
!     ---------------------------------------------------------
!     |*** STORE MATRIX IN W ARRAY AND SUBTRACT EIGENVALUE ***|
!     ---------------------------------------------------------
      X(1:M) = 0
      DO I = 1,M
         W(J) = U(I)
         W(J+1) = D(I) - E
         W(J+2) = L(I)
         J = J + 4
      END DO
      W(J) = 0
      W(J+1) = D(N) - E
      O = 65536
      O = O**(-4)
      T = O/2
      S = T
      T = T/2
      P = S
      S = S + T
      IF (S .GE. O) GO TO 40
      DO WHILE(S + T > S)
         T = T/2
         P = S
         S = S + T
         IF (S .GE. O) GO TO 40
      END DO
   40 CONTINUE
      R = W(3)
      V = ABS(R) + ABS(W(4))
      F = 4
      M = 4*N - 3
      G = -2
!     ---------------------------
!     |*** FACTOR THE MATRIX ***|
!     ---------------------------
      G = G + 4
      G1 = G
      DO G = G1, M, 4
         H = G - 1
         I = G + 2
         J = G + 5
!     --------------------
!     |*** FIND PIVOT ***|
!     --------------------
         Q = W(I)
         Y = ABS(Q)
         Z = ABS(R)
         IF (Z < Y) THEN
!     -------------------
!     |*** SWAP ROWS ***|
!     -------------------
            IF (V < Y) THEN
               V = Y
               F = I
            ENDIF
            T = W(G)
            W(G) = W(J)
            W(J) = T
            T = R/Q
            K = G + 1
            W(K) = Q
            K = J - 1
            S = W(K)
            IF (S .NE. ZERO) THEN
               IF (S .EQ. O) S = P
               W(K) = -S*T
               W(H) = S
               GO TO 100
            ENDIF
            W(K) = S
            W(H) = O
            GO TO 100
         ENDIF
         W(H) = 0
         IF (V < Z) THEN
            V = Z
            F = I
         ENDIF
         IF (R .EQ. ZERO) GO TO 120
         T = Q/R
!     -------------------
!     |*** ELIMINATE ***|
!     -------------------
  100    CONTINUE
         R = W(J) - T*W(G)
         W(J) = R
         W(I) = T
      END DO
      IF (ABS(R) < V) THEN
         V = R
         F = J + 1
!     ---------------------------------------------------
!     |*** COMPUTE INITIAL EIGENVECTOR APPROXIMATION ***|
!     ---------------------------------------------------
      ENDIF
  120 CONTINUE
      J = F/4
      X(J) = ONE
      IF (J .NE. 1) THEN
         K = F - 5
         J = J - 1
         X(J) = (X(J)-W(K-1)*X(J+1))/W(K)
         J1 = J
         IF (J1 - 1 > 0) THEN
            DO J = J1, 2, -1
               K = K - 4
               T = W(K-2)
               IF (T .EQ. O) T = 0
               X(J-1) = (X(J-1)-W(K-1)*X(J)-T*X(J+1))/W(K)
            END DO
         ENDIF
         GO TO 140
      ENDIF
  140 CONTINUE
      IF (V .NE. ZERO) THEN
         S = MAXVAL(ABS(X(:N)))
         S = ONE/S
         X(:N) = S*X(:N)
         R = SUM(X(:N)*X(:N))
         K = 0
         J = 1
         Y = X(1)
!     -----------------------------------------------------
!     |*** APPLY ONE ITERATION OF INVERSE POWER METHOD ***|
!     -----------------------------------------------------
         J2 = J
         J4 = MAX(N - 1,J2)
         DO J = J2, J4
            K = K + 4
            I = J
            S = W(K)
            W(K) = Y
            Y = X(J+1)
            IF (W(K-3) .NE. ZERO) THEN
               T = X(J+1)
               X(J+1) = X(I)
               X(I) = T
            ENDIF
            X(J+1) = X(J+1) - S*X(I)
         END DO
!     ---------------------------
!     |*** BACK SUBSTITUTION ***|
!     ---------------------------
         S = X(J)/W(K+3)
         X(J) = S
         T = ABS(S)
         V = S*Y
         J = J - 1
         K = K - 1
         S = (X(J)-W(K-1)*S)/W(K)
         X(J) = S
         T = MAX(ABS(S),T)
         V = V + S*W(K+1)
         J3 = J
         IF (J3 - 1 > 0) THEN
            DO J = J3, 2, -1
               K = K - 4
               Z = W(K-2)
               IF (Z .EQ. O) Z = 0
               S = (X(J-1)-W(K-1)*S-Z*X(J+1))/W(K)
               X(J-1) = S
               T = MAX(ABS(S),T)
               V = V + S*W(K+1)
            END DO
         ENDIF
         IF (V .NE. ZERO) V = R/V
         T = ONE/T
         S = SUM((X(:N)*V)**2)
         Z = SUM((T*X(:N))**2)
         T = T/SQRT(Z)
         X(:N) = T*X(:N)
!     --------------------------------------------------------------
!     |*** USE RAYLEIGH QUOTIENT TO IMPROVE EIGENVALUE ESTIMATE ***|
!     --------------------------------------------------------------
         IF (R + R .GE. S) E = E + V
         RETURN
      ENDIF
      T = MAXVAL(ABS(X(:N)))
      T = ONE/T
      Z = SUM((T*X(:N))**2)
      T = T/SQRT(Z)
      X(:N) = T*X(:N)

      END SUBROUTINE TVECT
