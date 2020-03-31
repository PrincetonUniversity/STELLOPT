      SUBROUTINE svdcmp(A, M, N, MP, NP, W, V)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: M, N, MP, NP
      REAL(rprec), DIMENSION(MP,NP), INTENT(inout), TARGET :: A
      REAL(rprec), DIMENSION(MIN(MP,NP)) :: W
      REAL(rprec), DIMENSION(NP,NP)      :: V
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: ZERO = 0, ONE = 1, TWO = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, L, K, J, ITS, NM
      INTEGER :: M1, N1
      REAL(rprec), DIMENSION(MP*NP) :: RV1
      REAL(rprec), POINTER          :: A1(:,:)
      REAL(rprec) :: G, SCALE, ANORM, S, F, H, C, Y, Z, X
      LOGICAL :: LTRANS
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PYTHAG
C-----------------------------------------------
!
!     COMPUTES THE SVD DECOMPOSITION OF THE MATRIX A SUCH
!     THAT A = U * W * VT, WHERE U=A ON OUTPUT. W IS A DIAGONAL
!     MATRIX OF SINGULAR VALUES.
!
!     TAKEN FROM NUMERICAL RECIPES. GENERALIZED BY SPH TO ACCEPT M < N CASE,
!     BY ANALYZING TRANSPOSE AND SWITCHING U,V AT END
!
      IF (M .GT. MP) STOP 'M > MP IN SVDCMP'
      IF (N .GT. NP) STOP 'N > NP IN SVDCMP'
!      IF (M .LT. N)  STOP 'M < N IN SVDCMP'
      LTRANS = (M .LT. N)
      IF (LTRANS) THEN
         M1 = N;  N1 = M
         ALLOCATE (A1(M1,N1))
         A1 = TRANSPOSE(A)
      ELSE
         M1 = M;  N1 = N
         A1 => A
      END IF

!     1. Householder reduction to bidiagonal form

      G = zero
      SCALE = zero
      ANORM = zero
      DO I = 1, N1
         L = I + 1
         RV1(I) = SCALE*G
         G = zero
         S = zero
         SCALE = zero
         IF (I .le. M1) THEN
            SCALE = SUM(ABS(A1(I:M1,I)))
            IF (SCALE .ne. zero) THEN
               A1(I:M1,I) = A1(I:M1,I)/SCALE
               S = SUM(A1(I:M1,I)*A1(I:M1,I))
               F = A1(I,I)
               G = -SIGN(SQRT(S),F)
               H = F*G - S
               A1(I,I) = F - G
               IF (I .ne. N1) THEN
                  DO J = L, N1
                     S = SUM(A1(I:M1,I)*A1(I:M1,J))
                     F = S/H
                     A1(I:M1,J) = A1(I:M1,J) + F*A1(I:M1,I)
                  END DO
               ENDIF
               A1(I:M1,I) = SCALE*A1(I:M1,I)
            ENDIF
         ENDIF
         W(I) = SCALE*G
         G = zero
         S = zero
         SCALE = zero
         IF (I.le.M1 .AND. I.ne.N1) THEN
            SCALE = SUM(ABS(A1(I,L:N1)))
            IF (SCALE .ne. zero) THEN
               A1(I,L:N1) = A1(I,L:N1)/SCALE
               S = SUM(A1(I,L:N1)*A1(I,L:N1))
               F = A1(I,L)
               G = -SIGN(SQRT(S),F)
               H = F*G - S
               A1(I,L) = F - G
               RV1(L:N1) = A1(I,L:N1)/H
               IF (I .ne. M1) THEN
                  DO J = L, M1
                     S = SUM(A1(J,L:N1)*A1(I,L:N1))
                     A1(J,L:N1) = A1(J,L:N1) + S*RV1(L:N1)
                  END DO
               ENDIF
               A1(I,L:N1) = SCALE*A1(I,L:N1)
            ENDIF
         ENDIF
         ANORM = MAX(ANORM,ABS(W(I))+ABS(RV1(I)))
      END DO

!     2. Accumulation of right-hand (V) transformations
      DO I = N1, 1, -1
         IF (I .lt. N1) THEN
            IF (G .ne. zero) THEN
               V(L:N1,I) = (A1(I,L:N1)/A1(I,L))/G
               DO J = L, N1
                  S = SUM(A1(I,L:N1)*V(L:N1,J))
                  V(L:N1,J) = V(L:N1,J) + S*V(L:N1,I)
               END DO
            ENDIF
            V(I,L:N1) = zero
            V(L:N1,I) = zero
         ENDIF
         V(I,I) = one
         G = RV1(I)
         L = I
      END DO

!     3. Accumulation of left-hand (A) transformations
      DO I = N1, 1, -1
         L = I + 1
         G = W(I)
         IF (I .lt. N1) THEN
            A1(I,L:N1) = zero
         ENDIF
         IF (G .ne. zero) THEN
            G = one/G
            IF (I .ne. N1) THEN
               DO J = L, N1
                  S = SUM(A1(L:M1,I)*A1(L:M1,J))
                  F = (S/A1(I,I))*G
                  A1(I:M1,J) = A1(I:M1,J) + F*A1(I:M1,I)
               END DO
            ENDIF
            A1(I:M1,I) = A1(I:M1,I)*G
         ELSE
            A1(I:M1,I) = zero
         ENDIF
         A1(I,I) = A1(I,I) + one
      END DO

!     4. Diagonalization of the bidiagonal form
!     4a.Loop over singular values
      L49: DO K = N1, 1, -1
!        4b. Loop over allowed iterations
         DO ITS = 1, 50
!           4c. Test for splitting
            DO L = K, 1, -1
               NM = L - 1
               IF (ABS(RV1(L)) + ANORM .eq. ANORM) GO TO 2
               IF (ABS(W(NM)) + ANORM .eq. ANORM) EXIT
            END DO
!           4d. Cancellation of RV1(L), if L>1
            C = zero
            S = one
            DO I = L, K
               F = S*RV1(I)
               IF (ABS(F) + ANORM .ne. ANORM) THEN
                  G = W(I)
                  H = PYTHAG(F,G)
c              H=SQRT(F*F+G*G)
                  W(I) = H
                  H = one/H
                  C = G*H
                  S = -F*H
 
                  CALL ROTATE_SVD (A1(1,NM), A1(1,I), C, S, M1)
!                  DO J = 1, M
!                     Y = A(J,NM)
!                    Z = A(J,I)
!                    A(J,NM) = Y*C + Z*S
!                    A(J,I)  =-Y*S + Z*C
!                  END DO

               ENDIF
            END DO
    2       CONTINUE
            Z = W(K)

!          4e. Convergence, make singular value non-negative
            IF (L .eq. K) THEN
               IF (Z .lt. zero) THEN
                  W(K) = -Z
                  V(:N1,K) = -V(:N1,K)
               ENDIF
               CYCLE  L49
            ENDIF
            IF (ITS .eq. 50) THEN
               WRITE (*, '(2A)') 'PAUSE ',
     1            'No convergence in 50 iterations'
               READ *
            ENDIF
!           4f. Shift from bottom 2 X 2 minor
            X = W(L)
            NM = K - 1
            Y = W(NM)
            G = RV1(NM)
            H = RV1(K)
            F = ((Y - Z)*(Y + Z) + (G - H)*(G + H))/(TWO*H*Y)
            G = PYTHAG(F,one)
c          G=SQRT(F*F + one)
            F = ((X - Z)*(X + Z) + H*(Y/(F + SIGN(G,F)) - H))/X

!           4f. Next QR transformation
            C = one
            S = one
            DO J = L, NM
               I = J + 1
               G = RV1(I)
               Y = W(I)
               H = S*G
               G = C*G
               Z = PYTHAG(F,H)
c            Z=SQRT(F*F+H*H)
               RV1(J) = Z
               C = F/Z
               S = H/Z
               F = X*C + G*S
               G =-X*S + G*C
               H = Y*S
               Y = Y*C

               CALL ROTATE_SVD (V(1,J), V(1,I), C, S, N1)
!               DO JJ = 1, N
!                  X = V(JJ,J)
!                  Z = V(JJ,I)
!                  V(JJ,J) = X*C + Z*S
!                  V(JJ,I) =-X*S + Z*C
!               END DO
               Z = PYTHAG(F,H)
c            Z=SQRT(F*F+H*H)
               W(J) = Z

!              4g. Rotation can be arbitrary if Z=0
               IF (Z .ne. zero) THEN
                  Z = one/Z
                  C = F*Z
                  S = H*Z
               ENDIF
               F = C*G + S*Y
               X =-S*G + C*Y
               CALL ROTATE_SVD(A1(1,J), A1(1,I), C, S, M1)
!               DO JJ = 1, M
!                  Y = A(JJ,J)
!                  Z = A(JJ,I)
!                  A(JJ,J) = Y*C + Z*S
!                  A(JJ,I) =-Y*S + Z*C
!               END DO
            END DO
            RV1(L) = zero
            RV1(K) = F
            W(K) = X
         END DO
      END DO L49

      IF (LTRANS) THEN
         A(1:M,1:M) = V(1:N1,1:N1)  !U is now square: M X M
         V(1:N,1:M) = A1(1:N,1:M)   !V not square anymore, but N X M
         DEALLOCATE (A1)
      END IF

      END SUBROUTINE SVDCMP


      SUBROUTINE ROTATE_SVD (A, B, C, S, M)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: M
      REAL(rprec), INTENT(inout), DIMENSION(M) :: A, B
      REAL(rprec), INTENT(in)                  :: C, S
C-----------------------------------------------
      REAL(rprec) :: YTEMP(M)
C-----------------------------------------------
      YTEMP = A
      A = YTEMP*C + B*S
      B =-YTEMP*S + B*C

      END SUBROUTINE ROTATE_SVD
