      SUBROUTINE SORTSVD(M, N, MP, NP, W, U, V)
c------------------------------
c     Sorts the weights and U, V matrices from SVD decomposition
c     Sorts weights in DECREASING order
c     Based on routine SORT2 from Numerical Recipes, pg 231
c     Created by Prashant Valanju, (Sept 1998)
c------------------------------
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(NP) :: W
      REAL(rprec), DIMENSION(MP,NP) :: U
      REAL(rprec), DIMENSION(NP,NP) :: V
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: L, IR, I, J
      REAL(rprec), DIMENSION(MP) :: RU
      REAL(rprec), DIMENSION(NP) :: RV
      REAL(rprec) :: RW
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SORTSVD'
      IF (N .GT. NP) STOP 'N > NP IN SORTSVD'

      L = N/2 + 1
      IR = N
   10 CONTINUE
      IF (L .gt. 1) THEN
         L = L - 1
         RW = W(L)
         RU(:M) = U(:M,L)
         RV(:N) = V(:N,L)
      ELSE
         RW = W(IR)
         RU(:M) = U(:M,IR)
         RV(:N) = V(:N,IR)
         W(IR) = W(1)
         U(:M,IR) = U(:M,1)
         V(:N,IR) = V(:N,1)
         IR = IR - 1
         IF (IR .eq. 1) THEN
            W(1) = RW
            U(:M,1) = RU(:M)
            V(:N,1) = RV(:N)
            RETURN
         ENDIF
      ENDIF
      I = L
      J = L + L
   20 CONTINUE
      IF (J .le. IR) THEN
         IF (J .lt. IR) THEN
C                                    !Just change these 2 for increasing
            IF (W(J) .gt. W(J+1)) J = J + 1
         ENDIF
         IF (RW .gt. W(J)) THEN         !Just change these 2 for increasing
            W(I) = W(J)
            U(:M,I) = U(:M,J)
            V(:N,I) = V(:N,J)
            I = J
            J = J + J
         ELSE
            J = IR + 1
         ENDIF
         GO TO 20
      ENDIF
      W(I) = RW
      U(:M,I) = RU(:M)
      V(:N,I) = RV(:N)
      GO TO 10

      END SUBROUTINE SORTSVD
