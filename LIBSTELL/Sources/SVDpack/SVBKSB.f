      SUBROUTINE SVBKSB(U, W, V, M, N, MP, NP, B, X, LNOTRANS)
      USE stel_kinds
      IMPLICIT NONE
!
!     Solves A*X = B for a vector X, whre A is specified by the arrays
!     U, W, V returned from an SVD decomposition.
!
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER M, N, MP, NP
      REAL(rprec), DIMENSION(MP,NP), INTENT(in) :: U
      REAL(rprec), DIMENSION(NP), INTENT(in)    :: W
      REAL(rprec), DIMENSION(NP,NP), INTENT(in) :: V
      REAL(rprec), DIMENSION(MP), INTENT(in)    :: B
      REAL(rprec), DIMENSION(NP), INTENT(out)   :: X
      LOGICAL, INTENT(in)                       :: LNOTRANS
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J
      REAL(rprec), DIMENSION(NP) :: TMP
      REAL(rprec) :: S
C-----------------------------------------------
      IF (M .GT. MP) STOP 'M > MP IN SVBKSB'
      IF (N .GT. NP) STOP 'N > NP IN SVBKSB'

      DO J = 1, N
         S = zero
         IF (W(J) .ne. zero) THEN
            S = SUM(U(:M,J)*B(:M))
            S = S/W(J)
         ENDIF
         TMP(J) = S
      END DO

      DO J = 1, N
         IF (LNOTRANS) THEN
            S = SUM(V(J,:N)*TMP(:N))
         ELSE
            S = SUM(V(:N,J)*TMP(:N))
         END IF
         X(J) = S
      END DO

      END SUBROUTINE SVBKSB
