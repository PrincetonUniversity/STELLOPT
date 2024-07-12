      SUBROUTINE TVAL(E, K, L, D, U, N, W)

!
! ===================================
! NIST Guide to Available Math Software.
! Fullsource for MODULE TVAL from package NAPACK.
! Retrieved from NETLIB on Fri Aug 7 17:41:59 1998.
! ===================================
!
!      ________________________________________________________
!     |                                                        |
!     |   FIND THE K-TH SMALLEST EIGENVALUE OF A TRIDIAGONAL   |
!     |  MATRIX WHOSE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE  |
!     |USING BOTH THE BISECTION METHOD AND A NEWTON-LIKE METHOD|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         K     --INDEX OF DESIRED EIGENVALUE (REPLACE K |
!     |                 BY -K TO GET K-TH LARGEST EIGENVALUE)  |
!     |                                                        |
!     |         L     --SUBDIAGONAL (CAN BE IDENTIFIED WITH U) |
!     |                                                        |
!     |         D     --DIAGONAL                               |
!     |                                                        |
!     |         U     --SUPERDIAGONAL                          |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |         W     --WORK ARRAY (LENGTH AT LEAST N)         |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MAX                        |
!     |    PACKAGE SUBROUTINES: CP,EQL,INP,STM                 |
!     |________________________________________________________|
!
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: K, N
      REAL(rprec):: E
      REAL(rprec), DIMENSION(N) :: L, D, U, W
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, M
      REAL(rprec) :: S, T, Y, Z
!-----------------------------------------------
      IF (N > 1) THEN
!     -------------------------------------------------------
!     |*** BOUND EIGENVALUES USING GERSCHGORIN'S THEOREM ***|
!     -------------------------------------------------------
         M = N - 1
         Y = D(1)
         Z = D(1)
         S = 0
         DO I = 1, M
            W(I) = L(I)*U(I)
            IF (W(I) < ZERO) GO TO 50
            S = S + ABS(U(I))
            T = D(I) + S
            Z = MAX(T,Z)
            T = D(I) - S
            Y = MIN(T,Y)
            S = ABS(L(I))
         END DO
         T = D(N) + S
         Z = MAX(T,Z)
         T = D(N) - S
         Y = MIN(T,Y)
         M = K
         IF (K .LE. 0) M = N + K + 1
         T = ONE
         T = T/2
         S = ONE + T
 1003    CONTINUE
         IF (S .LE. ONE) GO TO 1002
         T = T/2
         S = ONE + T
         IF (S .LE. ONE) GO TO 1002
         T = T/2
         S = ONE + T
         IF (S .LE. ONE) GO TO 1002
         T = T/2
         S = ONE + T
         IF (S .LE. ONE) GO TO 1002
         T = T/2
         S = ONE + T
         GO TO 1003
 1002    CONTINUE
         T = 8*N*T*MAX(ABS(Y),ABS(Z))
         IF (T .NE. ZERO) THEN
            CALL STM (E, Y, Z, M, T, D, W, N)
            RETURN
         ENDIF
      ENDIF
      E = D(1)
      RETURN

   50 CONTINUE
      WRITE (6, *) 'ERROR: SUBROUTINE TVAL CAN ONLY BE APPLIED'
      WRITE (6, *) 'WHEN THE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE'
      STOP

      END SUBROUTINE TVAL
