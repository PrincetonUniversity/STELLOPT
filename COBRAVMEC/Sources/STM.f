      SUBROUTINE STM(E, Y, Z, K, T, D, P, N)
!      ________________________________________________________
!     |                                                        |
!     |   FIND THE K-TH SMALLEST EIGENVALUE OF A TRIDIAGONAL   |
!     |  MATRIX WHOSE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE  |
!     |USING BOTH THE BISECTION METHOD AND A NEWTON-LIKE METHOD|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         Y,Z   --ENDS OF AN INTERVAL THAT CONTAINS THE  |
!     |                 DESIRED EIGENVALUE                     |
!     |                                                        |
!     |         K     --INDEX OF DESIRED EIGENVALUE            |
!     |                                                        |
!     |         T     --TOLERANCE (ITERATIONS CONTINUE UNTIL   |
!     |                 THE ERROR IN THE EIGENVALUE .LE. T)    |
!     |                                                        |
!     |         D     --DIAGONAL OF THE COEFFICIENT MATRIX A   |
!     |                                                        |
!     |         P     --CROSS-DIAGONAL PRODUCTS A SUB I+1,I    |
!     |                 TIMES A SUB I,I+1                      |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MAX,SIGN,SQRT              |
!     |    PACKAGE SUBROUTINES: CP,EQL,INP                     |
!     |________________________________________________________|
!
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: K, N
      REAL(rprec):: E, Y, Z, T
      REAL(rprec), DIMENSION(N) :: D, P
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, M
      REAL(rprec) :: A,B,FL,FR,L,Q,R,S,U,V,
     1   W,D2,D3,D4,P2,P3,P4,D34,D42,D23
      REAL(rprec), DIMENSION(4) :: F, H, X
      REAL(rprec), DIMENSION(4) :: G
      REAL(rprec) :: C, GL, GR
!-----------------------------------------------
      L = Y
      R = Z
      IF ( L .LT. R ) GOTO 10
      L = Z
      R = Y
10    IF ( N .EQ. 1 ) GOTO 210
      S = T + T
      U = ONE
20    U = U/2
      A = ONE + U
      IF (A .GT. ONE) GOTO 20
      U = 5*U
      V = U/2
      CALL CP(L,FL,GL,I,D,P,N)
      CALL CP(R,FR,GR,J,D,P,N)
      IF ( I .GE. K ) GOTO 190
      IF ( J .LT. K ) GOTO 200
C     --------------------------------
C     |*** ISOLATE THE EIGENVALUE ***|
C     --------------------------------
30    E = R - L
      IF ( E .LE. S ) GOTO 180
      IF ( E .LE. U*(ABS(L)+ABS(R)) ) GOTO 180
      IF ( J .EQ. I+1 ) GOTO 70
40    A = .5*(L+R)
      CALL CP(A,B,C,M,D,P,N)
      IF ( K .LE. M ) GOTO 50
      L = A
      I = M
      FL = B
      GL = C
      GOTO 30
50    R = A
      J = M
      FR = B
      GR = C
      GOTO 30
60    E = R - L
      IF ( E .LE. S ) GOTO 180
      IF ( E .LE. U*(ABS(L)+ABS(R)) ) GOTO 180
70    X(1) = L
      F(1) = FL
      G(1) = GL
      X(2) = R
      F(2) = FR
      G(2) = GR
      CALL EQL(X,F,G,H,2)
      IF ( H(1) .EQ. H(2) ) GOTO 160
C     ---------------------
C     |*** SECANT STEP ***|
C     ---------------------
      A = X(1) - H(1)*(X(1)-X(2))/(H(1)-H(2))
      Q = A
      W = MAX(T,V*(ABS(L)+ABS(R)))
      IF ( ABS(A-L) .LT. W ) A = L + W
      IF ( ABS(A-R) .LT. W ) A = R - W
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 80
      R = A
      FR = B
      GR = C
      GOTO 90
80    L = A
      FL = B
      GL = C
90    X(3) = A
      F(3) = B
      G(3) = C
      W = R - L
      IF ( W .LE. S ) GOTO 220
      IF ( W .LE. U*(ABS(L)+ABS(R)) ) GOTO 220
      CALL EQL(X,F,G,H,3)
C     --------------------------------------
C     |*** QUADRATIC INTERPOLATION STEP ***|
C     --------------------------------------
      CALL INP(A,X(1),X(2),X(3),H(1),H(2),H(3),L,R)
      B = L
      IF ( ABS(A-L) .GT. ABS(A-R) ) B = R
C     ------------------------------------
C     |*** APPLY PSEUDO-NEWTON METHOD ***|
C     ------------------------------------
100   Q = A
      W = MAX(T,V*(ABS(L)+ABS(R)))
      IF ( ABS(A-L) .LT. W ) GOTO 110
      IF ( ABS(A-R) .GT. W ) GOTO 130
110   IF ( A+A .GT. L+R ) GOTO 120
      A = L + W
      GOTO 130
120   A = R - W
130   IF ( A .LE. L ) GOTO 160
      IF ( A .GE. R ) GOTO 160
      E = .5*E
      IF ( E .LT. ABS(B-A) ) GOTO 160
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 140
      R = A
      FR = B
      GR = C
      GOTO 150
140   L = A
      FL = B
      GL = C
150   W = R - L
      IF ( W .LE. S ) GOTO 220
      IF ( W .LE. U*(ABS(L)+ABS(R)) ) GOTO 220
      X(4) = A
      F(4) = B
      G(4) = C
      CALL EQL(X,F,G,H,4)
      IF ( X(1) .LT. L ) GOTO 160
      IF ( X(1) .GT. R ) GOTO 160
      B = X(1)
      D4 = X(4) - B
      D3 = X(3) - B
      D2 = X(2) - B
      D34 = X(3) - X(4)
      D42 = X(4) - X(2)
      D23 = X(2) - X(3)
      P2 = D2*(ONE+((D2/D3)*D42+(D2/D4)*D23)/D34)
      P3 = D3*(ONE+((D3/D2)*D34+(D3/D4)*D23)/D42)
      P4 = D4*(ONE+((D4/D2)*D34+(D4/D3)*D42)/D23)
      IF (P2 .NE. ZERO) P2 = (H(2)-H(1))/P2
      IF (P3 .NE. ZERO) P3 = (H(3)-H(1))/P3
      IF (P4 .NE. ZERO) P4 = (H(4)-H(1))/P4
      P2 = P2 + P3 + P4
      IF (P2 .EQ. ZERO) GOTO 160
      A = B - H(1)/P2
      GOTO 100
C     --------------------------
C     |*** BISECTION METHOD ***|
C     --------------------------
160   A = (L+R)/2
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 170
      R = A
      FR = B
      GR = C
      GOTO 60
170   L = A
      FL = B
      GL = C
      GOTO 60
180   E = (L+R)/2
      RETURN
190   E = L
      RETURN
200   E = R
      RETURN
210   E = D(1)
      IF ( L .GT. E ) E = L
      IF ( R .LT. E ) E = R
      RETURN
220   E = Q
      RETURN
      END SUBROUTINE STM
