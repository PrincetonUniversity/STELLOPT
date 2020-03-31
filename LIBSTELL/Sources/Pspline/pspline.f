C  PSPLINE -- modified from SPLINE.FOR -- dmc 3 Oct 1995
C
C	THE CODES (SPLINE & SEVAL) ARE TAKEN FROM:
C	FORSYTHE,MALCOLM AND MOLER, "COMPUTER METHODS FOR
C	MATHEMATICAL COMPUTATIONS",PRENTICE-HALL, 1977.
c
c  PSPLINE -- periodic spline -- adaptation of SPLINE.FOR by D. McCune
c  October 1995:  The function to be splined is taken as having periodic
C  derivatives (it need not be periodic itself), so that
c  the interpolating function satisfies
c            y'(x(1))=y'(x(N))
c           y''(x(1))=y''(x(N))
c  this results in a fully specified system (no addl BC assumptions
c  required); it is not a tridiagonal system but a modified NM1xNM1
c  tridiagonal system with non-zero offdiagonal corner (1,NM1), (NM1,1)
c  elements.  It remains symmetric & diagonally dominant, and is
c  solved by Gaussian elimination w/o pivoting.  NM1=N-1.
c
      SUBROUTINE PSPLINE (N, X, Y, B, C, D, WK)
      INTEGER N
      REAL X(N), Y(N), B(N), C(N), D(N), WK(N)
C
C  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
C  FOR A periodic CUBIC INTERPOLATING SPLINE
C
C    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C    FOR  X(I) .LE. X .LE. X(I+1)
C
C    WK(...) is used as a workspace during the Guassian elimination.
C    it represents the rightmost column of the array
C
C    note B(N),C(N),D(N) give coeffs for a periodic continuation into
C    the next period, up to X(N)+(X(2)-X(1)) only.
C
C    SEVAL can be used to evaluate the spline, but, it is up to the
C    SEVAL caller to bring the argument into the normal range X(1) to X(N).
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
C        includes a complete period of the data, the last point being
C        a repeat of the first point.
C    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
C        X(N)-X(1) is the period of the periodic function represented.
C    Y = THE ORDINATES OF THE KNOTS
C
C  OUTPUT..
C
C    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
C
C  USING  P  TO DENOTE DIFFERENTIATION,
C
C    Y(I) = S(X(I))
C    B(I) = SP(X(I))
C    C(I) = SPP(X(I))/2
C    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
C
CCCCCCCCCCCCCCC
C  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
C  TO EVALUATE THE SPLINE.
C
C
      INTEGER NM1, NM2, IB, I
      REAL T
C
C-------------------------------------
C
      NM1 = N-1
      NM2 = NM1-1
C
      IF ( N .LT. 4 ) THEN
         write(6,9901)
 9901    format(/
     >   ' ?PSPLINE -- at least 4 pts required for periodic spline.')
         return
      ENDIF
C
C  SET UP MODIFIED NM1 x NM1 TRIDIAGONAL SYSTEM:
C  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
C  WK(1:NM2) = rightmost column above diagonal
C  WK(NM1)   = lower left corner element
C
      D(1) = X(2) - X(1)
      D(NM1) = X(N) - X(NM1)
      B(1) = 2.0*(D(1)+D(NM1))
      C(1) = (Y(N) - Y(NM1))/D(NM1)
      C(2) = (Y(2) - Y(1))/D(1)
      C(1) = C(2) - C(1)
      WK(1) = D(NM1)
      DO 10 I = 2, NM1
         D(I) = X(I+1) - X(I)
         B(I) = 2.*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
         WK(I) = 0.0
   10 CONTINUE
      WK(NM2) = D(NM2)
      WK(NM1) = D(NM1)
C
C  END CONDITIONS -- implied by periodicity
C     C(1)=C(N)  B(1)=B(N)  D(1)=D(N) -- no additional assumption needed.
C
C  FORWARD ELIMINATION
C   WK(1)--WK(NM2) represent the rightmost column above the
C   diagonal; WK(NM1) represents the non-zero lower left corner element
C   which in each step is moved one column to the right.
C
   15 DO 20 I = 2, NM2
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
         WK(I) = WK(I) - T*WK(I-1)
         Q = WK(NM1)/B(I-1)
         WK(NM1) = -Q*D(I-1)
         B(NM1) = B(NM1) - Q*WK(I-1)
         C(NM1) = C(NM1) - Q*C(I-1)
   20 CONTINUE
C
C  correct the (NM1,NM2) element
C
      WK(NM1) = WK(NM1) + D(NM2)
C
C  complete the forward elimination:  now WK(NM1) and WK(NM2) are
C  the off diagonal elements of the 2x2 at the lower right hand corner
C
      T = WK(NM1)/B(NM2)
      B(NM1) = B(NM1) - T*WK(NM2)
      C(NM1) = C(NM1) - T*C(NM2)
C
C  BACK SUBSTITUTION
C
      C(NM1) = C(NM1)/B(NM1)
      C(NM2) = (C(NM2) - WK(NM2)*C(NM1))/B(NM2)
C
      DO 30 IB = 3, NM1
         I = N-IB
         C(I) = (C(I) - D(I)*C(I+1) - WK(I)*C(NM1))/B(I)
   30 CONTINUE
C
C  PERIODICITY:
C
      C(N)=C(1)
C
C  C(I) IS NOW THE SIGMA(I) OF THE TEXT
C
C  COMPUTE POLYNOMIAL COEFFICIENTS
C
      DO 40 I = 1, NM1
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
   40 CONTINUE
C periodic continuation...
      B(N) = B(1)
      C(N) = C(1)
      D(N) = D(1)
C
      RETURN
      END
