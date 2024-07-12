C******************** START FILE SPLINE.FOR ; GROUP TRKRLIB ******************
CC
CCCCCCCCCCCCCCCC
C	THE CODES (SPLINE & SEVAL) ARE TAKEN FROM:
C	FORSYTHE,MALCOLM AND MOLER, "COMPUTER METHODS FOR
C	MATHEMATICAL COMPUTATIONS",PRENTICE-HALL, 1977.
C
C	THE CODES (SPLEEN,SPLAAN & SPEVAL) ARE ADAPTATIONS
C	BY R.M. WIELAND FOR SPECIAL CASES ... SEE COMMENTS
c
c  this is nspline -- natural cubic spline
c    second derivative of interpolant set to zero at endpts
c
 
      SUBROUTINE NSPLINE (N, X, Y, B, C, D)
      INTEGER N
      REAL X(N), Y(N), B(N), C(N), D(N)
C
C  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
C  FOR A CUBIC INTERPOLATING SPLINE
C
C    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C    FOR  X(I) .LE. X .LE. X(I+1)
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
C    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
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
      NM1 = N-1
      IF ( N .LT. 2 ) RETURN
C
C
C DMC - CHECK ORDINATES
      DO 5 I=1,NM1
        IP1=I+1
        IF(X(IP1).LE.X(I)) then
             write(6,9000)
             write(6,9001) (x(ii),ii=1,N)
             zbomb=sqrt(-1.0/(x(1)*x(1)))
             write(6,*) zbomb
             write(6,*) ' nspline -- code stop.'
             call exit(1)
 9000        FORMAT(/
     >' ? UNORDERED ORDINATE ARRAY PASSED TO SPLINE SUBROUTINE:')
 9001        FORMAT(1X,5(1X,1PE12.5))
        endif
 5    CONTINUE
C
      IF ( N .LT. 3 ) GO TO 50
C
C  SET UP TRIDIAGONAL SYSTEM
C
C  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
C
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      DO 10 I = 2, NM1
         D(I) = X(I+1) - X(I)
         B(I) = 2.*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
   10 CONTINUE
C
C  END CONDITIONS:  2nd deriv = 0 at endpts; solve reduced n-2 x n-2
C  system.
      C(1) = 0.
      C(N) = 0.
C
C  FORWARD ELIMINATION
C
   15 DO 20 I = 3, NM1
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
   20 CONTINUE
C
C  BACK SUBSTITUTION
C
      NM2=NM1-1
      C(NM1) = C(NM1)/B(NM1)
      DO 30 IB = 2, NM2
         I = N-IB
         C(I) = (C(I) - D(I)*C(I+1))/B(I)
   30 CONTINUE
C
C  C(I) IS NOW THE SIGMA(I) OF THE TEXT
C
C  COMPUTE POLYNOMIAL COEFFICIENTS
C
      B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2.*C(N))
      DO 40 I = 1, NM1
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
   40 CONTINUE
      C(N) = 3.*C(N)
      D(N) = D(N-1)
      RETURN
C
   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.
      RETURN
      END
C******************** END FILE SPLINE.FOR ; GROUP TRKRLIB ******************
