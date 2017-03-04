  MODULE QSHEP2D_MOD
  PRIVATE
  PUBLIC QSHEP2, QS2VAL, QS2GRD
  CONTAINS
  SUBROUTINE QSHEP2 (N, X, Y, F, NQ, NW, NR, LCELL, LNEXT, XMIN, &
      YMIN, DX, DY, RMAX, RSQ, A, IER)
!
! QSHEP2 computes an interpolant to scattered data in the plane.
!
! QSHEP2 computes a set of parameters, A and RSQ, defining a smooth, 
! once continuously differentiable, bi-variate function Q(X,Y) which 
! interpolates data values, F, at scattered nodes (X,Y).  
!
! The interpolant function Q(X,Y) may be evaluated at an arbitrary point 
! by passing the parameters A and RSQ to the function QS2VAL. The
! first derivatives dQ/dX(X,Y) and dQ/dY(X,Y) may be evaluated by 
! subroutine QS2GRD.
!
! The interpolation scheme is a modified quadratic Shepard method:
!
!   Q = (W(1) * Q(1) + W(2) * Q(2) + ... + W(N) * Q(N)) 
!       / (W(1)      + W(2)        + ... + W(N))
!
!   for bivariate functions W(K) and Q(K). The nodal functions are given by
!
!   Q(K)(X,Y) =   F(K)
!               + A(1,K) * (X - X(K))**2 
!               + A(2,K) * (X - X(K)) * (Y - Y(K))
!               + A(3,K) * (Y - Y(K))**2 
!               + A(4,K) * (X - X(K))
!               + A(5,K) * (Y - Y(K)).
!
!    Thus, Q(K) is a quadratic function which interpolates the
!    data value at node K. Its coefficients A(*,K) are obtained
!    by a weighted least squares fit to the closest NQ data
!    points with weights similar to W(K). Note that the radius
!    of influence for the least squares fit is fixed for each
!    K, but varies with K.
!
!    The weights are taken to be
!
!        W(K)(X,Y) = ((R(K)-D(K))+ / R(K) * D(K))**2
!
!    where (R(K)-D(K))+ = 0 if R(K) <= D(K) and D(K)(X,Y) is
!    the euclidean distance between (X,Y) and (X(K),Y(K)).  The
!    radius of influence R(K) varies with K and is chosen so
!    that NW nodes are within the radius.  Note that W(K) is
!    not defined at node (X(K),Y(K)), but Q(X,Y) has limit F(K)
!    as (X,Y) approaches (X(K),Y(K)).
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   N, integer, the number of nodes (X,Y) at which data values
!   are given.  N must be at least 6.
!
!   X(N), Y(N), real, the coordinates of the nodes at which
!   data has been supplied.
!
!   F(N), real, the data values.
!
!   NQ, integer, the number of data points to be used in the least
!   squares fit for coefficients defining the nodal functions Q(K).  
!   A highly recommended value is NQ = 13.  
!   NQ must be at least 5, and no greater than MIN(40,N-1).
!
!   NW, integer, the number of nodes within (and defining) the radii
!   of influence R(K) which enter into the weights W(K). 
!   For N sufficiently large, a recommended value is NW = 19.   
!   NW must be at least 1, and no greater than MIN(40,N-1).
!
!   NR, integer, the number of rows and columns in the cell grid 
!   defined in subroutine STORE2. A rectangle containing the nodes 
!   is partitioned into cells in order to increase search efficiency.  
!   NR = SQRT(N/3) is recommended. NR must be at least 1.
!
!  OUTPUT 
!   LCELL(NR,NR), integer, array of nodal indices associated
!   with cells.
!
!   LNEXT(N), integer, contains next-node indices (or their 
!   negatives).
!
!   XMIN, YMIN, DX, DY, real, the minimum nodal X, Y coordinates,
!   and the X, Y dimensions of a cell.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius of influence.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining the interpolant Q.
!
!   A(5,N), real, the coefficients for the nodal functions 
!   defining the interpolant Q.
!
!   IER, integer, error indicator.
!       0, if no errors were encountered.
!       1, if N, NQ, NW, or NR is out of range.
!       2, if duplicate nodes were encountered.
!       3, if all nodes are collinear.
!
! Local variables:
!
!   AV =        Root-mean-square distance between node K and the
!               nodes in the least squares system (unless
!               additional nodes are introduced for stability). 
!               The first 3 columns of the matrix are scaled by 
!               1/AVSQ, the last 2 by 1/AV.
!   
!   AVSQ =      AV * AV.
!
!   B =         Transpose of the augmented regression matrix.
!
!   C =         First component of the plane rotation used to
!               zero the lower triangle of B**T.
!               Computed by the subroutine GIVENS.
!
!   DDX, DDY =  Local variables for DX and DY.
!
!   DMIN =      Minimum of the magnitudes of the diagonal
!               elements of the regression matrix after
!               zeros are introduced below the diagonal.
!
!   DTOL =      Tolerance for detecting an ill-conditioned system. 
!               The system is accepted when DMIN >= DTOL.
!
!   FK =        Data value at node K (same as F(K)).
!
!   I =         Index for A, B, and NPTS.
!
!   IB =        Do-loop index for back solve.
!
!   IERR =      Error flag for the call to STORE2.
!
!   IROW =      Row index for B.
!
!   J =         Index for A and B.
!
!   JP1 =       J+1.
!
!   K =         Nodal function index and column index for A.
!
!   LMAX =      Maximum number of NPTS elements 
!               (must be consistent with the dimension of NPTS in the 
!               variable declaration section).
!
!   LNP =       Current length of NPTS.
!
!   NEQ =       Number of equations in the least squares fit.
!
!   NN, NNR =   Local copies of N and NR.
!
!   NP =        NPTS element.
!
!   NPTS =      Array containing the indices of a sequence of
!               nodes to be used in the least squares fit
!               or to compute RSQ. The nodes are ordered
!               by distance from node K and the last element
!               (usually indexed by LNP) is used only to
!               determine RQ, or RSQ(K) if NW > NQ.
!
!   NQWMAX =    MAX(NQ,NW).
!
!   RQ =        Radius of influence which enters into the
!               weights for Q(K) (see subroutine SETUP2).
!
!   RS =        Squared distance between K and NPTS(LNP).
!               Used to compute RQ and RSQ(K).
!
!   RSMX =      Maximum RSQ element encountered.
!
!   RSOLD =     Squared distance between K and NPTS(LNP-1).
!               Used to compute a relative change in RS
!               between succeeding NPTS elements.
!
!   RTOL =      Tolerance for detecting a sufficiently large
!               relative change in RS. If the change is
!               not greater than RTOL, the nodes are
!               treated as being the same distance from K.
!
!   RWS =       Current value of RSQ(K).
!
!   S =         Second component of the plane Givens rotation.
!
!   SF =        Marquardt stabilization factor used to damp
!               out the first 3 solution components (second
!               partials of the quadratic) when the system
!               is ill-conditioned. As SF increases, the
!               fitting function becomes more linear.
!
!   SUM2 =      Sum of squared Euclidean distances between
!               node K and the nodes used in the least
!               squares fit (unless additional nodes are
!               added for stability).
!
!   T =         Temporary variable for accumulating a scalar
!               product in the back solve.
!
!   XK, YK =    Coordinates of node K.
!
!   XMN, YMN =  Local variables for XMIN and YMIN.
!
    USE REAL_PRECISION
    IMPLICIT NONE
    
    INTEGER  :: I, IER, IERR, IROW, J, JP1, K, LMAX, LNP, N, NEQ, NN, &
          NNR, NP, NQ, NQWMAX, NR, NW
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION(40) :: NPTS
    INTEGER, DIMENSION(NR,NR) :: LCELL

    REAL(KIND=R8):: AV, AVSQ, C, DDX, DDY, DMIN
    REAL(KIND=R8):: DTOL
    REAL(KIND=R8):: DX, DY, FK, RMAX, RQ, RS, RSMX, RSOLD
    REAL(KIND=R8), PARAMETER :: RTOL = 1.0E-05_R8
    REAL(KIND=R8):: RWS, S
    REAL(KIND=R8), PARAMETER :: SF = 1.0E+00_R8
    REAL(KIND=R8):: SUM2, T, XK, XMIN, XMN, YK, YMIN, YMN
    REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X, Y
    REAL(KIND=R8), DIMENSION(5,N) :: A
    REAL(KIND=R8), DIMENSION(6,6) :: B

    DTOL = SQRT(EPSILON(1.0_R8))
    NN = N
    NNR = NR
    NQWMAX = MAX (NQ,NW)
    LMAX = MIN (40,N-1)

    IF (5 > NQ) THEN
        IER = 1
        RETURN
    ELSE IF (1 > NW) THEN
        IER = 1
        RETURN
    ELSE IF (NQWMAX > LMAX) THEN
        IER = 1
        RETURN
    ELSE IF (NR < 1) THEN
        IER = 1
        RETURN
    END IF
!
! Create the cell data structure, and initialize RSMX.
!
    CALL STORE2 (NN, X, Y, NNR, LCELL, LNEXT, XMN, YMN, DDX, DDY, IERR)

    IF (IERR /= 0) THEN
        XMIN = XMN
        YMIN = YMN
        DX = DDX
        DY = DDY
        IER = 3
        RETURN
    END IF

    RSMX = 0.0E+00_R8
!
! Outer loop on node K.
!
    DO K = 1, NN
        XK = X(K)
        YK = Y(K)
        FK = F(K)
!
! Mark node K to exclude it from the search for nearest neighbors.
!
        LNEXT(K) = -LNEXT(K)
!
! Initialize for loop on NPTS.
!
        RS = 0.0E+00_R8
        SUM2 = 0.0E+00_R8
        RWS = 0.0E+00_R8
        RQ = 0.0E+00_R8
        LNP = 0
!
! Compute NPTS, LNP, RWS, NEQ, RQ, and AVSQ.
!

        MAIN: DO 
            SUM2 = SUM2 + RS
            IF (LNP == LMAX) THEN
!
! All LMAX nodes are included in NPTS. RWS and/or RQ**2 are
! (arbitrarily) taken to be 10 percent larger than the
! distance rs to the last node included.
!
                IF (RWS == 0.0E+00_R8) THEN
                    RWS = 1.1E+00_R8 * RS
                END IF

                IF (RQ == 0.0E+00_R8) THEN
                    NEQ = LMAX
                    RQ = SQRT (1.1E+00_R8 * RS)
                    AVSQ = SUM2 / REAL(NEQ, KIND=R8)
                END IF
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                RSQ(K) = RWS
                RSMX = MAX (RSMX, RWS)
                AV = SQRT ( AVSQ )
                EXIT
            END IF
            LNP = LNP + 1
            RSOLD = RS

            CALL GETNP2 (XK, YK, X, Y, N, NNR, LCELL, LNEXT, XMN, YMN, &
                  DDX, DDY, NP, RS)

            IF (RS == 0.0E+00_R8) THEN
                IER = 2
                RETURN
            END IF
            NPTS(LNP) = NP

            IF ((RS - RSOLD) / RS < RTOL) THEN
                CYCLE MAIN
            END IF

            IF (RWS == 0.0E+00_R8 .AND. LNP > NW) THEN
                RWS = RS
            END IF
!
! RQ = 0 (not yet computed) and LNP > NQ.     
! RQ = SQRT(RS) is sufficiently large to (strictly) include NQ nodes.  
!
! The least squares fit will include NEQ = LNP - 1 equations for 
! 5 <= NQ <= NEQ < LMAX <= N-1.
!
            IF (RQ == 0.0E+00_R8 .AND. LNP > NQ) THEN
                NEQ = LNP - 1
                RQ = SQRT (RS)
                AVSQ = SUM2 / REAL (NEQ, KIND=R8)
            END IF

            IF ( LNP > NQWMAX ) THEN
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                RSQ(K) = RWS
                RSMX = MAX (RSMX, RWS)
                AV = SQRT ( AVSQ )
            ELSE
                CYCLE MAIN
            END IF
 
        END DO MAIN
!
! Set up the augmented regression matrix (transposed) as the
! columns of B, and zero out the lower triangle (upper
! triangle of B) with Givens rotations. QR decomposition
! with orthogonal matrix Q is not stored.
!
        I = 0
        PTS: DO
            I = I + 1
            NP = NPTS(I)
            IROW = MIN ( I, 6 )
            CALL SETUP2 (XK, YK, FK, X(NP), Y(NP), F(NP), AV, AVSQ, &
                    RQ, B(1,IROW))
            IF ( I == 1 ) THEN
                CYCLE PTS
            END IF

            DO J = 1, IROW-1
                JP1 = J + 1
                CALL GIVENS (B(J,J), B(J,IROW), C, S)
                CALL ROTATE (6-J, C, S, B(JP1,J), B(JP1,IROW))
            END DO

            IF (I < NEQ) THEN
                CYCLE PTS
            END IF
!
! Test the system for ill-conditioning.
!
            DMIN =  MIN (ABS (B(1,1)), ABS (B(2,2)), ABS (B(3,3)), &
                  ABS (B(4,4)), ABS (B(5,5)))
        
            IF ( DMIN * RQ >= DTOL ) THEN
                EXIT
            END IF

            IF (NEQ == LMAX) THEN
                EXIT
            END IF
!
! Increase RQ and add another equation to the system to improve conditioning.  
! The number of NPTS elements is also increased if necessary.
!
            TOL: DO
                RSOLD = RS
                NEQ = NEQ + 1

                IF (NEQ == LMAX) THEN
                    RQ = SQRT ( 1.1E+00_R8 * RS )
                    CYCLE PTS
                END IF
!
! NEQ < LNP.
!
                IF (NEQ /= LNP) THEN
                    NP = NPTS(NEQ+1)
                    RS = (X(NP) - XK)**2 + (Y(NP) - YK)**2
                    IF ((RS - RSOLD) / RS < RTOL) THEN
                        CYCLE TOL
                    END IF
                    RQ = SQRT(RS)
                    CYCLE PTS
                END IF
!
! Add an element to NPTS.
!
                LNP = LNP + 1
                CALL GETNP2 (XK, YK, X, Y, N, NNR, LCELL, LNEXT, XMN, &
                      YMN, DDX, DDY, NP, RS)
                IF ( NP == 0 ) THEN
                    IER = 2
                    RETURN
                END IF

                NPTS(LNP) = NP

                IF (( RS - RSOLD) / RS < RTOL) THEN
                    CYCLE TOL
                END IF

                RQ = SQRT (RS)
                CYCLE PTS
            END DO TOL
        END DO PTS

!
! Stabilize the system by damping the second partials. Add multiples of the 
! first three unit vectors to the first three equations.
!
        IF (NEQ == LMAX) THEN
            DO I = 1, 3
                B(I,6) = SF
                DO J = I+1, 6
                    B(J,6) = 0.0E+00_R8
                END DO
                DO J = I, 5
                    JP1 = J + 1
                    CALL GIVENS (B(J,J), B(J,6), C, S)
                    CALL ROTATE (6-J, C, S, B(JP1,J), B(JP1,6))
                END DO
            END DO
!
! Test the stabilized system for ill-conditioning.
!
            DMIN = MIN (ABS(B(1,1)), ABS(B(2,2)), ABS(B(3,3)), &
                    ABS(B(4,4)), ABS(B(5,5)))

            IF (DMIN * RQ < DTOL) THEN
                XMIN = XMN
                YMIN = YMN
                DX = DDX
                DY = DDY
                IER = 3
                RETURN
            END IF
        END IF
!
! Solve the 5 by 5 triangular system for the coefficients.
!
        DO I = 5, 1, -1
            T = 0.0E+00_R8
            DO J = I+1, 5
                T = T + B(J,I) * A(J,K)
            END DO
            A(I,K) = (B(6,I) - T) / B(I,I)
        END DO
!
! Scale the coefficients to adjust for the column scaling.
!
        DO I = 1, 3
            A(I,K) = A(I,K) / AVSQ
        END DO

        A(4,K) = A(4,K) / AV
        A(5,K) = A(5,K) / AV
!
! Unmark K and the elements of NPTS.
!
        LNEXT(K) = - LNEXT(K)
        DO I = 1, LNP
            NP = NPTS(I)
            LNEXT(NP) = - LNEXT(NP)
        END DO

    END DO  !End outer loop on node K
!
! No errors encountered.
!
    XMIN = XMN
    YMIN = YMN
    DX = DDX
    DY = DDY
    RMAX = SQRT ( RSMX )
    IER = 0
    RETURN
  END SUBROUTINE QSHEP2
!  
!=========================================================================
!  
  FUNCTION QS2VAL (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
      YMIN, DX, DY, RMAX, RSQ, A)
!
! QS2VAL evaluates the interpolant function at the point (PX,PY).
!
! QS2VAL returns the value Q(PX,PY) where Q is the weighted sum of 
! quadratic nodal functions defined by QSHEP2. If the spatial 
! derivatives of Q are also desired, call QS2GRD instead.
!
! Input parameters are not altered by this function.  The
! parameters other than PX and PY should be input unaltered
! from their values on output from QSHEP2.  This function
! should not be called if a nonzero error flag was returned
! by QSHEP2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, real, the (X,Y) coordinates of the point P at
!   which Q is to be evaluated.
!
!   N, integer, the number of nodes. N must be at least 6.
!
!   X(N), Y(N), real, the coordinates of the nodes at which
!   data has been supplied.
!
!   F(N), real, the data values at the nodes.
!
!   NR, integer, the number of rows and columns in the cell grid.
!   Refer to the subroutine STORE2. NR must be at least 1.
!
!   LCELL(NR,NR), integer, the array of nodal indices associated
!   with cells. Refer to STORE2.
!
!   LNEXT(N), integer, the next-node indices. Refer to STORE2.
!
!   XMIN, YMIN, DX, DY, real, the minimum nodal X, Y coordinates,
!   and the X, Y dimensions of a cell. Computed by QSHEP2.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius of influence. Computed by QSHEP2.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining the interpolant Q. Computed by QSHEP2.
!
!   A(5,N), real, the coefficients for the nodal functions 
!   defining the interpolant Q. Computed by QSHEP2.
!
!   
!  OUTPUT
!   QS2VAL, real, the interpolated function value at (PX,PY).
!
    USE REAL_PRECISION
    IMPLICIT NONE

    INTEGER :: I, IMAX, IMIN, J, JMAX, JMIN, K, KP, N, NR
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION(NR, NR) :: LCELL
    
    REAL(KIND=R8):: DELX, DELY, DS, DX, DY, PX, PY, QS2VAL, RD, RDS, &
          RMAX, RS, SW, SWQ, W, XMIN, XP, YMIN, YP
    REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X, Y
    REAL(KIND=R8), DIMENSION(5,N) :: A
    
    XP = PX
    YP = PY
    QS2VAL = 0.0E+00_R8
    IF (N < 6) THEN
        RETURN  
    ELSE IF ( NR < 1  ) THEN
        RETURN
    ELSE IF ( DX <= 0.0E+00_R8 ) THEN
        RETURN
    ELSE IF ( DY <= 0.0E+00_R8 ) THEN
        RETURN
    ELSE IF ( RMAX < 0.0E+00_R8 ) THEN
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indices defining
! the range of the search for nodes whose radii include P.  
! The cells which must be searched are those intersected 
! by (or contained in) a circle of radius RMAX centered at P.
!
    IMIN = INT (( XP - XMIN - RMAX) / DX) + 1
    IMIN = MAX (IMIN, 1)

    IMAX = INT (( XP - XMIN + RMAX) / DX) + 1
    IMAX = MIN (IMAX, NR)

    JMIN = INT (( YP - YMIN - RMAX) / DY) + 1
    JMIN = MAX (JMIN, 1)

    JMAX = INT (( YP - YMIN + RMAX) / DY) + 1
    JMAX = MIN (JMAX, NR)
!
! Test for no cells within the circle of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX) THEN
        QS2VAL = 0.0E+00_R8
        RETURN
    END IF
!
!  Accumulate weight values in SW and weighted nodal function
!  values in SWQ.  The weights are W(K) = ((r-d)+/(r*d))**2
!  for r**2 = RSQ(K) and d = distance between node p and node k.
!
    SW = 0.0E+00_R8
    SWQ = 0.0E+00_R8

    DO J = JMIN, JMAX
        DO I = IMIN, IMAX
            K = LCELL(I,J)
            IF (K /= 0) THEN
                EVAL: DO
                    DELX = XP - X(K)
                    DELY = YP - Y(K)
                    DS = DELX * DELX + DELY * DELY
                    RS = RSQ(K)
                    IF (DS < RS) THEN
                        IF (DS == 0.0E+00_R8) THEN
                            QS2VAL = F(K)
                            RETURN
                        END IF
                        RDS = RS * DS
                        RD = SQRT ( RDS )
                        W = (RS + DS - RD - RD) / RDS
                        SW = SW + W
                        SWQ = SWQ + W * ( F(K) + A(1,K) * DELX * DELX &
                              + A(2,K) * DELX * DELY + A(3,K) * DELY * DELY &
                              + A(4,K) * DELX + A(5,K) * DELY )
                    END IF
                    KP = K
                    K = LNEXT(KP)
                    IF (K /= KP) THEN
                        CYCLE EVAL
                    END IF
                    EXIT 
                END DO EVAL
            END IF
        END DO
    END DO
!
! SW = 0 if and only if P is not within the radius R(K) for any node K.
!
    IF ( SW == 0.0E+00_R8 ) THEN
        QS2VAL = 0.0E+00_R8
    ELSE
        QS2VAL = SWQ / SW
    END IF
    RETURN
  END FUNCTION QS2VAL
!  
!=========================================================================
!
  SUBROUTINE QS2GRD (PX, PY, N, X, Y, F, NR, LCELL, LNEXT, XMIN, &
    YMIN, DX, DY, RMAX, RSQ, A, Q, QX, QY, IER)
!
! QS2GRD evaluates the interpolant and its first spatial derivatives.
!
! QS2GRD computes the value and the gradient at the point (PX,PY) 
! of the interpolatory function Q, defined by QSHEP2 for a given set
! of scattered data.  Q(X,Y) is a weighted sum of quadratic
! nodal functions.
!
! Input parameters are not altered by this subroutine.  The parameters 
! other than PX and PY should be input unaltered from their values 
! on output from QSHEP2.  This subroutine should not be called if a 
! nonzero error flag was returned by QSHEP2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
!  Parameters:
!
!   INPUT
!    PX, PY, real, the coordinates of the point at which the
!    interpolant and its derivatives are to be evaluated.
!
!    N, integer, the number of nodes at which data values are
!    given to define the interpolating function. N must be at 
!    least 6. 
!
!    X(N), Y(N), real, the coordinates of the nodes at which
!    data has been supplied.
!
!    F(N), real, the data values at the nodes.
!
!    NR, integer, the number of rows and columns in the cell 
!    grid. Refer to the subroutine STORE2 for details. 
!    NR must be at least 1.
!
!    LCELL(NR,NR), integer, the array of nodal indices associated
!    with cells.
!
!    LNEXT(N), integer, contains next-node indices.
!
!    XMIN, YMIN, DX, DY, real, the minimum nodal X, Y coordinates,
!    and the X, Y dimensions of a cell.  Computed by QSHEP2.
!
!    RMAX, real, the square root of the largest element in RSQ,
!    the maximum radius of influence. Computed by QSHEP2.
!
!    RSQ(N), real, the squared radii which enter into the weights 
!    defining the interpolant Q. Computed by QSHEP2.
!
!    A(5,N), real, the coefficients for the nodal functions 
!    defining the interpolant Q. Computed by QSHEP2.
!
!   OUTPUT
!    Q, QX, QY, real, the value of the interpolant, and its 
!    derivatives with respect to X and Y, at (PX,PY).
!
!    IER, integer, error indicator.
!      0, if no errors were encountered.
!      1, if N, NR, DX, DY or RMAX is invalid.
!      2, if no errors were encountered but (PX,PY) is not within the 
!         radius R(K) for any node K and thus Q = QX = QY = 0.
!
    USE REAL_PRECISION
    IMPLICIT NONE
  
    INTEGER :: I, IER, IMAX, IMIN, J, JMAX, JMIN, K, KP, N, NR
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION(NR, NR) :: LCELL

    REAL(KIND=R8) :: DELX, DELY, DS, DX, DY, PX, PY, Q, QK, QKX, &
      QKY, QX, QY, RD, RDS, RMAX, RS, SW, SWQ, SWQX, SWQY, SWS, &
      SWX, SWY, T, W, WX, WY, XMIN, XP, YMIN, YP
    REAL(KIND=R8), DIMENSION(N) :: F, RSQ, X, Y
    REAL(KIND=R8), DIMENSION(5,N) :: A

    XP = PX
    YP = PY

    IF (N < 6) THEN
        IER = 1
        RETURN
    ELSE IF (NR < 1) THEN
        IER = 1
        RETURN
    ELSE IF (DX <= 0.0E+00_R8) THEN
        IER = 1
        RETURN
    ELSE IF (DY <= 0.0E+00_R8) THEN
        IER = 1
        RETURN
    ELSE IF (RMAX < 0.0E+00_R8) THEN
        IER = 1
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indices defining
! the range of the search for nodes whose radii include P.
! The cells which must be searched are those intersected 
! by (or contained in) a circle of radius RMAX centered at P.
!
    IMIN = INT (( XP - XMIN - RMAX) / DX) + 1
    IMIN = MAX (IMIN, 1)

    IMAX = INT (( XP - XMIN + RMAX) / DX) + 1
    IMAX = MIN (IMAX, NR)

    JMIN = INT (( YP - YMIN - RMAX) / DY) + 1
    JMIN = MAX (JMIN, 1)

    JMAX = INT (( YP - YMIN + RMAX) / DY) + 1
    JMAX = MIN (JMAX, NR)
!
! Test for no cells within the circle of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX) THEN
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        IER = 2
        RETURN
    END IF
!
! Q = SWQ/SW = SUM(W(K)*Q(K))/SUM(W(K)) where the SUM is
! from K = 1 to N, Q(K) is the quadratic nodal function,
! and W(K) = ((R-D)+/(R*D))**2 for radius R(K) and distance D(K).  
! Thus, QX = (SWQX*SW - SWQ*SWX)/SW**2  and
!       QY = (SWQY*SW - SWQ*SWY)/SW**2
! where SWQX and SWX are partial derivatives with respect
! to X of SWQ and SW, respectively. SWQY and SWY are 
! defined similarly.
!
    SW = 0.0E+00_R8
    SWX = 0.0E+00_R8
    SWY = 0.0E+00_R8
    SWQ = 0.0E+00_R8
    SWQX = 0.0E+00_R8
    SWQY = 0.0E+00_R8
!
! Outer loop on cells (I,J).
!
    DO J = JMIN, JMAX
        DO I = IMIN, IMAX
            K = LCELL(I,J)
!
!  Inner loop on nodes K.
!
            IF (K /= 0) THEN
                EVAL: DO
                    DELX = XP - X(K)
                    DELY = YP - Y(K)
                    DS = DELX * DELX + DELY * DELY
                    RS = RSQ(K)
                    IF (DS == 0.0E+00_R8) THEN
                        Q = F(K)
                        QX = A(4,K)
                        QY = A(5,K)
                        IER = 0
                        RETURN
                    END IF
                    IF (DS < RS) THEN
                        RDS = RS * DS
                        RD = SQRT (RDS)
                        W = (RS + DS - RD - RD) / RDS
                        T = 2.0E+00_R8 * (RD - RS) / (DS * RDS)
                        WX = DELX * T
                        WY = DELY * T
                        QKX = 2.0E+00_R8 * A(1,K) * DELX + A(2,K) * DELY
                        QKY = A(2,K) * DELX + 2.0E+00_R8 * A(3,K) * DELY
                        QK = (QKX * DELX + QKY * DELY) / 2.0E+00_R8
                        QKX = QKX + A(4,K)
                        QKY = QKY + A(5,K)
                        QK = QK + A(4,K) * DELX + A(5,K) * DELY + F(K)
                        SW = SW + W
                        SWX = SWX + WX
                        SWY = SWY + WY
                        SWQ = SWQ + W * QK
                        SWQX = SWQX + WX * QK + W * QKX
                        SWQY = SWQY + WY * QK + W * QKY
                    END IF
                    KP = K
                    K = LNEXT(KP)
                    IF (K /= KP) THEN
                        CYCLE EVAL
                    END IF
                    EXIT
                END DO EVAL
            END IF
        END DO
    END DO
!
! SW = 0 if and only if P is not within the radius R(K) for any node K.
!
    IF (SW /= 0.0E+00_R8) THEN
        Q = SWQ / SW
        SWS = SW * SW
        QX = (SWQX * SW - SWQ * SWX) / SWS
        QY = (SWQY * SW - SWQ * SWY) / SWS
        IER = 0
    ELSE
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        IER = 2
    END IF
    RETURN
  END SUBROUTINE QS2GRD
!  
!=========================================================================
!
  SUBROUTINE GETNP2 (PX, PY, X, Y, N, NR, LCELL, LNEXT, XMIN, YMIN, &
       DX, DY, NP, DSQ)
!
! GETNP2 finds the closest unmarked node to a point.
!
! GETNP2 uses the cell method to find the closest unmarked node NP
! to a specified point P, given a set of N nodes and the data structure 
! defined by the subroutine STORE2.
!
! NP is then marked by negating LNEXT(NP). Thus, the closest M nodes to
! P may be determined by a sequence of M calls to this routine.  
!
! If the point P is itself actually a node K, and you want to find the
! nearest point to P that is not node K, then you must be sure to mark
! node K before calling GETNP2.
!
! The search is begun in the cell containing or closest to P and proceeds 
! outward in rectangular layers until all cells which contain points 
! within distance R of P have been searched. R is the distance from P to 
! the first unmarked node encountered, or infinite if no unmarked nodes
! are present.
!
! Input parameters other than LNEXT are not altered by this routine.  
! With the exception of (PX,PY) and the signs of LNEXT elements, 
! these parameters should not be altered from their values on output 
! from the subroutine STORE2.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, real, the (X,Y) coordinates of the point P whose
!   nearest unmarked neighbor is to be found.
!
!   X(N), Y(N), real, the coordinates of the nodes at which
!   data has been supplied.
!
!   NR, integer, the number of rows and columns in the cell grid.
!   NR must be at least 1.
!
!   LCELL(NR,NR), integer, array of nodal indices associated
!   with cells.
!
!   LNEXT(N), integer, contains next-node indices (or their 
!   negatives). On return, if the output value of NP is nonzero, then
!   LNEXT(NP) will be negated (changed to a negative value).
!
!   XMIN, YMIN, DX, DY, real, the minimum nodal X,Y coordinates,
!   and the X,Y dimensions of a cell. DX and DY must be positive.
!
!  OUTPUT
!   NP, integer, the index into the vectors X and Y of the nearest
!   unmarked node to the point P.  NP will be 0 if all nodes are marked 
!   or if the values of NR, DX, DY are illegal.  LNEXT(NP) will be less
!   than 0 if NP is nonzero (this marks node NP as being used now).
!
!   DSQ, real, if NP is nonzero, then DSQ is the squared distance
!   between P and node NP.
!
!   Local variables:
!
!   FIRST = TRUE iff the first unmarked node has yet to be encountered.
!
!   IMIN, IMAX, JMIN, JMAX = cell indices defining the range of the search.
!
!   DELX, DELY = PX-XMIN AND PY-YMIN.
!
!   I0, J0 = cell containing or closest to P.
!
!   I1, I2, J1, J2 = cell indices of the layer whose intersection with 
!   the range defined by IMIN,...,JMAX is currently being searched.
!
    USE REAL_PRECISION
    IMPLICIT NONE
    
    INTEGER :: I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2, JMAX, JMIN, &
          L, LMIN, LN, N, NP, NR
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION(NR,NR) :: LCELL
  
    REAL(KIND=R8) :: DELX, DELY, DSQ, DX, DY, PX, PY, R, RSMIN, RSQ, &
          XMIN, XP, YMIN, YP
    REAL(KIND=R8), DIMENSION(N) :: X, Y
  
    LOGICAL FIRST

    XP = PX
    YP = PY
!
!  Test for invalid input parameters.
!
    IF (NR < 1 .OR. DX <= 0.0E+00_R8 .OR. DY <= 0.0E+00_R8) THEN
        NP = 0
        DSQ = 0.0E+00_R8
    END IF
!
!  Initialize parameters.
!
    FIRST = .TRUE.
    IMIN = 1
    IMAX = NR
    JMIN = 1
    JMAX = NR
    DELX = XP - XMIN
    DELY = YP - YMIN

    I0 = INT (DELX / DX) + 1
    I0 = MAX (I0, 1)
    I0 = MIN (I0, NR)

    J0 = INT (DELY / DY) + 1
    J0 = MAX (J0, 1)
    J0 = MIN (J0, NR)

    I1 = I0
    I2 = I0
    J1 = J0
    J2 = J0
!
! Outer loop on layers, inner loop on layer cells, excluding
! those outside the range (IMIN,IMAX) x (JMIN,JMAX).
!
    MAIN: DO

        J_LOOP: DO J = J1, J2
            IF (J > JMAX) EXIT
            IF (J < JMIN) CYCLE J_LOOP

            I_LOOP: DO I = I1, I2
                IF (I > IMAX) EXIT
                IF (I < IMIN) CYCLE I_LOOP

                IF (J /= J1 .AND. J /= J2 .AND. I /= I1 .AND. I /= I2) THEN
                    CYCLE I_LOOP
                END IF
!
! Search cell (I,J) for unmarked nodes L.
!
                L = LCELL(I,J)
                IF (L > 0) THEN
!
! Loop on nodes in cell (i,j).
!
                    CELL_LOOP: DO
                        LN = LNEXT(L)
!
! Node L is the first unmarked neighbor of P encountered.
!
! Initialize LMIN to the current candidate for NP, and
! RSMIN to the squared distance from P to LMIN. IMIN,
! IMAX, JMIN, and JMAX are updated to define the smallest 
! rectangle containing a circle of radius R = SQRT(RSMIN) 
! centered at P, and contained in (1,NR) x (1,NR) 
! (except that, if P is outside the rectangle
! defined by the nodes, it is possible that 
! IMIN > NR, IMAX < 1, JMIN > NR, or JMAX < 1).
!
                        IF (LN >= 0) THEN
                            RSQ = (X(L) - XP)**2 + (Y(L) - YP)**2
                            IF (FIRST) THEN
                                LMIN = L
                                RSMIN = RSQ
                                R = SQRT (RSMIN)
                                IMIN = INT ((DELX - R) / DX) + 1
                                IMIN = MAX (IMIN, 1 )
                                IMAX = INT (( DELX + R) / DX) + 1
                                IMAX = MIN (IMAX, NR)
                                JMIN = INT (( DELY - R) / DY) + 1
                                JMIN = MAX (JMIN, 1)
                                JMAX = INT (( DELY + R) / DY) + 1
                                JMAX = MIN (JMAX, NR)
                                FIRST = .FALSE.
                            ELSE
                                IF ( RSQ < RSMIN ) THEN
                                    LMIN = L
                                    RSMIN = RSQ
                                END IF
                            END IF
                        END IF
                        IF (ABS (LN) /= L) THEN
                            L = ABS (LN)
                            CYCLE CELL_LOOP
                        END IF
                        EXIT 
                    END DO CELL_LOOP
                END IF  
            END DO I_LOOP
          END DO J_LOOP
!
!  Test for termination of loop on cell layers.
!
          IF (I1 > IMIN .OR. I2 < IMAX .OR. J1 > JMIN .OR. J2 < JMAX) THEN
              I1 = I1 - 1
              I2 = I2 + 1
              J1 = J1 - 1
              J2 = J2 + 1
              CYCLE MAIN
          END IF
          EXIT
      END DO MAIN
  
      IF (FIRST) THEN
          NP = 0
          DSQ = 0.0E+00_R8
      ELSE
          NP = LMIN
          DSQ = RSMIN
          LNEXT(LMIN) = -LNEXT(LMIN)
      END IF
      RETURN
  END SUBROUTINE GETNP2
!  
!=========================================================================
!  
  SUBROUTINE SETUP2 (XK, YK, FK, XI, YI, FI, S1, S2, R, ROW)
!
! SETUP2 sets up a row of the least squares regression matrix.
!
! SETUP2 sets up the I-th row of an augmented regression matrix for 
! a weighted least-squares fit of a quadratic function Q(X,Y) to a set 
! of data values F, where Q(XK,YK) = FK.  
!
! The first 3 columns are quadratic terms, and are scaled by 1/S2.
! The fourth and fifth columns represent linear terms, and are scaled 
! by 1/S1.  
!
! If D = 0, or D >= R, the weight is 0,
! else if D < R, the weight is (R-D)/(R*D), 
! where D is the distance between nodes I and K, and R is a maximum
! influence distance.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!  INPUT
!   XK, YK, FK, real, the coordinates and value of the data
!   at data node K.
!
!   XI, YI, FI, real, the coorindates and value of the data
!   at data node I.
!
!   S1, S2, real, reciprocals of the scale factors.
!
!   R, real, the maximum radius of influence about node K.
!
!  OUTPUT
!   ROW(6), real, a row of the augmented regression matrix.
!
    USE REAL_PRECISION
    IMPLICIT NONE

    REAL(KIND=R8):: D, DX, DY, FI, FK, R, S1, S2, W, XI, YI, XK, YK
    REAL(KIND=R8), DIMENSION(6):: ROW
    
    DX = XI - XK
    DY = YI - YK

    D = SQRT (DX * DX + DY * DY)

    IF (D <= 0.0E+00_R8 .OR. D >= R) THEN
        ROW(1:6) = 0.0E+00_R8
    ELSE
        W = ( R - D ) / R / D
        ROW(1) = DX * DX * W / S2
        ROW(2) = DX * DY * W / S2
        ROW(3) = DY * DY * W / S2
        ROW(4) = DX * W / S1
        ROW(5) = DY * W / S1
        ROW(6) = ( FI - FK ) * W
    END IF
    RETURN
  END SUBROUTINE SETUP2
!  
!=========================================================================
!
  SUBROUTINE STORE2 (N, X, Y, NR, LCELL, LNEXT, XMIN, YMIN, DX, DY, IER)
!
! STORE2 creates a cell data structure for the scattered data.
!
! STORE2 is given a set of N arbitrarily distributed nodes in the 
! plane and creates a data structure for a cell-based method of 
! solving closest-point problems. The smallest rectangle containing 
! all the nodes is partitioned into an NR by NR uniform grid of cells, 
! and nodes are associated with cells.      
!
! In particular, the data structure stores the indices of the nodes 
! contained in each cell. For a uniform random distribution of nodes, 
! the nearest node to an arbitrary point can be determined in constant
! expected time.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to Fortran 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
! 
! Parameters:
!  
!  INPUT
!   N, integer, the number of data nodes. N must be at least 2.
!
!   X(N), Y(N), real, the coordinates of the data nodes.
!
!   NR, integer, the number of rows and columns in the grid. The
!   cell density, or average number of data nodes per cell, is
!     D = N / (NR * NR).
!   A recommended value, based on empirical evidence, is  D = 3. 
!   Hence, the corresponding value of NR is recommended to be about
!     NR = SQRT (N / 3).  NR must be at least 1.
!
!  OUTPUT
!
!   LCELL(NR,NR), integer, an array set up so that LCELL(I,J)
!   contains the index (for X and Y) of the first data node (that is, the
!   data node with smallest index) in the (I,J) cell. LCELL(I,J) will be 0 
!   if no data nodes are contained in the (I,J) cell. The upper right 
!   corner of the (I,J) cell has coordinates (XMIN + I * DX, YMIN + J * DY).
!
!   LNEXT(N), integer, an array of next-node indices. LNEXT(K)
!   contains the index of the next node in the cell which contains node K, 
!   or LNEXT(K) = K if K is the last node in the cell.
!
!   The data nodes contained in a cell are ordered by their indices.
!   If, for example, cell (I,J) contains the nodes 2, 3, and 5 and no others, 
!   then:
!     LCELL(I,J) = 2, (index of the first data node)
!     LNEXT(2) = 3, 
!     LNEXT(3) = 5,
!     LNEXT(5) = 5.
!
!   XMIN, YMIN, real, the X, Y coordinates of the lower left
!   corner of the rectangle defined by the data nodes. The upper right 
!   corner is (XMAX,YMAX), where
!     XMAX = XMIN + NR * DX,
!     YMAX = YMIN + NR * DY.
!
!   DX, DY, real, the X and Y dimensions of the individual cells.
!     DX = (XMAX - XMIN) / NR
!     DY = (YMAX - YMIN) / NR,
!   where XMIN, XMAX, YMIN and YMAX are the extrema of X and Y.
!
!    IER, integer, an error indicator.
!        0, if no errors were encountered.
!        1, if N < 2 or NR < 1.
!        2, if DX = 0 or DY = 0.
!
    USE REAL_PRECISION
    IMPLICIT NONE
  
    INTEGER :: I, IER, J, K, L, N, NR
    INTEGER, DIMENSION(N) :: LNEXT
    INTEGER, DIMENSION (NR,NR) :: LCELL

    REAL(KIND=R8) :: DX, DY, XMAX, XMIN, YMAX, YMIN
    REAL(KIND=R8), DIMENSION(N) :: X, Y
  
    IER = 0
    IF (N < 2) THEN
        IER = 1
        RETURN
    END IF
    IF (NR < 1) THEN
        IER = 1
        RETURN
    END IF
!
! Compute the dimensions of the (X,Y) rectangle containing all the data nodes.
!
    XMIN = MINVAL (X(1:N))
    XMAX = MAXVAL (X(1:N))
    YMIN = MINVAL (Y(1:N))
    YMAX = MAXVAL (Y(1:N))
!
! Compute the dimensions of a single cell.
!
    DX = (XMAX - XMIN) / REAL (NR, KIND=R8)
    DY = (YMAX - YMIN) / REAL (NR, KIND=R8)
!
! Test for zero area.
!
    IF ( DX == 0.0E+00_R8 .OR. DY == 0.0E+00_R8 ) THEN
        IER = 2
        RETURN
    END IF
!
! Initialize LCELL.
!
    DO J = 1, NR
        DO I = 1, NR
            LCELL(I,J) = 0
        END DO
    END DO
!
! Loop on nodes, storing indices in LCELL and LNEXT.
!
    DO K = N, 1, -1
        I = INT (( X(K) - XMIN) / DX) + 1
        I = MIN (I, NR)
        J = INT (( Y(K) - YMIN) / DY) + 1
        J = MIN (J, NR)
        L = LCELL(I,J)
        IF (L /= 0) THEN
            LNEXT(K) = L
        ELSE
            LNEXT(K) = K
        END IF
        LCELL(I,J) = K
    END DO
    RETURN
  END SUBROUTINE STORE2
!  
!=========================================================================
!
END MODULE QSHEP2D_MOD
