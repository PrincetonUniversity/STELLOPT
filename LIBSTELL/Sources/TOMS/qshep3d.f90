  MODULE QSHEP3D_MOD
  PRIVATE
  PUBLIC QSHEP3, QS3VAL, QS3GRD
  CONTAINS
  SUBROUTINE QSHEP3 (N, X, Y, Z, F, NQ, NW, NR, LCELL, LNEXT, &
    XYZMIN, XYZDEL, RMAX, RSQ, A, IER)
!
! QSHEP3 defines a smooth trivariate interpolant of scattered 3D data.
!
! Discussion:
!
!   This subroutine computes a set of parameters, A and RSQ,
!   defining a smooth (once continuously differentiable) trivariate 
!   function Q(X,Y,Z) which interpolates data values
!   F at scattered nodes (X,Y,Z). The interpolant Q may be
!   evaluated at an arbitrary point by the function QS3VAL, and
!   its first derivatives are computed by the subroutine QS3GRD.
!
!   The interpolation scheme is a modified quadratic Shepard
!   method. 
!
!   Q = (W(1)*Q(1)+W(2)*Q(2)+...+W(N)*Q(N)) / (W(1)+W(2)+...+W(N))
!   for trivariate functions W(K) and Q(K).  
!
!   The nodal functions are given by
!     Q(K)(X,Y,Z) =  A(1,K) * DX**2 
!                  + A(2,K) * DX * DY 
!                  + A(3,K) * DY**2
!                  + A(4,K) * DX * DZ
!                  + A(5,K) * DY * DZ 
!                  + A(6,K) * DZ**2
!                  + A(7,K) * DX 
!                  + A(8,K) * DY 
!                  + A(9,K) * DZ 
!                  + F(K), 
!
!   where DX = (X-X(K)), DY = (Y-Y(K)), and DZ = (Z-Z(K)).
!
!   Thus, Q(K) is a quadratic function which interpolates the
!   data value at node K. Its coefficients A(*,K) are obtained
!   by a weighted least squares fit to the closest NQ data
!   points with weights similar to W(K). Note that the radius
!   of influence for the least squares fit is fixed for each
!   K, but varies with K.
!
!   The weights are taken to be
!     W(K)(X,Y,Z) = ( (R(K)-D(K))+ / R(K)*D(K) )**2, 
!
!   where (R(K)-D(K))+ = 0 if R(K) <= D(K), and D(K)(X,Y,Z)
!   is the Euclidean distance between (X,Y,Z) and node K. The
!   radius of influence R(K) varies with K and is chosen so
!   that NW nodes are within the radius. Note that W(K) is
!   not defined at node (X(K),Y(K),Z(K)), but Q(X,Y,Z) has
!   limit F(K) as (X,Y,Z) approaches (X(K),Y(K),Z(K)).
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT 
!   N, integer, the number of nodes and associated data 
!   values. N >= 10.
!
!   X(N), Y(N), Z(N), real, the coordinates of the nodes.
!
!   F(N), real, the data values at the nodes.
!
!   NQ, integer, the number of data points to be used in the least
!   squares fit for coefficients defining the nodal functions Q(K).  
!   A recommended value is NQ = 17. 9 <= NQ <= MIN (40, N-1).
!
!   NW, integer, the number of nodes within (and defining) the radii
!   of influence R(K), which enter into the weights W(K). For N sufficiently
!   large, a recommended value is NW = 32. 1 <= NW <= MIN(40,N-1).
!
!   NR, integer, the number of rows, columns, and planes in the cell
!   grid defined in the subroutine STORE3. A box containing the nodes is
!   partitioned into cells in order to increase search efficiency.  
!   NR = (N/3)**(1/3) is recommended. NR >= 1.
!
!  OUTPUT
!   LCELL(NR,NR,NR), integer, nodal indices associated with cells.  
!   Refer to STORE3.
!
!   LNEXT(N), integer, next-node indices. Refer to STORE3.
!
!   XYZMIN(3), XYZDEL(3) real,  arrays of length 3 containing 
!   minimum nodal coordinates and cell dimensions, respectively.  
!   Refer to STORE3.
!
!   RMAX real,  square root of the largest element in RSQ,
!   the maximum radius of influence.
!
!   RSQ(N) real,  array containing the squares of the radii R(K),
!   which enter into the weights W(K).
!
!   A(9,N), real, the coefficients for the quadratic nodal 
!   function Q(K) in column K.
!
!   integer IER, error indicator.
!     0, if no errors were encountered.
!     1, if N, NQ, NW, or NR is out of range.
!     2, if duplicate nodes were encountered.
!     3, if all nodes are coplanar.
!
! Local variables:
!
!   AV =        Root-mean-square distance between node K and the nodes in the 
!               least squares system (unless additional nodes are introduced
!               for stability). The first 6 columns of the matrix are 
!               scaled by 1/AVSQ, the last 3 by 1/AV.
!
!   AVSQ =      AV*AV.
!
!   B =         Transpose of the augmented regression matrix.
!
!   C =         First component of the plane rotation used to zero 
!               the lower triangle of B**T, computed by the subroutine GIVENS.
!
!   DMIN =      Minimum of the magnitudes of the diagonal elements of 
!               the regression matrix after zeros are introduced below 
!               the diagonal.
!
!   DTOL =      Tolerance for detecting an ill-conditioned system. 
!               The system is accepted when DMIN >= DTOL.
!
!   FK =        Data value F(K) at node K. 
!
!   I =         Index for A, B, NPTS, XYZMIN, XYZMN, XYZDEL, and XYZDL.
!
!   IB =        Do-loop index for back solve.
!
!   IERR =      Error flag for the call to STORE3.
!
!   IP1 =       I+1.
!
!   IRM1 =      IROW-1.
!
!   IROW =      Row index for B.
!
!   J =         Index for A and B.
!
!   JP1 =       J+1.
!
!   K =         Nodal function index and column index for A.
!
!   LMAX =      Maximum number of NPTS elements (must be consistent 
!               with the dimension of NPTS in the variable declaration section).
!
!   LNP =       Current length of NPTS.
!
!   NEQ =       Number of equations in the least squares fit.
!
!   NN =        Local copy of N.
!
!   NNQ =       Local copy of NQ.
!
!   NNR =       Local copy of NR.
!
!   NNW =       Local copy of NW.
!
!   NP =        NPTS element.
!
!   NPTS =      Array containing the indices of a sequence of nodes to be
!               used in the least squares fit or to compute RSQ. The nodes 
!               are ordered by distance from K and the last element
!               (usually indexed by LNP) is used only to determine RQ, 
!               or RSQ(K) if NW > NQ.
!
!   NQWMAX =    MAX(NQ,NW).
!
!   RQ =        Radius of influence, which enters into the weights for Q(K).
!               (See subroutine SETUP3).
!
!   RS =        Squared distance between node K and NPTS(LNP). 
!               Used to compute RQ and RSQ(K).
!
!   RSMX =      Maximum RSQ element encountered.
!
!   RSOLD =     Squared distance between node K and NPTS(LNP-1). Used to 
!               compute a relative change in RS between succeeding 
!               NPTS elements.
!
!   RTOL =      Tolerance for detecting a sufficiently large relative 
!               change in RS. If the change is not greater than RTOL,
!               the nodes are treated as being the same distance from K.
!
!   RWS =       Current value of RSQ(K).
!
!   S =         Second component of the plane Givens rotation.
!
!   SF =        Marquardt stabilization factor used to damp out the 
!               first 6 solution components (second partials of the 
!               quadratic) when the system is ill-conditioned. As 
!               SF increases, the fitting function becomes more linear.
!
!   SUM2 =      Sum of squared Euclidean distances between node K and 
!               the nodes used in the least squares fit 
!               (unless additional nodes are added for stability).
!
!   T =         Temporary variable for accumulating a scalar product 
!               in the back solve.
!
!   XK,YK,ZK =  Coordinates X(K), Y(K), Z(K) of node K.
!
!   XYZDL =     Local variables for XYZDEL.
!
!   XYZMN =     Local variables for XYZMIN.
!
    USE REAL_PRECISION
    IMPLICIT NONE
    
    INTEGER :: I, IB, IER, IERR, IP1, IRM1, IROW, J, JP1, K, LMAX, LNP, &
        N, NEQ, NN, NNQ, NNR, NNW, NP, NQ, NQWMAX, NR, NW
    INTEGER, DIMENSION (N) :: LNEXT
    INTEGER, DIMENSION (5000) :: NPTS
    !INTEGER, DIMENSION (40) :: NPTS
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL

    REAL(KIND=R8):: AV, AVSQ, C, DMIN
    REAL(KIND=R8):: DTOL
    REAL(KIND=R8):: FK, RMAX, RQ, RS, RSMX, RSOLD
    REAL(KIND=R8), PARAMETER :: RTOL = 1.0E-05
    REAL(KIND=R8):: RWS, S
    REAL(KIND=R8), PARAMETER :: SF = 1.0E+00
    REAL(KIND=R8):: SUM2, T, XK, YK, ZK
    REAL(KIND=R8), DIMENSION (N) :: F, RSQ, X
    REAL(KIND=R8), DIMENSION (3) :: XYZDEL, XYZDL, XYZMIN, XYZMN
    REAL(KIND=R8), DIMENSION (N) :: Y, Z
    REAL(KIND=R8), DIMENSION (9,N) :: A
    REAL(KIND=R8), DIMENSION (10,10) :: B

    DTOL = SQRT(EPSILON(1.0_R8))
    NN = N
    NNQ = NQ
    NNW = NW
    NNR = NR
    NQWMAX = MAX(NNQ,NNW)
    !LMAX = MIN(40,NN-1)
    LMAX = MIN(5000,NN-1)

    IF (9 > NNQ .OR.  1 > NNW  .OR.  NQWMAX > &
          LMAX  .OR.  NNR < 1) THEN
!
! N, NQ, NW, or NR is out of range.
!
        IER = 1
        RETURN
    END IF
!
! Create the cell data structure, and initialize RSMX.
!
    CALL STORE3 (NN, X, Y, Z, NNR, LCELL, LNEXT, XYZMN, XYZDL, IERR)

    IF (IERR /= 0) THEN
        XYZMIN(1:3) = XYZMN(1:3)
        XYZDEL(1:3) = XYZDL(1:3)
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
        ZK = Z(K)
        FK = F(K)
!
!  Mark node K to exclude it from the search for nearest neighbors.
!
        LNEXT(K) = -LNEXT(K)
!
!  Initialize for loop on NPTS.
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
! distance RS to the last node included.
!
                IF (RWS == 0.0E+00_R8) RWS = 1.1E+00_R8 * RS
                IF (RQ == 0.0E+00_R8) THEN
                    NEQ = LMAX
                    RQ = SQRT (1.1_R8 * RS)
                    AVSQ = SUM2 / REAL (NEQ, KIND=R8)
                END IF
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                RSQ(K) = RWS
                IF (RWS > RSMX) RSMX = RWS
                AV = SQRT(AVSQ)
                EXIT
            END IF    
            LNP = LNP + 1
            RSOLD = RS
            CALL GETNP3 (N, XK, YK, ZK, X, Y, Z, NNR, LCELL, LNEXT, XYZMN, &
                  XYZDL, NP, RS)
            IF (RS == 0.0E+00_R8) THEN
!
! Duplicate nodes were encountered by GETNP3.
!
                IER = 2
                RETURN
            END IF
            NPTS(LNP) = NP
            IF ((RS-RSOLD)/RS <  RTOL) CYCLE MAIN
            IF (RWS == 0.0E+00_R8 .AND.  LNP > NNW) RWS = RS
            IF (.NOT.( RQ /= 0.0E+00_R8  .OR.  LNP <= NNQ )) THEN
!
! RQ = 0 (not yet computed) and LNP > NQ.  
! RQ = SQRT(RS) is sufficiently large to (strictly) include
! NQ nodes. The least squares fit will include 
! NEQ = LNP-1 equations for 9 <= NQ <= NEQ < LMAX <= N-1.
!
                NEQ = LNP - 1
                RQ = SQRT(RS)
                AVSQ = SUM2 / REAL (NEQ, KIND=R8)
!
! Bottom of loop, test for termination.
!
            ELSE
                IF (LNP > NQWMAX) THEN
!
! Store RSQ(K), update RSMX if necessary, and compute AV.
!
                    RSQ(K) = RWS
                    IF (RWS > RSMX) RSMX = RWS
                    AV = SQRT(AVSQ)
                    EXIT
                END IF
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
            IROW = MIN0(I,10)
            CALL  SETUP3 (XK, YK, ZK, FK, X(NP), Y(NP), Z(NP), F(NP),&
                    AV, AVSQ, RQ, B(1,IROW))
            IF (I == 1) CYCLE PTS
            IRM1 = IROW-1
            DO J = 1, IROW-1
                JP1 = J + 1
                CALL GIVENS (B(J,J),B(J,IROW),C,S)
                CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,IROW))
            END DO
            IF (I < NEQ) CYCLE PTS
!
! Test the system for ill-conditioning.
!
            DMIN = MIN (ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)), &
                    ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)), &
                    ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)))

            IF (DMIN * RQ >= DTOL) EXIT
            IF (NEQ == LMAX) EXIT
!
! Increase RQ and add another equation to the system to
! improve the conditioning. The number of NPTS elements
! is also increased if necessary.
!
            TOL: DO
                RSOLD = RS
                NEQ = NEQ + 1
                IF (NEQ == LMAX) THEN
                    RQ = SQRT (1.1E+00_R8 * RS)
                    CYCLE PTS
                END IF
                IF (NEQ == LNP) THEN
!
! Add an element to NPTS.
!
                    LNP = LNP + 1
                    CALL GETNP3 (N, XK, YK, ZK, X, Y, Z, NNR, LCELL, LNEXT, &
                          XYZMN, XYZDL, NP, RS)
                    IF (NP == 0) THEN
!
! Duplicate nodes were encountered by GETNP3.
!
                        IER = 2
                        RETURN
                    END IF
                    NPTS(LNP) = NP
                    IF ((RS-RSOLD)/RS < RTOL) CYCLE TOL
                    RQ = SQRT (RS)
                    CYCLE PTS
                END IF
!
! NEQ < LNP.
!
                NP = NPTS(NEQ+1)
                RS = (X(NP)-XK)**2 + (Y(NP)-YK)**2 + (Z(NP)-ZK)**2
                IF ((RS-RSOLD)/RS < RTOL) CYCLE TOL
                RQ = SQRT(RS)
                CYCLE PTS
            END DO TOL
        END DO PTS
!
! Stabilize the system by damping second partials. Add
! multiples of the first six unit vectors to the first
! six equations.
!
        IF (NEQ == LMAX) THEN
            DO I = 1, 6
                B(I,10) = SF
                IP1 = I + 1
                DO J = IP1, 10
                    B(J,10) = 0.0E+00_R8
                END DO
                DO J = I, 9
                    JP1 = J + 1
                    CALL GIVENS (B(J,J),B(J,10),C,S)
                    CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,10))
                END DO
            END DO
!
! Test the stabilized system for ill-conditioning.
!
            DMIN = MIN ( ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)), &
                  ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)), &
                  ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)) )
            IF (DMIN * RQ < DTOL) THEN
!
! No unique solution due to collinear nodes.
!
                XYZMIN(1:3) = XYZMN(1:3)
                XYZDEL(1:3) = XYZDL(1:3)
                PRINT *,' ERROR HERE'
                IER = 3
                RETURN
            END IF
        END IF
!
! Solve the 9 by 9 triangular system for the coefficients.
!

        DO IB = 1, 9
            I = 10-IB
            T = 0.0E+00_R8
            DO J = I+1, 9
                T = T + B(J,I)*A(J,K)
            END DO
            A(I,K) = (B(10,I)-T)/B(I,I)
        END DO
!
! Scale the coefficients to adjust for the column scaling.
!
        A(1:6,K) = A(1:6,K) / AVSQ
        A(7,K) = A(7,K)/AV
        A(8,K) = A(8,K)/AV
        A(9,K) = A(9,K)/AV
!
! Unmark K and the elements of NPTS.
!
        LNEXT(K) = -LNEXT(K)
        DO I = 1, LNP
            NP = NPTS(I)
            LNEXT(NP) = -LNEXT(NP)
        END DO
    END DO
!
! No errors encountered.
!
    XYZMIN(1:3) = XYZMN(1:3)
    XYZDEL(1:3) = XYZDL(1:3)
    RMAX = SQRT (RSMX)
    IER = 0
    RETURN
  END SUBROUTINE QSHEP3
!  
!=========================================================================
!  
  SUBROUTINE STORE3 (N, X, Y, Z, NR, LCELL, LNEXT, XYZMIN, XYZDEL, IER)
!
! STORE3 sets up a data structure for N scattered nodes in 3D.
!
! Discussion:
!
!   Given a set of N arbitrarily distributed nodes in three-space, 
!   this subroutine creates a data structure for a cell-based method of 
!   solving closest-point problems. The smallest box containing the nodes 
!   is partitioned into an NR by NR by NR uniform grid of cells, and 
!   nodes are associated with cells. In particular, the data structure
!   stores the indices of the nodes contained in each cell. For a 
!   uniform random distribution of nodes, the nearest node to an 
!   arbitrary point can be determined in constant expected time.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   N, integer, the number of nodes.  N >= 2.
!
!   X(N), Y(N), Z(N), real, the coordinates of the nodes.
!
!   NR, integer, the number of rows, columns, and planes in the
!   grid. The cell density (average number of nodes per cell) is 
!   D = N/(NR**3). A recommended value, based on empirical evidence, 
!   is D = 3, so NR = (N/3)**(1/3). NR >= 1.
!
!  OUTPUT
!   LCELL(NR,NR,NR), integer, a cell array such that 
!   LCELL(I,J,K) contains the index (for X, Y, and Z) of the 
!   first node (node with smallest index) in cell (I,J,K), or 
!   LCELL(I,J,K) = 0 if no nodes are contained in the cell. 
!   The corner of cell (I,J,K) which is farthest from the box 
!   corner defined by XYZMIN has coordinates XMIN+I*DX, 
!   YMIN+J*DY, ZMIN+K*DZ), where (XMIN, YMIN, ZMIN) are 
!   the elements of XYZMIN. LCELL is not defined if IER /= 0.
!
!   LNEXT(N), integer, next-node indices such that
!   LNEXT(L) contains the index of the next node in the cell
!   which contains node L, or LNEXT(L) = L if L is the last 
!   node in the cell for L = 1,...,N. (The nodes contained
!   in a cell are ordered by their indices.)
!
!     If, for example, cell (I,J,K) contains nodes
!     2, 3, and 5 (and no others), then LCELL(I,J,K) = 2, 
!     LNEXT(2) = 3, LNEXT(3) = 5, and LNEXT(5) = 5.  
!     LNEXT is not defined if IER /= 0.
!
!   XYZMIN(3), real, the minimum nodal coordinates
!   XMIN, YMIN, and ZMIN (in that order) unless IER = 1.  
!   The opposite corner of the box defined by the nodes is
!   (XMIN+NR*DX, YMIN+NR*DY, ZMIN+NR*DZ).
!
!   XYZDEL(3), real, the dimensions of the cells 
!   unless IER = 1. XYZDEL(1) = (XMAX-XMIN)/NR, 
!                   XYZDEL(2) = (YMAX-YMIN)/NR,
!               and XYZDEL(3) = (ZMAX-ZMIN)/NR, 
!   where XMIN, XMAX, YMIN, YMAX, ZMIN, and ZMAX are the
!   extrema of X, Y, and Z.
!
!   IER, integer, error indicator.
!     0, if no errors were encountered.
!     1, if N < 2 or NR < 1.
!     2, if a component of XYZDEL is not positive.
!
    USE REAL_PRECISION    
    IMPLICIT NONE         

    INTEGER :: I, IER, J, K, L, LB, LL, N, NN, NR, NNR, NP1
    INTEGER, DIMENSION (N) :: LNEXT            
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL     

    REAL(KIND=R8) :: DELX, DELY, DELZ, XMN, XMX, YMN, YMX, ZMN, ZMX
    REAL(KIND=R8), DIMENSION (N) :: X                    
    REAL(KIND=R8), DIMENSION (3) :: XYZDEL, XYZMIN       
    REAL(KIND=R8), DIMENSION (N) :: Y, Z                 

    NN = N
    NNR = NR
    IF (NN < 2 .OR. NNR < 1) THEN
        IER = 1
        RETURN
    END IF
!
! Compute the dimensions of the box containing the nodes.
!
    XMN = MINVAL (X(1:NN))
    XMX = MAXVAL (X(1:NN))
    YMN = MINVAL (Y(1:NN))
    YMX = MAXVAL (Y(1:NN))
    ZMN = MINVAL (Z(1:NN))
    ZMX = MAXVAL (Z(1:NN))
    XYZMIN(1) = XMN
    XYZMIN(2) = YMN
    XYZMIN(3) = ZMN
!
! Compute cell dimensions and test for zero area.
!
    DELX = (XMX-XMN)/REAL(NNR,KIND=R8)
    DELY = (YMX-YMN)/REAL(NNR,KIND=R8)
    DELZ = (ZMX-ZMN)/REAL(NNR,KIND=R8)
    XYZDEL(1) = DELX
    XYZDEL(2) = DELY
    XYZDEL(3) = DELZ
    IF (DELX == 0.0E+00_R8 .OR. DELY == 0.0E+00_R8 .OR. &
          DELZ == 0.0E+00_R8) THEN
        IER = 2
        RETURN
    END IF
!
! Initialize LCELL.
!
    LCELL(1:NNR,1:NNR,1:NNR) = 0
!
! Loop on nodes, storing indices in LCELL and LNEXT.
!
    NP1 = NN + 1
    DO LL = 1, NN
        LB = NP1 - LL
        I = INT((X(LB)-XMN)/DELX) + 1
        IF (I > NNR) I = NNR
        J = INT((Y(LB)-YMN)/DELY) + 1
        IF (J > NNR ) J = NNR
        K = INT((Z(LB)-ZMN)/DELZ) + 1
        IF (K > NNR) K = NNR
        L = LCELL(I,J,K)
        LNEXT(LB) = L
        IF (L == 0) LNEXT(LB) = LB
        LCELL(I,J,K) = LB
    END DO
    IER = 0
    RETURN
  END SUBROUTINE STORE3
!  
!=========================================================================
!  
  SUBROUTINE GETNP3 (N, PX, PY, PZ, X, Y, Z, NR, LCELL, LNEXT, XYZMIN, &
    XYZDEL, NP, DSQ)
!
! GETNP3 finds the closest node to a given point.
!
! Discussion:
!
!   Given a set of N nodes and the data structure defined in
!   subroutine STORE3, this subroutine uses the cell method to
!   find the closest unmarked node, NP, to a specified point P
!   (PX,PY,PZ).
!
!   Node NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  
!   (A node is marked if and only if the corresponding LNEXT element 
!   is negative.  The absolute values of LNEXT elements,
!   however, must be preserved.)  
!   Thus, the closest M nodes to P may be determined by a sequence
!   of M calls to this routine. Note that if the nearest neighbor to
!   node K is to be determined (PX = X(K), PY = Y(K), and PZ = Z(K)), 
!   then node K should be marked before the call to this routine.
!
!   The search is begun in the cell containing (or closest
!   to) node P and proceeds outward in box-shaped layers until all
!   cells which contain points within distance R of node P have
!   been searched, where R is the distance from node P to the first
!   unmarked node encountered (infinite if no unmarked nodes
!   are present).
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, PZ, real, the coordinates of the point P whose 
!   nearest unmarked neighbor is to be found.
!
!   X(N), Y(N), Z(N), real, the coordinates of the nodes.
!
!   NR, integer, the number of rows, columns, and planes in the cell
!   grid. NR >= 1.
!
!   LCELL(NR,NR,NR), integer, nodal indices associated with cells.
!
!   LNEXT(N), integer, next-node indices (or their negatives). May have
!   a value negated upon return.
!
!   XYZMIN(3), XYZDEL(3), real, minimum nodal coordinates and cell 
!   dimensions, respectively. XYZDEL elements must be positive.
!
!  OUTPUT
!   NP, integer, index of the nearest unmarked node to P, or 0 
!   if all nodes are marked or NR < 1 or an element of XYZDEL is not 
!   positive. LNEXT(NP) < 0 if NP /= 0.
!
!   DSQ, real, squared Euclidean distance between P and node
!   NP, or 0 if NP = 0.
!
    USE REAL_PRECISION
    IMPLICIT NONE
    
    INTEGER :: I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2, JMAX, &
          JMIN, K, K0, K1, K2, KMAX, KMIN, L, LMIN, LN, N, NP, NR
    INTEGER, DIMENSION (N) :: LNEXT
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL
  
    REAL(KIND=R8):: DELX, DELY, DELZ, DSQ, DX, DY, DZ, PX, PY, &
          PZ, R, RSMIN, RSQ, XP, YP, ZP
    REAL(KIND=R8), DIMENSION (N) :: X
    REAL(KIND=R8), DIMENSION (3) :: XYZDEL, XYZMIN 
    REAL(KIND=R8), DIMENSION (N) :: Y, Z    
    
    LOGICAL FIRST

    XP = PX
    YP = PY
    ZP = PZ
    DX = XYZDEL(1)
    DY = XYZDEL(2)
    DZ = XYZDEL(3)
!
! Test for invalid input parameters.
!
    IF (NR < 1 .OR. DX <= 0.0E+00_R8 .OR. DY <= 0.0E+00_R8 & 
          .OR. DZ <= 0.0E+00_R8) THEN
        NP = 0
        DSQ = 0.0E+00_R8
        RETURN
    END IF
!
! Initialize parameters 
!
! FIRST = TRUE iff the first unmarked node has yet to be encountered,
!   IMIN,...,KMAX = cell indices defining the range of the search,
!   DELX, DELY, DELZ = PX-XYZMIN(1), PY-XYZMIN(2), and PZ-XYZMIN(3),
!   I0,J0,K0 = cell containing or closest to P,
!   I1,...,K2 = cell indices of the layer whose intersection
!   with the range defined by IMIN,...,KMAX is
!   currently being searched.
!
    FIRST = .TRUE.
    IMIN = 1
    IMAX = NR
    JMIN = 1
    JMAX = NR
    KMIN = 1
    KMAX = NR
    DELX = XP - XYZMIN(1)
    DELY = YP - XYZMIN(2)
    DELZ = ZP - XYZMIN(3)
    I0 = INT(DELX/DX) + 1
    IF (I0 < 1) I0 = 1
    IF (I0 > NR) I0 = NR
    J0 = INT (DELY / DY) + 1
    IF (J0 < 1) J0 = 1
    IF (J0 > NR) J0 = NR
    K0 = INT(DELZ/DZ) + 1
    IF (K0 < 1) K0 = 1
    IF (K0 > NR) K0 = NR
    I1 = I0
    I2 = I0
    J1 = J0
    J2 = J0
    K1 = K0
    K2 = K0
!
! Outer loop on layers, inner loop on layer cells, excluding
! those outside the range (IMIN,IMAX) X (JMIN,JMAX) X (KMIN,KMAX).
!
    MAIN: DO
        LOOP_K: DO K = K1, K2
            IF (K > KMAX) EXIT 
            IF (K < KMIN) CYCLE LOOP_K
            LOOP_J: DO J = J1, J2
                IF (J > JMAX) EXIT
                IF (J < JMIN) CYCLE LOOP_J
                LOOP_I: DO I = I1, I2
                    IF (I > IMAX) EXIT
                    IF (I < IMIN) CYCLE LOOP_I
                    IF (K /= K1  .AND.  K /= K2  .AND. J /= J1  .AND.  &
                          J /= J2  .AND. I /= I1  .AND.  I /= I2) CYCLE LOOP_I
!
! Search cell (I,J,K) for unmarked nodes L.
!
                    L = LCELL(I,J,K)
                    IF (L == 0) CYCLE LOOP_I
!
! Loop on nodes in cell (I,J,K).
!
                    CELL_LOOP: DO
                        LN = LNEXT(L)
                        IF (LN < 0) THEN 
                            IF (ABS(LN) == L) CYCLE LOOP_I
                            L = ABS(LN)
                            CYCLE CELL_LOOP
                        END IF
!
! Node L is not marked.
!
                        RSQ = (X(L)-XP)**2 + (Y(L)-YP)**2 + (Z(L)-ZP)**2
                        IF (.NOT. FIRST) THEN
!
! Update LMIN and RSMIN.
!
                            IF (RSQ < RSMIN) THEN
                                LMIN = L
                                RSMIN = RSQ
                            END IF
!
! Test for termination of loop on nodes in cell (i,j,k).
!
                            IF ( ABS(LN) == L ) CYCLE LOOP_I
                            L = ABS(LN)
                            CYCLE CELL_LOOP
                        END IF
!
! Node L is the first unmarked neighbor of node P encountered.
! Initialize LMIN to the current candidate for NP, and
! RSMIN to the squared distance from node P to node LMIN. 
! IMIN, IMAX, JMIN, JMAX, KMIN, and KMAX are updated to define
! the smallest rectangle containing a sphere of radius
! R = SQRT(RSMIN) centered at P, and contained in 
!(1,NR) X (1,NR) X (1,NR) (except that, if P is outside the
! box defined by the nodes, it is possible that IMIN > NR or 
! IMAX < 1, etc.). FIRST is reset to false.
!
                        LMIN = L
                        RSMIN = RSQ
                        R = SQRT(RSMIN)
                        IMIN = INT((DELX-R)/DX) + 1
                        IF (IMIN < 1) IMIN = 1
                        IMAX = INT((DELX+R)/DX) + 1
                        IF (IMAX > NR) IMAX = NR
                        JMIN = INT((DELY-R)/DY) + 1
                        IF (JMIN < 1) JMIN = 1
                        JMAX = INT((DELY+R)/DY) + 1
                        IF (JMAX > NR) JMAX = NR
                        KMIN = INT((DELZ-R)/DZ) + 1
                        IF (KMIN < 1) KMIN = 1
                        KMAX = INT((DELZ+R)/DZ) + 1
                        IF (KMAX > NR) KMAX = NR
                        FIRST = .FALSE.
                        IF (ABS(LN) == L) CYCLE LOOP_I
                        L = ABS(LN)
                        CYCLE CELL_LOOP
!
!  Test for node L closer than node LMIN to node P.
!
                    END DO CELL_LOOP
                END DO LOOP_I
            END DO LOOP_J
        END DO LOOP_K
!
! Test for termination of loop on cell layers.
!
        IF (I1 <= IMIN  .AND.  I2 >= IMAX  .AND. &
              J1 <= JMIN  .AND.  J2 >= JMAX  .AND. &
              K1 <= KMIN  .AND.  K2 >= KMAX) EXIT  
        I1 = I1 - 1
        I2 = I2 + 1
        J1 = J1 - 1
        J2 = J2 + 1
        K1 = K1 - 1
        K2 = K2 + 1
        CYCLE MAIN
    END DO MAIN
!
! Unless no unmarked nodes were encountered, LMIN is the
! closest unmarked node to node P.
!
    IF (.NOT. FIRST) THEN
        NP = LMIN
        DSQ = RSMIN
        LNEXT(LMIN) = -LNEXT(LMIN)
    ELSE
        NP = 0
        DSQ = 0.0E+00_R8
    END IF
    RETURN
  END SUBROUTINE GETNP3
!  
!=========================================================================
!  
  SUBROUTINE SETUP3 (XK, YK, ZK, FK, XI, YI, ZI, FI, S1, S2, R, ROW)
!
! SETUP3 sets up the weighted least-squares fit of the data.
!
! Discussion:
!   This routine sets up the I-th row of an augmented regression matrix 
!   for a weighted least-squares fit of a quadratic function Q(X,Y,Z) 
!   to a set of data values F, where Q(XK,YK,ZK) = FK.  
!
!   The first 6 columns (quadratic terms) are scaled by 1/S2, and 
!   columns 7, 8, and 9 (linear terms) are scaled by 1/S1.  
!   The weight is (R-D)/(R*D) if R > D, and 0 if R <= D, 
!   where D is the distance between nodes I and K.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   XK, YK, ZK, FK, real, coordinates and data value at node K
!    (interpolated by Q).
!
!   XI, YI, ZI, FI, real, coordinates and data value at node I.
!
!   S1, S2, real, reciprocals of the scale factors.
!
!   R, real, radius of influence about node K defining the weight.
!
!  OUTPUT
!   real ROW(10), a row of the augmented regression matrix.
!
! Local variables:
!
!    D =  Distance between nodes K and I.
!
!    W =  Weight associated with a row.
!
!    W1 = W/S1.
!
!    W2 = W/S2.
!
    USE REAL_PRECISION
    IMPLICIT NONE
    
    REAL(KIND=R8):: D, DX, DXSQ, DY, DYSQ, DZ, DZSQ, FI, FK, R, S1, S2, &
          W, W1, W2, XI, XK, YI, YK, ZI, ZK
    REAL(KIND=R8), DIMENSION(10) :: ROW
    
    DX = XI - XK
    DY = YI - YK
    DZ = ZI - ZK
    DXSQ = DX*DX
    DYSQ = DY*DY
    DZSQ = DZ*DZ
    D = SQRT (DXSQ + DYSQ + DZSQ)

    IF (D <= 0.0E+00_R8 .OR. D >= R) THEN
        ROW(1:10) = 0.0E+00_R8
        RETURN
    END IF
    W = (R - D) / R / D
    W1 = W/S1
    W2 = W/S2
    ROW(1) = DXSQ*W2
    ROW(2) = DX*DY*W2
    ROW(3) = DYSQ*W2
    ROW(4) = DX*DZ*W2
    ROW(5) = DY*DZ*W2
    ROW(6) = DZSQ*W2
    ROW(7) = DX*W1
    ROW(8) = DY*W1
    ROW(9) = DZ*W1
    ROW(10) = (FI - FK)*W
    RETURN
  END SUBROUTINE SETUP3
!  
!=========================================================================
!  
  SUBROUTINE QS3GRD (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, &
        LNEXT, XYZMIN, XYZDEL, RMAX, RSQ, A, Q, QX, QY, QZ, IER)
!
! QS3GRD computes the value and gradient of the interpolant function.
!
! Discussion:
!
! This subroutine computes the value and gradient at (PX,PY,PZ) of 
! the interpolatory function Q defined in the subroutine QSHEP3.  
!
! Q(X,Y,Z) is a weighted sum of quadratic nodal functions.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, PZ, real, the point P at which Q and its partials are 
!   to be evaluated.
!
!   N, integer, the number of nodes and data values defining Q.
!   N >= 10.
!
!   X(N), Y(N), Z(N), F(N), real, the node coordinates and
!   data values interpolated by Q.
!
!   NR, integer, the number of rows, columns and planes in the cell
!   grid. Refer to STORE3. NR >= 1.
!
!   LCELL(NR,NR,NR), integer, nodal indices associated with cells.  
!   Refer to STORE3.
!
!   LNEXT(N), integer, the next-node indices. Refer to STORE3.
!
!   XYZMIN(3), XYZDEL(3), real, the minimum nodal coordinates and 
!   cell dimensions, respectively. XYZDEL elements must be positive.  
!   Refer to STORE3.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining Q.
!
!   A(9,N), real, the coefficients for the nodal functions defining Q.
!
!  OUTPUT
!   Q, real, the value of Q at (PX,PY,PZ) unless IER == 1, in
!   which case no values are returned.
!
!   QX, QY, QZ, real, the first partial derivatives of Q at
!   (PX,PY,PZ) unless IER == 1.
!
!   IER, integer, error indicator
!      0, if no errors were encountered.
!      1, if N, NR, XYZDEL, or RMAX are invalid.
!      2, if no errors were encountered but (PX.PY.PZ) is not within the 
!         radius R(K) for any node K (and thus Q = QX = QY = QZ = 0).
!
    USE REAL_PRECISION
    IMPLICIT NONE
  
    INTEGER :: I, IER, IMAX, IMIN, J, JMAX, JMIN, K, KMAX, &
          KMIN, L, LP, N, NR
    INTEGER, DIMENSION (N) :: LNEXT
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL
    
    REAL(KIND=R8):: DELX, DELY, DELZ, DS, DX, DXSQ, DY, DYSQ, DZ, DZSQ, &
          PX, PY, PZ, Q, QL, QLX, QLY, QLZ, QX, QY, QZ, RD, RDS, RMAX, RS, &
          SW, SWQ, SWQX, SWQY, SWQZ, SWS, SWX, SWY, SWZ, T, W, WX, WY, WZ, &
          XMIN, XP, YMIN, YP, ZMIN, ZP
    REAL(KIND=R8), DIMENSION (N) :: F, RSQ, X
    REAL(KIND=R8), DIMENSION (3) :: XYZDEL, XYZMIN 
    REAL(KIND=R8), DIMENSION (N) :: Y, Z
    REAL(KIND=R8), DIMENSION (9,N) :: A
 
    XP = PX
    YP = PY
    ZP = PZ
    XMIN = XYZMIN(1) 
    YMIN = XYZMIN(2)
    ZMIN = XYZMIN(3)
    DX = XYZDEL(1)
    DY = XYZDEL(2)
    DZ = XYZDEL(3)

    IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0E+00_R8 &
          .OR.  DY <= 0.0E+00_R8  .OR.  DZ <= 0.0E+00_R8  .OR. &
          RMAX < 0.0E+00_R8) THEN
        IER = 1
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, JMAX, KMIN, and KMAX to cell indices
! defining the range of the search for nodes whose radii
! include P.  The cells which must be searched are those
! intersected by (or contained in) a sphere of radius RMAX
! centered at P.
!
    IMIN = INT((XP-XMIN-RMAX)/DX) + 1
    IMIN = MAX (IMIN, 1)
    IMAX = INT((XP-XMIN+RMAX)/DX) + 1
    IMAX = MIN (IMAX, NR)
    JMIN = INT((YP-YMIN-RMAX)/DY) + 1
    JMIN = MAX (JMIN, 1)
    JMAX = INT((YP-YMIN+RMAX)/DY) + 1
    JMAX = MIN (JMAX, NR)
    KMIN = INT((ZP-ZMIN-RMAX)/DZ) + 1
    KMIN = MAX (KMIN, 1)
    KMAX = INT((ZP-ZMIN+RMAX)/DZ) + 1
    KMAX = MIN (KMAX, NR)
!
! Test for no cells within the sphere of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX .OR. KMIN > KMAX) THEN
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        QZ = 0.0E+00_R8
        IER = 2
        RETURN
    END IF
!
! Q = SWQ/SW = sum(W(L)*Q(L))/sum(W(L)) where the sum is
! from L = 1 to N, Q(L) is the quadratic nodal function,
! and W(L) = ((R-D)+/(R*D))**2 for radius R(L) and distance D(L).  
! Thus
!   QX = (SWQX*SW - SWQ*SWX)/SW**2
!   QY = (SWQY*SW - SWQ*SWY)/SW**2
!   QZ = (SWQZ*SW - SWQ*SWZ)/SW**2,
! where SWQX and SWX are partial derivatives with respect
! to X of SWQ and SW, respectively. SWQY, SWY, SWQZ, and
! SWZ are defined similarly.
!
    SW = 0.0E+00_R8
    SWX = 0.0E+00_R8
    SWY = 0.0E+00_R8
    SWZ = 0.0E+00_R8
    SWQ = 0.0E+00_R8
    SWQX = 0.0E+00_R8
    SWQY = 0.0E+00_R8
    SWQZ = 0.0E+00_R8
!
! Outer loop on cells (I,J,K).
!
    DO K = KMIN, KMAX
        DO J = JMIN, JMAX
            DO I = IMIN, IMAX
                L = LCELL(I,J,K)
                IF (L == 0) THEN
                    CYCLE
                END IF
!
! Inner loop on nodes L.
!
                DO
                    DELX = XP - X(L)
                    DELY = YP - Y(L)
                    DELZ = ZP - Z(L)
                    DXSQ = DELX*DELX
                    DYSQ = DELY*DELY
                    DZSQ = DELZ*DELZ
                    DS = DXSQ + DYSQ + DZSQ
                    RS = RSQ(L)

                    IF (DS < RS) THEN
                        IF (DS == 0.0E+00_R8) THEN
                            Q = F(L)
                            QX = A(7,L)
                            QY = A(8,L)
                            QZ = A(9,L)
                            IER = 0
                            RETURN
                        END IF
                        RDS = RS*DS
                        RD = SQRT(RDS)
                        W = (RS+DS-RD-RD)/RDS
                        T = 2.0E+00_R8 *(RD-RS)/(DS*RDS)
                        WX = DELX*T
                        WY = DELY*T
                        WZ = DELZ*T
                        QLX = 2.0E+00_R8 *A(1,L)*DELX + A(2,L)*DELY + A(4,L)*DELZ
                        QLY = A(2,L)*DELX + 2.0E+00_R8 * A(3,L)*DELY + A(5,L)*DELZ
                        QLZ = A(4,L)*DELX + A(5,L)*DELY + 2.0E+00_R8 * A(6,L)*DELZ
                        QL = (QLX*DELX + QLY*DELY + QLZ*DELZ) / 2.0E+00_R8 + &
                        A(7,L)*DELX + A(8,L)*DELY + A(9,L)*DELZ + F(L)
                        QLX = QLX + A(7,L)
                        QLY = QLY + A(8,L)
                        QLZ = QLZ + A(9,L)
                        SW = SW + W
                        SWX = SWX + WX
                        SWY = SWY + WY
                        SWZ = SWZ + WZ
                        SWQ = SWQ + W*QL
                        SWQX = SWQX + WX*QL + W*QLX
                        SWQY = SWQY + WY*QL + W*QLY
                        SWQZ = SWQZ + WZ*QL + W*QLZ
                    END IF
                    LP = L
                    L = LNEXT(LP)
                    IF ( L == LP ) THEN
                        EXIT
                    END IF
                END DO
            END DO
        END DO
    END DO
!
! SW = 0 iff node P is not within the radius R(L) for any node L.
!
    IF (SW /= 0.0E+00_R8) THEN
        Q = SWQ/SW
        SWS = SW*SW
        QX = (SWQX*SW - SWQ*SWX)/SWS
        QY = (SWQY*SW - SWQ*SWY)/SWS
        QZ = (SWQZ*SW - SWQ*SWZ)/SWS
        IER = 0
!
! No cells contain a point within RMAX of node P, or
! SW = 0 and thus DS >= RSQ(L) for all L.
!
    ELSE
        Q = 0.0E+00_R8
        QX = 0.0E+00_R8
        QY = 0.0E+00_R8
        QZ = 0.0E+00_R8
        IER = 2
    END IF
    RETURN
  END SUBROUTINE QS3GRD
!  
!=========================================================================
!  
  FUNCTION QS3VAL (PX, PY, PZ, N, X, Y, Z, F, NR, LCELL, LNEXT, &
    XYZMIN, XYZDEL, RMAX, RSQ, A)
!
! QS3VAL evaluates the interpolant function Q(X,Y,Z) 
! created by QSHEP3.
!
! Discussion:
!
!  This function returns the value Q(PX,PY,PZ) where Q is
!  the weighted sum of quadratic nodal functions defined in
!  subroutine QSHEP3. QS3GRD may be called to compute a
!  gradient of Q along with the value, or to test for errors.
!
!  This function should not be called if a nonzero error flag was
!  returned by QSHEP3.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Updated to FORTRAN 95 by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 December 2005  
!
! Parameters:
!
!  INPUT
!   PX, PY, PZ, real, the point P at which Q is to be evaluated.
!
!   N, integer, the number of nodes and data values defining Q.
!   N >= 10.
!
!   X(N), Y(N), Z(N), F(N), real, the node coordinates
!   and data values interpolated by Q.
!
!   NR, integer, the number of rows, columns and planes in the cell
!   grid. Refer to STORE3. NR >= 1.
!
!   LCELL(NR,NR,NR), integer, the nodal indices associated with cells.  
!   Refer to STORE3.
!
!   LNEXT(N), integer, the next-node indices. Refer to STORE3.
!
!   XYZMIN(3), XYZDEL(3), real, the minimum nodal coordinates and 
!   cell dimensions, respectively. XYZDEL elements must be positive.  
!   Refer to STORE3.
!
!   RMAX, real, the square root of the largest element in RSQ,
!   the maximum radius.
!
!   RSQ(N), real, the squared radii which enter into the weights 
!   defining Q.
!
!   A(9,N), real, the coefficients for the nodal functions defining Q.
!
!  OUTPUT
!   QS3VAL, real, the function value Q(PX,PY,PZ) unless N, NR,
!   XYZDEL, or RMAX are invalid, in which case the value 0 is returned.
!
    USE REAL_PRECISION
    IMPLICIT NONE
    
    INTEGER :: I, IMAX, IMIN, J, JMAX, JMIN, K, KMAX, KMIN, &
          L, LP, N, NR
    INTEGER, DIMENSION (N) :: LNEXT
    INTEGER, DIMENSION (NR,NR,NR) :: LCELL

   REAL(KIND=R8):: DELX, DELY, DELZ, DS, DX, DXSQ, DY, DYSQ, DZ, &
         DZSQ, PX, PY, PZ, QS3VAL, RD, RDS, RMAX, RS, SW, SWQ, W,  &
         XMIN, XP, YMIN, YP, ZMIN, ZP
    REAL(KIND=R8), DIMENSION (N) :: F, RSQ, X
    REAL(KIND=R8), DIMENSION (3) :: XYZDEL, XYZMIN 
    REAL(KIND=R8), DIMENSION (N) :: Y, Z
    REAL(KIND=R8), DIMENSION (9,N) :: A
  
    XP = PX
    YP = PY
    ZP = PZ
    XMIN = XYZMIN(1)
    YMIN = XYZMIN(2)
    ZMIN = XYZMIN(3)
    DX = XYZDEL(1)
    DY = XYZDEL(2)
    DZ = XYZDEL(3)

    IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0 &
           .OR.  DY <= 0.0  .OR.  DZ <= 0.0  .OR. &
           RMAX < 0.0 ) THEN
        QS3VAL = 0.0E+00
        RETURN
    END IF
!
! Set IMIN, IMAX, JMIN, JMAX, KMIN, and KMAX to cell indices
! defining the range of the search for nodes whose radii
! include node P. The cells which must be searched are those
! intersected by (or contained in) a sphere of radius RMAX
! centered at node P.
!
    IMIN = INT((XP-XMIN-RMAX)/DX) + 1
    IMIN = MAX (IMIN, 1)
    IMAX = INT((XP-XMIN+RMAX)/DX) + 1
    IMAX = MIN (IMAX, NR)
    JMIN = INT((YP-YMIN-RMAX)/DY) + 1
    JMIN = MAX (JMIN, 1)
    JMAX = INT((YP-YMIN+RMAX)/DY) + 1
    JMAX = MIN (JMAX, NR)
    KMIN = INT((ZP-ZMIN-RMAX)/DZ) + 1
    KMIN = MAX (KMIN, 1)
    KMAX = INT((ZP-ZMIN+RMAX)/DZ) + 1
    KMAX = MIN (KMAX, NR)
!
! Test for no cells within the sphere of radius RMAX.
!
    IF (IMIN > IMAX .OR. JMIN > JMAX .OR. KMIN > KMAX) THEN
        QS3VAL = 0.0E+00_R8
        RETURN
    END IF
!
! Accumulate weight values in SW and weighted nodal function
! values in SWQ.  The weights are W(L) = ((R-D)+/(R*D))**2
! for R**2 = RSQ(L) and D = distance between node P and node L.
!
    SW = 0.0E+00_R8
    SWQ = 0.0E+00_R8
!
! Outer loop on cells (I,J,K).
!
    DO K = KMIN, KMAX
        DO J = JMIN, JMAX
            DO I = IMIN, IMAX
                L = LCELL(I,J,K)
                IF (L == 0) THEN
                    CYCLE
                END IF
!
! Inner loop on nodes L.
!
                DO
                    DELX = XP - X(L)
                    DELY = YP - Y(L)
                    DELZ = ZP - Z(L)
                    DXSQ = DELX*DELX
                    DYSQ = DELY*DELY
                    DZSQ = DELZ*DELZ
                    DS = DXSQ + DYSQ + DZSQ
                    RS = RSQ(L)
                    IF (DS < RS) THEN
                        IF (DS == 0.0E+00_R8) THEN
                            QS3VAL = F(L)
                            RETURN
                        END IF
                        RDS = RS*DS
                        RD = SQRT(RDS)
                        W = (RS+DS-RD-RD)/RDS
                        SW = SW + W
                        SWQ = SWQ + W *( A(1,L)*DXSQ + A(2,L)*DELX*DELY + &
                                   A(3,L)*DYSQ + A(4,L)*DELX*DELZ + &
                                   A(5,L)*DELY*DELZ + A(6,L)*DZSQ + &
                                   A(7,L)*DELX + A(8,L)*DELY + &
                                   A(9,L)*DELZ + F(L) )

                    END IF
                    LP = L
                    L = LNEXT(LP)
                    IF (L == LP) THEN
                        EXIT
                    END IF
                END DO
            END DO
        END DO
    END DO
!
! SW = 0 iff P is not within the radius R(L) for any node L.
!
    IF (SW == 0.0E+00_R8) THEN
        QS3VAL = 0.0E+00_R8
    ELSE
        QS3VAL = SWQ / SW
    END IF
    RETURN
  END FUNCTION QS3VAL
!  
!=========================================================================
!  
END MODULE QSHEP3D_MOD
