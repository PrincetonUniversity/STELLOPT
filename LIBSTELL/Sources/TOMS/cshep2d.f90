      MODULE CSHEP2D_MOD
      PRIVATE
      PUBLIC CSHEP2, CS2VAL, CS2HES, CS2GRD
      CONTAINS
      SUBROUTINE CSHEP2 (N,X,Y,F,NC,NW,NR, LCELL,LNEXT,XMIN,&
                        YMIN,DX,DY,RMAX,RW,A,IER)
      USE REAL_PRECISION
      IMPLICIT NONE
      INTEGER:: N, NC, NW, NR, IER
      INTEGER, DIMENSION (N)::LNEXT
      INTEGER, DIMENSION(NR,NR)::LCELL
      REAL (KIND=R8):: XMIN, YMIN, DX, DY, RMAX
      REAL (KIND=R8), DIMENSION(N):: X,Y,F,RW
      REAL (KIND=R8), DIMENSION(9,N)::A
!
!***********************************************************
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/13/97
! Updated to Fortran 95 by:
!  William Thacker
!  Winthrop University
!  6/10/2008
!
!   This subroutine computes a set of parameters defining a
! C2 (twice continuously differentiable) bivariate function
! C(X,Y) which interpolates data values F at a set of N
! arbitrarily distributed points (X,Y) in the plane (nodes).
! The interpolant C may be evaluated at an arbitrary point
! by function CS2VAL, and its first partial derivatives are
! computed by the subroutine CS2GRD.
!
!   The interpolation scheme is a modified Cubic Shepard
! method:
!
! C = [W(1)*C(1)+W(2)*C(2)+..+W(N)*C(N)]/[W(1)+W(2)+..+W(N)]
!
! for bivariate functions W(k) and C(k).  The nodal func-
! tions are given by
!
!  C(k)(x,y) = A(1,k)*(x-X(k))**3 +
!              A(2,k)*(x-X(k))**2*(y-Y(k)) +
!              A(3,k)*(x-X(k))*(y-Y(k))**2 +
!              A(4,k)*(y-Y(k))**3 + A(5,k)*(x-X(k))**2 +
!              A(6,k)*(x-X(k))*(y-Y(k)) + A(7,k)*(y-Y(k))**2
!              + A(8,k)*(x-X(k)) + A(9,k)*(y-Y(k)) + F(k) .
!
! Thus, C(k) is a cubic function which interpolates the data
! value at node k.  Its coefficients A(,k) are obtained by a
! weighted least squares fit to the closest NC data points
! with weights similar to W(k).  Note that the radius of
! influence for the least squares fit is fixed for each k,
! but varies with k.
!
! The weights are taken to be
!
!   W(k)(x,y) = ( (R(k)-D(k))+ / R(k)*D(k) )**3 ,
!
! where (R(k)-D(k))+ = 0 if R(k) < D(k), and D(k)(x,y) is
! the Euclidean distance between (x,y) and (X(k),Y(k)).  The
! radius of influence R(k) varies with k and is chosen so
! that NW nodes are within the radius.  Note that W(k) is
! not defined at node (X(k),Y(k)), but C(x,y) has limit F(k)
! as (x,y) approaches (X(k),Y(k)).
!
! On input:
!
!       N = Number of nodes and data values.  N >= 10.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       F = Array of length N containing the data values
!           in one-to-one correspondence with the nodes.
!
!       NC = Number of data points to be used in the least
!            squares fit for coefficients defining the nodal
!            functions C(k).  Values found to be optimal for
!            test data sets ranged from 11 to 25.  A recom-
!            mended value for general data sets is NC = 17.
!            For nodes lying on (or close to) a rectangular
!            grid, the recommended value is NC = 11.  In any
!            case, NC must be in the range 9 to Min(40,N-1).
!
!       NW = Number of nodes within (and defining) the radii
!            of influence R(k) which enter into the weights,
!            W(k).  For N sufficiently large, a recommended
!            value is NW = 30.  In general, NW should be
!            about 1.5*NC.  1 <= NW <= Min(40,N-1).
!
!       NR = Number of rows and columns in the cell grid de-
!            fined in the subroutine STORE2.  A rectangle con-
!            taining the nodes is partitioned into cells in
!            order to increase search efficiency.  NR =
!            Sqrt(N/3) is recommended.  NR >= 1.
!
! The above parameters are not altered by this routine.
!
!       LCELL = Array of length >= NR**2.
!
!       LNEXT = Array of length >= N.
!
!       RW = Array of length >= N.
!
!       A = Array of length >= 9N.
!
! On output:
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to the subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to the subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  Refer to the subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW (maximum of the radii, R(k)).
!
!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k).
!
!       A = 9 by N array containing the coefficients for the
!           cubic nodal function C(k) in column k.
!
!   Note that the output parameters described above are not
! defined unless IER = 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NC, NW, or NR is outside its
!                     valid range.
!             IER = 2 if duplicate nodes were encountered.
!             IER = 3 if all nodes are collinear.
!
!***********************************************************
!
      INTEGER, PARAMETER:: LMX=40
      INTEGER:: I, IERR, IROW, J, JP1, K, LMAX,&
          LNP, NEQ, NN, NNC, NNR, NNW, NP, NCWMAX
      INTEGER, DIMENSION(LMX):: NPTS
      REAL(KIND=R8):: C, DDX, DDY, DMIN, FK, RC, RS, RSMX, &
          RSOLD, RWS, S, SF, SFC, SFS, STF, SUM, T, XK, XMN, &
          YK, YMN
      REAL(KIND=R8), DIMENSION (10,10):: B
      REAL(KIND=R8), PARAMETER::RTOL=1.0E-05_R8
      REAL(KIND=R8)::DTOL
!
! Local variables:
!
! B =          Transpose of the augmented regression matrix
! C =          First component of the plane rotation used to
!                zero the lower triangle of B**T -- computed
!                by the subroutine GIVENS
! DDX,DDY =    Local variables for DX and DY
! DMIN =       Minimum of the magnitudes of the diagonal
!                elements of the regression matrix after
!                zeros are introduced below the diagonal
! DTOL =       Tolerance for detecting an ill-conditioned
!                system.  The system is accepted when
!                DMIN*RC >= DTOL.
! FK =         Data value at node K (F(K))
! I =          Index for A, B, and NPTS
! IERR =       Error flag for the call to the subroutine STORE2
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       Row index for B
! J =          Index for A and B
! JP1 =        J+1
! K =          Nodal function index and column index for A
! LMAX =       Maximum number of NPTS elements
! LMX =        Maximum value of LMAX
! LNP =        Current length of NPTS
! NEQ =        Number of equations in the least squares fit
! NN,NNC,NNR = Local copies of N, NC, and NR
! NNW =        Local copy of NW
! NP =         NPTS element
! NPTS =       Array containing the indexes of a sequence of
!                nodes to be used in the least squares fit
!                or to compute RW.  The nodes are ordered
!                by distance from K, and the last element
!                (usually indexed by LNP) is used only to
!                determine RC, or RW(K) if NW > NC.
! NCWMAX =     Max(NC,NW)
! RC =         Radius of influence which enters into the
!                weights for C(K) (see the subroutine SETUPC2)
! RS =         Squared distance between K and NPTS(LNP).
!                Used to compute RC and RW(K)
! RSMX =       Maximum squared RW element encountered
! RSOLD =      Squared distance between K and NPTS(LNP-1).
!                Used to compute a relative change in RS
!                between succeeding NPTS elements
! RTOL =       Tolerance for detecting a sufficiently large
!                relative change in RS.  If the change is
!                not greater than RTOL, the nodes are
!                treated as being the same distance from K
! RWS =        Current squared value of RW(K)
! S =          Second component of the plane rotation deter-
!                mined by the subroutine GIVENS
! SF =        Scale factor for the linear terms (columns 8
!               and 9) in the least squares fit. It is the inverse
!               of the root-mean-square distance between K
!               and the nodes (other than K) in the least
!               squares fit
! SFS =       Scale factor for the quadratic terms (columns
!               5, 6, and 7) in the least squares fit 
!               (SF*SF)
! SFC =       Scale factor for the cubic terms (first 4
!               columns) in the least squares fit (SF**3)
! STF =        Marquardt stabilization factor used to damp
!                out the first 4 solution components (third
!                partials of the cubic) when the system is
!                ill-conditioned.  As STF increases, the
!                fitting function approaches a quadratic
!                polynomial.
! SUM =        Sum of squared Euclidean distances between
!                node K and the nodes used in the least
!                squares fit (unless additional nodes are
!                added for stability)
! T =          Temporary variable for accumulating a scalar
!                product in the back solve
! XK,YK =      Coordinates of node K (X(K), Y(K))
! XMN,YMN =    Local variables for XMIN and YMIN
!
      DTOL=SQRT(EPSILON(1.0_R8))
      NN = N
      NNC = NC
      NNW = NW
      NNR = NR
      NCWMAX = MAX(NNC,NNW)
      LMAX = MIN(LMX,NN-1)
      IF (NNC < 9  .OR.  NNW < 1  .OR.  NCWMAX >  &
          LMAX  .OR.  NNR < 1) THEN
          IER=1
          RETURN
      END IF
!
! Create the cell data structure, and initialize RSMX.
!
      CALL STORE2 (NN,X,Y,NNR, LCELL,LNEXT,XMN,YMN,DDX,DDY, &
          IERR)
      IF (IERR /= 0) THEN
          XMIN=XMN
          YMIN=YMN
          DX=DDX
          DY=DDY
          IER=3
          RETURN
      END IF

      RSMX = 0.0E+00_R8
!
! Outer loop on node K:
!
      DO K = 1,NN
          XK = X(K)
          YK = Y(K)
          FK = F(K)
!
! Mark node K to exclude it from the search for nearest
!   neighbors.
!
          LNEXT(K) = -LNEXT(K)
!
! Initialize for loop on NPTS.
!
          RS = 0.0E+00_R8
          SUM = 0.0E+00_R8
          RWS = 0.0E+00_R8
          RC = 0.0E+00_R8
          LNP = 0
!
! Compute NPTS, LNP, RWS, NEQ, RC, and SFS.
!
          MAIN: DO
              SUM = SUM + RS
              IF (LNP == LMAX) THEN
! 
! All LMAX nodes are included in NPTS.  RWS and/or RC**2 is
!   (arbitrarily) taken to be 10 percent larger than the
!   distance RS to the last node included.
!

                  IF (RWS == 0.0E+00_R8) THEN
                      RWS = 1.1E+00_R8*RS
                  END IF
                  IF (RC == 0.0E+00_R8) THEN
                      NEQ = LMAX
                      RC = SQRT(1.1E+00_R8*RS)
                      SFS = REAL(NEQ,KIND=R8)/SUM
                  END IF
             
!
! Store RW(K), update RSMX if necessary, and compute SF
!   and SFC.
!
                  RW(K) = SQRT(RWS)
                  IF (RWS > RSMX) THEN
                      RSMX = RWS
                  END IF
                  SF = SQRT(SFS)
                  SFC = SF*SFS
                  EXIT
              END IF

              LNP = LNP + 1
              RSOLD = RS
              CALL GETNP2 (XK,YK,X,Y,NN,NNR,LCELL,LNEXT,XMN,YMN, &
                  DDX,DDY, NP,RS)
              IF (RS == 0.0E+00_R8) THEN
                  IER=2
                  RETURN
              END IF

              NPTS(LNP) = NP
              IF ( (RS-RSOLD)/RS < RTOL ) THEN
                  CYCLE MAIN
              END IF
              IF (RWS == 0.0E+00_R8  .AND.  LNP > NNW) THEN
                  RWS = RS
              END IF
              IF (RC == 0.0E+00_R8  .AND.  LNP > NNC) THEN
!
!   RC = 0 (not yet computed) and LNP > NC.  RC = Sqrt(RS)
!     is sufficiently large to (strictly) include NC nodes.
!     The least squares fit will include NEQ = LNP - 1
!     equations for 9 <= NC <= NEQ < LMAX <= N-1.
!
                  NEQ = LNP - 1
                  RC = SQRT(RS)
                  SFS = REAL(NEQ,KIND=R8)/SUM
              END IF
!
!   Bottom of loop -- test for termination.
!
              IF (LNP > NCWMAX) THEN
! 
! Store RW(K), update RSMX if necessary, and compute SF
!   and SFC.
!

                  RW(K) = SQRT(RWS)
                  IF (RWS > RSMX) THEN
                      RSMX = RWS
                  END IF
                  SF = SQRT(SFS)
                  SFC = SF*SFS
                  EXIT 
              END IF

          END DO MAIN
!
! A Q-R decomposition is used to solve the least squares
!   system.  The transpose of the augmented regression
!   matrix is stored in B with columns (rows of B) defined
!   as follows:  1-4 are the cubic terms, 5-7 are the quad-
!   ratic terms, 8 and 9 are the linear terms, and the last
!   column is the right hand side.
!
! Set up the equations and zero out the lower triangle with
!   Givens rotations.
!
          I = 0
          PTS: DO
              I = I + 1
              NP = NPTS(I)
              IROW = MIN(I,10)
              CALL SETUPC2 (XK,YK,FK,X(NP),Y(NP),F(NP),SF,SFS, &
                  SFC,RC, B(1,IROW))
              IF (I == 1) THEN
                  CYCLE PTS
              END IF

              DO  J = 1,IROW-1
                  JP1 = J + 1
                  CALL GIVENS (B(J,J),B(J,IROW),C,S)
                  CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,IROW))
              END DO
              IF (I < NEQ) THEN
                  CYCLE PTS
              END IF
!
! Test the system for ill-conditioning.
!
              DMIN = MIN( ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)), &
                  ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)), &
                  ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)) )
              IF (DMIN*RC >= DTOL) THEN
                  EXIT
              END IF
              IF (NEQ == LMAX) THEN
                  EXIT
              END IF
!
! Increase RC and add another equation to the system to
!   improve the conditioning.  The number of NPTS elements
!   is also increased if necessary.
!
              TOL: DO
                  RSOLD = RS
                  NEQ = NEQ + 1
                  IF (NEQ == LMAX) THEN
                      RC = SQRT(1.1E+00_R8*RS)
                      CYCLE PTS
                  END IF
                  IF (NEQ < LNP) THEN
!
!   NEQ < LNP.
!
                      NP = NPTS(NEQ+1)
                      RS = (X(NP)-XK)**2 + (Y(NP)-YK)**2
                      IF ( (RS-RSOLD)/RS < RTOL ) THEN
                         CYCLE TOL
                      END IF

                      RC = SQRT(RS)
                      CYCLE PTS
                  END IF
!
!   NEQ = LNP.  Add an element to NPTS.
!
                  LNP = LNP + 1
                  CALL GETNP2 (XK,YK,X,Y,NN,NNR,LCELL,LNEXT,XMN,YMN, &
                      DDX,DDY, NP,RS)
                  IF (NP == 0) THEN
                      IER=2
                      RETURN
                  END IF

                  NPTS(LNP) = NP
                  IF ( (RS-RSOLD)/RS < RTOL ) THEN
                      CYCLE TOL
                  END IF

                  RC = SQRT(RS)
              END DO TOL
          END DO PTS
!
! Stabilize the system by damping third partials -- add
!   multiples of the first four unit vectors to the first
!   four equations.
!
          IF (NEQ == LMAX) THEN
              STF = 1.0E+00_R8/RC
              DO  I = 1,4
                  B(I,10) = STF
                  DO  J = I+1,10
                      B(J,10) = 0.0E+00_R8
                  END DO
                  DO  J = I,9
                      JP1 = J + 1
                      CALL GIVENS (B(J,J),B(J,10),C,S)
                      CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,10))
                  END DO
              END DO
!
! Test the damped system for ill-conditioning.
!
              DMIN = MIN( ABS(B(5,5)),ABS(B(6,6)),ABS(B(7,7)), &
                   ABS(B(8,8)),ABS(B(9,9)) )
              IF (DMIN*RC < DTOL) THEN
!
! No unique solution due to collinear nodes.
!
                  XMIN = XMN
                  YMIN = YMN
                  DX = DDX
                  DY = DDY
                  IER = 3
                  RETURN
              END IF

          END IF
!
! Solve the 9 by 9 triangular system for the coefficients.
!
          DO  I = 9,1,-1
              T = 0.0E+00_R8
              IF (I /= 9) THEN
                  DO  J = I+1,9
                      T = T + B(J,I)*A(J,K)
                  END DO
              END IF
              A(I,K) = (B(10,I)-T)/B(I,I)
          END DO
!
! Scale the coefficients to adjust for the column scaling.
!
          DO  I = 1,4
              A(I,K) = A(I,K)*SFC
          END DO
          A(5,K) = A(5,K)*SFS
          A(6,K) = A(6,K)*SFS
          A(7,K) = A(7,K)*SFS
          A(8,K) = A(8,K)*SF
          A(9,K) = A(9,K)*SF
!
! Unmark K and the elements of NPTS.
!
          LNEXT(K) = -LNEXT(K)
          DO  I = 1,LNP
              NP = NPTS(I)
              LNEXT(NP) = -LNEXT(NP)
          END DO
      END DO
!
! No errors encountered.
!
      XMIN = XMN
      YMIN = YMN
      DX = DDX
      DY = DDY
      RMAX = SQRT(RSMX)
      IER = 0
      RETURN
      END SUBROUTINE CSHEP2





      FUNCTION CS2VAL (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN, &
          YMIN,DX,DY,RMAX,RW,A)
      USE REAL_PRECISION
      IMPLICIT NONE
      INTEGER:: N, NR
      INTEGER, DIMENSION(NR, NR)::LCELL
      INTEGER, DIMENSION(N)::LNEXT
      REAL(KIND=R8):: PX, PY, XMIN, YMIN, DX, DY, RMAX,CS2VAL
      REAL(KIND=R8), DIMENSION(N)::X, Y, F, RW
      REAL(KIND=R8), DIMENSION(9,N)::A
!
!***********************************************************
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
! Updated to Fortran 95 by:
!  William Thacker
!  Winthrop University
!  6/10/2008
!
!   This function returns the value C(PX,PY), where C is the
! weighted sum of cubic nodal functions defined in the subrou-
! tine CSHEP2.  CS2GRD may be called to compute a gradient
! of C along with the value, and/or to test for errors.
! CS2HES may be called to compute a value, first partial
! derivatives, and second partial derivatives at a point.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P at
!               which C is to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N >= 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to the subroutine STORE2.  NR >= 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to the subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to the subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to the subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW (maximum of the radii R(k)).
!
!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for the
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this function.  The
! parameters other than PX and PY should be input unaltered
! from their values on output from CSHEP2.  This function
! should not be called if a nonzero error flag was returned
! by CSHEP2.
!
! On output:
!
!       CS2VAL = Function value C(PX,PY) unless N, NR, DX,
!                DY, or RMAX is invalid, in which case no
!                value is returned.
!
!
!***********************************************************
!
      INTEGER:: I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      REAL(KIND=R8) D, DELX, DELY, R, SW, SWC, W, XP, YP
!
! Local variables:
!
! D =         Distance between node P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at node P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at node P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of the weights W(K)
! SWC =       Sum of the weighted nodal function values at node P
! W =         Weight W(K) value at node P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! XP,YP =     Local copies of PX and PY (coordinates of P)
!
      XP = PX
      YP = PY
      IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0E+00_R8  .OR. &
          DY <= 0.0E+00_R8  .OR.  RMAX < 0.0E+00_R8) RETURN
!
! Set IMIN, IMAX, JMIN, and JMAX to the cell indexes defining
!   the range of the search for nodes whose radii include
!   node P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at node P.
!
      IMIN = INT((XP-XMIN-RMAX)/DX) + 1
      IMAX = INT((XP-XMIN+RMAX)/DX) + 1
      IMIN=MAX(IMIN,1)
      IMAX=MIN(IMAX,NR)
      JMIN = INT((YP-YMIN-RMAX)/DY) + 1
      JMAX = INT((YP-YMIN+RMAX)/DY) + 1
      JMIN=MAX(JMIN,1)
      JMAX=MIN(JMAX,NR)
!
! The following is a test for no cells within the circle
!   of radius RMAX.
!
      IF (IMIN > IMAX  .OR.  JMIN > JMAX) THEN
!
! All weights are 0 at P.
!
          CS2VAL=0.0E+00_R8
          RETURN
      END IF
!
! Accumulate weight values in SW and weighted nodal function
!   values in SWC.  The weights are W(K) = ((R-D)+/(R*D))**3
!   for R = RW(K) and D = distance between P and node K.
!
      SW = 0.0E+00_R8
      SWC = 0.0E+00_R8
!
! Outer loop on cells (I,J).
!
      DO  J = JMIN,JMAX
          DO  I = IMIN,IMAX
              K = LCELL(I,J)
              IF (K /= 0) THEN
!
! Inner loop on nodes K.
!
                  EVAL: DO
                      DELX = XP - X(K)
                      DELY = YP - Y(K)
                      D = SQRT(DELX*DELX + DELY*DELY)
                      R = RW(K)
                      IF (D < R) THEN
                          IF (D == 0.0E+00_R8) THEN
!
! (PX,PY) = (X(K),Y(K)).
!
                              CS2VAL=F(K)
                              RETURN
                          END IF
                          W = (1.0E+00_R8/D - 1.0E+00_R8/R)**3
                          SW = SW + W
                          SWC = SWC + W*( ( (A(1,K)*DELX+A(2,K)*DELY+ &
                          A(5,K))*DELX + (A(3,K)*DELY+ &
                          A(6,K))*DELY + A(8,K) )*DELX + &
                          ( (A(4,K)*DELY+A(7,K))*DELY + &
                          A(9,K) )*DELY + F(K) )
                      END IF
!
! Bottom of loop on nodes in cell (I,J).
!
                      KP = K
                      K = LNEXT(KP)
                      IF (K == KP) THEN
                          EXIT
                      END IF
                  END DO EVAL
              END IF 
          END DO
      END DO
!
! SW = 0 iff P is not within the radius R(K) for any node K.
!
      IF (SW == 0.0E+00_R8) THEN
          CS2VAL=0.0E+00_R8
      ELSE
          CS2VAL = SWC/SW
      END IF
      RETURN
      END FUNCTION CS2VAL



      SUBROUTINE CS2GRD (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN, &
                        YMIN,DX,DY,RMAX,RW,A, C,CX,CY,IER)
      USE REAL_PRECISION
      IMPLICIT NONE
      INTEGER:: N, NR, IER
      INTEGER, DIMENSION(NR,NR)::LCELL
      INTEGER, DIMENSION(N)::LNEXT
      REAL(KIND=R8) PX, PY, XMIN, YMIN, DX, DY, RMAX, C, CX, CY
      REAL(KIND=R8), DIMENSION(N)::X,Y,F,RW
      REAL(KIND=R8), DIMENSION(9,N):: A
!
!***********************************************************
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
! Updated to Fortran 95 by:
!  William Thacker
!  Winthrop University
!  6/10/2008
!
!   This subroutine computes the value and gradient at P =
! (PX,PY) of the interpolatory function C defined in the sub-
! routine CSHEP2.  C is a weighted sum of cubic nodal
! functions.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the node P at
!               which C and its partial derivatives are
!               to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N >= 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to the subroutine STORE2.  NR >= 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to the subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to the subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to the subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW (maximum of the radii R(k)).
!
!       RW = Array of length N containing the radii R(k)
!            which enter into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for the
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2.
!
! On output:
!
!       C = Value of C at (PX,PY) unless IER = 1, in
!           which case no values are returned.
!
!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER = 1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C=CX=CY=0).
!
!***********************************************************
!
      INTEGER:: I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      REAL(KIND=R8):: CK, CKX, CKY, D, DELX, DELY, R, SW, SWC,&
          SWCX, SWCY, SWS, SWX, SWY, T, W, WX, WY, XP, YP
!
! Local variables:
!
! CK =        Value of the cubic nodal function C(K) at node P
! CKX,CKY =   Partial derivatives of C(K) with respect to X
!               and Y, respectively
! D =         Distance between node P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at node P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at node P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at node P
! SWCX,SWCY = Partial derivatives of SWC with respect to X
!               and Y, respectively
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X
!               and Y, respectively
! T =         Temporary variable
! W =         Weight W(K) value at node P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X
!               and Y, respectively
! XP,YP =     Local copies of PX and PY (coordinates of node P)
!
      XP = PX
      YP = PY
      IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0E+00_R8  .OR. &
          DY <= 0.0E+00_R8  .OR.  RMAX < 0.0E+00_R8) THEN
!
! Invalid input parameter.
!
          IER=1
          RETURN
      END IF
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   node P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at node P.
!
      IMIN = INT((XP-XMIN-RMAX)/DX) + 1
      IMAX = INT((XP-XMIN+RMAX)/DX) + 1
      IMIN=MAX(IMIN,1)
      IMAX=MIN(IMAX,NR)
      JMIN = INT((YP-YMIN-RMAX)/DY) + 1
      JMAX = INT((YP-YMIN+RMAX)/DY) + 1
      JMIN=MAX(JMIN,1)
      JMAX=MIN(JMAX,NR)
!
! The following is a test for no cells within the circle
!   of radius RMAX.
!
      IF (IMIN > IMAX  .OR.  JMIN > JMAX) THEN
!
! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.
          C=0.0E+00_R8
          CX=0.0E+00_R8
          CY=0.0E+00_R8
          IER=2
          RETURN
      END IF
!
! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
!   from K = 1 to N, C(K) is the cubic nodal function value,
!   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
!   ance D(K).  Thus
!
!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2
!
!   where SWCX and SWX are partial derivatives with respect
!   to X of SWC and SW, respectively.  SWCY and SWY are de-
!   fined similarly.
!
      SW = 0.0E+00_R8
      SWX = 0.0E+00_R8
      SWY = 0.0E+00_R8
      SWC = 0.0E+00_R8
      SWCX = 0.0E+00_R8
      SWCY = 0.0E+00_R8
!
! Outer loop on cells (I,J).
!
      DO  J = JMIN,JMAX
          DO  I = IMIN,IMAX
              K = LCELL(I,J)
              IF (K /= 0) THEN
!
! Inner loop on nodes K.
!
                  EVAL: DO
                      DELX = XP - X(K)
                      DELY = YP - Y(K)
                      D = SQRT(DELX*DELX + DELY*DELY)
                      R = RW(K)
                      IF (D < R) THEN
                          IF (D == 0.0E+00_R8) THEN
!
! (PX,PY) = (X(K),Y(K)).
!
                              C=F(K)
                              CX=A(8,K)
                              CY=A(9,K)
                              IER=0
                              RETURN
                          END IF
                          T = (1.0E+00_R8/D - 1.0E+00_R8/R)
                          W = T**3
                          T = -3.0E+00_R8*T*T/(D**3)
                          WX = DELX*T
                          WY = DELY*T
                          T = A(2,K)*DELX + A(3,K)*DELY + A(6,K)
                          CKY = ( 3.0E+00_R8*A(4,K)*DELY + A(3,K)*DELX + &
                              2.0E+00_R8*A(7,K) )*DELY + T*DELX + A(9,K)
                          T = T*DELY + A(8,K)
                          CKX = ( 3.0E+00_R8*A(1,K)*DELX + A(2,K)*DELY + &
                              2.0E+00_R8*A(5,K) )*DELX + T
                          CK = ( (A(1,K)*DELX+A(5,K))*DELX + T )*DELX + &
                              ( (A(4,K)*DELY+A(7,K))*DELY + A(9,K) )*DELY + &
                              F(K)
                          SW = SW + W
                          SWX = SWX + WX
                          SWY = SWY + WY
                          SWC = SWC + W*CK
                          SWCX = SWCX + WX*CK + W*CKX
                          SWCY = SWCY + WY*CK + W*CKY
                      END IF
!
! Bottom of loop on nodes in cell (I,J).
!
                      KP = K
                      K = LNEXT(KP)
                      IF (K == KP) THEN
                          EXIT
                      END IF
                  END DO EVAL
              END IF
          END DO
      END DO
!
! SW = 0 iff P is not within the radius R(K) for any node K.
!
      IF (SW == 0.0E+00_R8) THEN
          C=0.0E+00_R8
          CX=0.0E+00_R8
          CY=0.0E+00_R8
          IER=2
          RETURN
      END IF
      C = SWC/SW
      SWS = SW*SW
      CX = (SWCX*SW - SWC*SWX)/SWS
      CY = (SWCY*SW - SWC*SWY)/SWS
      IER = 0
      RETURN
      END SUBROUTINE CS2GRD




      SUBROUTINE CS2HES (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN, &
          YMIN,DX,DY,RMAX,RW,A,C,CX,CY,CXX,CXY,CYY,IER)
      USE REAL_PRECISION
      IMPLICIT NONE
      INTEGER:: N, NR, IER
      INTEGER, DIMENSION(NR,NR):: LCELL
      INTEGER, DIMENSION(N):: LNEXT
      REAL(KIND=R8):: PX, PY, XMIN, YMIN, DX, DY, RMAX, C, CX, &
          CY, CXX, CXY, CYY
      REAL(KIND=R8), DIMENSION(N):: X, Y, F, RW
      REAL(KIND=R8), DIMENSION(9,N):: A
!
!***********************************************************
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
! Updated to Fortran 95 by:
!  William Thacker
!  Winthrop University
!  6/10/2008
!
!   This subroutine computes the value, gradient, and
! Hessian at P = (PX,PY) of the interpolatory function C
! defined in the subroutine CSHEP2.  C is a weighted sum of
! cubic nodal functions.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the node P at
!               which C and its partial derivatives are
!               to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N >= 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to the subroutine STORE2.  NR >= 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to the subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to the subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW (maximum of the radii R(k)).
!
!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for the
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2.
!
! On output:
!
!       C = Value of C at (PX,PY) unless IER = 1, in
!           which case no values are returned.
!
!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER = 1.
!
!       CXX,CXY,CYY = Second partial derivatives of C at
!                     (PX,PY) unless IER = 1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C = 0).
!
!***********************************************************
!
      INTEGER:: I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      REAL(KIND=R8) CK, CKX, CKXX, CKXY, CKY, CKYY, D, DELX,&
          DELY, DXSQ, DYSQ, R, SW, SWC, SWCX, SWCXX, SWCXY, &
          SWCY, SWCYY, SWS, SWX, SWXX, SWXY, SWY, SWYY, T1, &
          T2, T3, T4, W, WX, WXX, WXY, WY, WYY, XP, YP
!
! Local variables:
!
! CK =        Value of the cubic nodal function C(K) at node P
! CKX,CKY =   Partial derivatives of C(K) with respect to X
!               and Y, respectively
! CKXX,CKXY,CKYY = Second partial derivatives of CK
! D =         Distance between node P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! DXSQ,DYSQ = DELX**2, DELY**2
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at node P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at node P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at node P
! SWCX,SWCY = Partial derivatives of SWC with respect to X
!               and Y, respectively
! SWCXX,SWCXY,SWCYY = Second partial derivatives of SWC
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X
!               and Y, respectively
! SWXX,SWXY,SWYY = Second partial derivatives of SW
! T1,T2,T3,T4 = Temporary variables
! W =         Weight W(K) value at node P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X
!               and Y, respectively
! WXX,WXY,WYY = Second partial derivatives of W
! XP,YP =     Local copies of PX and PY -- coordinates of node P
!
      XP = PX
      YP = PY
      IF (N < 10  .OR.  NR < 1  .OR.  DX <= 0.0E+00_R8  .OR. &
          DY <= 0.0E+00_R8  .OR.  RMAX < 0.0E+00_R8) THEN
!
! Invalid input parameter.
!
          IER = 1
          RETURN
      END IF
!
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   node P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at node P.
!
      IMIN = INT((XP-XMIN-RMAX)/DX) + 1
      IMAX = INT((XP-XMIN+RMAX)/DX) + 1
      IMIN=MAX(IMIN,1)
      IMAX=MIN(IMAX,NR)
      JMIN = INT((YP-YMIN-RMAX)/DY) + 1
      JMAX = INT((YP-YMIN+RMAX)/DY) + 1
      JMIN=MAX(JMIN,1)
      JMAX=MIN(JMAX,NR)
!
! The following is a test for no cells within the circle
!   of radius RMAX.
!
      IF (IMIN > IMAX  .OR.  JMIN > JMAX) THEN
! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.
!
          C = 0.0E+00_R8
          CX = 0.0E+00_R8
          CY = 0.0E+00_R8
          CXX = 0.0E+00_R8
          CXY = 0.0E+00_R8
          CYY = 0.0E+00_R8
          IER = 2
          RETURN
      END IF
!
! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
!   from K = 1 to N, C(K) is the cubic nodal function value,
!   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
!   ance D(K).  Thus
!
!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2
!
!   where SWCX and SWX are partial derivatives with respect
!   to X of SWC and SW, respectively.  SWCY and SWY are de-
!   fined similarly.  The second partials are
!
!        CXX = ( SW*(SWCXX -    2*SWX*CX) - SWC*SWXX )/SW**2
!        CXY = ( SW*(SWCXY-SWX*CY-SWY*CX) - SWC*SWXY )/SW**2
!        CYY = ( SW*(SWCYY -    2*SWY*CY) - SWC*SWYY )/SW**2
!
!   where SWCXX and SWXX are second partials with respect
!   to X, SWCXY and SWXY are mixed partials, and SWCYY and
!   SWYY are second partials with respect to y.
!
      SW = 0.0E+00_R8
      SWX = 0.0E+00_R8
      SWY = 0.0E+00_R8
      SWXX = 0.0E+00_R8
      SWXY = 0.0E+00_R8
      SWYY = 0.0E+00_R8
      SWC = 0.0E+00_R8
      SWCX = 0.0E+00_R8
      SWCY = 0.0E+00_R8
      SWCXX = 0.0E+00_R8
      SWCXY = 0.0E+00_R8
      SWCYY = 0.0E+00_R8
!
! Outer loop on cells (I,J).
!
      DO  J = JMIN,JMAX
          DO  I = IMIN,IMAX
              K = LCELL(I,J)
              IF (K /= 0) THEN
!
! Inner loop on nodes K.
!
                  DO
                      DELX = XP - X(K)
                      DELY = YP - Y(K)
                      DXSQ = DELX*DELX
                      DYSQ = DELY*DELY
                      D = SQRT(DXSQ + DYSQ)
                      R = RW(K)
                      IF (D < R) THEN
                          IF (D == 0.0E+00_R8) THEN
!
! (PX,PY) = (X(K),Y(K)).
!
                              C = F(K)
                              CX = A(8,K)
                              CY = A(9,K)
                              CXX = 2.0E+00_R8*A(5,K)
                              CXY = A(6,K)
                              CYY = 2.0E+00_R8*A(7,K)
                              IER = 0
                              RETURN
                          END IF
                          T1 = (1.00E+00_R8/D - 1.00E+00_R8/R)
                          W = T1**3
                          T2 = -3.0E+00_R8*T1*T1/(D**3)
                          WX = DELX*T2
                          WY = DELY*T2
                          T1 = 3.0E+00_R8*T1*(2.0E+00_R8+3.0E+00_R8*D*T1)/(D**6)
                          WXX = T1*DXSQ + T2
                          WXY = T1*DELX*DELY
                          WYY = T1*DYSQ + T2
                          T1 = A(1,K)*DELX + A(2,K)*DELY + A(5,K)
                          T2 = T1 + T1 + A(1,K)*DELX
                          T3 = A(4,K)*DELY + A(3,K)*DELX + A(7,K)
                          T4 = T3 + T3 + A(4,K)*DELY
                          CK = (T1*DELX + A(6,K)*DELY + A(8,K))*DELX + &
                             (T3*DELY + A(9,K))*DELY + F(K)
                          CKX = T2*DELX + (A(3,K)*DELY+A(6,K))*DELY + A(8,K)
                          CKY = T4*DELY + (A(2,K)*DELX+A(6,K))*DELX + A(9,K)
                          CKXX = T2 + 3.0E+00_R8*A(1,K)*DELX
                          CKXY = 2.0E+00_R8*(A(2,K)*DELX + A(3,K)*DELY) + A(6,K)
                          CKYY = T4 + 3.0E+00_R8*A(4,K)*DELY
                          SW = SW + W
                          SWX = SWX + WX
                          SWY = SWY + WY
                          SWXX = SWXX + WXX
                          SWXY = SWXY + WXY
                          SWYY = SWYY + WYY
                          SWC = SWC + W*CK
                          SWCX = SWCX + WX*CK + W*CKX
                          SWCY = SWCY + WY*CK + W*CKY
                          SWCXX = SWCXX + W*CKXX + 2.0E+00_R8*WX*CKX + CK*WXX
                          SWCXY = SWCXY + W*CKXY + WX*CKY + WY*CKX + CK*WXY
                          SWCYY = SWCYY + W*CKYY + 2.0E+00_R8*WY*CKY + CK*WYY
!
! Bottom of loop on nodes in cell (I,J).
!
                      END IF
                      KP = K
                      K = LNEXT(KP)
                      IF (K == KP) THEN
                          EXIT
                      END IF
                  END DO 
              END IF
          END DO
      END DO
!
! SW = 0 iff P is not within the radius R(K) for any node K.
!
      IF (SW /= 0.0E+00_R8) THEN
          C = SWC/SW
          SWS = SW*SW
          CX = (SWCX*SW - SWC*SWX)/SWS
          CY = (SWCY*SW - SWC*SWY)/SWS
          CXX = (SW*(SWCXX-2.0E+00_R8*SWX*CX) - SWC*SWXX)/SWS
          CXY = (SW*(SWCXY-SWY*CX-SWX*CY) - SWC*SWXY)/SWS
          CYY = (SW*(SWCYY-2.0E+00_R8*SWY*CY) - SWC*SWYY)/SWS
          IER = 0
          RETURN
      END IF
! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.
!
      C = 0.0E+00_R8
      CX = 0.0E+00_R8
      CY = 0.0E+00_R8
      CXX = 0.0E+00_R8
      CXY = 0.0E+00_R8
      CYY = 0.0E+00_R8
      IER = 2
      RETURN
      END SUBROUTINE CS2HES



      SUBROUTINE SETUPC2 (XK,YK,ZK,XI,YI,ZI,S1,S2,S3,R, ROW)
      USE REAL_PRECISION
      IMPLICIT NONE
      REAL(KIND=R8):: XK, YK, ZK, XI, YI, ZI, S1, S2, S3, R
      REAL(KIND=R8), DIMENSION(10)::ROW
!
!***********************************************************
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
! Updated to Fortran 95 by:
!  William Thacker
!  Winthrop University
!  6/10/2008
!
!   This subroutine sets up the I-th row of an augmented re-
! gression matrix for a weighted least squares fit of a
! cubic function f(x,y) to a set of data values z, where
! f(XK,YK) = ZK.  The first four columns (cubic terms) are
! scaled by S3, the next three columns (quadratic terms)
! are scaled by S2, and the eighth and ninth columns (lin-
! ear terms) are scaled by S1.
!
! On input:
!
!       XK,YK = Coordinates of node K.
!
!       ZK = Data value at node K to be interpolated by f.
!
!       XI,YI,ZI = Coordinates and data value at node I.
!
!       S1,S2,S3 = Scale factors.
!
!       R = Radius of influence about node K defining the
!           weight.
!
! The above parameters are not altered by this routine.
!
!       ROW = Array of length 10.
!
! On output:
!
!       ROW = Array containing a row of the augmented re-
!             gression matrix.
!
!
!***********************************************************
!
      INTEGER:: I
      REAL(KIND=R8):: D, DX, DXSQ, DY, DYSQ, W, W1, W2, W3
!
! Local variables:
!
! D =    Distance between nodes K and I
! DX =   XI - XK
! DXSQ = DX*DX
! DY =   YI - YK
! DYSQ = DY*DY
! I =    DO-loop index
! W =    Weight associated with the row:  (R-D)/(R*D)
!          (0 if D = 0 or D > R)
! W1 =   S1*W
! W2 =   S2*W
! W3 =   W3*W
!
      DX = XI - XK
      DY = YI - YK
      DXSQ = DX*DX
      DYSQ = DY*DY
      D = SQRT(DXSQ + DYSQ)
      IF (D > 0.0E+00_R8  .AND.  D < R) THEN
          W = (R-D)/R/D
          W1 = S1*W
          W2 = S2*W
          W3 = S3*W
          ROW(1) = DXSQ*DX*W3
          ROW(2) = DXSQ*DY*W3
          ROW(3) = DX*DYSQ*W3
          ROW(4) = DYSQ*DY*W3
          ROW(5) = DXSQ*W2
          ROW(6) = DX*DY*W2
          ROW(7) = DYSQ*W2
          ROW(8) = DX*W1
          ROW(9) = DY*W1
          ROW(10) = (ZI - ZK)*W
          RETURN
      END IF
!
! Nodes K and I coincide or node I is outside of the radius
!   of influence.  Set ROW to the zero vector.
!
      DO  I = 1,10
          ROW(I) = 0.0E+00_R8
      END DO
      RETURN
      END SUBROUTINE SETUPC2

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
        RETURN
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
                IF (L /= 0) THEN
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

END MODULE CSHEP2D_MOD
