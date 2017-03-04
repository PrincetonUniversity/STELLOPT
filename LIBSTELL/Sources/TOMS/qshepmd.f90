MODULE QSHEPMD_MOD
PRIVATE
PUBLIC QSHEPM, QSMVAL
CONTAINS

 SUBROUTINE QSHEPM(M, N, X, F, NQ, NW, NR, RMAX, IER)
!*************************************************************
!                         QSHEPM
!                                              ROBERT RENKA
!                                      UNIV. OF NORTH TEXAS
!                                            (817) 565-2767
!                                                  10/25/87
!
!                                       Translated into C++
!                                           by Karen Minser
!                                        Univ. of Tennessee
!                                                    7/6/98
!
!                                       Translated into F95
!                                        by William Thacker
!                                       Winthrop University
!                                                 1/10/2009
! 
! This subroutine computes a set of parameters defining a
! smooth (once continuously differentiable) function Q(P)
! which interpolates data values F at scattered nodes X in
! Euclidiean M-space (P and X are vectors of length M).  The
! interpolant Q may be evaluated at an arbitrary point by
! function QSMVAL, and its gradient may be computed by the
! subroutine QSMGRD.
!   The interpolation scheme is a modified quadratic Shepard
! method:
! 
! Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))
! 
! For functions W(K) and Q(K) of M variables.  The nodal 
! functions are given by
! 
!   Q(K)(P) = SUM[A(K,L)*(P(I)-X(I,K))*(P(J)-X(J,K))] +
!             SUM[A(K,MT+I)*(P(I)-X(I,K)] + F(K)
! 
! where the sums are over all I and J such that 
! 1<=I<=J<=M, with L=(J-1)*J/2+I and MT=M*(M+1)/2.
! Thus, Q(K) is a quadratic function which interpolates the
! data value at node K.  Its coefficients A(K, ) are ordered
! by ordering the upper triangular array of MT quadratic
! terms by columns, followed by the M linear terms.  These
! coefficients are obtained by a weighted least squares fit 
! to the closest NQ data points with weights similar to W(K)
! (defined below).  The radius of influence for the least 
! squares fit is fixed for each K, but varies with K.
!   The weights are taken to be
! 
!       W(K)(P) = ( (R(K)-D(K))+ / (R(K)*D(K)) )**2
! 
! where (R(K)-D(K))+ = 0 if R(K)<=D(K), and D(K)(P) is
! the Euclidean distance between P and X( ,K).  THe radius
! of influence R(K) varies with K and is chosen so that NW 
! nodes are within the radius.  Note that W(K) is not 
! defined at node K, but Q(P) has limit F(K) as P approaches
! X( ,K)
!   A cell-based search method, described in subroutines
! STOREM and GETNPM, is used to improve efficiency in both
! the preprocessing qhase (QSHEPM) and the evaluation phase
! (QSMVAL and QSMGRD).  On a swquential processor, with a
! reasonably uniform distribution of nodes, operation counts
! are O(N) for preprocessing, and constant for each evaluation.
! However, computation time increases exponentially 
! with the dimension M.  Also, preprocessing time increases
! linearly with NQ, and evaluation time is proportional to
! NW.
! 
!   Input Parameters:
!        M, integer, dimension of the problem
! 
!        N, integer, the number of nodes and associated 
!           data values.
! 
!        X(M,N), real, the M-coordinates of the N nodes.
! 
!        F(N), real, the data values at the nodes
! 
!        NQ, integer, the number of data points to be used
!            in the least squares fit for coefficients
!            defining the nodal functions Q(K).
! 
!        NW, integer, the number of nodes within (and
!            defining) the radii of influence R(K), which
!            enter into the weights W(K).
! 
!        NR, integer, the number of divisions in each dimension
!            for the cell grid defined in the subroutine
!            STOREM.  A hyperbox containing the nodes is
!            partitioned into cells in order to increase
!            search efficiency.  NR=(N/3)**(1/M) is
!            recommended.  NR>=1.
! 
!   Output Parameters:
!       RMAX, real, square root of the largest element in
!           RSQ, the maximum radius of influence.
! 
!       IER, integer, error indicator.
!          0, no errors encountered.
!          1, if M, N, NQ, NW, or NR is out of range.
!          2, if duplicate nodes were encountered.
!          3, if all the nodes lie in an affine subspace
!             of dimension M-1.
!          4, if space could not be allocated for temporary
!             arrays.
! 
!   Local Variables:
!          I, IB, integer, loop and array indices.
!          IERR, integer, error indicator for STOREM.
!          II,IJ, integer, indices for WS -- elements (I,I)
!             and (I,J) of the augmented regression matrix.
!          IRM1, integer, IROW-1
!          IROW, integer, number of rows currently included in 
!                the regression matrix.
!         J, integer, Row index for X, column indes for A, and
!            loop index.
!         JJ, integer, Index for WS -- element (J,J).
!         K, integer, nodal index for F, Lnext, columns of X,
!            and rows of A
!         LMAX, integer, Maximum number of NPTS elements
!            (must be consistent with the dimension of NPTS).
!         LNP, integer, current length of NTPS.
!         MT, integer, number of quadratic terms in Q(K).
!         NEQ, integer, number of equations in the least squares fit.
!         NI, NJ, integer, indices for WS - elements (NT,I) and
!            (NT,J).
!         NP, integer, NPTS element.
!         NPTS(100), integer, array containing the indices of a 
!           sequence of nodes to be used in the least squares fit
!           or to compute RSQ.  The nodes are ordered by distance
!           from node K, and the last element (usually indexed by
!           LNP) is used only to determine RQ, or RSQ(K) if NW>NQ.
!         NQWMAX, integer, max(NQ, NW).
!         NTT, integer, Number of terms in each quadratic nodal
!            function Q(K).
!         NTM1, integer, NTT-1, the number of coefficients which
!             along with F(K), define Q(K).
!         NTM1NT, integer, (NTT-1)*NTT.
!         NTP1, integer, NTT+1.
!         
!         AV, real, root-mean-square distance between node K and the
!           nodes in the least squares system (unless additional nodes
!           are introduced for stability).  The first MT columns of the
!           matrix (quadratic terms) are scaled by 1/AVSQ, the last
!           M (linear terms) by 1/AV.
!         AVSQ, real, AV*AV.
!         C, real, first component of the plane rotation use to zero
!           the lower triangle of the regression matrix.  Computed
!           by GIVENS.
!         DMIN, real, minimum of the magnitudes of the diagonal elements
!           of the regression matrix after zeros are introduces below
!           the diagonal.
!         DTOL, real, tolerance for detecting an ill-conditioned system.
!         FK, real, data value at node K (F(K)).
!         RQ, real, radius of influence which enters into the weights
!           for Q(K) (see the subroutine SETUPM).
!         RS, real, squared distance between node K and NPTS(LNP-1).
!           Used to compute RQ and RSQ(K).
!         RSMX, real, maximum RSQ element encountered.
!         RSOLD, real, squared distance between node K and NPTS(LNP-1).
!            Used to compute a relative change in RS between succeeding
!            NPTS elements. 
!         RTOL, real, tolerance for detecting a sufficiently large 
!           relative change in RS.  If the changes is not greater than
!           RTOL, the nodes are treated as being the same distance from
!           node K
!         RWS, real current value of RSQ(K).
!         S, real, second component of the plane Givens rotation.
!         SF, real, Marquardt stabilization factor used to damp out
!            the first MT solution components (second partials of the
!            quadratic) when the system is ill-conditioned.  As SF 
!            increases, the fitting function approaches a linear.
!         SUM, real, sum of squared Euclidean distances between node
!            K and the nodes used in the least squares fit (unless
!            additional nodes are added for stability)
!         T, real, temporary variable for accumulating a scalar 
!            product in the back substitution loop.
!
!         NECESSARY, NECESSARY4, NECESSARY7, BRK, logical, conditionals
!            for loops.


 USE REAL_PRECISION
 USE QSHEPMDATA
 IMPLICIT NONE

 INTEGER :: M, N, NQ, NW, NR, IER
 REAL(KIND=R8), DIMENSION(M,N):: X
 REAL(KIND=R8), DIMENSION(N) :: F
 REAL(KIND=R8) :: RMAX

 INTEGER :: I, IB, IERR, II, IJ, IRM1, IROW, J, JJ, K, LMAX, LNP, MT
 INTEGER :: NEQ, NI, NJ, NP, NQWMAX, NTT, NTM1, NTM1NT, NTP1
 INTEGER, DIMENSION(100) :: NPTS

 REAL(KIND=R8) :: AV, AVSQ, C, DMIN, DTOL, FK, RQ, RS, RSMX, RSOLD
 REAL(KIND=R8) :: RTOL, RWS, S, SF, SUM, T

 LOGICAL :: NECESSARY, NECESSARY4, NECESSARY7, BRK

 DTOL=0.01_R8
 RTOL=1.0E-05_R8
 SF=1.0_R8
 
 NQWMAX=MAX(NQ,NW)
 LMAX=MIN(100,N-1)
 NTM1=(M*(M+3))/2
 NTT=NTM1+1
 NTP1=NTT+1
 MT=NTM1-M
 NTM1NT=NTM1*NTT
 C=0.0_R8
 S=0.0_R8
 ALLOCATE(A(N,NTT), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(DX(M), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(XMIN(M), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(RSQ(N), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(LCELL(NR**M), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(LNEXT(N), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(WS(NTT*NTT), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

 ALLOCATE(IW(M,5), STAT=IER)
    IF (IER /= 0) THEN
       IER=4
       RETURN
    END IF

! Make sure the variables are legitimate.  If not, exit the program.
 IF ((M<1) .OR. (NQ<NTM1) .OR. (NW<1) .OR. (NQWMAX>LMAX) .OR. &
     (NR<1)) THEN
 write(6,*)m,nq,ntm1,nw,nqwmax,lmax,nr
       IER=1
       RETURN
 END IF

 CALL STOREM(M, N, X, NR, IERR)
 IF (IERR /= 0) THEN
    IER=3
    RETURN
 END IF


 RSMX=0.0_R8

! Outer loop on node K
 DO K=1,N
    FK=F(K)
! Mark node K to exclude it from the search for nearest neighbors.

    LNEXT(K) = -LNEXT(K)

! Get the index of the nearest neighbor, test for duplicate nodes,
!    and initialize the loop on the remaining NPTS elements.

    CALL GETNPM(M, N, X, K, NR, NP, RS)
    IF (RS == 0.0_R8) THEN
       IER=2
       RETURN
    END IF

    LNP=1
    NPTS(1)=NP
    SUM=RS
    RWS=0.0_R8
    RQ=0.0_R8

! Initialize conditional variable for while-loop
    NECESSARY=.TRUE.

! Initialize break variable
    BRK=.FALSE.

! Compute NPTS, LNP, RWS, RQ, NEQ, and AVSQ
    DO WHILE (NECESSARY)
       LNP=LNP+1
       RSOLD=RS
       CALL GETNPM(M, N, X, K, NR, NP, RS)
       IF (NP == 0) THEN
!  Duplicate node encountered
          IER=2
          RETURN
       END IF

       NPTS(LNP)=NP
       IF (((RS-RSOLD)/RSOLD) < RTOL) THEN
          SUM=SUM+RS
       ELSE
          IF ((RWS == 0.0_R8) .AND. (LNP > NW)) THEN
             RWS=RS
          END IF
          IF ((RQ == 0.0_R8) .AND. (LNP > NQ)) THEN
! RQ=0 (Not yet computed) and LNP>NQ.  RQ=sqrt(RS) is 
! sufficiently large to (strictly include NQ nodes.  The
! least squares fit will include NEQ=LNP-1 equations for
! NT-1 <= NQ < LMAX <= N-1

             RQ=SQRT(RS)
             NEQ=LNP-1
             AVSQ=SUM/NEQ
          END IF

! Test for termination on LNP>NQ and LNP>NW
          IF (LNP>NQWMAX) THEN
             BRK=.TRUE.
             NECESSARY=.FALSE.
          ELSE
             SUM=SUM+RS
          END IF
      END IF
      IF (LNP == LMAX) THEN
         NECESSARY=.FALSE.
      ENDIF
    END DO 

! If BRK=0 all LMAX nodes are included in NPTS.  RWS
! and/or RQ**2 is (arbitrarily) taken to be 10% larger
! than the distance RS to the last node included.

    IF (.NOT. BRK) THEN
       IF (RWS == 0.0_R8) THEN
          RWS=RS*1.1_R8
       END IF
       IF (RQ == 0) THEN
          RQ=SQRT(RS*1.1_R8)
          NEQ=LMAX
          AVSQ=SUM/NEQ
       END IF
    END IF

! Store RSQ(K), update RSMX if necessary, and compute AV.
    RSQ(K)=RWS
    IF (RWS>RSMX) THEN
       RSMX=RWS
    END IF
    AV=SQRT(AVSQ)
   
! Set up the augmented regression matrix (stored by rows in WS,
! and zero out the lower triangle with Givens rotations - QR
! decomposition with orthogonal matrix Q not stored.

    BRK=.FALSE.
    NECESSARY4=.TRUE.
    I=0
    DO WHILE (NECESSARY4)
       I=I+1
       NP=NPTS(I)
       IROW=MIN(I,NTT)
       IRM1=IROW-1
       IJ=IRM1*NTT+1
       CALL SETUPM(M, N, X, F, K, NP, AV, AVSQ, RQ, IJ)
       IF (I /= 1) THEN
          JJ=1
          DO J=1,IRM1
             CALL GIVENS(WS(JJ), WS(IJ), C, S)
             IJ=IJ+1
             JJ=JJ+1
             CALL ROTATE(NTT-J, C, S, JJ, IJ)
             JJ=JJ+NTT
          END DO
          IF (I >= NEQ) THEN
! Test the system for ill-conditioning

             DMIN=ABS(WS(1))
             JJ=1
             DO J=2,NTM1
                JJ=JJ+NTP1
                DMIN=MIN(DMIN,ABS(WS(JJ)))
             END DO
             IF (DMIN*RQ < DTOL) THEN
                IF (NEQ /= LMAX) THEN
! Increase RQ and add another equation to the system to 
! improve the conditioning.  The number of NPTS elements
! is also increased if necessary.
                   NECESSARY7=.TRUE.
                   DO WHILE (NECESSARY7)
                      RSOLD=RS
                      NEQ=NEQ+1
                      IF (NEQ /= LMAX) THEN
                         IF (NEQ < LNP) THEN
                            NP=NPTS(NEQ+1)
                            RS=0.0_R8
                            DO J=1,M
                               RS=RS+(X(J,NP)-X(J,K))*(X(J,NP)-X(J,K))
                            END DO
                         ELSE
                            NEQ=LNP
                            LNP=LNP+1
                            CALL GETNPM(M, N, X, K, NR, NP, RS)
                            IF (NP == 0) THEN
                               IER=2
                               RETURN
                            END IF
                            NPTS(LNP)=NP
                         END IF
                         IF ((RS-RSOLD)/RSOLD >= RTOL) THEN
                            RQ=SQRT(RS)
                            NECESSARY7=.FALSE.
                         END IF
                      ELSE
                         RQ=SQRT(RS*1.1_R8)
                         NECESSARY7=.FALSE.
                      END IF
                   END DO 

                ELSE
                   NECESSARY4=.FALSE.
                   BRK=.TRUE.
                END IF
             ELSE
                NECESSARY4=.FALSE.
             END IF
          END IF
       END IF

    END DO

! Stabilize the system by damping second partials - add multples of
! the first MT unit vectors to the first MT equation, where
! MT=M*(M+1)/2 is the number of quadratic terms.

      IF (BRK) THEN
         NI=NTM1NT
         DO I=1,MT
            NI=NI+1
            WS(NI)=SF
            NJ=NI
            DO J=I+1, NTT
               NJ=NJ+1
               WS(NJ)=0.0_R8
            END DO
            JJ=NTP1*I-NTT
            NJ=NI
            DO J=I,NTM1
               CALL GIVENS(WS(JJ), WS(NJ), C, S)
               JJ=JJ+1
               NJ=NJ+1
               CALL ROTATE(NTT-J, C, S, JJ, NJ)
               JJ=JJ+NTT
            END DO
         END DO
! Test the stabilized system for ill-conditioning.
         DMIN=ABS(WS(1))
         JJ=1
         DO J=2,NTM1
            JJ=JJ+NTP1
            DMIN=MIN(DMIN,ABS(WS(JJ)))
         END DO

         IF (DMIN*RQ <DTOL) THEN
            IER=3
            RETURN
         END IF
      END IF

! Solve the order NT-1 upper triangular system for the coefficients.

     II=NTM1NT-1
     DO IB=1,NTM1
        I=NTT-IB
        T=0.0_R8
        IF (I /= NTM1) THEN
           IJ=II
           DO J=I+1,NTM1
              IJ=IJ+1
              T=T+WS(IJ)*A(K,J)
           END DO
        END IF
        A(K,I)=(WS(II+IB)-T)/WS(II)
        II=II-NTP1
     END DO

! Scale the coefficients to adjust for the column scaling.

     DO I=1,MT
        A(K,I)=A(K,I)/AVSQ
     END DO

     DO I=MT+1, NTM1
        A(K,I)=A(K,I)/AV
     END DO

! Unmark K and the elements of NPTS.

     LNEXT(K)= -LNEXT(K)
     DO I=1,LNP
        NP=NPTS(I)
        LNEXT(NP)= -LNEXT(NP)
     END DO
 END DO
 
! No errors encountered

 RMAX=SQRT(RSMX)
 IER=0
 RETURN
 END SUBROUTINE QSHEPM

!***********************************************************************
!                                                   Robert Renka
!                                           Univ. of North Texas
!                                                 (817) 565-2767
!                                                       10/25/87
! 
!                                            Translated into C++
!                                                by Karen Minser
!                                             Univ. of Tennessee
!                                                         7/6/98
!  
!                                            Translated into F95
!                                             by William Thacker
!                                            Winthrop University
!                                                      1/10/2009
! 
!                      QSMVAL
!   This function returns the value Q(P) where Q is the weighted sum
! of the quadratic nodal functions defined in the function QSHEPM.
!   Upon completion, QSMVAL returns the function value Q(P) unless
! M, N, NR, an element of DX or RMAX is invalid.   
!***********************************************************************

 FUNCTION QSMVAL(M, N, P, X, F, NR, RMAX)

! Input Parameters:
!    M, integer, dimension of the problem.
!   
!    N, integer, number of nodes and data values.
! 
!    P(M), real, vector of length m containing the Cartesian 
!      coordinates of the point at which Q is to be evaluated.
! 
!    X(M,N), real, the node coordinates of the N, M-dimensional points
! 
!    F(N), real, the data values at the N points specified by X.
! 
!    NR, integer, the number of divisions used in the cell method.
!
!    RMAX, real, square root of the largest element in RSQ, the maximum
!          radius of influence.
! 
! Local Variables:
!    I, J, K, integer, loop and array indices.
!    KP, integer, temporary variable for K and index into LNEXT.
!    IMIN, integer, minimum cell index at I.
!    IMAX, integer, maximum cell index at I.
!    L, integer, counter and index.
!    LCI, integer, index into LCELL.
!    MT, integer, number of quadratic terms in the nodal functions Q(K).
!    NNR, integer, local variable for NR.
!    NTT, integer, total number of terms in Q(K).
!    NRP, integer, affects LCI index in proportion to NR.
!    DS, real, used in determining weights.
!    QK, real, value of Q(K)(P).
!    RD, real, used in determining weights.
!    RDS, real, used in determining weights.
!    RM, real, local variable for RMAX.
!    RS, real, RSQ(k).
!    SW, real, accumulated weight values.
!    SWQ, real, accumulated weighted nodal function values.
!    T, real, variable used in determining QK.
!    WW, real, W(K).
!    QSMVAL, real, final value Q(P) that is returned.
!    CELLS, logical, conditional for do-while that loops over cells.
!    NODES, logical, conditional for do-while that loops over nodes.
!    BRK, logical, signals if there was as break out of a do-while loop.

 USE REAL_PRECISION
 USE QSHEPMDATA
 IMPLICIT NONE
 
 INTEGER :: M, N, NR
 REAL(KIND=R8), DIMENSION(M) :: P
 REAL(KIND=R8), DIMENSION(M,N) :: X
 REAL(KIND=R8), DIMENSION(N) :: F
 REAL(KIND=R8) :: RMAX

 INTEGER :: I, J, K, KP, IMIN, IMAX, L, LCI, MT, NNR, NTT, NRP
 REAL(KIND=R8) :: DS, QK, RD, RDS, RM, RS, SW, SWQ, T, WW, QSMVAL
 LOGICAL :: CELLS, NODES, BRK

 NNR=NR
 RM=RMAX
 MT=(M*(M+1))/2
 NTT=MT+M+1

 IF ((M<1) .OR. (N<NTT) .OR. (NNR<1) .OR. (RM<0.0_R8)) THEN
    QSMVAL=0.0_R8
    RETURN
 END IF

! Cell (L(1), L(2),..., L(M)) has index LCI=SUM[(L(I)-1)*
!   NR**(I-1)]+1 where the sum is over I=1 to M.  Set
!   IW( ,1) and IW( ,2) to the minimum and maximum cell indices
!   (M-tuples) defining the range of the search for nodes
!   whose radii include P.  The cells which must be searched
!   are those intersected by an M-ball of radius RMAX centered
!   at P.  No cells are within RMAX of P if IMIN=
!   IW(I,1)>IMAX=IW(I,2) for any I.  IW( ,3) is initialized
!   to IW( ,1) and LCI is initialized to the index
!   associated with IW( ,3) by Horner's method.

 LCI=0
 DO I=M,1,-1

    IF (DX(I)<=0.0_R8) THEN
       QSMVAL=0.0_R8
       RETURN
    END IF

    IMIN=INT((P(I)-XMIN(I)-RM)/DX(I))+1
    IMIN=MAX(IMIN,1)
    IMAX=INT((P(I)-XMIN(I)+RM)/DX(I))+1
    IMAX=MIN(IMAX,NNR)

    IF (IMIN>IMAX) THEN
       QSMVAL=0.0_R8
       RETURN
    ENDIF

    LCI=LCI*NNR+IMIN-1
    IW(I,1)=IMIN
    IW(I,2)=IMAX
    IW(I,3)=IMIN
 END DO

 LCI=LCI+1

! Accumulate weight values in SW and weighted nodal function
! values in SWQ.  The weights are W(K)=( (R-D)+ / (R*D))**2
! for R**2=RSQ(k) and D=distance between P and node K.

 SW=0.0_R8
 SWQ=0.0_R8

! Outer loop on cells LCI

 CELLS=.TRUE.
 DO WHILE (CELLS)
    BRK=.FALSE.
    K=LCELL(LCI)
    IF (K /= 0) THEN
       NODES=.TRUE.
! Inner loop on nodes K.  Compute WW=W(K) and update SW.
       DO WHILE (NODES)
          DS=0.0_R8
          DO I=1,M
             DS=DS+(P(I)-X(I,K))*(P(I)-X(I,K))
          END DO
          RS=RSQ(K)
          IF (RS>DS) THEN
             IF (DS /= 0.0_R8) THEN
                RDS=RS*DS
                RD=SQRT(RDS)
                WW=(RS+DS-RD-RD)/RDS
                SW=SW+WW
! compute QK=Q(K)(P) and update SWQ
                QK=F(K)
                L=0
                DO J=1,M
                   T=A(K,MT+J)
                   DO I=1,J
                      L=L+1
                      T=T+A(K,L)*(P(I)-X(I,K))
                   END DO
                   QK=QK+T*(P(J)-X(J,K))
                END DO
                SWQ=SWQ+WW*QK
             ELSE
                QSMVAL=F(K)
                RETURN
             END IF
          END IF
          KP=K
          K=LNEXT(KP)
          IF (K == KP) THEN
             NODES=.FALSE.
          END IF
       END DO
    END IF

    NRP=1
  CHECK:  DO I=1,M
       IF (IW(I,3) < IW(I,2)) THEN
          IW(I,3)=IW(I,3)+1
          LCI=LCI+NRP
          BRK=.TRUE.
          EXIT CHECK
       ELSE
          IF (I<M) THEN
             IW(I,3)=IW(I,1)
             LCI=LCI-(IW(I,2)-IW(I,1))*NRP
             NRP=NRP*NNR
          END IF
       ENDIF
    END DO CHECK
    IF (.NOT. BRK) THEN
       CELLS=.FALSE.
    ENDIF
 END DO
 
! SW=0 iff P is not within the radius R(K) for any node K

 IF (SW == 0.0_R8) THEN
    QSMVAL=0.0_R8
 ELSE
    QSMVAL=SWQ/SW
 ENDIF

 RETURN
 END FUNCTION QSMVAL


!***********************************************************************
!                                                   Robert Renka
!                                           Univ. of North Texas
!                                                 (817) 565-2767
!                                                       10/25/87
! 
!                                            Translated into C++
!                                                by Karen Minser
!                                             Univ. of Tennessee
!                                                         7/6/98
!  
!                                            Translated into F95
!                                             by William Thacker
!                                            Winthrop University
!                                                      1/10/2009
! 
!                      GETNPM
!   Given a set of N nodes and the data structure defined in
! the subroutine STOREM, this subroutine uses the cell method to
! find the closest unmarked node NP to a specified point P.  
! NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
! node is marked if and only if the corresponding LNEXT element
! is negative.  The absolute values of LNEXT elements, 
! however, must be preserved.)  Thus, the closest j nodes to 
! P may be determined by a sequence of j calls to this routine.
! Note that if the nearest neighbor to node K is to 
! be determined (P=X( ,K)), then K should be marked before
! the call to this routine.
!   The search is begun in the cell containing (or closest
! to) P and proceeds outward in (M-dimensional) rectangular
! layers until all cells which contain points within distance
! R of P have been searched, where R is the distance
! from P to the first unmarked node encountered (infinite if
! no unmarked nodes are present).
!***********************************************************************

 SUBROUTINE GETNPM(M, N, X, K, NR, NP, DSQ)

! Input Parameters:
!        M, integer, dimension of the problem
! 
!        N, integer, the number of nodes and associated 
!           data values.
! 
!        X(M,N), real, the M-coordinates of the N nodes.
! 
!        K, integer, current node whose nearest unmarked neighbor
!            is to be found.

!        NR, integer, the number of divisions in each dimension
!            for the cell grid defined in the subroutine
!            STOREM.  A hyperbox containing the nodes is
!            partitioned into cells in order to increase
!            search efficiency.  NR=(N/3)**(1/M) is
!            recommended.  NR>=1.
! 
!   Output Parameters:
!       NP, integer, node index (column index of X) of the
!           nearest unmarked node to P=(X( ,K)) or 0 if all
!           nodes are marked, M<1, NR<1, or an element of DX is
!           not positive. LNEXT(NP)<0 unless NP=0.
!
!       DSQ, real, squared Euclidean distance between P and node
!          NP, or 0 if NP=0.
! 
! Local Variables:
! I, J, integer, loop indices.
! LCI, integer, index for LCELL.
!      Sum[(IW(I,3)-1)*NR**(I-1)]+1 where the sum is over I=1 to M
! L, integer values of LCELL used as indices as well.
! LL, integer, temporary variable for determining minimum and maximum
!     values.
! LMIN, integer, current node index of the nearest unmarked node.
! LN, integer, value of LNEXT used as indices.
! RRSQ, real, local RSQ.
! RSMIN, real, RSQ of LMIN.
! R, real, square root of RSMIN used to find min and max nodal 
!    indices of the search range of the nearest neighbor.
! TMP, real, temporary for updating IW.
! FIRST, logical, signals the first unmarked neighbor of K.  Is equal
!    to TRUE iff the first unmarked node has not yet been encountered.
! NODES, logical, conditional for do-while loop over nodes.
! CELLS, logical, conditional for do-while loop over cells.
! LAYERS, logical, conditional for do-while loop over layers.
! BRK, logical, signal for breaking out of loop.

!Description of use of IW.
! IW(I,4), IW(I,5), Minimum and maximum cell indices, respectively, 
!        defining the range of the search.
! IW(I,1), IW(I,2), Minimum and maximum cell indices, respectively,
!        defining the layer whose intersection with the range is 
!        currently being searched.  These are initialized to the 
!        cell containing or closest to P.
! IW(I,3), Indices of the current cell in the range
!        max(IW(I,1),IW(I,4)) to min(IW(I,2),IW(I,5)).

 USE REAL_PRECISION
 USE QSHEPMDATA
 IMPLICIT NONE
 
 INTEGER :: M, N, K, NR, NP
 REAL(KIND=R8), DIMENSION(M,N) :: X
 REAL(KIND=R8) :: DSQ

 INTEGER :: I, J, LCI, L, LL, LMIN, LN
 REAL(KIND=R8) :: RRSQ, RSMIN, R, TMP
 LOGICAL :: FIRST, NODES, CELLS, LAYERS, BRK

 L=0
 LN=0
 RSMIN=0.0_R8
 IF ((M<1) .OR. (NR<1)) THEN
    NP=0
    DSQ=0.0_R8
    RETURN
 END IF

 FIRST=.TRUE.

 DO I=1,M
    IW(I,4)=1
    IW(I,5)=NR
    IF (DX(I)<=0.0_R8) THEN
       NP=0
       DSQ=0.0_R8
       RETURN
    END IF

    LL=MIN(NR, INT((X(I,K)-XMIN(I))/DX(I))+1)
    IW(I,1)=MAX(1,LL)
    IW(I,2)=IW(I,1)
 END DO

 LAYERS=.TRUE.

! Top of outer loop on layers: initialize IW(I,3)

 DO WHILE (LAYERS)
    BRK=.FALSE.
    DO I=1,M
       IW(I,3)=MAX(IW(I,1),IW(I,4))
    END DO

    CELLS=.TRUE.
! Top of loop on cells: bypass cells interior to the layer.

    DO WHILE (CELLS)
       BRK=.FALSE.
       DO I=1,M
          IF ((IW(I,3) == IW(I,1)) .OR. (IW(I,3) == IW(I,2))) THEN
             BRK=.TRUE.
             EXIT
          END IF
       END DO
       
       IF (BRK) THEN
          LCI=0
! Compute LCI by Horner's method and test for an empty cell.
          DO J=M,1,-1
             LCI=LCI*NR+IW(J,3)-1
          END DO

          LCI=LCI+1
          L=LCELL(LCI)
          IF (L /= 0) THEN
             NODES=.TRUE.
! Loop on nodes in cell LCI
             DO WHILE (NODES)
                LN=LNEXT(L)
                IF (LN>=0) THEN
! Node L is not marked.  Set RSQ to its squared distance from P.
                   RRSQ=0.0_R8
                   DO J=1,M
                      RRSQ=RRSQ+(X(J,K)-X(J,L))*(X(J,K)-X(J,L)) 
                   END DO
                   IF (FIRST) THEN
! Node L is the first unmarked neighbor of P encountered
! and hence the first candidate for NP.  Initialize LMIN
! and RSMIN to L and RSQ, and update the search range
! IW(I,4) and IW(I,5) to the smallest M-box containing
! a hypersphere of radius R=sqrt(RSMIN) centered at P
! and contained in [1][NR]**M.  First is reset to FALSE.
                      LMIN=L
                      RSMIN=RRSQ
                      R=SQRT(RSMIN)
                      DO J=1,M
                         IW(J,4)=MAX(1,INT((X(J,K)-R-XMIN(J))/DX(J))+1)
                         IW(J,5)=MIN(NR,INT((X(J,K)+R-XMIN(J))/DX(J))+1)
                      END DO
                      FIRST=.FALSE.
                   ELSE
                      IF (RRSQ<RSMIN) THEN
                         LMIN=L
                         RSMIN=RRSQ
                      END IF
                   END IF
                END IF
! Test for termination of loop on nodes in cell LCI.
                IF (ABS(LN) /= L) THEN
                   L=ABS(LN)
                ELSE
                   NODES=.FALSE.
                END IF
             END DO
          END IF
! Bottom of loop on cells.  Update the current cell indices.
          BRK=.FALSE.
          DO J=1,M
             TMP=MIN(IW(J,2), IW(J,5))
             IF (IW(J,3)<TMP) THEN
                IW(J,3)=IW(J,3)+1
                BRK=.TRUE.
                EXIT
             ELSE
                IW(J,3)=MAX(IW(J,1),IW(J,4))
             END IF
          END DO

          ELSE
! Bottom of loop on cells.  Update the current cell items.
             DO J=1,M
                TMP=MIN(IW(J,2),IW(J,5))
                IF (IW(J,3)<TMP) THEN
                   IW(J,3)=IW(J,3)+1
                   BRK=.TRUE.
                   EXIT
                ELSE
                   IW(J,3)=MAX(IW(J,1),IW(J,4))
                ENDIF
             END DO
          END IF
    
! If no break-out above, then drop out of do-while on cells.

          IF (.NOT. BRK) THEN
             CELLS=.FALSE.
          ENDIF
       END DO

! Test for termination of loop on cell layers.

       BRK=.FALSE.
       DO I=1,M
          IF ((IW(I,1)>IW(I,4)) .OR. (IW(I,2)<IW(I,5))) THEN
! Update IW(I,1) and IW(I,2) to the next layer out.
             DO J=1,M
                IW(J,1)=IW(J,1)-1
                IW(J,2)=IW(J,2)+1
             END DO
             BRK=.TRUE.
             EXIT
          END IF
        END DO

        IF (.NOT. BRK) THEN
           LAYERS=.FALSE.
        END IF
     END DO

! Unless no unmarked nodes were encountered, LMIN is the closest
! unmarked node to P.

     IF (FIRST) THEN
        NP=0
        DSQ=0.0_R8
        RETURN
     ELSE
        NP=LMIN
        DSQ=RSMIN
        LNEXT(LMIN)= -LNEXT(LMIN)
        RETURN
     END IF
 END SUBROUTINE GETNPM


!   *************************************************************
!                         GIVENS
!                                              ROBERT RENKA
!                                      UNIV. OF NORTH TEXAS
!                                            (817) 565-2767
!                                                  10/25/87
!
!                                       Translated into C++
!                                           by Karen Minser
!                                        Univ. of Tennessee
!                                                    7/6/98
!
!                                       Translated into F95
!                                        by William Thacker
!                                       Winthrop University
!                                                 1/10/2009
! 
!              Givens
! This routine constructs the Givens plan rotation
!         ( C S)
!    G =  (    ), where C*C+S*S=1, which zeros the second
!         (-S C)
! entry of the 2-vector (A B)-Transpose.  A call to GIVENS
! is normally followed by a call to ROTATE, which applies
! the transformation to a 2 by N matrix.  This routine was
! taken from LINPACK.

 SUBROUTINE GIVENS(A, B, C, S)

! Input parameters:
!  A, B, real, components of the 2-vector to be rotated
!
! Output parameters:
!  A, real, overwritten by R= +|- sqrt(A*A+B*B)
!  B, real, overwritten by a value Z which allows C and
!    S to be recovered as follows:
!      C=sqrt(1-Z*Z), S=Z  if abs(z)<=1.
!      C=1/Z, S=sqrt(1-C*C) if abs(z)>1.
!  C, real, + | - (A/R).
!  S, real, + | - (B/R).
!
! Local variables:
!  AA, real, local copy of A.
!  BB, real, local copy of B.
!  R, real, C*A+S*B = +|- sqrt(A*A+B*B).
!  U, V, real, variables used to scale A and B for computing R.


 USE REAL_PRECISION
 REAL(KIND=R8) :: A, B, C, S
 REAL(KIND=R8) :: AA, BB, R, U, V

 AA=A
 BB=B
 IF (ABS(AA) > ABS(BB)) THEN
    U=AA+AA
    V=BB/U
    R=SQRT(V*V+0.25_R8)*U
    C=AA/R
    S=V*(C+C)
! Note: R has the sign of A, C>0, and S has Sign(A)*Sign(B).
    B=S
    A=R
    RETURN
 ELSE
! abs(AA)<=abs(BB).
    IF (BB == 0.0_R8) THEN
       C=1.0_R8
       S=0.0_R8
       RETURN
    ELSE
       U=BB+BB
       V=AA/U
       A=SQRT(V*V+0.25_R8)*U
       S=BB/A
       C=V*(S+S)
       B=1
       IF (C /= 0.0_R8) THEN
          B=1.0_R8/C
          RETURN
       END IF
    END IF
 END IF
 RETURN
 END SUBROUTINE GIVENS

!   *************************************************************
!              ROTATE
!                                              ROBERT RENKA
!                                      UNIV. OF NORTH TEXAS
!                                            (817) 565-2767
!                                                  10/25/87
!
!                                       Translated into C++
!                                           by Karen Minser
!                                        Univ. of Tennessee
!                                                    7/6/98
!
!                                       Translated into F95
!                                        by William Thacker
!                                       Winthrop University
!                                                 1/10/2009
! 
!                                          ( C S)
! This routine applies the Givens rotation (    ) to the
!                                          (-S C) 
!
!                (X(1) ... X(N))
! 2 by N matrix  (             ).
!                (Y(1) ... Y(N))
!
!
! Input Parameters:
!   N, integer, the number of columns to be rotated.
!
!   C, S, real, elements of the Givens rotation.  These
!        may be determined by the subroutine GIVENS.
!
!   II, JJ, integer, indices into arrays of length >= N
!        containing the vectors to be rotated.
!
! Local Variables:
!  I, integer, loop index.
!  XI, YI, real, values from the X and Y vectors
! 


! Note: instead of passing vectors which are part of a larger 
! vector as was done in the original FORTRAN code, the indices
! marking the smaller vectors were passed instead.  The 
! larger vector, WS is in the module qshepmdata.   

 SUBROUTINE ROTATE(N, C, S, IJ, JJ)
 USE REAL_PRECISION
 USE QSHEPMDATA
 IMPLICIT NONE

 INTEGER :: N, IJ, JJ
 REAL(KIND=R8) :: C, S

 INTEGER :: I
 REAL(KIND=R8) :: XI, YI

 IF ((N <= 0) .OR. (C == 1.0_R8) .OR. (S == 0.0_R8)) THEN
    RETURN
 END IF

 DO I=1,N
    XI=WS(I+IJ-1)
    YI=WS(I+JJ-1)
    WS(I+IJ-1)=C*XI+S*YI
    WS(I+JJ-1)=-S*XI+C*YI
 END DO
 RETURN
 END SUBROUTINE ROTATE

!   *************************************************************
!                         SETUPM
!                                              ROBERT RENKA
!                                      UNIV. OF NORTH TEXAS
!                                            (817) 565-2767
!                                                  10/25/87
!
!                                       Translated into C++
!                                           by Karen Minser
!                                        Univ. of Tennessee
!                                                    7/6/98
!
!                                       Translated into F95
!                                        by William Thacker
!                                       Winthrop University
!                                                 1/10/2009
! 
! This routine sets up the L-th row of an augmented regression
! matrix for a weighted least squares fit of a quadratic
! function Q(X) to a set of data values F, where Q(XK)=FK and
! Q is a function of M variables (X is an element of 
! Euclidean M-space).  The upper triangular array of
! M(M+1)/2 quadratic terms (XL(I)-XK(I))*(XL(J)-XK(J)), I<J,
! are ordered by columns and scaled by W/S2.  These are 
! followed by the M linear terms, scaled by W/S1, and the 
! last element of ROW is the right hand side W*(FL-FK).  The
! weight is W=(R-D)/(R*D) if R>D, or W=0 if R<=D, where D is the
! distance between nodes L and K.
! 
!   Input Parameters:
!        M, integer, dimension of the problem
! 
!        N, integer, the number of nodes and associated 
!           data values.
! 
!        X(M,N), real, the M-coordinates of the N nodes.
! 
!        F(N), real, the data values at the nodes.
!
!        K, integer, current node (column of X) which 
!           contains the Cartesian coordinates of node K.
!
!        NP, integer, index morking the column of X which
!           contains the coordinates of node L.  Also, the
!           index into array F which is the data value at 
!           node L.
!
!        S1, S2, real, reciprocals of the scale factors for
!           the linear and quadratic terms, respectively.
! 
!        R, real, radius of influence about node K defining
!           the weight W.
!
!        IROW, integer, index into global array WS to store
!           values.
! 
!   Local Variables:
!       I, J, integer, loop indices.
!       IJ, integer, counter and index.
!       MT, integer, number of quadratic terms in the nodal functions
!           Q(K).
!       NNT, integer, total number of terms in Q(K).
!       ISTART, offset into WS to start storing.
!       D, real, distance between nodes L and K.
!       W0, W1, W2, real, weight values.
!       DXJ, real, Difference between respective coordinates of L and K.
 

  SUBROUTINE SETUPM(M, N, X, F, K, NP, S1, S2, R, IROW)
  USE REAL_PRECISION
  USE QSHEPMDATA
  IMPLICIT NONE

  INTEGER :: M, N, K, NP, IROW
  REAL(KIND=R8) :: S1, S2, R
  REAL(KIND=R8), DIMENSION(M,N) :: X
  REAL(KIND=R8), DIMENSION(N) :: F

  INTEGER :: I, J, IJ, MT, NNT, ISTART
  REAL(KIND=R8) :: D, W0, W1, W2, DXJ, FK

  ISTART=IROW-1
  MT=(M*(M+1))/2
  NNT=MT+M+1
  FK=F(K)
  IF (M >= 1) THEN
! Compute distance D, weights W=(R-D)/(R*D), W1=W/S1
     D=0.0_R8
     DO I=1,M
        D=D+(X(I,NP)-X(I,K))*(X(I,NP)-X(I,K))
     END DO

     D=SQRT(D)
     IF ((D>0.0_R8) .AND. (D < R)) THEN
        W0=(R-D)/R/D
        W1=W0/S1
        W2=W0/S2
! Store the row
        IJ=0
        DO J=1,M
           DXJ=X(J,NP)-X(J,K)
           DO I=1,J
              IJ=IJ+1
              WS(IJ+ISTART)=(X(I,NP)-X(I,K))*DXJ*W2
           END DO
           WS(MT+J+ISTART)=DXJ*W1
        END DO

        WS(NNT+ISTART)=(F(NP)-FK)*W0
        RETURN
     END IF
  END IF
  
! M<1, D=0 (nodes K and L coincide), or node L is outside the
! radius of influence.  Set WS to the zero vector. 

  DO I=1,NNT
     WS(I+ISTART)=0.0_R8
  END DO
  RETURN
  END SUBROUTINE SETUPM


!   *************************************************************
!                         STOREM
!                                              ROBERT RENKA
!                                      UNIV. OF NORTH TEXAS
!                                            (817) 565-2767
!                                                  10/25/87
!
!                                       Translated into C++
!                                           by Karen Minser
!                                        Univ. of Tennessee
!                                                    7/6/98
!
!                                       Translated into F95
!                                        by William Thacker
!                                       Winthrop University
!                                                 1/10/2009
! 
! Given a set of N nodes arbitrarily distributed in
! euclidean M-space, this function creates a data structure
! for a cell-based method for solving closest-point problems.
! The smallest M-dimensional box containing the nodes
! is partitioned into an NR**M uniform grid of cells, and the
! indices of the nodes contained in each cell are stored in the 
! data structure.  For a uniform random distribution of nodes,
! the nearest node to an arbitrary point can be determined in
! constant expected time.  Refer to the function GETNPM.
! 
!   Input Parameters:
!        M, integer, dimension of the problem
! 
!        N, integer, the number of nodes and associated 
!           data values.
! 
!        X(M,N), real, the M-coordinates of the N nodes.
! 
!        NR, integer, the number of divisions in each dimension
!            for the cell grid defined in the subroutine
!            STOREM.  A hyperbox containing the nodes is
!            partitioned into cells in order to increase
!            search efficiency.  NR=(N/3)**(1/M) is
!            recommended.  NR>=1.
! 
!   Output Parameters:
!       IER, integer, error indicator.
!          0, no errors encountered.
!          1, if N<2 or NR<1.
!          2, if a component of DX is not positive.
! 
! 
! Upon completion:
! LCELL - Array of length NR**M containing nodal indices of the
! first node in each cell, with 0 entries corresponding to empty
! cells.  THe nodes in each cell are ordered by their indices,
! so that the first node refers to the one with smallest index.
! Cells are ordered in the manner of FORTRAN array storage, with
! the first subscript varying fastest.  Thus, LCELL(LCI) contains
! the index of the first node in cell (L(1), L(2),..., L(M)) for
! LCI=SUM((L(I)-1*NR**(I-1)))+1 (with the sum over I=1 to M), 
! where the cell is defined by intervals [XMIN(I)+(L(I)-1)*DX(I),
! XMIN(I)+L(I)*DX(I)], for each L(I) in the range 1 to NR.  LCELL
! is not defined if IER /= 0 on output.
! 
! LNEXT - Array of lenght N containing next-node indices such that
! LNEXT(K) is the index of the next node in the cell which contains
! node K, or LNEXT(K)=K if K is the last node in the cell.  If, for 
! example, cell LCI contains nodes 2, 3, and 5 (and no others), 
! then LCELL(LCI)=2, LNEXT(2)=3, LNEXT(3)=5, and LNEXT(5)=5.  LNEXT
! is not defined if IER /= 0 on output.
! 
! XMIN - Array of length M containing the minimum nodal coordinates
! unless IER=1.
! 
! DX - Array of length M containing the cell dimensions (interval
! lengths) unless IER=1.  DX(I)=(XMAX(I)-XMIN(I))/NR where XMIN
! and XMAX contain the extrema defining the smallest M-box which
! contains the nodes.

!   Local Variables:
!   I, K, integer, loop indices.
!   L, integer index for LNEXT.
!   LCI, integer, index into LCELL.
!   LI, integer, used to help determine LCI.
!   LN, integer, used with L.
!   NC, integer, loop index = total number of cells.
!   RNR, real, local variable for NR.


  SUBROUTINE STOREM(M, N, X, NR, IER)
  USE REAL_PRECISION
  USE QSHEPMDATA
  IMPLICIT NONE
 
  INTEGER :: M, N, NR, IER
  REAL(KIND=R8), DIMENSION(M,N) :: X

  INTEGER :: I, K, L, LCI, LI, LN, NC
  REAL(KIND=R8) :: RNR 

  IF ((M<1) .OR. (N<2) .OR. (NR<1)) THEN
     IER=1
     RETURN
  END IF

! Compute the corner coordinates XMIN and XMAX (stored temporarily
! in DX) of the smallest box containing the nodes.

  DO I=1,M
     XMIN(I)=X(I,1)
     DX(I)=X(I,1)
  END DO

  DO K=2, N
     DO I=1,M
        XMIN(I)=MIN(XMIN(I),X(I,K))
        DX(I)=MAX(DX(I),X(I,K))
     END DO
  END DO

! Compute the interval widths and test for a zero width.

  RNR=NR
  DO I=1,M
     DX(I)=(DX(I)-XMIN(I))/RNR
     IF (DX(I) == 0.0_R8) THEN
        IER=2
        RETURN
     END IF
  END DO

! Initialize LCELL to zeros.

  NC=NR**M
  DO LCI=1,NC
     LCELL(LCI)=0
  END DO

! Outer loop on nodes.

  DO K=1,N
! Compute the index LCI of the cell containing node K
! using Horner's method.

     LCI=0
     DO I=M,1,-1
        LI=((X(I,K)-XMIN(I))/DX(I))+1 
        LI=MIN(LI,NR)
        LCI=LCI*NR+LI-1
     END DO
     LCI=LCI+1

! Add K to the data structure as the last node in cell LCI.
     LN=LCELL(LCI)
     IF (LN == 0) THEN
! K is the first node in cell LCI.
        LCELL(LCI)=K
     ELSE
! Find the index L such that LNEXT(L)=L
        L=0
        DO WHILE (LN /= L) 
           L=LN
           LN=LNEXT(L)
        END DO
        LNEXT(L)=K
     END IF
     
! Bottom of loop: set LNEXT(K)=K to indicate that K is the 
! last node in the cell.

     LNEXT(K)=K
  END DO
  IER=0
  RETURN
  END SUBROUTINE STOREM
END MODULE QSHEPMD_MOD
