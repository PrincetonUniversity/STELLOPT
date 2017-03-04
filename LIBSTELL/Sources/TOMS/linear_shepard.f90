MODULE LINEAR_SHEPARD_MOD
INTERFACE 
   SUBROUTINE DGELSS(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
        WORK, LWORK, INFO)
     USE REAL_PRECISION
     INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
     INTEGER, INTENT(OUT) :: RANK, INFO
     REAL(KIND = R8), INTENT(IN) :: RCOND
     REAL(KIND = R8), DIMENSION(LDA, N), INTENT(INOUT) :: A
     REAL(KIND = R8), DIMENSION(LDB), INTENT(INOUT) :: B
     REAL(KIND = R8), DIMENSION(N), INTENT(OUT) :: JPVT
     REAL(KIND = R8), DIMENSION(LWORK), INTENT(OUT) :: WORK
   END SUBROUTINE DGELSS
END INTERFACE
PRIVATE
PUBLIC LSHEP, LSHEPVAL, RIPPLE
CONTAINS
    
SUBROUTINE LSHEP ( M, N, X, F, A, RW, IER, RLSHEP )
! This subroutine computes a set of parameters defining a function that 
! interpolates N data values F(i) at scattered nodes X(i) in M dimensions. The 
! interpolant may be evaluated at an arbitrary point by the function LSHEPVAL.
! The interpolation scheme is a modified linear Shepard method, and 
! can use either the Shepard weighted least squares fit or robust 
! M-estimation for each local approximation.
!
! Input parameters:  
!   M is the dimension of the data.
!   N is the number of nodes.
!   X(M, N) contains the coordinates of the nodes.
!   F(N) contains the function values of the nodes.
! Output parameters:
!   A(M, N) contains the coefficients of the local fits.
!   RW(N) is an array containing the radius of influence for each point.
!   IER 
!   = 0, if no errors were encountered.
!   = 1, if N is too small relative to M.
!   = 2, if any least squares problem is rank deficient, in which
!        case an SVD based minimum norm solution is returned; the
!        local fit should still be reasonable.
!   = 3, if the IRLS subroutine returns an error, meaning the IRLS
!        iteration has failed, and the returned local fit may 
!        be bad.
! Optional input:
!   RLSHEP specifies by its presence that robust M-estimation is to be used to 
!          compute the local approximation.
USE REAL_PRECISION
IMPLICIT NONE
    
INTEGER, INTENT(IN) :: M, N
REAL(KIND=R8), DIMENSION(M,N), INTENT(IN) :: X
REAL(KIND=R8), DIMENSION(N), INTENT(IN) :: F
REAL(KIND=R8), DIMENSION(M,N), INTENT(OUT) :: A
REAL(KIND=R8), DIMENSION(N), INTENT(OUT) :: RW
INTEGER, INTENT(OUT) :: IER
LOGICAL, OPTIONAL :: RLSHEP

! Local variables.         
INTEGER :: FUNC  ! IRLS influence function switch, 1 = Huber, 2 = Tukey's 
                 ! bisquare.
INTEGER :: IER1  ! error indicator for IRLS subroutine.
INTEGER :: INFO  ! error indicator for DGELSS.
INTEGER :: ITER  ! number of iterations the IRLS algorithm will perform. 
INTEGER :: I, J, K  ! loop control variables.
INTEGER :: LWORK ! length of the work array WORK for DGELSS.
INTEGER :: NP    ! number of nodes used for the local fit.
INTEGER :: RANK  ! rank of the array (used for DGELSS).
INTEGER, DIMENSION(N-1) :: IDIST ! index array for DIST. 
REAL(KIND=R8) :: DIAM  ! diameter of the whole date set.
REAL(KIND=R8) :: MAR   ! the median absolute residual.
REAL(KIND=R8) :: RCOND ! reciprocal condition number used in DGELSS.
REAL(KIND=R8) :: RP    ! least squares fit radius R_p.
REAL(KIND=R8), DIMENSION(M) :: JPVT   ! array used for DGELSS.
REAL(KIND=R8), DIMENSION(N-1) :: DIST ! array that stores all the distances. 
REAL(KIND=R8), DIMENSION(MIN(N-1, (3*M+1)/2)) :: BETA, BETA_IN 
               ! variables used in DGELSS.
REAL(KIND=R8), DIMENSION(MIN(N-1, (3*M+1)/2)) :: OMEGA   
               ! weights used for linear least squares fit.
REAL(KIND=R8), DIMENSION(MIN(N-1, (3*M+1)/2)) :: RES   
               ! the residuals from the linear fit.
REAL(KIND=R8), DIMENSION(3*M+2*MIN(N-1, (3*M+1)/2)) :: WORK  
               ! work array of length LWORK for DGELSS.
REAL(KIND=R8), DIMENSION(MIN(N-1, (3*M+1)/2), M) :: ALPHA, ALPHA_IN 
               ! variables used in DGELSS. 

IER = 0
! Check the validity of M and N.    
IF ( ( M < 1 ) .OR. ( N <= (M + 1) ) ) THEN
   IER = 1
   RETURN
END IF
! Set the number of nodes used in local linear fit.
NP = MIN( N - 1,  (3 * M + 1) / 2 ) 
! Set the RCOND parameter used in DGELSS.
RCOND = NP * EPSILON( 1.0_R8 )
! Set the LWORK parameter used in DGELSS.
LWORK = 3 * M + 2 * NP
! Calculate RW and A. 
DIAM = 0.0_R8
DO K = 1, N
   J = 0
   DO I = 1, N
      IF (I /= K) THEN
         J = J + 1
         DIST(J) = DOT_PRODUCT( X(:, I) - X(:, K), X(:, I) - X(:, K) )
         IDIST(J) = I
      END IF
   END DO
   DIAM = MAX( DIAM, MAXVAL( DIST(1:N-1) ) )
   CALL SORT( DIST, IDIST, NP, N - 1 )
   DIST(1:NP) = SQRT( DIST(1:NP) ) 
   BETA_IN(1:NP) = F(IDIST(1:NP)) - F(K)
   DO I = 1, NP
      ALPHA_IN(I, 1:M) = X(1:M, IDIST(I)) - X(1:M, K)
   END DO
   RW(K) = DIST(NP)
   IF ( PRESENT(RLSHEP) ) THEN
      OMEGA(1:NP) = 1.0_R8
   ELSE
      RP = 1.1_R8 * RW(K)
      OMEGA(1:NP) =  (RP - DIST(1:NP)) / (RP * DIST(1:NP)) 
   END IF
   BETA(1:NP) = OMEGA(1:NP) * BETA_IN(1:NP)
   DO I = 1, NP
      ALPHA(I, 1:M) = OMEGA(I) * ALPHA_IN(I, 1:M)
   END DO
   CALL DGELSS( NP, M, 1, ALPHA(1:NP, 1:M), NP, BETA(1:NP), NP, JPVT, &
        RCOND, RANK, WORK, LWORK, INFO )
   A(1:M, K) = BETA(1:M)
   IF ( RANK < M )  IER = 2
   ROBUST: IF ( PRESENT(RLSHEP) ) THEN
      RES(1:NP) = MATMUL( ALPHA_IN(1:NP, 1:M), BETA(1:M) ) - BETA_IN(1:NP)
      MAR = MAD( RES(1:NP), NP )
      IF ( MAR > EPSILON( 1.0_R8 ) ) THEN
         FUNC = 1
         ITER = 5
         CALL IRLS(M, N, MAR, ITER, FUNC, OMEGA, BETA, BETA_IN, ALPHA,&
              ALPHA_IN, RES, NP, IER1 )
         IF ( IER1 /= 0 ) THEN
            IER = 3
            CYCLE
         ELSE
            A(1:M, K) = BETA(1:M)
         END IF
         MAR = MAD( RES(1:NP), NP )
         IF ( MAR > EPSILON( 1.0_R8 ) ) THEN
            FUNC = 2
            ITER = 5
            CALL IRLS(M, N, MAR, ITER, FUNC, OMEGA, BETA, BETA_IN, &
                 ALPHA, ALPHA_IN, RES, NP, IER1 )
            IF ( IER1 /= 0 ) THEN
               IER = 3
               CYCLE
            ELSE
               A(1:M, K) = BETA(1:M)
            END IF
         END IF
      END IF
! Adjust RW(K) according to weights.
      DO I = 1, NP
         IF ( OMEGA(I) < 0.8_R8 ) THEN
            IF ( I == 1 ) THEN
               RW(K) = DIST(I) / 2.0_R8
            ELSE
               RW(K) = ( DIST(I - 1) + DIST(I) ) / 2.0_R8
            END IF
            EXIT
         END IF
      END DO
   END IF ROBUST
END DO
RW(1:N) = MIN( SQRT( DIAM ) / 2.0_R8, RW(1:N) )
RETURN
END SUBROUTINE LSHEP

SUBROUTINE IRLS( M, N,  S, ITER, FUNC, OMEGA, BETA, BETA_IN, ALPHA, &
           ALPHA_IN, RES, NP, IER )
! IRLS returns the coefficients of the linear least squares fit at the current 
! node using robust M-estimation. For Tukey's bisquare influence function, 
! there is an additional convergence test.
!
! Input parameters:
!  M    is the dimension of the problem.
!  N    is the number of points.
!  S    is the scale constant.
!  ITER is the number of iterations to be used in the algorithm.
!  FUNC identifies which IRLS influence function to use here, 1 = Huber,
!       2 = Tukey's bisquare.
!  RES is the residuals from the linear fit.
!  BETA_IN holds distances.
!  ALPHA, ALPHA_IN are used for DGELSS calls.
! Output parameters:
!  OMEGA stores the final weights.
!  BETA stores the coefficients of the linear fit using robust M-estimation.
!  NP is the size of OMEGA and BETA.
!  IER 
!   = 0, if no errors were encountered.
!   = 1, if any least squares problem is rank deficient.
!   = 2, if the bisquare IRLS did not converge.
!   = 3, if the input parameters are invalid.
USE REAL_PRECISION
IMPLICIT NONE
REAL(KIND=R8), INTENT(INOUT) :: S 
INTEGER, INTENT(IN) :: M, N, ITER, FUNC, NP 
REAL(KIND=R8), INTENT(OUT), DIMENSION(:) :: OMEGA, BETA
REAL(KIND=R8), INTENT(IN), DIMENSION(:) :: BETA_IN
REAL(KIND=R8), INTENT(INOUT), DIMENSION(:) :: RES
REAL(KIND=R8), INTENT(IN), DIMENSION(:,:) :: ALPHA_IN
REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: ALPHA
INTEGER, INTENT(OUT) :: IER

! Local variables.  
INTEGER :: I,J     ! loop control variables.
INTEGER:: INFO   ! error code from DGELSS.
INTEGER :: LWORK ! size of work array for DGELSS.
INTEGER :: RANK  ! rank of array (used for DGELSS).
REAL(KIND=R8) :: RCOND ! reciprocal condition number used in DGELSS.
REAL(KIND=R8), DIMENSION(M) :: JPVT ! array used for DGELSS.
REAL(KIND=R8) :: ERR_NEW      ! objective function value after bisquare IRLS.
REAL(KIND=R8) :: ERR_OLD      ! objective function value before bisquare IRLS.
REAL(KIND=R8), DIMENSION(3*M+2*MIN(N-1, (5*M+1)/2)) :: WORK 
REAL(KIND=R8), PARAMETER, DIMENSION(2) :: TUNING = (/1.0_R8, 3.0_R8/)       
               ! tuning parameters for influence functions.
IF ((FUNC /= 1) .AND. (FUNC /= 2)) THEN
    IER=3
    RETURN
END IF
S = TUNING(FUNC) * S
! Set the RCOND parameter used in DGELSS.
RCOND = NP * EPSILON( 1.0_R8 )
! Set the LWORK parameter used in DGELSS.
LWORK = 3 * M + 2 * NP 
IF ( FUNC ==  2 ) THEN
   ERR_OLD = SUM( 1.0_R8 - (MAX( 0.0_R8, 1.0_R8 - (RES(1:NP) /  S)**2 ))**3 )
END IF
DO I = 1, ITER
! Update the weight from the previous linear fit residual.
   SELECT CASE ( FUNC )
! Huber function
   CASE (1)
      DO J = 1, NP
         IF ( RES(J) == 0.0_R8 ) THEN
            OMEGA(J) = 1.0_R8
         ELSE
            OMEGA(J) = SQRT( MIN( 1.0_R8, S / ABS( RES(J) ) ) )
         END IF
      END DO
! Bisquare function
   CASE(2)
      OMEGA(1:NP) = MAX( 0.0_R8, 1.0_R8 - (RES(1:NP) / S)**2 )
   END SELECT
! Compute the linear fit with a new weight.
   BETA(1:NP) = OMEGA(1:NP) * BETA_IN(1:NP)
   DO J = 1, NP
      ALPHA(J, 1:M) = OMEGA(J) * ALPHA_IN(J, 1:M)
   END DO
   CALL DGELSS( NP, M, 1, ALPHA(1:NP,1:M), NP, BETA(1:NP), NP, JPVT, & 
        RCOND, RANK, WORK, LWORK, INFO )
   IF ( RANK < M ) THEN
      IER = 1
      RETURN
   END IF
   RES(1:NP) = MATMUL( ALPHA_IN(1:NP, 1:M), BETA(1:M) ) - BETA_IN(1:NP)
END DO
IF ( FUNC ==  2 ) THEN
   ERR_NEW = SUM( 1.0_R8 - (MAX( 0.0_R8, 1.0_R8 - (RES(1:NP) / S)**2 ))**3 )
   IF ( ERR_NEW > ERR_OLD ) THEN
      IER = 2
      RETURN
   END IF
END IF
IER = 0
RETURN
END SUBROUTINE IRLS

FUNCTION MAD( RES, NP )
! The function MAD computes the median absolute deviation (MAD) scaled 
! estimate for robust M-estimation from the residual array RES.
USE REAL_PRECISION
IMPLICIT NONE
REAL(KIND=R8), DIMENSION(:), INTENT(IN) :: RES
INTEGER, INTENT(IN)::NP
REAL(KIND=R8) :: MAD

! Local variables.  
INTEGER :: NMID
REAL(KIND=R8) :: FACTOR, MED
  
FACTOR = 1.0_R8 / 0.6745_R8
NMID = NP / 2
IF ( MOD(NP, 2) == 0 ) THEN
   MED = (QUANTILE( RES(:), NMID ) + QUANTILE( RES(:), NMID + 1 )) / 2.0_R8
   MAD = FACTOR * (QUANTILE( ABS( RES(:) - MED ), NMID ) + & 
        QUANTILE( ABS( RES(:) - MED ), NMID + 1 )) / 2.0_R8
ELSE
   MED = QUANTILE( RES(:), NMID + 1 )
   MAD = FACTOR * QUANTILE( ABS( RES(:) - MED ), NMID + 1 )
END IF
RETURN
END FUNCTION MAD

RECURSIVE FUNCTION QUANTILE( A, K ) RESULT( VALUE )
! Recursive function QUANTILE returns the K-th smallest element in array A.
USE REAL_PRECISION
IMPLICIT NONE
REAL(KIND=R8), DIMENSION (:), INTENT(IN) :: A
INTEGER, INTENT(IN) :: K
REAL(KIND=R8) :: VALUE

! Local variables.
INTEGER :: J
REAL(KIND=R8) :: AK

AK = A(K)
J = COUNT( A(:) < AK )
IF ( J >= K ) THEN
   VALUE = QUANTILE( PACK( A(:), A(:) < AK ), K )
ELSE
   J = COUNT( A(:) > AK ) + K - SIZE( A(:) )
   IF ( J > 0 ) THEN
      VALUE = QUANTILE( PACK( A(:), A(:) > AK ), J )
   ELSE
      VALUE = AK
   END IF
END IF
RETURN
END FUNCTION QUANTILE

SUBROUTINE SORT( DIST, IDIST, NUM, LENGTH )
! The subroutine SORT sorts the real array DIST of length LENGTH in ascending 
! order for the smallest NUM elements. IDIST stores the original label for 
! each element in array DIST. 
! Local variables.
USE REAL_PRECISION
IMPLICIT NONE
REAL(KIND=R8), INTENT(INOUT), DIMENSION(:)::DIST
INTEGER, INTENT(INOUT), DIMENSION(:)::IDIST
INTEGER :: I, ITEMP, J, LENGTH, NUM
REAL(KIND=R8) :: TEMP
  
DO I = 2, LENGTH
   DO J = 1, MIN( I-1, NUM )
      IF ( DIST(I) < DIST(J) ) THEN
         TEMP = DIST(I)
         ITEMP = IDIST(I)
         DIST(J+1:I) = DIST(J:I-1) 
         IDIST(J+1:I) = IDIST(J:I-1)
         DIST(J) = TEMP
         IDIST(J) = ITEMP
         EXIT 
      END IF
   END DO
END DO
END SUBROUTINE SORT  


FUNCTION LSHEPVAL( XP, M, N, X, F, A, RW, IER )
! LSHEPVAL returns the linear Shepard approximation at the point XP, using the 
! local linear approximations computed by LSHEP.
!
! Input parameters:
!  XP is the point at which the linear Shepard interpolant function 
!     approximation function is to be evaluated.
!  M is the dimension of the data.
!  N is the number of interpolation points.
!  X contains the interpolation nodes, by column.
!  F contains the function values of the interpolation nodes.
!  A is an M by N matrix containing the coefficients for linear nodal functions
!    returned by LSHEP.
!  RW contains the radius of influence about each interpolation node returned
!    by LSHEP.
! Output parameter:
!  IER 
!   = 0, normal returns.
!   = 1, if the point XP is outside the radius of influence RW(i) for all 
!        nodes, in which case LSHEPVAL is computed using the original Shepard 
!        algorithm with the M+1 closest points.
!   = 2, if the hyperplane normals of the local approximations with positive
!        weights are significantly different. For a nonlinear underlying 
!        function f(x), e.g., quadratic f(x), very different normals are 
!        typical. For a piecewise linear underlying function f(x), IER = 2 
!        signals a potentially large error in LSHEPVAL, since local 
!        approximations from different facets of f(x) have been used.
USE REAL_PRECISION
IMPLICIT NONE

INTEGER, INTENT(IN)  :: M, N
REAL(KIND=R8), DIMENSION(M), INTENT(IN) :: XP 
REAL(KIND=R8), DIMENSION(M, N), INTENT(IN) :: X
REAL(KIND=R8), DIMENSION(N), INTENT(IN) :: F
REAL(KIND=R8), DIMENSION(M, N), INTENT(IN) :: A
REAL(KIND=R8), DIMENSION(N), INTENT(IN) :: RW
INTEGER, INTENT(OUT):: IER
REAL(KIND=R8) :: LSHEPVAL 

! Local variables.
INTEGER :: I   ! number of nodal functions used in Shepard approximation.
INTEGER :: J   ! temporary integer variable.
INTEGER :: K   ! loop control variable.
INTEGER, DIMENSION(M+1) :: ID ! data indices in the original Shepard scheme.
REAL(KIND=R8) :: DIST ! distance between XP and the interpolation nodes.
REAL(KIND=R8) :: TEMP ! temporary real variable.
REAL(KIND=R8) :: W    ! weight value, depends on D.    
REAL(KIND=R8), DIMENSION(M+1) :: D ! weights in the original Shepard scheme.
REAL(KIND=R8), DIMENSION(M+1) :: SLOPE ! the hyperplane normal of the Shepard
               ! approximation function.
REAL(KIND=R8), DIMENSION(M+1, N) :: SSLOPE ! the hyperplane normals of the 
               ! local approximations with positive weights.  
REAL(KIND=R8), DIMENSION(N) :: SW   ! weight values. 
REAL(KIND=R8), DIMENSION(N) :: SWC  ! nodal (linear) function values.

IER = 0
I = 0
J = 1
D(1:M+1) = 0.0_R8
TEMP = 0.0_R8
DO K = 1, N
   DIST = SQRT( DOT_PRODUCT( XP(:) - X(:, K), XP(:) - X(:, K) ) )
   IF ( DIST < RW(K) ) THEN
      IF ( RW(K) - DIST == RW(K) ) THEN
         LSHEPVAL = F(K)
         RETURN
      END IF
      I = I + 1
      W = ((RW(K) - DIST) / (RW(K) * DIST))**2
      SW(I) = W
      SWC(I) = DOT_PRODUCT( A(:, K), (XP(:) - X(:, K)) ) + F(K)
      SSLOPE(1:M, I) = A(1:M, K)
      SSLOPE(M+1, I) = -1.0_R8
   ELSE
      W = 1.0_R8 / DIST
      IF ( W > TEMP ) THEN
         D(J) = W
         ID(J) = K
         J = MINLOC( D(1:M+1), DIM = 1 )
         TEMP = D(J)
      END IF
   END IF
END DO
! I = 0 iff the point XP is not within the radius RW(K) of X(:,K) for all K; 
! IER = 1.
IF ( I == 0 ) THEN
   D = D**2
   LSHEPVAL = DOT_PRODUCT( D / SUM( D ), F(ID) )
   IER = 1 
   RETURN
ELSE
   SW(1:I) = SW(1:I) / SUM( SW(1:I) )
   LSHEPVAL = DOT_PRODUCT( SW(1:I), SWC(1:I) )
! Return IER = 2 iff the angle between two local approximation hyperplanes is 
! too large.
   SLOPE(1:M+1) = MATMUL( SSLOPE(1:M+1, 1:I), SW(1:I) )
   SLOPE(:) = SLOPE(:) / SQRT( DOT_PRODUCT( SLOPE(:), SLOPE(:) ) )
   DO K = 1, I
      SSLOPE(:, K) = SSLOPE(:, K) / SQRT( DOT_PRODUCT( SSLOPE(:, K), &
           SSLOPE(:, K) ) )
      IF ( DOT_PRODUCT( SLOPE(1:M+1), SSLOPE(1:M+1, K) ) < 0.9_R8 ) THEN
         IER = 2
         EXIT
      END IF
   END DO
END IF
RETURN
END FUNCTION LSHEPVAL

SUBROUTINE RIPPLE( M, N, X, F, A, RW, IER )
! The subroutine RIPPLE, residual initiated polynomial-time piecewise linear
! estimation, computes a set of parameters defining a function that 
! interpolates N data values F(i) at scattered nodes X(i) in M dimensions. 
! The interpolant may be evaluated at an arbitrary point by the function 
! LSHEPVAL.
! The interpolation scheme is a modified linear Shepard method.
! The time complexity of the algorithm RIPPLE is between the robust 
! M-estimation and LMS estimation.
!
! Input parameters:  
!   M is the dimension of the data.
!   N is the number of nodes.
!   X(M, N) contains the coordinates of the nodes.
!   F(N) contains the function values of the nodes.
! Output parameters:
!   A(M, N) contains the coefficients of the local fits.
!   RW(N) is an array containing the radius of influence for each point.
!   IER 
!   = 0, if no errors were encountered.
!   = 1, if N is too small relative to M.
!   = 2, if the IRLS subroutine returns an error.
USE REAL_PRECISION
IMPLICIT NONE
    
INTEGER, INTENT(IN) :: M, N
REAL(KIND=R8), DIMENSION(M,N), INTENT(IN) :: X
REAL(KIND=R8), DIMENSION(N), INTENT(IN) :: F
REAL(KIND=R8), DIMENSION(M,N), INTENT(OUT) :: A
REAL(KIND=R8), DIMENSION(N), INTENT(OUT) :: RW
INTEGER, INTENT(OUT) :: IER

! Local variables
INTEGER :: FUNC ! IRLS influence function switch, 1 = Huber, 2 = Tukey's 
                ! bisquare.
INTEGER :: IER1 ! error indicator for IRLS subroutine.
INTEGER :: INFO ! error indicator for DGELSS.
INTEGER :: ITER ! number of iterations for IRLS algorithm.
INTEGER :: I, J, K, L  ! loop control variables.
INTEGER :: LWORK ! length of the work array WORK for DGELSS.
INTEGER :: NP    ! number of nodes used for the local fit.
INTEGER :: RANK  ! rank of array (used for DGELSS).
INTEGER, DIMENSION(N-1) :: IDIST    ! index array for DIST.
INTEGER, DIMENSION(N, MIN(N-1, (3*M+5)/2)) :: S   
                    ! set of indices close to X(:, K).
INTEGER, DIMENSION(MIN(N-1, (3*M+1)/2), M+3) :: D 
                    ! distance matrix of indices.
REAL(KIND=R8) :: DIAM  ! diameter of the whole date set.
REAL(KIND=R8) :: LSE   ! the least squares error.
REAL(KIND=R8) :: MAR   ! the median absolute residual.
REAL(KIND=R8) :: RCOND ! reciprocal condition number used in DGELSS.
REAL(KIND=R8), DIMENSION(M) :: JPVT ! array used for DGELSS.
REAL(KIND=R8), DIMENSION(MIN(N-1, (5*M+1)/2)) :: BETA, BETA_IN
               ! variables used for DGELSS.
REAL(KIND=R8), DIMENSION(MIN(N-1, (5*M+1)/2)) :: OMEGA 
               ! weights used for least squares fit.  
REAL(KIND=R8), DIMENSION(MIN(N-1, (5*M+1)/2)) :: RES   
               ! the residuals from the linear fit.
REAL(KIND=R8), DIMENSION(N-1) :: DIST  ! array that stores all distances.
REAL(KIND=R8), DIMENSION(3*M+2*MIN(N-1, (5*M+1)/2)) :: WORK 
               ! work array for DGELSS.
REAL(KIND=R8), DIMENSION(MIN(N-1, (5*M+1)/2), M):: ALPHA, ALPHA_IN  
               ! variables used for DGELSS.

IER = 0
! Check the validity of M and N.    
IF ( ( M < 1 ) .OR. ( N <= (M + 3) ) ) THEN
   IER = 1
   RETURN
END IF
! Compute the set of indices S.
DIAM = 0.0_R8
NP = MIN( N - 1, (3 * M + 5) / 2 )
DO K = 1, N
   J = 0
   DO I = 1, N
      IF (I /= K) THEN
         J = J + 1
         DIST(J) = DOT_PRODUCT( X(:, I) - X(:, K), X(:, I) - X(:, K) )
         IDIST(J) = I
      END IF
   END DO 
   DIAM = MAX( DIAM, MAXVAL( DIST(1:N-1) ) )
   CALL SORT( DIST, IDIST, NP, N - 1 )
   S(K, 1:NP) = IDIST(1:NP)
END DO
! Compute the local linear approximation.
DO K = 1, N
   NP = MIN( N - 1,  (3 * M + 1) / 2 ) 
   D(1:NP, 1) = S(K, 1:NP)
   DO J = 2, M + 3
      DO I = 1, NP 
         DO L = 1, J
            IF ( ALL( S(D(I, J-1), L) /= D(I, 1:J-1) ) .AND. &
                 S(D(I, J-1), L) /= K )  D(I, J) = S(D(I, J-1), L)
         END DO
      END DO
   END DO
   LSE = HUGE(1.0_R8)
! Set the number of nodes used in searching for the set R.
   DO J = 1, NP
      BETA_IN(1:M+3) = F(D(J, 1:M+3)) - F(K)
      DO I = 1, M + 3
         ALPHA_IN(I, 1:M) = X(1:M, D(J, I)) - X(1:M, K)
      END DO
! Eliminate 2 of the M+3 points.
      RCOND = (M + 1) * EPSILON( 1.0_R8 )
      LWORK = 5 * M + 2
      DO I = 3, M + 3
         DO L = 2, I - 1
            ALPHA(1:L-1, 1:M) = ALPHA_IN(1:L-1, 1:M)
            ALPHA(L:I-2, 1:M) = ALPHA_IN(L+1:I-1, 1:M)
            ALPHA(I-1:M+1, 1:M) = ALPHA_IN(I+1:M+3, 1:M)
            BETA(1:L-1) = BETA_IN(1:L-1)
            BETA(L:I-2) = BETA_IN(L+1:I-1)
            BETA(I-1:M+1) = BETA_IN(I+1:M+3)
            CALL DGELSS( M+1, M, 1, ALPHA(1:M+1,1:M), M+1, BETA(1:M+1), M+1, & 
                 JPVT, RCOND, RANK, WORK, LWORK, INFO ) 
            IF ( ABS( BETA(M+1) ) < LSE )  THEN
               LSE = ABS( BETA(M+1) )
               A(1:M, K) = BETA(1:M)   
               IDIST(1:L-1) = D(J, 1:L-1)
               IDIST(L:I-2) = D(J, L+1:I-1)
               IDIST(I-1:M+1) = D(J, I+1:M+3)
            END IF
         END DO
      END DO
   END DO
! Compute the scale estimate and distance vector from the optimal set R.
   BETA_IN(1:M+1) = F(IDIST(1:M+1)) - F(K)
   DO I = 1, M + 1
      ALPHA_IN(I, 1:M) = X(1:M, IDIST(I)) - X(1:M, K)
      DIST(I) = DOT_PRODUCT( ALPHA_IN(I, 1:M), ALPHA_IN(I, 1:M) )
   END DO
   RES(1:M+1) = MATMUL( ALPHA_IN(1:M+1, 1:M), A(1:M, K) ) - BETA_IN(1:M+1)
   MAR = MAX( MAD( RES(1:M+1),NP ), 1.0E-6_R8 )
! Determine the set T for the final local approx.
   BETA_IN(1:NP) = F(S(K, 1:NP)) - F(K)
   DO I = 1, NP
      ALPHA_IN(I, 1:M) = X(1:M, S(K, I)) - X(1:M, K)
   END DO
   RW(K) = DOT_PRODUCT( ALPHA_IN(NP, :), ALPHA_IN(NP, :) ) 
   CALL SORT( DIST, IDIST, M + 1, M + 1 )
   DO I = 1, M + 1
      IF ( DIST(I) > RW(K) .OR. ( DIST(I) == RW(K) .AND. ALL( IDIST(I) & 
           /= S(K, 1: MIN( N - 1, (3 * M + 1) / 2 )) ) ) ) THEN
         NP = NP + 1
         BETA_IN(NP) = F(IDIST(I)) - F(K)
         ALPHA_IN(NP, 1:M) = X(1:M, IDIST(I)) - X(1:M, K)
      END IF
   END DO
   RW(K) = SQRT( RW(K) )
   RES(1:NP) = MATMUL( ALPHA_IN(1:NP, 1:M), A(1:M, K) ) - BETA_IN(1:NP)
   FUNC = 1
   ITER = 5
   CALL IRLS(M, N, MAR, ITER, FUNC, OMEGA, BETA, BETA_IN, ALPHA, ALPHA_IN, &
        RES, NP, IER1 )
   IF ( IER1  /= 0 ) THEN
      IER = 2
      CYCLE
   ELSE
      A(1:M, K) = BETA(1:M)
   END IF
   FUNC = 2
   ITER = 5
   CALL IRLS(M, N, MAR, ITER, FUNC, OMEGA, BETA, BETA_IN, ALPHA, ALPHA_IN, &
        RES, NP, IER1 )
   IF ( IER1 /= 0 ) THEN
      IER = 2
      CYCLE
   ELSE 
      A(1:M, K) = BETA(1:M)
   END IF
! Adjust RW(K) according to the weights.
   DO I = 1, NP
      IF ( OMEGA(I) < 0.8_R8 ) THEN
         IF ( I == 1 ) THEN
            RW(K) = SQRT( DOT_PRODUCT( ALPHA_IN(I, :), ALPHA_IN(I, :) ) ) &
                 / 2.0_R8
         ELSE
            RW(K) = (SQRT( DOT_PRODUCT( ALPHA_IN(I-1, :), ALPHA_IN(I-1, :) ) )&
                 + SQRT( DOT_PRODUCT( ALPHA_IN(I, :), ALPHA_IN(I, :) ) )) & 
                 / 2.0_R8
         END IF
         EXIT
      END IF
   END DO
END DO
RW(1:N) = MIN( SQRT( DIAM ) / 2.0_R8, RW(1:N) )    
RETURN
END SUBROUTINE RIPPLE

END MODULE LINEAR_SHEPARD_MOD
