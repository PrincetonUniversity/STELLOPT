      SUBROUTINE PCHSP(IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N, INCFD, NWK, IERR
      INTEGER, DIMENSION(2) :: IC
      REAL(rprec), DIMENSION(2) :: VC
      REAL(rprec), DIMENSION(N) :: X
      REAL(rprec), DIMENSION(INCFD,N) :: F, D
      REAL(rprec), DIMENSION(2,N) :: WK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IBEG, IEND, INDEX, J, NM1
      REAL(rprec) :: G, HALF, ONE
      REAL(rprec), DIMENSION(3) :: STEMP
      REAL(rprec) :: THREE, TWO
      REAL(rprec), DIMENSION(4) :: XTEMP
      REAL(rprec) :: ZERO
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PCHDF
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHSP
!***DATE WRITTEN   820503   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(PCHSP-S DPCHSP-D),
!             CUBIC HERMITE INTERPOLATION,PIECEWISE CUBIC INTERPOLATION,
!             SPLINE INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine the Hermite represen-
!            tation of the cubic spline interpolant to given data, with
!            specified boundary conditions.
!***DESCRIPTION
!
!          PCHSP:   Piecewise Cubic Hermite Spline
!
!     Computes the Hermite representation of the cubic spline inter-
!     polant to the data given in X and F satisfying the boundary
!     conditions specified by IC and VC.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite FUNCTION may be evaluated
!     by PCHFE or PCHFD.
!
!     NOTE:  This is a modified version of C. de Boor''S cubic spline
!            routine CUBSPL.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        REAL  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) INTEGER array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at END of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor''S cubic spline routine CUBSPL.
!              < This is the DEFAULT boundary condition. >
!           IBEG = 1  IF first derivative at X(1) is given in VC(1).
!           IBEG = 2  IF second derivative at X(1) is given in VC(1).
!           IBEG = 3  to USE the 3-point difference formula for D(1).
!                     (Reverts to the DEFAULT b.c. IF N.LT.3 .)
!           IBEG = 4  to USE the 4-point difference formula for D(1).
!                     (Reverts to the DEFAULT b.c. IF N.LT.4 .)
!          NOTES:
!           1. An error RETURN is taken IF IBEG is out of range.
!           2. For the "natural" boundary condition, USE IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In CASE IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error RETURN is taken IF IEND is out of range.
!           2. For the "natural" boundary condition, USE IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) REAL array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set ONLY IF IC(1) = 1 or 2 .
!           VC(2) need be set ONLY IF IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error RETURN IF N.LT.2 .)
!
!     X -- (input) REAL array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error RETURN IF not.)
!
!     F -- (input) REAL array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!
!     D -- (output) REAL array of derivative values at the data points.
!           These values will determine the cubic spline interpolant
!           with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error RETURN IF  INCFD.LT.1 .)
!
!     WK -- (scratch) REAL array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error RETURN IF NWK.LT.2*N .)
!
!     IERR -- (output) error flag.
!           Normal RETURN:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  IF N.LT.2 .
!              IERR = -2  IF INCFD.LT.1 .
!              IERR = -3  IF the X-array is not strictly increasing.
!              IERR = -4  IF IBEG.LT.0 or IBEG.GT.4 .
!              IERR = -5  IF IEND.LT.0 of IEND.GT.4 .
!              IERR = -6  IF both of the above are true.
!              IERR = -7  IF NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (the d-array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D-array may have been changed in this case.)
!             (             Do **NOT** USE it!                )
!
!***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
!                 VERLAG (NEW YORK, 1978), PP. 53-59.
!***ROUTINES CALLED  PCHDF
!***END PROLOGUE  PCHSP
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-04   Converted to SLATEC library version.
!     87-07-07   Minor cosmetic changes to prologue.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a Double precision version, simply:
!        a. Change PCHSP to DPCHSP wherever it occurs,
!        b. Change the real declarations to double precision, and
!        c. Change the constants ZERO, HALF, ... to Double precision.
!
!  DECLARE ARGUMENTS.
!
!
!  DECLARE LOCAL VARIABLES.
!
!
      DATA ZERO/0/
      DATA HALF/0.5_DP/
      DATA ONE/1/
      DATA TWO/2/
      DATA THREE/3/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHSP
      IF (N >= 2) THEN
         IF (INCFD < 1) GO TO 5002
         DO J = 2, N
            IF (X(J) <= X(J-1)) GO TO 5003
         END DO
!
         IBEG = IC(1)
         IEND = IC(2)
         IERR = 0
         IF (IBEG<0 .OR. IBEG>4) IERR = IERR - 1
         IF (IEND<0 .OR. IEND>4) IERR = IERR - 2
         IF (IERR < 0) GO TO 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
         IF (NWK < 2*N) GO TO 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
         WK(1,2:N) = X(2:N) - X(:N-1)
         WK(2,2:N) = (F(1,2:N)-F(1,:N-1))/WK(1,2:N)
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
         IF (IBEG > N) IBEG = 0
         IF (IEND > N) IEND = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
         IF (IBEG==1 .OR. IBEG==2) THEN
            D(1,1) = VC(1)
         ELSE IF (IBEG > 2) THEN
!        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
            DO J = 1, IBEG
               INDEX = IBEG - J + 1
!           INDEX RUNS FROM IBEG DOWN TO 1.
               XTEMP(J) = X(INDEX)
               IF (J < IBEG) STEMP(J) = WK(2,INDEX)
            END DO
!                 --------------------------------
            D(1,1) = PCHDF(IBEG,XTEMP,STEMP,IERR)
!                 --------------------------------
            IF (IERR /= 0) GO TO 5009
            IBEG = 1
         ENDIF
!
         IF (IEND==1 .OR. IEND==2) THEN
            D(1,N) = VC(2)
         ELSE IF (IEND > 2) THEN
!        PICK UP LAST IEND POINTS.
            DO J = 1, IEND
               INDEX = N - IEND + J
!           INDEX RUNS FROM N+1-IEND UP TO N.
               XTEMP(J) = X(INDEX)
               IF (J < IEND) STEMP(J) = WK(2,INDEX+1)
            END DO
!                 --------------------------------
            D(1,N) = PCHDF(IEND,XTEMP,STEMP,IERR)
!                 --------------------------------
            IF (IERR /= 0) GO TO 5009
            IEND = 1
         ENDIF
!
! --------------------( BEGIN CODING FROM CUBSPL )--------------------
!
!  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
!  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
         IF (IBEG == 0) THEN
            IF (N == 2) THEN
!           NO CONDITION AT LEFT END AND N = 2.
               WK(2,1) = ONE
               WK(1,1) = ONE
               D(1,1) = TWO*WK(2,2)
            ELSE
!           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
               WK(2,1) = WK(1,3)
               WK(1,1) = WK(1,2) + WK(1,3)
               D(1,1) = ((WK(1,2)+TWO*WK(1,1))*WK(2,2)
     1            *WK(1,3)+WK(1,2)**2*WK(2,3))/WK(1,1)
            ENDIF
         ELSE IF (IBEG == 1) THEN
!        SLOPE PRESCRIBED AT LEFT END.
            WK(2,1) = ONE
            WK(1,1) = ZERO
         ELSE
!        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
            WK(2,1) = TWO
            WK(1,1) = ONE
            D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
         ENDIF
!
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
         NM1 = N - 1
         IF (NM1 > 1) THEN
            DO J = 2, NM1
               IF (WK(2,J-1) == ZERO) GO TO 5008
               G = -WK(1,J+1)/WK(2,J-1)
               D(1,J) = G*D(1,J-1) + THREE*(WK(1,J)
     1            *WK(2,J+1)+WK(1,J+1)*WK(2,J))
               WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J)+WK(1,J+1))
            END DO
         ENDIF
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!     AT THIS POINT.
         IF (IEND /= 1) THEN
!
            IF (IEND == 0) THEN
               IF (N==2 .AND. IBEG==0) THEN
!           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
                  D(1,2) = WK(2,2)
                  GO TO 30
               ELSE IF (N==2 .OR. N==3 .AND. IBEG==0) THEN
!           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!           NOT-A-KNOT AT LEFT END POINT).
                  D(1,N) = TWO*WK(2,N)
                  WK(2,N) = ONE
                  IF (WK(2,N-1) == ZERO) GO TO 5008
                  G = -ONE/WK(2,N-1)
               ELSE
!           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
!           KNOT AT LEFT END POINT.
                  G = WK(1,N-1) + WK(1,N)
!           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
                  D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)
     1            +WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
                  IF (WK(2,N-1) == ZERO) GO TO 5008
                  G = -G/WK(2,N-1)
                  WK(2,N) = WK(1,N-1)
               ENDIF
            ELSE
!        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
               D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
               WK(2,N) = TWO
               IF (WK(2,N-1) == ZERO) GO TO 5008
               G = -ONE/WK(2,N-1)
            ENDIF
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
            WK(2,N) = G*WK(1,N-1) + WK(2,N)
            IF (WK(2,N) == ZERO) GO TO 5008
            D(1,N) = (G*D(1,N-1)+D(1,N))/WK(2,N)
!
!  CARRY OUT BACK SUBSTITUTION
!
         ENDIF
   30    CONTINUE
         DO J = NM1, 1, -1
            IF (WK(2,J) == ZERO) GO TO 5008
            D(1,J) = (D(1,J)-WK(1,J)*D(1,J+1))/WK(2,J)
         END DO
! --------------------(  END  CODING FROM CUBSPL )--------------------
!
!  NORMAL RETURN.
!
         RETURN
!
!  ERROR RETURNS.
!
      ENDIF
!     N.LT.2 RETURN.
      IERR = -1
      STOP 'PCHSP -- NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      STOP 'PCHSP -- INCREMENT LESS THAN ONE'
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      STOP 'PCHSP -- X-ARRAY NOT STRICTLY INCREASING'
      RETURN
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      STOP 'PCHSP -- IC OUT OF RANGE'
      RETURN
!
 5007 CONTINUE
!     NWK TOO SMALL RETURN.
      IERR = -7
      STOP 'PCHSP -- WORK ARRAY TOO SMALL'
      RETURN
!
 5008 CONTINUE
!     SINGULAR SYSTEM.
!   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
!   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
      STOP 'PCHSP -- SINGULAR LINEAR SYSTEM'
      RETURN
!
 5009 CONTINUE
!     ERROR RETURN FROM PCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      STOP 'PCHSP -- ERROR RETURN FROM PCHDF'
      RETURN
!------------- LAST LINE OF PCHSP FOLLOWS ------------------------------
      END SUBROUTINE PCHSP
