      SUBROUTINE pchfd(n, x, f, d, incfd, skip, ne, xe, fe,
     1                 de, ierr)
      USE stel_kinds
      IMPLICIT NONE
!- ----------------------------------------------
!   D u m m y   A r g u m e n t s
!- ----------------------------------------------
      INTEGER :: n, incfd, ne, ierr
      LOGICAL :: skip
      REAL(rprec), DIMENSION(n), INTENT(in) :: x
      REAL(rprec), DIMENSION(incfd,n), INTENT(in) :: f, d
      REAL(rprec), DIMENSION(ne), INTENT(in) :: xe
      REAL(rprec), DIMENSION(ne), INTENT(out) :: fe, de
!- ----------------------------------------------
!   L o c a l   V a r i a b l e s
!- ----------------------------------------------
      INTEGER :: i, ierc, ir, j, jfirst
      INTEGER :: next(2)
      INTEGER :: nj
!- ----------------------------------------------
!*  **BEGIN PROLOGUE  PCHFD
!*  **DATE WRITTEN   811020   (YYMMDD)
!*  **REVISION DATE  870707   (YYMMDD)
!*  **CATEGORY NO.  E3,H1
!*  **KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(PCHFD-S DPCHFD-D),
!             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
!             HERMITE INTERPOLATION,PIECEWISE CUBIC EVALUATION
!*  **AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!*  **PURPOSE  Evaluate a piecewise cubic Hermite function and its first
!            derivative at an array of points.  May be used by itself
!            for Hermite interpolation, or as an evaluator for PCHIM
!            or PCHI!.   If only function values are required, use
!            PCHFE instead.
!*  **DESCRIPTION
!
!          PCHFD:  Piecewise Cubic Hermite Function and Derivative
!                  evaluator
!
!     Evaluates the cubic Hermite FUNCTION defined by  N, X, F, D,  to-
!     gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, Use PCHFE, instead.
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
C ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        REAL X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
!        LOGICAL  SKIP
!
!        CALL  PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error RETURN IF N.LT.2 .)
!
!     X -- (input) REAL array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error RETURN IF not.)
!
!     F -- (input) REAL array of FUNCTION values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) REAL array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error RETURN IF  INCFD.LT.1 .)
!
!     SKIP -- (input/output) LOGICAL variable which should be set to
!           .TRUE. IF the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will SAVE time in CASE these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on normal RETURN.
!
!     NE -- (input) number of evaluation points.  (Error RETURN IF
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient If the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  ALL K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) REAL array of values of the cubic Hermite FUNCTION
!           defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) REAL array of values of the first derivative of
!           the same FUNCTION at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal RETURN:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  IF N.LT.2 .
!              IERR = -2  IF INCFD.LT.1 .
!              IERR = -3  IF the X-array is not strictly increasing.
!              IERR = -4  IF NE.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT*! been validated.
!              IERR = -5  IF an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY*! if it does.
!
!*  **REFERENCES  (NONE)
!*  **ROUTINES CALLED  CHFDV
!*  **END PROLOGUE  PCHFD
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-03   Minor cosmetic changes for release 1.
!     87-07-07   Minor cosmetic changes to prologue.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     1. To produce a Double precision version, simply:
!        a. Change PCHFD to DPCHFD, and CHFDV to DCHFDV, wherever they
!           occur,
!        b. Change the real declaration to Double precision,
!
!     2. Most of the coding between the call to CHFDV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. CHFDv does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of PCHFD that assumes a strict-
!        ly decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which i have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points .LT.
!        X(N-1), followed by points .GT.X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.
!
!  DECLARE ARGUMENTS.
!
!
!  DECLARE LOCAL VARIABLES.
!
!
!  VALIDITY-CHECK ARGUMENTS.
!
!*  **FIRST EXECUTABLE STATEMENT  PCHFD
      IF (.NOT.SKIP) THEN
!
         IF (N < 2) GO TO 5001
         IF (INCFD < 1) GO TO 5002
         DO I = 2, N
            IF (X(I) <= X(I-1)) GO TO 5003
         END DO
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
      ENDIF
      IF (NE < 1) GO TO 5004
      IERR = 0
      SKIP = .TRUE.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
!
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
      IF (JFIRST > NE) GO TO 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
      DO J = JFIRST, NE
         IF (XE(J) >= X(IR)) GO TO 30
      END DO
      J = NE + 1
      GO TO 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30 CONTINUE
      IF (IR == N) J = NE + 1
!
   40 CONTINUE
      NJ = J - JFIRST
!
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
      IF (NJ .NE. 0) THEN
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
         CALL CHFDV (X(IR-1), X(IR), F(1,IR-1), F(1,IR),
     1      D(1,IR-1), D(1,IR), NJ, XE(JFIRST), FE(JFIRST),
     2      DE(JFIRST), NEXT, IERC)
!       ----------------------------------------------------------------
         IF (IERC < 0) GO TO 5005
!
         IF (NEXT(2) .NE. 0) THEN
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
            IF (IR .GE. N) THEN
!           IF (IR .EQ. N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
            ELSE
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
!           ENDIF
!        ENDIF
            ENDIF
         ENDIF
!
         IF (NEXT(1) .NE. 0) THEN
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
            IF (IR .LE. 2) THEN
!           IF (IR .EQ. 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
            ELSE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO I = JFIRST, J - 1
                  IF (XE(I) < X(IR-1)) GO TO 45
               END DO
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN CHFDV.
               GO TO 5005
!
   45          CONTINUE
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO I = 1, IR - 1
                  IF (XE(J) < X(I)) EXIT
               END DO
!              AT THIS POINT, EITHER  XE(J) .LT. X(1)
!                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
               IR = MAX(1,I - 1)
!           ENDIF
!        ENDIF
            ENDIF
         ENDIF
!
         JFIRST = J
!
!     END OF IR-LOOP.
!
      ENDIF
      IR = IR + 1
      IF (IR <= N) GO TO 10
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      STOP 'PCHFD -- NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      STOP 'PCHFD -- INCREMENT LESS THAN ONE'
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      STOP 'PCHFD -- X-ARRAY NOT STRICTLY INCREASING'
      RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
      IERR = -4
      STOP 'PCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
      RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM CHFDV.
!   **! THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      STOP 'PCHFD -- ERROR RETURN FROM CHFDV -- FATAL'
      RETURN

      END SUBROUTINE PCHFD
