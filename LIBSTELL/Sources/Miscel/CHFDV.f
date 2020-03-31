      SUBROUTINE CHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)
      USE stel_kinds
      IMPLICIT NONE
!- ----------------------------------------------
!   D u m m y   A r g u m e n t s
!- ----------------------------------------------
      INTEGER :: NE, IERR
      REAL(rprec), INTENT(IN) :: X1, X2, F1, F2, D1, D2
      INTEGER, INTENT(OUT) :: NEXT(2)
      REAL(rprec), DIMENSION(NE), INTENT(IN) :: XE
      REAL(rprec), DIMENSION(NE), INTENT(OUT) :: FE, DE
!- ----------------------------------------------
!   L o c a l   V a r i a b l e s
!- ----------------------------------------------
      REAL(rprec), PARAMETER :: ZERO = 0
      INTEGER :: I
      REAL(rprec) :: C2,C2T2,C3,C3T3,DEL1,DEL2,DELTA,H,
     1 X,XMI,XMA
!- ----------------------------------------------
!*  **BEGIN PROLOGUE  CHFDV
!*  **DATE WRITTEN   811019   (YYMMDD)
!*  **REVISION DATE  870707   (YYMMDD)
!*  **CATEGORY NO.  E3,H1
!*  **KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(CHFDV-S DCHFDV-D),
!             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
!             CUBIC POLYNOMIAL EVALUATION
!*  **AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!*  **PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
!            first derivative at an array of points.  While designed for
!            use by PCHFD, it may be useful directly as an evaluator for
!            a piecewise cubic Hermite functioN in applications, such as
!            graphing, where the interval is known in advance.
!            If only function values are required, use CHFEV instead.
!*  **DESCRIPTION
!
!        CHFDV:  Cubic Hermite Function and Derivative Evaluator
!
!     Evaluates the cubic polynomial determined by FUNCTION values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
!     its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use CHFEV, instead.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        REAL X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
!
!        CALL  CHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error RETURN IF  X1.EQ.X2 .)
!
!     F1,F2 -- (input) values of FUNCTION at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error RETURN IF
!           NE.LT.1 .)
!
!     XE -- (input) REAL array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) REAL array of values of the cubic function defined
!           by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) REAL array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) INTEGER array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal RETURN:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  IF NE.LT.1 .
!              IERR = -2  IF X1.EQ.X2 .
!                (Output arrays have not been changed in either CASE.)
!
!*  **REFERENCES  (NONE)
!*  **END PROLOGUE  CHFDV
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-08-03   Minor cosmetic changes for release 1.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a Double precision version, simply:
!        a. Change CHFDV to DCHFDV wherever it occurs,
!        b. Change the Real declaration to Double precision,
!        c. Change the constant ZERO to Double precision, and
!        d. Change the names of the Fortran functions:  MAX, MIN.
!
!  DECLARE ARGUMENTS.
!
!
!  VALIDITY-CHECK ARGUMENTS.
!
!*  **FIRST EXECUTABLE STATEMENT  CHFDV
      IF (NE >= 1) THEN
         H = X2 - X1
         IF (H == ZERO) GO TO 5002
!
!  INITIALIZE.
!
         IERR = 0
         NEXT(1) = 0
         NEXT(2) = 0
         XMI = MIN(ZERO,H)
         XMA = MAX(ZERO,H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
         DELTA = (F2 - F1)/H
         DEL1 = (D1 - DELTA)/H
         DEL2 = (D2 - DELTA)/H
!                                           (DELTA IS NO LONGER NEEDED.)
         C2 = -(DEL1 + DEL1 + DEL2)
         C2T2 = C2 + C2
         C3 = (DEL1 + DEL2)/H
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
         C3T3 = C3 + C3 + C3
!
!  EVALUATION LOOP.
!
         DO I = 1, NE
            X = XE(I) - X1
            FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
            DE(I) = D1 + X*(C2T2 + X*C3T3)
!          COUNT EXTRAPOLATION POINTS.
            IF (X < XMI) NEXT(1) = NEXT(1) + 1
            IF (X > XMA) NEXT(2) = NEXT(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
         END DO
!
!  NORMAL RETURN.
!
         RETURN
!
!  ERROR RETURNS.
!
      ENDIF

!     NE.LT.1 RETURN.
      IERR = -1
      STOP 'CHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
      RETURN
!
 5002 CONTINUE
!     X1.EQ.X2 RETURN.
      IERR = -2
      STOP 'CHFDV -- INTERVAL ENDPOINTS EQUAL'

      END SUBROUTINE CHFDV
