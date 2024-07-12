      SUBROUTINE pchev(n, x, f, d, nval, xval, fval, dval, ierr)
      USE stel_kinds
      IMPLICIT NONE
!- ----------------------------------------------
!   D u m m y   A r g u m e n t s
!- ----------------------------------------------
      INTEGER :: n, nval, ierr
      REAL(rprec), DIMENSION(n), INTENT(in) :: x, f, d
      REAL(rprec), DIMENSION(nval), INTENT(in)  :: xval
      REAL(rprec), DIMENSION(nval), INTENT(out) :: fval, dval
!- ----------------------------------------------
!   L o c a l   V a r i a b l e s
!- ----------------------------------------------
      INTEGER :: incfd
      LOGICAL :: skip
!- ----------------------------------------------
!*  **BEGIN PROLOGUE  PCHEV
!*  **DATE WRITTEN   870828   (YYMMDD)
!*  **REVISION DATE  870828   (YYMMDD)
!*  **CATEGORY NO.  E3,H1
!*  **KEYWORDS  CUBIC HERMITE OR SPLINE DIFFERENTIATION,CUBIC HERMITE
!             EVALUATION,EASY TO USE SPLINE OR CUBIC HERMITE EVALUATOR
!*  **AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             ROOM A161, TECHNOLOGY BUILDING
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!*  **PURPOSE  Evaluates the function and first derivative of a piecewise
!            cubic Hermite or spline function at an array of points XVAL,
!            easy to use.
!*  **DESCRIPTION
!
!          PCHEV:  Piecewise Cubic Hermite or Spline Derivative Evaluator,
!                  Easy to Use.
!
!     From the book "Numerical Methods and Software"
!          by  D. Kahaner, C.  Moler, S. Nash
!                 Prentice HALL 1988
!
!     Evaluates the function and first derivative of the cubic Hermite
!     or spline function defined by  N, X, F, D, at the array of points
!
!     This is an easy to use driver for the routines by F.N. Fritsch
!     described in reference (2) below. Those also have other capabilities.
!
! ----------------------------------------------------------------------
!
!  Calling sequence: CALL  PCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
!
!     INTEGER N, NVAL, IERR
!     REAL X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error RETURN IF N.LT.2 .)
!
!     X -- (input) REAL array of independent variable values.  The
!           elements of X must be strictly increasing:
!             X(I-1) .LT. X(I),  I = 2(1)N. (Error RETURN IF not.)
!
!     F -- (input) REAL array of FUNCTION values.  F(I) is
!           the value corresponding to X(I).
!
!     D -- (input) REAL array of derivative values.  D(I) is
!           the value corresponding to X(I).
!
!  NVAL -- (input) number of points at which the functions are to be
!           evaluated. ( error return if nval.lt.1 )
!
!  XVAL -- (input) real array of points at which the functions are to
!           be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XVAL are increasing relative to X;
!              that is,   XVAL(J) .GE. X(I)
!              implies    XVAL(K) .GE. X(I),  ALL K.GE.J .
!           2. If any of the XVAL are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!  FVAL -- (output) REAL array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XVAL.
!
!  DVAL -- (output) REAL array of values of the first derivative of
!           the same function at the points  XVAL.
!
!  IERR -- (output) error flag.
!           Normal RETURN:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  IF N.LT.2 .
!              IERR = -3  IF the X-array is not strictly increasing.
!              IERR = -4  IF NVAL.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT*! been validated.
!              IERR = -5  IF an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY*! If it does.
!
! ----------------------------------------------------------------------
!*  **REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC INTERPOLATION,'
!                   SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980), 238-246.
!                2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION PACKAGE, FINAL SPECIFICATIONS',
!                   LAWRENCE LIVERMORE NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194, AUGUST 1982.
!*  **ROUTINES CALLED  PCHFD
!*  **END PROLOGUE  PCHEV
!
!  DECLARE LOCAL VARIABLES.
!
      skip = .true.;  incfd = 1
!
!
!***FIRST EXECUTABLE STATEMENT  PCHEV
!
      CALL pchfd (n, x, f, d, incfd, skip, nval, xval, fval,
     1            dval, ierr)

      END SUBROUTINE pchev
