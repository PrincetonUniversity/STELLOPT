      SUBROUTINE PCHEZ(N, X, F, D, SPLINE, WK, LWK, IERR)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N, LWK, IERR
      LOGICAL SPLINE
      REAL(rprec), DIMENSION(N) :: X, F, D
      REAL(rprec), DIMENSION(LWK) :: WK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(2) :: IC
      INTEGER :: INCFD
      REAL(rprec), DIMENSION(2) :: VC
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHEZ
!***DATE WRITTEN   870821   (YYMMDD)
!***REVISION DATE  870908   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  CUBIC HERMITE MONOTONE INTERPOLATION, SPLINE
!             INTERPOLATION, EASY TO USE PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  KAHANER, D.K., (NBS)
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             GAITHERSBURG, MARYLAND 20899
!             (301) 975-3808
!***PURPOSE  Easy to USE spline or cubic Hermite interpolation.
!***DESCRIPTION
!
!          PCHEZ:  Piecewise Cubic Interpolation, Easy to Use.
!
!     From the book "Numerical Methods and Software"
!          by  D. Kahaner, C. Moler, S. Nash
!               Prentice HALL 1988
!
!     Sets derivatives for spline (two continuous derivatives) or
!     Hermite cubic (one continuous derivative) interpolation.
!     Spline interpolation is smoother, but may not "look" right IF the
!     data CONTAINS both "steep" and "flat" sections.  Hermite cubics
!     can produce a "visually pleasing" and monotone interpolant to
!     monotone data. This is an easy to USE driver for the routines
!     by F. N. Fritsch in reference (4) below. Various boundary
!     conditions are set to DEFAULT values by PCHEZ. MANY other choices
!     are available in the SUBROUTINEs PCHIC, PCHIM and PCHSP.
!
!     Use PCHEV to evaluate the resulting function and its derivative.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:   CALL  PCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
!
!     INTEGER  N, IERR,  LWK
!     REAL  X(N), F(N), D(N), WK(*)
!     LOGICAL SPLINE
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error RETURN IF N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) REAL array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error RETURN IF not.)
!
!     F -- (input) REAL array of dependent variable values to be inter-
!           polated.  F(I) is value corresponding to X(I).
!
!     D -- (output) REAL array of derivative values at the data points.
!
!     SPLINE -- (input) LOGICAL variable to specify IF the interpolant
!           is to be a spline with two continuous derivaties
!           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
!           continuous derivative (set SPLINE=.FALSE.).
!        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
!           DEFAULT "not-a-knot" boundary condition, with a continuous
!           third derivative at X(2) and X(N-1). See reference (3).
!              If SPLINE=.FALSE. the interpolating Hermite cubic will be
!           monotone IF the input data is monotone. Boundary conditions
!           computed from the derivative of a local quadratic unless this
!           alters monotonicity.
!
!     WK -- (scratch) REAL work array, which must be declared by the calling
!           program to be at least 2*N IF SPLINE is .TRUE. and not used
!           otherwise.
!
!     LWK -- (input) length of work array WK. (Error RETURN IF
!           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
!
!     IERR -- (output) error flag.
!           Normal RETURN:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  (can ONLY occur when SPLINE=.FALSE.) means that
!                 IERR switches in the direction of monotonicity were detected.
!                 When SPLINE=.FALSE.,  PCHEZ guarantees that IF the input
!                 data is monotone, the interpolant will be too. This warning
!                 is to alert you to the fact that the input data was not
!                 monotone.
!           "Recoverable" errors:
!              IERR = -1  IF N.LT.2 .
!              IERR = -3  IF the X-array is not strictly increasing.
!              IERR = -7  IF LWK is less than 2*N and SPLINE is .TRUE.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!----------------------------------------------------------------------
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC INTERPOLATION,'
!                  SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980), 238-246.
!               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,'
!                  LLNL PREPRINT UCRL-87559 (APRIL 1982).
!               3. CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
!                 VERLAG (NEW YORK, 1978).  (ESP. CHAPTER IV, PP.49-62.)
!               4. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION PACKAGE, FINAL SPECIFICATIONS',
!                  LAWRENCE LIVERMORE NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194, AUGUST 1982.
!***ROUTINES CALLED  PCHIM,PCHSP
!***END PROLOGUE  PCHEZ
!
!  DECLARE LOCAL VARIABLES.
!
      ic = 0;  incfd = 1
!
!***FIRST EXECUTABLE STATEMENT  PCHEZ
!
      IF (SPLINE) THEN
         CALL PCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
         CALL PCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
!
!  ERROR CONDITIONS ALREADY CHECKED IN PCHSP OR PCHIM

      RETURN
!------------- LAST LINE OF PCHEZ FOLLOWS ------------------------------
      END SUBROUTINE PCHEZ
