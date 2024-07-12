      SUBROUTINE PCHIM(N, X, F, D, INCFD, IERR)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: N, INCFD, IERR
      REAL(rprec), DIMENSION(N) :: X
      REAL(rprec), DIMENSION(INCFD,N) :: F, D
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: ZERO = 0, THREE = 3
      INTEGER :: I, NLESS1
      REAL(rprec) :: R1, DEL1, DEL2, D_max, D_min
      REAL(rprec) :: DRAT1, DRAT2, DSAVE, H1, H2
      REAL(rprec) :: HSUM, HSUMT3, W1, W2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec) , EXTERNAL :: PCHST
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHIM
!***DATE WRITTEN   811103   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!             TYPE=SINGLE PRECISION(PCHIM-S DPCHIM-D),
!             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
!             PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!             MATHEMATICS AND STATISTICS DIVISION
!             LAWRENCE LIVERMORE NATIONAL LABORATORY
!             P.O. BOX 808  (L-316)
!             LIVERMORE, CA  94550
!             FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine a monotone piecewise
!            cubic Hermite interpolant to given data.  Boundary values
!            are provided which are compatible with monotonicity.  The
!            interpolant will have an extremum at each point where mono-
!            tonicity switches direction.  (See PCHIC if user control is
!            desired over boundary or switch conditions.)
!***DESCRIPTION
!
!          PCHIM:  Piecewise Cubic Hermite Interpolation to
!                  Monotone data.
!
!     Sets derivatives needed to determine a monotone piecewise cubi!
!     Hermite interpolant to the data given in X and F.
!
!     Default boundary conditions are provided which are compatible
!     with monotonicity.  (See PCHIC IF user control of boundary con-
!     ditions is desired.)
!
!     If the data are ONLY piecewise monotonic, the interpolant will
!     have an extremum at each point where monotonicity switches direc-
!     tion.  (See PCHIC IF user control is desired in such cases.)
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite FUNCTION may be evaluated
!     by PCHFE or PCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  PCHIM (N, X, F, D, INCFD, IERR)
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
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!           PCHIM is designed for monotonic data, but it will work for
!           ANY F-array.  It will force extrema at points WHERE mono-
!           tonicity switches direction.  If some other treatment of
!           switch points is desired, PCHIC should be used instead.
!                                     -----
!     D -- (output) REAL array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error RETURN IF  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal RETURN:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  IF N.LT.2 .
!              IERR = -2  IF INCFD.LT.1 .
!              IERR = -3  IF the X-array is not strictly increasing.
!             (The D-array has not been changed in ANY of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC INTERPOLATION,'
!                 SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980), 238-246.
!              2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,'
!                 LLNL PREPRINT UCRL-87559 (APRIL 1982).
!***ROUTINES CALLED  PCHST
!***END PROLOGUE  PCHIM
!
! ----------------------------------------------------------------------
!
!  Change record:
!     82-02-01   1. Introduced  PCHST  to reduce possible over/under-
!                   flow problems.
!                2. Rearranged derivative formula for same reason.
!     82-06-02   1. Modified END conditions to be continuous functions
!                   of data when monotonicity switches in next interval.
!                2. Modified formulas so END conditions are less prone
!                   of over/underflow problems.
!     82-08-03   Minor cosmetic changes for release 1.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 IF they are of opposite sign.
!     2. To produce a Double precision version, simply:
!        a. Change PCHIM to DPCHIM wherever it occurs,
!        b. Change PCHST to DPCHST wherever it occurs,
!        c. Change all references to the fortran intrinsics to their
!           Double precision equivalents,
!        d. Change the Real declarations to Double precision, and
!        e. Change the constants ZERO and THREE to Double precision.
!
!  DECLARE ARGUMENTS.
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHIM
      IF (N >= 2) THEN
         IF (INCFD < 1) GO TO 5002
         DO I = 2, N
            IF (X(I) <= X(I-1)) GO TO 5003
         END DO
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
         IERR = 0
         NLESS1 = N - 1
         H1 = X(2) - X(1)
         DEL1 = (F(1,2)-F(1,1))/H1
         DSAVE = DEL1
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
         IF (NLESS1 <= 1) THEN
            D(1,1) = DEL1
            D(1,N) = DEL1
         ELSE
!
!  NORMAL CASE  (N .GE. 3).
!
            H2 = X(3) - X(2)
            DEL2 = (F(1,3)-F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
            HSUM = H1 + H2
            W1 = (H1 + HSUM)/HSUM
            W2 = -H1/HSUM
            D(1,1) = W1*DEL1 + W2*DEL2
            IF (PCHST(D(1,1),DEL1) <= ZERO) THEN
               D(1,1) = ZERO
            ELSE IF (PCHST(DEL1,DEL2) < ZERO) THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
               D_max = THREE*DEL1
               IF (ABS(D(1,1)) > ABS(D_max)) D(1,1) = D_max
            ENDIF
!
!  LOOP THROUGH INTERIOR POINTS.
!
            DO I = 2, NLESS1
               IF (I /= 2) THEN
!
                  H1 = H2
                  H2 = X(I+1) - X(I)
                  HSUM = H1 + H2
                  DEL1 = DEL2
                  DEL2 = (F(1,I+1)-F(1,I))/H2
               ENDIF
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
               D(1,I) = ZERO
               R1 = PCHST(DEL1,DEL2)
               IF (R1 >= ZERO) THEN
                  IF (R1 > ZERO) GO TO 45
                  IF (DEL2 == ZERO) CYCLE
                  IF (PCHST(DSAVE,DEL2) < ZERO) IERR = IERR + 1
                  DSAVE = DEL2
                  CYCLE
!
               ENDIF
               IERR = IERR + 1
               DSAVE = DEL2
               CYCLE
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45          CONTINUE
               HSUMT3 = HSUM + HSUM + HSUM
               W1 = (HSUM + H1)/HSUMT3
               W2 = (HSUM + H2)/HSUMT3
               D_max = MAX(ABS(DEL1),ABS(DEL2))
               D_min = MIN(ABS(DEL1),ABS(DEL2))
               DRAT1 = DEL1/D_max
               DRAT2 = DEL2/D_max
               D(1,I) = D_min/(W1*DRAT1 + W2*DRAT2)
!
            END DO
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
            W1 = -H2/HSUM
            W2 = (H2 + HSUM)/HSUM
            D(1,N) = W1*DEL1 + W2*DEL2
            IF (PCHST(D(1,N),DEL2) <= ZERO) THEN
               D(1,N) = ZERO
            ELSE IF (PCHST(DEL1,DEL2) < ZERO) THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
               D_max = THREE*DEL2
               IF (ABS(D(1,N)) > ABS(D_max)) D(1,N) = D_max
            ENDIF
!
!  NORMAL RETURN.
!
         ENDIF
         RETURN
!
!  ERROR RETURNS.
!
      ENDIF
!     N.LT.2 RETURN.
      IERR = -1
      STOP 'PCHIM -- NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      STOP 'PCHIM -- INCREMENT LESS THAN ONE'
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      STOP 'PCHIM -- X-ARRAY NOT STRICTLY INCREASING'
      RETURN
!------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
      END SUBROUTINE PCHIM
