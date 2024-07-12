      FUNCTION PCHST (ARG1, ARG2)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) :: ARG1, ARG2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: ZERO = 0, ONE = 1
      REAL(rprec) :: PCHST
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHST
!***REFER TO  PCHCE,PCHCI,PCHCS,PCHIM
!***ROUTINES CALLED  (NONE)
!***DESCRIPTION
!
!         PCHST:  PCHIP Sign-Testing Routine.
!
!
!     Returns:
!        -1. IF ARG1 and ARG2 are of opposite sign.
!         0. IF either argument is zero.
!        +1. IF ARG1 and ARG2 are of the same sign.
!
!     The object is to DO this without multiplying ARG1*ARG2, to avoid
!     possible over/underflow problems.
!
!  Fortran INTRINSICs used:  SIGN.
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                  Mathematics and Statistics Division,
!                  Lawrence Livermore National Laboratory.
!
!  Change record:
!     82-08-05   Converted to SLATEC library version.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     To produce a Double precision version, simply:
!        a. Change PCHST to DPCHST wherever it occurs,
!        b. Change all references to the fortran intrinsics to their
!           Double presision equivalents,
!        c. Change the Real declarations to Double precision, and
!        d. Change the constants  ZERO  and  ONE  to Double precision.
!***END PROLOGUE  PCHST
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = SIGN(ONE,ARG1)*SIGN(ONE,ARG2)
      IF (ARG1==ZERO .OR. ARG2==ZERO) PCHST = ZERO
!
      RETURN
!------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END FUNCTION PCHST
