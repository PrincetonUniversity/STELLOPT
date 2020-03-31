      FUNCTION PCHDF (K, X, S, IERR)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: K, IERR
      REAL(rprec), DIMENSION(K) :: X, S
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J
      REAL(rprec) :: VALUE, ZERO, PCHDF
!-----------------------------------------------
!***BEGIN PROLOGUE  PCHDF
!***REFER TO  PCHCE,PCHSP
!***DESCRIPTION
!
!          PCHDF:   PCHIP Finite Difference Formula
!
!     Uses a divided difference formulation to compute a K-point approx-
!     imation to the derivative at X(K) based on the data in X and S.
!
!     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary
!     derivative approximations.
!
! ----------------------------------------------------------------------
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error RETURN IF not).
!        X      CONTAINS the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      CONTAINS the ASSOCIATED slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need ONLY be of length K-1.)
!
!     On RETURN:
!        S      will be destroyed.
!        IERR   will be set to -1 IF K.LT.2 .
!        PCHDF  will be set to the desired derivative approximation IF
!               IERR=0 or to zero IF IERR=-1.
!
! ----------------------------------------------------------------------
!
!  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
!              Verlag (New York, 1978), pp. 10-16.
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
!        a. Change PCHDF to DPCHDF wherever it occurs,
!        b. Change the Real declarations to Double precision, and
!        c. Change the constant ZERO to Double precision.
!***END PROLOGUE  PCHDF
!
!  DECLARE LOCAL VARIABLES.
!
      DATA ZERO/0/
!
!  CHECK FOR LEGAL VALUE OF K.
!
!***FIRST EXECUTABLE STATEMENT  PCHDF
      IF (K >= 3) THEN
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
         DO J = 2, K - 1
            S(:K-J) = (S(2:K-J+1)-S(:K-J))/(X(1+J:K)-X(:K-J))
         END DO
!
!  EVALUATE DERIVATIVE AT X(K).
!
         VALUE = S(1)
         DO I = 2, K - 1
            VALUE = S(I) + VALUE*(X(K)-X(I))
         END DO
!
!  NORMAL RETURN.
!
         IERR = 0
         PCHDF = VALUE
         RETURN
!
!  ERROR RETURN.
!
      ENDIF
!     K.LT.3 RETURN.
      IERR = -1
      STOP 'PCHDF -- K LESS THAN THREE'
      PCHDF = ZERO
      RETURN
!------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
      END FUNCTION PCHDF
