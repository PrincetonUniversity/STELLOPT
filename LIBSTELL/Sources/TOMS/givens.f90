  SUBROUTINE GIVENS (A, B, C, S)
!
! GIVENS constructs a Givens plane rotation.
!
! The transformation has the form of a 2 by 2 matrix G(C,S):
! (  C  S)
! (- S  C)
!
! where C*C + S*S = 1, which zeroes the second entry of the
! the column vector (A, B) when C and S are properly chosen.
! A call to GIVENS is normally followed by a call to ROTATE
! which computes the product of G(C,S) with a 2 by N matrix.
!
! Author:
!   Robert Renka
!   University of North Texas
!
! Modified by:
!   Manjula Iyer
!   Virginia Polytechnic Institute and State University
! 
! Modified:   
!   15 November 2005  
!
! Parameters:
!
!   Input/output, real A, B.
!
!   On input, A and B define the 2-vector whose second entry (B) is
!   to be annihilated by a Givens rotation.
!
!   On output, A has been overwritten by a value
!     R = +/- SQRT ( A*A + B*B )
!     and B has been overwritten by a value Z which allows C
!     and S to be recovered as:
!
!       if | Z | <= 1, then
!         C = SQRT (1 - Z*Z), 
!         S = Z
!       else if | Z | > 1 then
!         C = 1 / Z, 
!         S = SQRT (1 - C*C).
!
!     Output, real C, S, the components of the Givens transformation, 
!     which may be computed by:
!       C = +/- A / SQRT (A*A + B*B)
!       S = +/- B / SQRT (A*A + B*B)
!
! Local parameters:
!
!   R = C*A + S*B = +/-SQRT(A*A+B*B)
!   U,V = variables used to scale A and B for computing R
!
!   ABS(A) > ABS(B)
!
!   Note that R has the sign of A, C > 0, and S has
!   SIGN(A)*SIGN(B).
!
    IMPLICIT NONE
    DOUBLE PRECISION :: A, B, C, R, S, U, V

    IF (ABS (A) > ABS (B)) THEN
      U = 2.0E+00 * A
      V = B / U
      R = SQRT ( 0.25E+00 + V * V ) * U
      C = A / R
      S = 2.0E+00 * V * C
      B = S
      A = R
!
! ABS(A) <= ABS(B)
! Store R in A.
! Note that R has the sign of B, S > 0, and C has SIGN(A)*SIGN(B).
!
      ELSE IF (B /= 0.0E+00) THEN
        U = 2.0E+00 * B
        V = A / U
        A = SQRT (0.25E+00 + V * V) * U
        S = B / A
        C = 2.0E+00 * V * S
        IF (C /= 0.0E+00) THEN
          B = 1.0E+00 / C
        ELSE
          B = 1.0E+00
        END IF
!
! A = B = 0.
!
      ELSE
        C = 1.0E+00
        S = 0.0E+00
    END IF
    RETURN
  END
