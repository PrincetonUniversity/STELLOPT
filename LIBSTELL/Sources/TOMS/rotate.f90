
  SUBROUTINE ROTATE ( N, C, S, X, Y )
!
! ROTATE applies a Givens rotation.
!
! The rotation has the form:
!
! (  C  S)
! (- S  C)
!
! and is essentially applied to a 2 by N matrix:
! (X(1) X(2) ... X(N))
! (Y(1) Y(2) ... Y(N))
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
!   Input, integer N, the dimension of the vectors.
!
!   Input, real C, S, the cosine and sine entries of the Givens
!   rotation matrix. These may be determined by subroutine GIVENS.
!
!   Input/output, real X(N), Y(N), the rotated vectors. 
!
    IMPLICIT NONE

    INTEGER I, N
    DOUBLE PRECISION :: C, S, XI, YI
    DOUBLE PRECISION, DIMENSION(N) :: X, Y
    
    IF (N <= 0) THEN
      RETURN
    ELSE IF (C == 1.0E+00 .AND. S == 0.0E+00) THEN
      RETURN
    END IF

    DO I = 1, N
      XI = X(I)
      YI = Y(I)
      X(I) =   C * XI + S * YI
      Y(I) = - S * XI + C * YI
    END DO
    RETURN
  END
