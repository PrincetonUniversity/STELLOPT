!============================================================================
      MODULE math_utilities
!--------------------------------------------------------------------------
! 
! FUNCTION: Module containing routines necessary to compute the line-integrated
!           Faraday rotation angle and the line-integrated polarization phase shift.
!
!------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE


      CONTAINS
!=============================================================================
      REAL(KIND=rprec) FUNCTION MAGNITUDE(vec)
!--------------------------------------------------------------------------
! FUNCTION: calculates the magnitude of 3D input vector.  
!
!---------------------------------------------------------------------------      
      USE stel_kinds
      IMPLICIT NONE


      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            vec           ! 3D vector

      MAGNITUDE = SQRT( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )

      RETURN
      END FUNCTION MAGNITUDE

!=============================================================================
      SUBROUTINE UNIT_VECTOR(v1, v1_unit)
!--------------------------------------------------------------------------
! FUNCTION: calculates the magnitude of 3D input vector.  
!
!---------------------------------------------------------------------------      
      USE stel_kinds
      IMPLICIT NONE

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            v1           ! 3D vector

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            v1_unit        ! 3D vector

      v1_unit = v1 / MAGNITUDE(v1)

      RETURN
      END SUBROUTINE UNIT_VECTOR


!=============================================================================
      REAL(KIND=rprec) FUNCTION DIST(v1, v2)
!--------------------------------------------------------------------------
! FUNCTION: calculates the absolute distance between 2 input vectors.  
!
!---------------------------------------------------------------------------      
      USE stel_kinds
      IMPLICIT NONE

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            v1, v2           ! 3D vectors

      REAL(rprec), DIMENSION(3) ::                                             &
     &            v3           ! difference vector btwn v1 & v2

      v3 = ABS( v2 - v1 )
      DIST = SQRT( v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3) )

      RETURN
      END FUNCTION DIST

!==================================================================================
      SUBROUTINE CROSS_PRODUCT(A, B, C)
!----------------------------------------------------------------------------------
!  FUNCTION:  computes the cross product (C) for input CARTESIAN vectors A and B.
!
!  LOGIC:  C = A .cross. B
!
!  created by J. Shields 2/15/06
!
!--------------------------------------------------------------------------------------
      USE stel_kinds
      IMPLICIT NONE

!........passed variables....................................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) ::  A   ! input Cartesian vector

      REAL(rprec), DIMENSION(3), INTENT(IN) ::  B   ! input Cartesian vector

      REAL(rprec), DIMENSION(3), INTENT(OUT) :: C   ! Cross product:  A .cross. B

      
!.................local variables...................................!
      INTEGER(iprec) :: i
      character(len=*), PARAMETER ::  subname = 'CROSS_PRODUCT'


!..................X component.....................................!
      C(1) = A(2) * B(3)  -  A(3) * B(2)  

!...................Y component.....................................!
      C(2) = A(3) * B(1)  -  A(1) * B(3) 
 
!...................Z component.....................................!
      C(3) = A(1) * B(2)  -  A(2) * B(1)  

      RETURN
      END SUBROUTINE CROSS_PRODUCT


      END MODULE math_utilities
