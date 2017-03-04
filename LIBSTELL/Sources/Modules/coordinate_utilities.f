!*******************************************************************************
!  File coordinate_utilities.f
!
!  Module is part of the LIBSTELL. This modules containes code to convert from
!  different coordinate systems.
!
!*******************************************************************************

      MODULE coordinate_utilities
      USE stel_kinds
      USE stel_constants

      PUBLIC :: cart_to_cyl, cyl_to_cart, cood_utils_test
      PRIVATE :: check

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Convert from cartesian coordinates to cylindical coordinates.
!*******************************************************************************
      FUNCTION cart_to_cyl(cart)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3), INTENT(in) :: cart
      REAL(rprec), DIMENSION(3)             :: cart_to_cyl

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      cart_to_cyl(1) = SQRT(DOT_PRODUCT(cart(1:2),cart(1:2)))
      cart_to_cyl(2) = ATAN2(cart(2), cart(1))
      cart_to_cyl(3) = cart(3)

      END FUNCTION

!*******************************************************************************
!  Convert from cylindical coordinates to cartesian coordinates.
!
!  NOTE: Phi must be in radians.
!
!*******************************************************************************
      FUNCTION cyl_to_cart(cyl)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3), INTENT(in) :: cyl
      REAL(rprec), DIMENSION(3)             :: cyl_to_cart

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      cyl_to_cart(1) = cyl(1)*COS(cyl(2))
      cyl_to_cart(2) = cyl(1)*SIN(cyl(2))
      cyl_to_cart(3) = cyl(3)

      END FUNCTION

!*******************************************************************************
!  Convert vector from cylindical coordinates to cartesian coordinates.
!
!  NOTE: Phi must be in radians.
!
!*******************************************************************************
      FUNCTION cyl_to_cart_vec(cyl, vec)
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(3), INTENT(in) :: cyl
      REAL(rprec), DIMENSION(3), INTENT(in) :: vec
      REAL(rprec), DIMENSION(3)             :: cyl_to_cart_vec

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Transformation matrix is
!
! [ v_x ]   [ cos[phi] -sin[phi] 0 ]   [  v_r  ]
! [ v_y ] = [ sin[phi]  cos[phi] 0 ] * [ v_phi ]
! [ v_z ]   [ 0         0        1 ]   [  v_z  ]

      cyl_to_cart_vec(1) = COS(cyl(2))*vec(1) - SIN(cyl(2))*vec(2)
      cyl_to_cart_vec(2) = SIN(cyl(2))*vec(1) + COS(cyl(2))*vec(2)
      cyl_to_cart_vec(3) = vec(3)

      END FUNCTION

!*******************************************************************************
!  Unit Test Functions
!*******************************************************************************
!*******************************************************************************
!  Main test function
!*******************************************************************************
      FUNCTION cood_utils_test()
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: cood_utils_test
      REAL(rprec), DIMENSION(3) :: result

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Test cart_to_cyl function. cart(0,0,1) = cyl(0,0,1)
      result = cart_to_cyl((/ 0.0d+0, 0.0d+0, 1.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 1, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 2, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(3), 3, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN

! Test cart_to_cyl function. cart(-1,0,0) = cyl(1,Pi,0)
      result = cart_to_cyl((/ -1.0d+0, 0.0d+0, 0.0d+0 /))
      cood_utils_test = check(1.0d+0, result(1), 4, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(pi, result(2), 5, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 6, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN

! Test cart_to_cyl function. cart(-1,0,-1) = cyl(1,Pi/2,-1)
      result = cart_to_cyl((/ 0.0d+0, 1.0d+0, -1.0d+0 /))
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(1), 7, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(pi/2.0d+0, result(2), 8, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(-1.0d+0, result(3), 9, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN

! Test cyl_to_cart function. cart(0,0,1) = cyl(0,0,1)
      result = cyl_to_cart((/ 0.0d+0, 0.0d+0, 1.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 1, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 2, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(3), 3, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN

! Test cyl_to_cart function. cart(-1,0,0) = cyl(1,Pi,0)
      result = cyl_to_cart((/ 1.0d+0, pi, 0.0d+0 /))
      cood_utils_test = check(-1.0d+0, result(1), 4, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0*SIN(pi), result(2), 5, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 6, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN

! Test cyl_to_cart function. cart(0,1,-1) = cyl(1,Pi/2,-1)
      result = cyl_to_cart((/ 1.0d+0, pi/2.0d+0, -1.0d+0 /))
      cood_utils_test =                                                        &
     &   check(COS(pi/2.0d+0), result(1), 7, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(2), 8, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(-1.0d+0, result(3), 9, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN

! Test cyl_to_cart_vec function. cart(1,0,0) = cyl(1,0,0) @ phi=0
      result = cyl_to_cart_vec((/ 1.0d+0, 0.0d+0, 0.0d+0 /),                   &
     &                         (/ 1.0d+0, 0.0d+0, 0.0d+0 /))
      cood_utils_test = check(1.0d+0, result(1), 1, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 2, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 3, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN

! Test cyl_to_cart_vec function. cart(-1,0,0) = cyl(0,1,0) @ phi=0
      result = cyl_to_cart_vec((/ 1.0d+0, 0.0d+0, 0.0d+0 /),                   &
     &                         (/ 0.0d+0, 1.0d+0, 0.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 4, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(2), 5, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 6, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN

! Test cyl_to_cart_vec function. cart(0,0,1) = cyl(0,0,1) @ phi=Pi/2
      result = cyl_to_cart_vec((/ 1.0d+0, pi/2.0d+0, 0.0d+0 /),                &
     &                         (/ 0.0d+0, 0.0d+0, 1.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 7, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 8, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(3), 9, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN

      END FUNCTION

!*******************************************************************************
!  Check Test result
!*******************************************************************************
      FUNCTION check(expected, recieved, testNum, name)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: check
      REAL(rprec), INTENT(in) :: expected, recieved
      INTEGER, INTENT(in) :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      check = expected .eq. recieved
      IF (.not.check) THEN
         write(*,*) "coordinate_utilities.f: ", name, " test", testNum,        &
     &              "failed."
         write(*,*) "Expected", expected, "Recieved", recieved
      END IF

      END FUNCTION

      END MODULE
