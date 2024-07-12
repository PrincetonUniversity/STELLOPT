!*******************************************************************************
!>  @file coordinate_utilities.f
!>  @brief Contains module @ref coordinate_utilities
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module is part of the LIBSTELL. This modules containes code to convert from
!>  different coordinate systems.
!*******************************************************************************

      MODULE coordinate_utilities
      USE stel_kinds
      USE stel_constants

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PUBLIC :: cart_to_cyl, cyl_to_cart, cyl_to_cart_vec,                     &
     &          cood_utils_test
      PRIVATE :: check

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Convert a point from cartes cartesian coordinates to cylindical
!>  coordinates.
!>
!>    r = Sqrt(x*x + y*y)
!>    phi = atan(y/x)
!>    z = z
!>
!>  @param[in] cart Point in cartesian coordinates.
!>  @returns The point in cyclindical coordinates.
!-------------------------------------------------------------------------------
      PURE FUNCTION cart_to_cyl(cart)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cart_to_cyl
      REAL (rprec), DIMENSION(3), INTENT(in) :: cart

!  Start of executable code
      cart_to_cyl(1) = SQRT(DOT_PRODUCT(cart(1:2), cart(1:2)))
      cart_to_cyl(2) = ATAN2(cart(2), cart(1))
      cart_to_cyl(3) = cart(3)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Convert a point from cylindical coordinates to cartesian coordinates.
!>
!>    x = r*cos(phi)
!>    y = r*sin(phi)
!>    z = z
!>
!>  @note Phi must be in radians.
!>
!>  @param[in] cart Point in cylindical coordinates.
!>  @returns The point in cartesian coordinates.
!-------------------------------------------------------------------------------
      PURE FUNCTION cyl_to_cart(cyl)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cyl_to_cart
      REAL (rprec), DIMENSION(3), INTENT(in) :: cyl

!  Start of executable code
      cyl_to_cart(1) = cyl(1)*COS(cyl(2))
      cyl_to_cart(2) = cyl(1)*SIN(cyl(2))
      cyl_to_cart(3) = cyl(3)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Convert vector from cylindical coordinates to cartesian coordinates.
!>
!>    [ v_x ]   [ cos[phi] -sin[phi] 0 ]   [  v_r  ]
!>    [ v_y ] = [ sin[phi]  cos[phi] 0 ] * [ v_phi ]
!>    [ v_z ]   [ 0         0        1 ]   [  v_z  ]
!>
!>  @note Phi must be in radians.
!>
!>  @param[in] cyl Point in cylindical coordinates.
!>  @param[in] vec Vector in cylindical coordinates.
!>  @returns The point in cartesian coordinates.
!-------------------------------------------------------------------------------
      PURE FUNCTION cyl_to_cart_vec(cyl, vec)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cyl_to_cart_vec
      REAL (rprec), DIMENSION(3), INTENT(in) :: cyl
      REAL (rprec), DIMENSION(3), INTENT(in) :: vec

!  Start of executable code
      cyl_to_cart_vec(1) = COS(cyl(2))*vec(1) - SIN(cyl(2))*vec(2)
      cyl_to_cart_vec(2) = SIN(cyl(2))*vec(1) + COS(cyl(2))*vec(2)
      cyl_to_cart_vec(3) = vec(3)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Convert vector from cartesian coordinates to cylindical coordinates.
!>
!>    [  v_r  ]   [  x/r y/r 0 ]   [ v_x ]
!>    [ v_phi ] = [ -y/r x/r 0 ] * [ v_y ]
!>    [  v_z  ]   [  0   0   1 ]   [ v_z ]
!>
!>  Where r = Sqrt(x*x + y*y)
!>
!>  @note Phi must be in radians.
!>
!>  @param[in] cart Point in cartesian coordinates.
!>  @param[in] vec  Point in cartesian coordinates.
!>  @returns The point in cylindical coordinates.
!-------------------------------------------------------------------------------
      PURE FUNCTION cart_to_cyl_vec(cart, vec)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cart_to_cyl_vec
      REAL (rprec), DIMENSION(3), INTENT(in) :: cart
      REAL (rprec), DIMENSION(3), INTENT(in) :: vec

!  local variables
      REAL (rprec)                           :: r

!  Start of executable code
      r = SQRT(DOT_PRODUCT(cart(1:2), cart(1:2)))

      cart_to_cyl_vec(1) = ( cart(1)*vec(1) + cart(2)*vec(2))/r
      cart_to_cyl_vec(2) = (-cart(2)*vec(1) + cart(1)*vec(2))/r
      cart_to_cyl_vec(3) = vec(3)

      END FUNCTION

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Coordinate utilities unit test function.
!>
!>  This runs the associated unit tests and returns the result.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION cood_utils_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                    :: cood_utils_test

!  local variables
      REAL (rprec), DIMENSION(3) :: result

!  Start of executable code
!  Test cart_to_cyl function. cart(0,0,1) = cyl(0,0,1)
      result = cart_to_cyl((/ 0.0d+0, 0.0d+0, 1.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 1, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 2, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(3), 3, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN

!  Test cart_to_cyl function. cart(-1,0,0) = cyl(1,Pi,0)
      result = cart_to_cyl((/ -1.0d+0, 0.0d+0, 0.0d+0 /))
      cood_utils_test = check(1.0d+0, result(1), 4, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(pi, result(2), 5, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 6, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN

!  Test cart_to_cyl function. cart(-1,0,-1) = cyl(1,Pi/2,-1)
      result = cart_to_cyl((/ 0.0d+0, 1.0d+0, -1.0d+0 /))
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(1), 7, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(pi/2.0d+0, result(2), 8, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(-1.0d+0, result(3), 9, "cart_to_cyl")
      IF (.not.cood_utils_test) RETURN

!  Test cyl_to_cart function. cart(0,0,1) = cyl(0,0,1)
      result = cyl_to_cart((/ 0.0d+0, 0.0d+0, 1.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 1, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 2, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(3), 3, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN

!  Test cyl_to_cart function. cart(-1,0,0) = cyl(1,Pi,0)
      result = cyl_to_cart((/ 1.0d+0, pi, 0.0d+0 /))
      cood_utils_test = check(-1.0d+0, result(1), 4, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0*SIN(pi), result(2), 5, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 6, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN

!  Test cyl_to_cart function. cart(0,1,-1) = cyl(1,Pi/2,-1)
      result = cyl_to_cart((/ 1.0d+0, pi/2.0d+0, -1.0d+0 /))
      cood_utils_test =                                                        &
     &   check(COS(pi/2.0d+0), result(1), 7, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(2), 8, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(-1.0d+0, result(3), 9, "cyl_to_cart")
      IF (.not.cood_utils_test) RETURN

!  Test cyl_to_cart_vec function. cart(1,0,0) = cyl(1,0,0) @ phi=0
      result = cyl_to_cart_vec((/ 1.0d+0, 0.0d+0, 0.0d+0 /),                   &
     &                         (/ 1.0d+0, 0.0d+0, 0.0d+0 /))
      cood_utils_test = check(1.0d+0, result(1), 1, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(2), 2, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 3, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN

!  Test cyl_to_cart_vec function. cart(-1,0,0) = cyl(0,1,0) @ phi=0
      result = cyl_to_cart_vec((/ 1.0d+0, 0.0d+0, 0.0d+0 /),                   &
     &                         (/ 0.0d+0, 1.0d+0, 0.0d+0 /))
      cood_utils_test = check(0.0d+0, result(1), 4, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(1.0d+0, result(2), 5, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN
      cood_utils_test = check(0.0d+0, result(3), 6, "cyl_to_cart_vec")
      IF (.not.cood_utils_test) RETURN

!  Test cyl_to_cart_vec function. cart(0,0,1) = cyl(0,0,1) @ phi=Pi/2
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
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check a value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check
      REAL (rprec), INTENT(in)      :: expected
      REAL (rprec), INTENT(in)      :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  local parameters
      REAL(rprec), PARAMETER :: range = 1.0E-15_dp

!  Start of executable code
      check = (expected .eq. received) .or.                                    &
     &        ((expected .lt. received + range) .and.                          &
     &         (expected .gt. received - range))
      IF (.not.check) THEN
         write(*,*) "coordinate_utilities.f: ", name, " test", testNum,        &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

      END MODULE
