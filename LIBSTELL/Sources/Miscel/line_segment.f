!*******************************************************************************
!  File line_segment.f
!
!  Module is part of the libstell.a. This module containes code to create a
!  profile constructed of line sigments. These line segments are assumed to be
!  specified such that xx(i) < xx(i + 1)
!
!*******************************************************************************

      MODULE line_segment
      USE stel_kinds
      USE stel_constants

      PUBLIC line_seg, line_seg_int, line_seg_test
      PRIVATE get_indices, y_value, y_value_int, slope, offset, check

      CONTAINS

!*******************************************************************************
!  Public Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Find the linearly interpolated value of y at position x from line sigments
!  specified by yy(xx).
!*******************************************************************************
      SUBROUTINE line_seg(x, y, xx, yy, n)
      IMPLICIT none

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: n

      REAL(rprec), INTENT(out) :: y
      REAL(rprec), INTENT(in)  :: x
      REAL(rprec), DIMENSION(n), INTENT(in)  :: xx, yy
!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: ilow, ihigh

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Check for valid line segments.
      IF (n .le. 1) THEN
         STOP "Line sigments require at least two points"
      ELSE IF (xx(1) .ge. xx(2)) THEN
         STOP "Line sigments must be specified in increasing order of x"
      END IF

! If x is outside of array, linearly interprolate from the last two points.
      IF (x .lt. xx(1)) THEN
         ilow = 1; ihigh = 2
      ELSE IF (x .gt. xx(n)) THEN
         ilow = n - 1; ihigh = n
      ELSE
         CALL get_indices(x, xx, 1, n, ilow, ihigh)
      END IF

      y = y_value(x, yy, xx, ilow, ihigh)

      END SUBROUTINE

!*******************************************************************************
!  Find the linearly integrated value of y at position x from line sigments
!  specified by yy(xx).
!*******************************************************************************
      SUBROUTINE line_seg_int(x, y, xx, yy, n)
      IMPLICIT none

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: n

      REAL(rprec), INTENT(out) :: y
      REAL(rprec), INTENT(in)  :: x
      REAL(rprec), DIMENSION(n), INTENT(in)  :: xx, yy
!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: ilow, ihigh, i0low, i0high, i

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Check for valid line segments.
      IF (n .le. 1) THEN
         STOP "Line sigments require at least two points"
      ELSE IF (xx(1) .ge. xx(2)) THEN
         STOP "Line sigments must be specified in increasing order of x"
      END IF

! If x is outside of array, linearly interprolate from the last two points.
      IF (x .lt. xx(1)) THEN
         ilow = 1; ihigh = 2
      ELSE IF (x .gt. xx(n)) THEN
         ilow = n - 1; ihigh = n
      ELSE
         CALL get_indices(x, xx, 1, n, ilow, ihigh)
      END IF

! Determine where zero is.
      IF (zero .lt. xx(1)) THEN
         i0low = 1; i0high = 2
      ELSE IF (zero .gt. xx(n)) THEN
         i0low = n - 1; i0high = n
      ELSE
         CALL get_indices(zero, xx, 1, n, i0low, i0high)
      END IF

! If 0 and x have the same low and high indices, then the integrate directly
! from 0 to x. Other wise Integrate from:
!   0 to xx(i0high); xx(i0high) to xx(ilow); xx(ilow) to x
      IF (ilow .eq. i0low) THEN
         y = y_value_int(zero, x, yy, xx, ilow, ihigh)
      ELSE
         y = y_value_int(zero, xx(i0high), yy, xx, i0low, i0high)
         DO i = i0high, ilow - 1
             y = y + y_value_int(xx(i), xx(i+1), yy, xx, i, i+1)
         END DO
         y = y + y_value_int(xx(ilow), x, yy, xx, ilow, ihigh)
      END IF

      END SUBROUTINE

!*******************************************************************************
!  Private Functions and Subroutines.
!*******************************************************************************
!*******************************************************************************
!  Find bounding indices of the xx array
!*******************************************************************************
      PURE RECURSIVE                                                           &
     &SUBROUTINE get_indices(x, xx, lBound, uBound, ilow, ihigh)
      IMPLICIT none

!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), INTENT(in)  :: x
      REAL(rprec), DIMENSION(:), INTENT(in)  :: xx
      INTEGER, INTENT(in) :: lBound, uBound ! Lower and upper bounds of the xx array
      INTEGER, INTENT(out) :: ilow, ihigh
!-------------------------------------------------------------------------------
!  Local Variable declarations.
!-------------------------------------------------------------------------------
      INTEGER :: hBound ! Half way index of the xx array.

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Perform a tree search to find the interval that x falls between. This
! algorthim works by dividing the array into two sub arrays then determining if
! x is in the lower or upper array. Then the algortrim recurses on the sub array
! until the size of the array equals 2. This assumes that the array xx is
! correctly sorted.
      IF (uBound - lBound .eq. 1) THEN
         ilow = lBound; ihigh = uBound
         RETURN
      END IF

      hBound = (uBound + lBound)/2
      IF (x .ge. xx(lBound) .and. x .lt. xx(hBound)) THEN
         CALL get_indices(x, xx, lBound, hBound, ilow, ihigh)
      ELSE
         CALL get_indices(x, xx, hBound, uBound, ilow, ihigh)
      END IF
      END SUBROUTINE

!*******************************************************************************
!  Find m of y = m*x + b
!*******************************************************************************
      PURE FUNCTION slope(yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: slope
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      slope = (yy(ihigh) - yy(ilow))/(xx(ihigh) - xx(ilow))

      END FUNCTION

!*******************************************************************************
!  Find b of y = m*x + b
!*******************************************************************************
      PURE FUNCTION offset(yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: offset
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      offset = (xx(ihigh)*yy(ilow) - xx(ilow)*yy(ihigh))                       &
     &       / (xx(ihigh) - xx(ilow))

      END FUNCTION

!*******************************************************************************
!  Find y(x)
!*******************************************************************************
      PURE FUNCTION y_value(x, yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      REAL(rprec), INTENT(in) :: x
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: y_value
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      y_value = slope(yy, xx, ilow, ihigh)*x                                   &
     &        + offset(yy, xx, ilow, ihigh)

      END FUNCTION

!*******************************************************************************
!  Find Int[y(x),x]
!*******************************************************************************
      PURE FUNCTION y_value_int(x0, x1, yy, xx, ilow, ihigh)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), INTENT(in) :: xx, yy
      REAL(rprec), INTENT(in) :: x0, x1
      INTEGER, INTENT(in) :: ilow, ihigh
      REAL(rprec) :: y_value_int
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      y_value_int = slope(yy, xx, ilow, ihigh)/2.0*(x1**2.0 - x0**2.0)         &
     &            + offset(yy, xx, ilow, ihigh)*(x1 - x0)

      END FUNCTION

!*******************************************************************************
!  Unit Test Functions
!*******************************************************************************
!*******************************************************************************
!  Main test function
!*******************************************************************************
      FUNCTION line_seg_test()
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: line_seg_test
      REAL(rprec) :: result

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
! Test slope function. yy(0.0,1.0) and xx(0.0,1.0) should have a slope of 1.0
      result = slope((/ 0.0d+0, 1.0d+0 /),                                     &
     &               (/ 0.0d+0, 1.0d+0 /),                                     &
     &               1, 2)
      line_seg_test = check(1.0d+0, result, 1, "slope()")
      IF (.not.line_seg_test) RETURN
! Test slope function. yy(1.0,0.0) and xx(0.0,1.0) should have a slope of -1.0
      result = slope((/ 1.0d+0, 0.0d+0 /),                                     &
     &               (/ 0.0d+0, 1.0d+0 /),                                     &
     &               1, 2)
      line_seg_test = check(-1.0d+0, result, 2, "slope()")
      IF (.not.line_seg_test) RETURN
! Test slope function. yy(1.0,1.0) and xx(0.0,1.0) should have a slope of 0.0
      result = slope((/ 1.0d+0, 1.0d+0 /),                                     &
     &               (/ 0.0d+0, 1.0d+0 /),                                     &
     &               1, 2)
      line_seg_test = check(0.0d+0, result, 3, "slope()")
      IF (.not.line_seg_test) RETURN

! Test slope function. yy(0.0,1.0) and xx(1.0,2.0) should have a slope of 1.0
      result = slope((/ 0.0d+0, 1.0d+0 /),                                     &
     &               (/ 1.0d+0, 2.0d+0 /),                                     &
     &               1, 2)
      line_seg_test = check(1.0d+0, result, 4, "slope()")

! Test offset function. yy(0.0,1.0) and xx(0.0,1.0) should have an offset of 0.0
      result = offset((/ 0.0d+0, 1.0d+0 /),                                    &
     &                (/ 0.0d+0, 1.0d+0 /),                                    &
     &                1, 2)
      line_seg_test = check(0.0d+0, result, 1, "offset()")
      IF (.not.line_seg_test) RETURN
! Test offset function. yy(1.0,0.0) and xx(0.0,1.0) should have an offset of 1.0
      result = offset((/ 1.0d+0, 0.0d+0 /),                                    &
     &                (/ 0.0d+0, 1.0d+0 /),                                    &
     &                1, 2)
      line_seg_test = check(1.0d+0, result, 2, "offset()")
      IF (.not.line_seg_test) RETURN
! Test offset function. yy(2.0,1.0) and xx(1.0,2.0) should have an offset of 3.0
      result = offset((/ 2.0d+0, 1.0d+0 /),                                    &
     &                (/ 1.0d+0, 2.0d+0 /),                                    &
     &                1, 2)
      line_seg_test = check(3.0d+0, result, 3, "offset()")
      IF (.not.line_seg_test) RETURN

! Test y_value function. yy(0.0,2.0) and xx(0.0,2.0) should have a value of 1.0
! at x = 1.0
      result = y_value(1.0d+0,                                                 &
     &                 (/ 0.0d+0, 2.0d+0 /),                                   &
     &                 (/ 0.0d+0, 2.0d+0 /),                                   &
     &                 1, 2)
      line_seg_test = check(1.0d+0, result, 1, "y_value()")
      IF (.not.line_seg_test) RETURN
! Test y_value function. yy(1.0,2.0) and xx(1.0,2.0) should have a value of 0.0
! at x = 0.0
      result = y_value(0.0d+0,                                                 &
     &                 (/ 1.0d+0, 2.0d+0 /),                                   &
     &                 (/ 1.0d+0, 2.0d+0 /),                                   &
     &                 1, 2)
      line_seg_test = check(0.0d+0, result, 2, "y_value()")
      IF (.not.line_seg_test) RETURN
! Test y_value function. yy(0.0,1.0) and xx(0.0,1.0) should have a value of 3.0
! at x = 3.0
      result = y_value(3.0d+0,                                                 &
     &                 (/ 0.0d+0, 1.0d+0 /),                                   &
     &                 (/ 0.0d+0, 1.0d+0 /),                                   &
     &                 1, 2)
      line_seg_test = check(3.0d+0, result, 3, "y_value()")
      IF (.not.line_seg_test) RETURN

! Test y_value_int function. yy(1.0,1.0) and xx(0.0,1.0) should have a value of
! 1.0 at x1 = 1.0 to x0 = 0.0
      result = y_value_int(0.0d+0, 1.0d+0,                                     &
     &                    (/ 1.0d+0, 1.0d+0 /),                                &
     &                    (/ 0.0d+0, 1.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(1.0d+0, result, 1, "y_value_int()")
      IF (.not.line_seg_test) RETURN

! Test y_value_int function. yy(1.0,1.0) and xx(1.0,2.0) should have a value of
! 1.0 at x1 = 1.0 to x0 = 0.0
      result = y_value_int(0.0d+0, 1.0d+0,                                     &
     &                    (/ 1.0d+0, 1.0d+0 /),                                &
     &                    (/ 1.0d+0, 2.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(1.0d+0, result, 2, "y_value_int()")
      IF (.not.line_seg_test) RETURN
! Test y_value_int function. yy(1.0,1.0) and xx(0.0,1.0) should have a value of
! 1.0 at x1 = 2.0 to x0 = 1.0
      result = y_value_int(1.0d+0, 2.0d+0,                                     &
     &                    (/ 1.0d+0, 1.0d+0 /),                                &
     &                    (/ 0.0d+0, 1.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(1.0d+0, result, 3, "y_value_int()")
      IF (.not.line_seg_test) RETURN
! Test y_value_int function. yy(0.0,1.0) and xx(0.0,1.0) should have a value of
! 0.5 at x1 = 1.0 to x0 = 0.0
      result = y_value_int(0.0d+0, 1.0d+0,                                     &
     &                    (/ 0.0d+0, 1.0d+0 /),                                &
     &                    (/ 0.0d+0, 1.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(0.5d+0, result, 4, "y_value_int()")
      IF (.not.line_seg_test) RETURN
! Test y_value_int function. yy(1.0,2.0) and xx(0.0,1.0) should have a value of
! 1.5 at x1 = 1.0 to x0 = 0.0
      result = y_value_int(0.0d+0, 1.0d+0,                                     &
     &                    (/ 1.0d+0, 2.0d+0 /),                                &
     &                    (/ 0.0d+0, 1.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(1.5d+0, result, 5, "y_value_int()")
      IF (.not.line_seg_test) RETURN
! Test y_value_int function. yy(1.0,2.0) and xx(1.0,2.0) should have a value of
! 0.5 at x1 = 1.0 to x0 = 0.0
      result = y_value_int(0.0d+0, 1.0d+0,                                     &
     &                    (/ 1.0d+0, 2.0d+0 /),                                &
     &                    (/ 1.0d+0, 2.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(0.5d+0, result, 6, "y_value_int()")
      IF (.not.line_seg_test) RETURN
! Test y_value_int function. yy(0.0,1.0) and xx(0.0,1.0) should have a value of
! 1.5 at x1 = 2.0 to x0 = 1.0
      result = y_value_int(1.0d+0, 2.0d+0,                                     &
     &                    (/ 1.0d+0, 2.0d+0 /),                                &
     &                    (/ 1.0d+0, 2.0d+0 /),                                &
     &                    1, 2)
      line_seg_test = check(1.5d+0, result, 7, "y_value_int()")
      IF (.not.line_seg_test) RETURN

! Test line_seg function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 1.0 at x = 0.5
      call line_seg(0.5d+0, result,                                            &
     &              (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                      &
     &              (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                      &
     &              4)
      line_seg_test = check(1.0d+0, result, 1, "line_seg()")
      IF (.not.line_seg_test) RETURN
! Test line_seg function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 0.0 at x = 2.5
      call line_seg(2.5d+0, result,                                            &
     &              (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                      &
     &              (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                      &
     &              4)
      line_seg_test = check(0.0d+0, result, 2, "line_seg()")
      IF (.not.line_seg_test) RETURN
! Test line_seg function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 0.5 at x = 1.5
      call line_seg(1.5d+0, result,                                            &
     &              (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                      &
     &              (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                      &
     &              4)
      line_seg_test = check(0.5d+0, result, 3, "line_seg()")
      IF (.not.line_seg_test) RETURN
! Test line_seg function. yy(1.0,1.0,0.0,0.0) and xx(1.0,2.0,3.0,4.0) should
! have a value of 1.0 at x = 0.0
      call line_seg(0.0d+0, result,                                            &
     &              (/ 1.0d+0, 2.0d+0, 3.0d+0, 4.0d+0 /),                      &
     &              (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                      &
     &              4)
      line_seg_test = check(1.0d+0, result, 4, "line_seg()")
      IF (.not.line_seg_test) RETURN
! Test line_seg function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 0.0 at x = 4.0
      call line_seg(4.0d+0, result,                                            &
     &              (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                      &
     &              (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                      &
     &              4)
      line_seg_test = check(0.0d+0, result, 5, "line_seg()")
      IF (.not.line_seg_test) RETURN

! Test line_seg_int function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 0.5 at x = 0.5
      call line_seg_int(0.5d+0, result,                                        &
     &                  (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                  &
     &                  (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                  &
     &                  4)
      line_seg_test = check(0.5d+0, result, 1, "line_seg_int()")
      IF (.not.line_seg_test) RETURN
! Test line_seg_int function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 1.375 at x = 1.5
      call line_seg_int(1.5d+0, result,                                        &
     &                  (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                  &
     &                  (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                  &
     &                  4)
      line_seg_test = check(1.375d+0, result, 2, "line_seg_int()")
      IF (.not.line_seg_test) RETURN
! Test line_seg_int function. yy(1.0,1.0,0.0,0.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 1.5 at x = 2.5
      call line_seg_int(2.5d+0, result,                                        &
     &                  (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                  &
     &                  (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                  &
     &                  4)
      line_seg_test = check(1.5d+0, result, 3, "line_seg_int()")
      IF (.not.line_seg_test) RETURN
! Test line_seg_int function. yy(1.0,1.0,0.0,0.0) and xx(1.0,2.0,3.0,4.0) should
! have a value of 0.5 at x = 0.5
      call line_seg_int(0.5d+0, result,                                        &
     &                  (/ 1.0d+0, 2.0d+0, 3.0d+0, 4.0d+0 /),                  &
     &                  (/ 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0 /),                  &
     &                  4)
      line_seg_test = check(0.5d+0, result, 4, "line_seg_int()")
      IF (.not.line_seg_test) RETURN
! Test line_seg_int function. yy(1.0,0.0,0.0,1.0) and xx(0.0,1.0,2.0,3.0) should
! have a value of 2.5 at x = 4.0
      call line_seg_int(4.0d+0, result,                                        &
     &                  (/ 0.0d+0, 1.0d+0, 2.0d+0, 3.0d+0 /),                  &
     &                  (/ 1.0d+0, 0.0d+0, 0.0d+0, 1.0d+0 /),                  &
     &                  4)
      line_seg_test = check(2.5d+0, result, 5, "line_seg_int()")
      IF (.not.line_seg_test) RETURN
! Test line_seg_int function. yy(1.0,0.0,0.0,1.0) and xx(1.0,2.0,3.0,4.0) should
! have a value of 4.0 at x = 5.0
      call line_seg_int(5.0d+0, result,                                        &
     &                  (/ 1.0d+0, 2.0d+0, 3.0d+0, 4.0d+0 /),                  &
     &                  (/ 1.0d+0, 0.0d+0, 0.0d+0, 1.0d+0 /),                  &
     &                  4)
      line_seg_test = check(4.0d+0, result, 6, "line_seg_int()")
      IF (.not.line_seg_test) RETURN

      END FUNCTION

!*******************************************************************************
!  Check Test result
!*******************************************************************************
      FUNCTION check(expected, received, testNum, name)
!-------------------------------------------------------------------------------
!  Variable declarations.
!-------------------------------------------------------------------------------
      LOGICAL :: check
      REAL(rprec), INTENT(in) :: expected, received
      INTEGER, INTENT(in) :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      check = expected .eq. received
      IF (.not.check) THEN
         write(*,*) "line_segment_mod.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

      END MODULE line_segment
