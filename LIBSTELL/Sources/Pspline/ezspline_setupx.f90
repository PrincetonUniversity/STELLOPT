subroutine ezspline_setup1_r8x(spline_o, f, ier)

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline1_r8) :: spline_o
  real(ezspline_r8) :: f(size(spline_o%fspl,2))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup1_r8x

subroutine ezspline_setup2_r8x(spline_o, f, ier)

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline2_r8) :: spline_o
  real(ezspline_r8) :: f(size(spline_o%fspl,2),size(spline_o%fspl,3))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup2_r8x

subroutine ezspline_setup3_r8x(spline_o, f, ier)

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline3_r8) :: spline_o
  real(ezspline_r8) :: f(size(spline_o%fspl,2),size(spline_o%fspl,3), &
       size(spline_o%fspl,4))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup3_r8x

subroutine ezspline_setup1_r4x(spline_o, f, ier)

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline1_r4) :: spline_o
  real(ezspline_r4) :: f(size(spline_o%fspl,2))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup1_r4x

subroutine ezspline_setup2_r4x(spline_o, f, ier)

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline2_r4) :: spline_o
  real(ezspline_r4) :: f(size(spline_o%fspl,2),size(spline_o%fspl,3))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup2_r4x

subroutine ezspline_setup3_r4x(spline_o, f, ier)

  !  setup call (set array size from object)

  use EZspline_obj
  use EZspline
  implicit NONE

  type(EZspline3_r4) :: spline_o
  real(ezspline_r4) :: f(size(spline_o%fspl,2),size(spline_o%fspl,3), &
       size(spline_o%fspl,4))
  integer, intent(out) :: ier

  call ezspline_setup(spline_o, f, ier)

end subroutine ezspline_setup3_r4x
