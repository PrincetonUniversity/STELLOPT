!/////
! R8 !
!/////
!
! map point into (xmin, xmax) cell when boundary conditions are periodic.

subroutine EZspline_modulo1_r8(spline_o, p1, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  real(ezspline_r8) :: p1 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  endif

end subroutine EZspline_modulo1_r8

subroutine EZspline_modulo_array1_r8(spline_o, k1, p1, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: k1 
  real(ezspline_r8) :: p1(k1) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  endif

end subroutine EZspline_modulo_array1_r8


subroutine EZspline_modulo2_r8(spline_o, p1, p2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  real(ezspline_r8) :: p1, p2 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2 = MOD( p2, spline_o%x2max - spline_o%x2min )
  endif

end subroutine EZspline_modulo2_r8

subroutine EZspline_modulo_array2_r8(spline_o, k1, k2, p1, p2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r8) :: p1(k1), p2(k2) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k2) = MOD( p2(1:k2), spline_o%x2max - spline_o%x2min )
  endif

end subroutine EZspline_modulo_array2_r8

subroutine EZspline_modulo_cloud2_r8(spline_o, k, p1, p2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8) :: p1(k), p2(k) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k) = MOD( p1(1:k), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k) = MOD( p2(1:k), spline_o%x2max - spline_o%x2min )
  endif

end subroutine EZspline_modulo_cloud2_r8


subroutine EZspline_modulo3_r8(spline_o, p1, p2, p3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  real(ezspline_r8) :: p1, p2, p3 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2 = MOD( p2, spline_o%x2max - spline_o%x2min )
  endif

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3 = MOD( p3, spline_o%x3max - spline_o%x3min )
  endif

end subroutine EZspline_modulo3_r8

subroutine EZspline_modulo_array3_r8(spline_o, k1, k2, k3, p1, p2, p3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k1, k2, k3
  real(ezspline_r8) :: p1(k1), p2(k2), p3(k3) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k2) = MOD( p2(1:k2), spline_o%x2max - spline_o%x2min )
  endif

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3(1:k3) = MOD( p3(1:k3), spline_o%x3max - spline_o%x3min )
  endif

end subroutine EZspline_modulo_array3_r8

subroutine EZspline_modulo_cloud3_r8(spline_o, k, p1, p2, p3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8) :: p1(k), p2(k), p3(k) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k) = MOD( p1(1:k), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k) = MOD( p2(1:k), spline_o%x2max - spline_o%x2min )
  endif

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3(1:k) = MOD( p3(1:k), spline_o%x3max - spline_o%x3min )
  endif

end subroutine EZspline_modulo_cloud3_r8
!/////
! R4 !
!/////
!
! map point into (xmin, xmax) cell when boundary conditions are periodic.

subroutine EZspline_modulo1_r4(spline_o, p1, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  real(ezspline_r4) :: p1 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  endif

end subroutine EZspline_modulo1_r4

subroutine EZspline_modulo_array1_r4(spline_o, k1, p1, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: k1 
  real(ezspline_r4) :: p1(k1) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  endif

end subroutine EZspline_modulo_array1_r4


subroutine EZspline_modulo2_r4(spline_o, p1, p2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  real(ezspline_r4) :: p1, p2 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2 = MOD( p2, spline_o%x2max - spline_o%x2min )
  endif

end subroutine EZspline_modulo2_r4

subroutine EZspline_modulo_array2_r4(spline_o, k1, k2, p1, p2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r4) :: p1(k1), p2(k2) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k2) = MOD( p2(1:k2), spline_o%x2max - spline_o%x2min )
  endif

end subroutine EZspline_modulo_array2_r4

subroutine EZspline_modulo_cloud2_r4(spline_o, k, p1, p2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4) :: p1(k), p2(k) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k) = MOD( p1(1:k), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k) = MOD( p2(1:k), spline_o%x2max - spline_o%x2min )
  endif

end subroutine EZspline_modulo_cloud2_r4


subroutine EZspline_modulo3_r4(spline_o, p1, p2, p3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  real(ezspline_r4) :: p1, p2, p3 ! the location
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1 = MOD( p1, spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2 = MOD( p2, spline_o%x2max - spline_o%x2min )
  endif

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3 = MOD( p3, spline_o%x3max - spline_o%x3min )
  endif

end subroutine EZspline_modulo3_r4

subroutine EZspline_modulo_array3_r4(spline_o, k1, k2, k3, p1, p2, p3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: k1, k2, k3
  real(ezspline_r4) :: p1(k1), p2(k2), p3(k3) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k1) = MOD( p1(1:k1), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k2) = MOD( p2(1:k2), spline_o%x2max - spline_o%x2min )
  endif

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3(1:k3) = MOD( p3(1:k3), spline_o%x3max - spline_o%x3min )
  endif

end subroutine EZspline_modulo_array3_r4

subroutine EZspline_modulo_cloud3_r4(spline_o, k, p1, p2, p3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4) :: p1(k), p2(k), p3(k) ! the points
  integer, intent(out) :: ier

  ier = 0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     p1(1:k) = MOD( p1(1:k), spline_o%x1max - spline_o%x1min )
  endif

  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     p2(1:k) = MOD( p2(1:k), spline_o%x2max - spline_o%x2min )
  endif

  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     p3(1:k) = MOD( p3(1:k), spline_o%x3max - spline_o%x3min )
  endif

end subroutine EZspline_modulo_cloud3_r4
