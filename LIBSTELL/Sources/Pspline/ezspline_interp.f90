!/////
! R8 !
!/////

!!!
!!! 1-d
!!!

subroutine EZspline_interp1_r8(spline_o, p1, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  real(ezspline_r8) p1  ! the location where the interpolation is sought
  real(ezspline_r8) f   ! the interpolation
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3)=(/1, 0, 0/)
  real(ezspline_r8) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isLinear == 1) then

     call r8pc1ev(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call r8evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, ansr, ifail)

  else

     call r8herm1ev(p1, &
          &         spline_o%x1(1), spline_o%n1, &
          &         spline_o%ilin1,&
          &         spline_o%fspl(1,1),  &
          &         ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_r8


subroutine EZspline_interp1_array_r8(spline_o, k, p1, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k) ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: i, ifail
  integer, parameter :: ict(3)=(/1,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isLinear == 1) then

     call r8vecpc1(ict, k, p1, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then

     call r8vecspline(ict, k, p1, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else

     call r8vecherm1(ict, k, p1, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_array_r8

!!!
!!! 2-d
!!!

subroutine EZspline_interp2_r8(spline_o, p1, p2, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  real(ezspline_r8) p1, p2  ! the location where the interpolation is sought
  real(ezspline_r8) f          ! the interpolation
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6)=(/1, 0, 0, 0, 0, 0 /)
  real(ezspline_r8) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3),  &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call r8pc2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call r8evbicub(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else

     call r8herm2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_r8

subroutine EZspline_interp2_array_r8(spline_o, k1, k2, p1, p2, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r8) :: p1(k1), p2(k2) ! location arrays
  real(ezspline_r8) :: f(k1,k2)  ! interpolated function array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8gridintrp2d( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8gridpc2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then

     call r8gridbicub( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call r8gridherm2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_array_r8

subroutine EZspline_interp2_cloud_r8(spline_o, k, p1, p2, f, ier)
  ! list of coordinate doublets
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k) ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  if (spline_o%isHybrid == 1) then

     call r8vecintrp2d(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8vecpc2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call r8vecbicub(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call r8vecherm2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_cloud_r8


!!!
!!! 3-d
!!!

subroutine EZspline_interp3_r8(spline_o, p1, p2, p3, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  real(ezspline_r8) p1, p2, p3 ! the location where the interpolation is sought
  real(ezspline_r8) f          ! the interpolation

  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  real(ezspline_r8) ansr(1)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8evintrp3d(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), size(spline_o%fspl,4), &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call r8pc3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call r8evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else

     call r8herm3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_r8

subroutine EZspline_interp3_array_r8(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer :: k1, k2, k3
  real(ezspline_r8) :: p1(k1), p2(k2), p3(k3)  ! location arrays
  real(ezspline_r8) :: f(k1,k2,k3)  ! interpolant array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8gridintrp3d( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8gridpc3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call r8gridtricub( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else

     call r8gridherm3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_array_r8

subroutine EZspline_interp3_cloud_r8(spline_o, k, p1, p2, p3, f, ier)
  ! list of coordinate triplets
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)  ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call r8vecpc3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if (spline_o%isHermite == 0) then
     !
     call r8vectricub(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call r8vecherm3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_cloud_r8
!/////
! R4 !
!/////

!!!
!!! 1-d
!!!

subroutine EZspline_interp1_r4(spline_o, p1, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  real(ezspline_r4) p1  ! the location where the interpolation is sought
  real(ezspline_r4) f   ! the interpolation
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3)=(/1, 0, 0/)
  real(ezspline_r4) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isLinear == 1) then

     call pc1ev(p1, &
          &         spline_o%x1(1), spline_o%n1, &
          &         spline_o%ilin1,&
          &         spline_o%fspl(1,1),  &
          &         ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, ansr, ifail)

  else

     call herm1ev(p1, &
          &         spline_o%x1(1), spline_o%n1, &
          &         spline_o%ilin1,&
          &         spline_o%fspl(1,1),  &
          &         ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_r4


subroutine EZspline_interp1_array_r4(spline_o, k, p1, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k) ! location arrays
  real(ezspline_r4), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: i, ifail
  integer, parameter :: ict(3)=(/1,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isLinear == 1) then

     call vecpc1(ict, k, p1, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn,ifail)

  else if (spline_o%isHermite == 0) then

     call vecspline(ict, k, p1, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else

     call vecherm1(ict, k, p1, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp1_array_r4

!!!
!!! 2-d
!!!

subroutine EZspline_interp2_r4(spline_o, p1, p2, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  real(ezspline_r4) p1, p2  ! the location where the interpolation is sought
  real(ezspline_r4) f          ! the interpolation
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6)=(/1, 0, 0, 0, 0, 0 /)
  real(ezspline_r4) ansr(1)

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call pc2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then
     call evbicub(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  else

     call herm2ev(p1, p2,  &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_r4

subroutine EZspline_interp2_array_r4(spline_o, k1, k2, p1, p2, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r4) :: p1(k1), p2(k2) ! location arrays
  real(ezspline_r4) :: f(k1,k2)  ! interpolated function array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call gridintrp2d( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call gridpc2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then

     call gridbicub( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call gridherm2( &
          & p1, k1, &
          & p2, k2, &
          & f, k1, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97


end subroutine EZspline_interp2_array_r4

subroutine EZspline_interp2_cloud_r4(spline_o, k, p1, p2, f, ier)
  ! list of coordinate doublets
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k), p2(k) ! location arrays
  real(ezspline_r4), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier
  integer :: ifail
  integer, parameter :: ict(6) = (/1,0,0,0,0,0/)
  integer:: iwarn = 0

  ier = 0
  ifail = 0

  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  if (spline_o%isHybrid == 1) then

     call vecintrp2d(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call vecpc2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call vecbicub(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call vecherm2(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp2_cloud_r4


!!!
!!! 3-d
!!!

subroutine EZspline_interp3_r4(spline_o, p1, p2, p3, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  real(ezspline_r4) p1, p2, p3 ! the location where the interpolation is sought
  real(ezspline_r4) f          ! the interpolation

  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  real(ezspline_r4) ansr(1)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call evintrp3d(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), size(spline_o%fspl,4), &
          &   ict, ansr, ifail)

  else if (spline_o%isLinear == 1) then

     call pc3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else if (spline_o%isHermite == 0) then

     call evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  else

     call herm3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, ansr, ifail)

  endif

  f=ansr(1)

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_r4

subroutine EZspline_interp3_array_r4(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer :: k1, k2, k3
  real(ezspline_r4) :: p1(k1), p2(k2), p3(k3)  ! location arrays
  real(ezspline_r4) :: f(k1,k2,k3)  ! interpolant array
  integer, intent(out) :: ier

  integer ifail
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier =94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call gridintrp3d( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call gridpc3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else if (spline_o%isHermite == 0) then
     !
     call gridtricub( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else

     call gridherm3( &
          & p1, k1, &
          & p2, k2, &
          & p3, k3, &
          & f, k1, k2, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_array_r4

subroutine EZspline_interp3_cloud_r4(spline_o, k, p1, p2, p3, f, ier)
  ! list of coordinate triplets
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)  ! location arrays
  real(ezspline_r4), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier

  integer :: ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if (spline_o%isLinear == 1) then

     call vecpc3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if (spline_o%isHermite == 0) then
     !
     call vectricub(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call vecherm3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 97

end subroutine EZspline_interp3_cloud_r4
