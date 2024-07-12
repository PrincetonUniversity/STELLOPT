!/////
! R8 !
!/////
!
! 1-D
!

subroutine EZspline_gradient1_r8(spline_o, &
     p1, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  real(ezspline_r8), intent(in) :: p1
  real(ezspline_r8), intent(out) :: df
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3) = (/0, 1, 0/)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if(spline_o%isLinear == 1) then

     call r8pc1ev(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, df, ifail)

  else if(spline_o%isHermite == 0) then

     call r8evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, df, ifail)

  else

     call r8herm1ev(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, df, ifail)

  endif
  if(ifail/=0) ier = 95

end subroutine EZspline_gradient1_r8


subroutine EZspline_gradient1_array_r8(spline_o, k1, &
     p1, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: k1
  real(ezspline_r8), intent(in) :: p1(k1)
  real(ezspline_r8), intent(out) :: df(k1) 
  !
  ! ier=95 some error occurred in EZspline_gradient
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(3) = (/0, 1, 0/)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if(spline_o%isLinear == 1) then

     call r8vecpc1(ict, k1, p1, k1, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else if(spline_o%isHermite == 0) then

     call r8vecspline(ict, k1, p1, k1, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)


  else

     call r8vecherm1(ict, k1, p1, k1, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 95
end subroutine EZspline_gradient1_array_r8


!!
!! 2-D
!!

subroutine EZspline_gradient2_r8(spline_o, &
     p1, p2, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  real(ezspline_r8), intent(in) :: p1, p2
  real(ezspline_r8), intent(out) :: df(2) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6) = (/0, 1, 1, 0, 0, 0/)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3),  &
          &   ict, df, ifail)

  else if(spline_o%isLinear == 1) then

     call r8pc2ev(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, df, ifail)

  else if(spline_o%isHermite == 0) then

     call r8evbicub(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, df, ifail)

  else

     call r8herm2ev(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, df, ifail)
  endif

  if(ifail/=0) ier = 95

end subroutine EZspline_gradient2_r8

subroutine EZspline_gradient2_cloud_r8(spline_o, k, &
     p1, p2, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out) :: df(k,2) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6) = (/0, 1, 1, 0, 0, 0/)
  integer::  iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8vecintrp2d(ict, k, p1, p2, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call r8vecpc2(ict, k, p1, p2, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn,ifail)

  else if(spline_o%isHermite == 0) then

     call r8vecbicub(ict, k, p1, p2, k, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call r8vecherm2(ict, k, p1, p2, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn,ifail)

  endif
  if(ifail /= 0) ier = 95

end subroutine EZspline_gradient2_cloud_r8


subroutine EZspline_gradient2_array_r8(spline_o, k1, k2, &
     p1, p2, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r8), intent(in) :: p1(k1), p2(k2)
  real(ezspline_r8), intent(out) :: df(k1, k2, 2) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6) = (/0, 1, 1, 0, 0, 0/)
  integer:: iwarn=0
  real(ezspline_r8), dimension(:), allocatable :: p1_cloud, p2_cloud
  integer k12

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  k12 = k1*k2
  allocate(p1_cloud(k12), p2_cloud(k12), stat=ifail)
  if(ifail /= 0) then
     ier = 30
     return
  endif

  p1_cloud = reshape( &
       & source=spread(source=p1, dim=2, ncopies=k2), shape=(/k12/))
  p2_cloud = reshape( &
       & source=spread(source=p2, dim=1, ncopies=k1), shape=(/k12/))

  if(spline_o%isHybrid == 1) then

     call r8vecintrp2d(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call r8vecpc2(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if(spline_o%isHermite == 0) then

     call r8vecbicub(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call r8vecherm2(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  deallocate(p1_cloud, p2_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 31
     return
  endif

  if(ifail /= 0) ier = 95

end subroutine EZspline_gradient2_array_r8


!!!
!!! 3-D
!!!

subroutine EZspline_gradient3_r8(spline_o, &
     p1, p2, p3, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  real(ezspline_r8), intent(in) :: p1, p2, p3
  real(ezspline_r8), intent(out) :: df(3) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10) = (/0, 1, 1, 1, 0, 0, 0, 0, 0, 0/)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
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
          &   ict, df, ifail)

  else if(spline_o%isLinear == 1) then

     call r8pc3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, df, ifail)

  else if(spline_o%isHermite == 0) then

     call r8evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, df, ifail)

  else

     call r8herm3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, df, ifail)

  endif

  if(ifail/=0) ier = 95

end subroutine EZspline_gradient3_r8

subroutine EZspline_gradient3_cloud_r8(spline_o, k, &
     p1, p2, p3, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r8), intent(out) :: df(k, 3) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10) = (/0, 1, 1, 1, 0, 0, 0, 0, 0, 0/)
  integer :: iwarn=0

  ier = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call r8vecintrp3d(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call r8vecpc3(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if(spline_o%isHermite == 0) then

     call r8vectricub(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call r8vecherm3(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 95


end subroutine EZspline_gradient3_cloud_r8


subroutine EZspline_gradient3_array_r8(spline_o, k1, k2, k3, &
     p1, p2, p3, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k1, k2, k3
  real(ezspline_r8), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(ezspline_r8), intent(out) :: df(k1, k2, k3, 3) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10) = (/0, 1, 1, 1, 0, 0, 0, 0, 0, 0/)
  integer:: iwarn=0
  real(ezspline_r8), dimension(:), allocatable :: p1_cloud, p2_cloud, p3_cloud
  integer k123

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  k123 = k1*k2*k3
  allocate(p1_cloud(k123), p2_cloud(k123), p3_cloud(k123), stat=ifail)
  if(ifail /= 0) then
     ier = 30
     return
  endif

  p1_cloud = reshape(source=spread( &
       & source=spread(source=p1, dim=2, ncopies=k2), &
       & dim=3, ncopies=k3), shape=(/k123/))
  p2_cloud = reshape(source=spread( &
       & source=spread(source=p2, dim=1, ncopies=k1), &
       & dim=3, ncopies=k3), shape=(/k123/))
  p3_cloud = reshape(source=spread( &
       & source=spread(source=p3, dim=1, ncopies=k1), &
       & dim=2, ncopies=k2), shape=(/k123/))

  if (spline_o%isHybrid == 1) then

     call r8vecintrp3d(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call r8vecpc3(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if(spline_o%isHermite == 0) then

     call r8vectricub(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call r8vecherm3(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 95

  deallocate(p1_cloud, p2_cloud, p3_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 31
     return
  endif

end subroutine EZspline_gradient3_array_r8
!/////
! R4 !
!/////
!
! 1-D
!


subroutine EZspline_gradient1_r4(spline_o, &
     p1, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  real(ezspline_r4), intent(in) :: p1
  real(ezspline_r4), intent(out) :: df
  integer, intent(out) :: ier

  integer ifail
  integer, parameter :: ict(3) = (/0, 1, 0/)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if(spline_o%isLinear == 1) then

     call pc1ev(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, df, ifail)

  else if(spline_o%isHermite == 0) then

     call evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, df, ifail)

  else

     call herm1ev(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, df, ifail)

  endif
  if(ifail/=0) ier = 95

end subroutine EZspline_gradient1_r4


subroutine EZspline_gradient1_array_r4(spline_o, k1, &
     p1, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: k1
  real(ezspline_r4), intent(in) :: p1(k1)
  real(ezspline_r4), intent(out) :: df(k1) 
  !
  ! ier=95 some error occurred in EZspline_gradient
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(3) = (/0, 1, 0/)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if(spline_o%isLinear == 1) then

     call vecpc1(ict, k1, p1, k1, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  else if(spline_o%isHermite == 0) then

     call vecspline(ict, k1, p1, k1, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)


  else

     call vecherm1(ict, k1, p1, k1, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 95
end subroutine EZspline_gradient1_array_r4


!!
!! 2-D
!!

subroutine EZspline_gradient2_r4(spline_o, &
     p1, p2, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  real(ezspline_r4), intent(in) :: p1, p2
  real(ezspline_r4), intent(out) :: df(2) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6) = (/0, 1, 1, 0, 0, 0/)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3),  &
          &   ict, df, ifail)

  else if(spline_o%isLinear == 1) then

     call pc2ev(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, df, ifail)

  else if(spline_o%isHermite == 0) then

     call evbicub(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, df, ifail)

  else

     call herm2ev(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, df, ifail)

  endif

  if(ifail/=0) ier = 95

end subroutine EZspline_gradient2_r4

subroutine EZspline_gradient2_cloud_r4(spline_o, k, &
     p1, p2, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k), p2(k)
  real(ezspline_r4), intent(out) :: df(k,2) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6) = (/0, 1, 1, 0, 0, 0/)
  integer::  iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  if (spline_o%isHybrid == 1) then

     call vecintrp2d(ict, k, p1, p2, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call vecpc2(ict, k, p1, p2, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn,ifail)

  else if(spline_o%isHermite == 0) then

     call vecbicub(ict, k, p1, p2, k, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)


  else

     call vecherm2(ict, k, p1, p2, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn,ifail)

  endif
  if(ifail /= 0) ier = 95

end subroutine EZspline_gradient2_cloud_r4


subroutine EZspline_gradient2_array_r4(spline_o, k1, k2, &
     p1, p2, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: k1, k2
  real(ezspline_r4), intent(in) :: p1(k1), p2(k2)
  real(ezspline_r4), intent(out) :: df(k1, k2, 2) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(6) = (/0, 1, 1, 0, 0, 0/)
  integer:: iwarn=0
  real(ezspline_r4), dimension(:), allocatable :: p1_cloud, p2_cloud
  integer k12

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif


  k12 = k1*k2
  allocate(p1_cloud(k12), p2_cloud(k12), stat=ifail)
  if(ifail /= 0) then
     ier = 30
     return
  endif

  p1_cloud = reshape( &
       & source=spread(source=p1, dim=2, ncopies=k2), shape=(/k12/))
  p2_cloud = reshape( &
       & source=spread(source=p2, dim=1, ncopies=k1), shape=(/k12/))

  if (spline_o%isHybrid == 1) then

     call vecintrp2d(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call vecpc2(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else if(spline_o%isHermite == 0) then

     call vecbicub(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  else

     call vecherm2(ict, k12, p1_cloud, p2_cloud, k12, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  deallocate(p1_cloud, p2_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 31
     return
  endif

  if(ifail /= 0) ier = 95

end subroutine EZspline_gradient2_array_r4


!!!
!!! 3-D
!!!

subroutine EZspline_gradient3_r4(spline_o, &
     p1, p2, p3, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  real(ezspline_r4), intent(in) :: p1, p2, p3
  real(ezspline_r4), intent(out) :: df(3) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10) = (/0, 1, 1, 1, 0, 0, 0, 0, 0, 0/)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
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
          &   ict, df, ifail)

  else if(spline_o%isLinear == 1) then

     call pc3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, df, ifail)

  else if(spline_o%isHermite == 0) then

     call evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, df, ifail)

  else

     call herm3ev(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, df, ifail)

  endif

  if(ifail/=0) ier = 95

end subroutine EZspline_gradient3_r4

subroutine EZspline_gradient3_cloud_r4(spline_o, k, &
     p1, p2, p3, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: k
  real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r4), intent(out) :: df(k, 3) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10) = (/0, 1, 1, 1, 0, 0, 0, 0, 0, 0/)
  integer :: iwarn=0

  ier = 0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if (spline_o%isHybrid == 1) then

     call vecintrp3d(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call vecpc3(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if(spline_o%isHermite == 0) then

     call vectricub(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call vecherm3(ict, k, p1, p2, p3, k, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 95


end subroutine EZspline_gradient3_cloud_r4


subroutine EZspline_gradient3_array_r4(spline_o, k1, k2, k3, &
     p1, p2, p3, df, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: k1, k2, k3
  real(ezspline_r4), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(ezspline_r4), intent(out) :: df(k1, k2, k3, 3) 
  !
  integer, intent(out) :: ier
  integer ifail
  integer, parameter :: ict(10) = (/0, 1, 1, 1, 0, 0, 0, 0, 0, 0/)
  integer:: iwarn=0
  real(ezspline_r4), dimension(:), allocatable :: p1_cloud, p2_cloud, p3_cloud
  integer k123

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  k123 = k1*k2*k3
  allocate(p1_cloud(k123), p2_cloud(k123), p3_cloud(k123), stat=ifail)
  if(ifail /= 0) then
     ier = 30
     return
  endif

  p1_cloud = reshape(source=spread( &
       & source=spread(source=p1, dim=2, ncopies=k2), &
       & dim=3, ncopies=k3), shape=(/k123/))
  p2_cloud = reshape(source=spread( &
       & source=spread(source=p2, dim=1, ncopies=k1), &
       & dim=3, ncopies=k3), shape=(/k123/))
  p3_cloud = reshape(source=spread( &
       & source=spread(source=p3, dim=1, ncopies=k1), &
       & dim=2, ncopies=k2), shape=(/k123/))

  if (spline_o%isHybrid == 1) then

     call vecintrp3d(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isLinear == 1) then

     call vecpc3(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if(spline_o%isHermite == 0) then

     call vectricub(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  else

     call vecherm3(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, df, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  endif

  if(ifail /= 0) ier = 95

  deallocate(p1_cloud, p2_cloud, p3_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 31
     return
  endif

end subroutine EZspline_gradient3_array_r4
