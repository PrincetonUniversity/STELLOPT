!/////
! R8 !
!/////
!
! 1-D
!

subroutine EZspline_derivative1_r8(spline_o, i1, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: i1
  real(ezspline_r8), intent(in) :: p1
  real(ezspline_r8), intent(out) :: f
  integer, intent(out) :: ier

  integer ifail

  integer ict_herm(2)
  integer ict(3)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1 < 0) then
     ier = 12
     return
  endif
  if(i1 > 3) then
     ier = 13
     return
  endif

  if(i1 >1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHermite==1 .or. spline_o%isLinear==1) then

     if(i1 == 1) then
        ict_herm = (/0, 1/)
     else
        ict_herm = (/1, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8herm1ev(p1, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%ilin1, &
             &   spline_o%fspl(1,1), &
             &   ict_herm, f, ifail)
     else
        call r8pc1ev(p1, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%ilin1, &
             &   spline_o%fspl(1,1), &
             &   ict_herm, f, ifail)
     endif

  else   

     call ezmake_ict1(i1,ict)

     call r8evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, f, ifail)

  endif
  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative1_r8

subroutine EZspline_derivative1_array_r8(spline_o, i1, k1, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: i1
  integer, intent(in) :: k1
  real(ezspline_r8), intent(in) :: p1(k1)
  real(ezspline_r8), intent(out) :: f(k1)
  integer, intent(out) :: ier

  integer ifail

  integer ict(3)
  integer ict_herm(2)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1 < 0) then
     ier = 12
     return
  endif
  if(i1 > 3) then
     ier = 13
     return
  endif

  if(i1 >1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHermite==1 .or. spline_o%isLinear==1) then

     if(i1 == 1) then
        ict_herm = (/0, 1/)
     else
        ict_herm = (/1, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8vecherm1(ict_herm, k1, p1, k1, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%fspl(1,1), &
             & iwarn, ifail)
     else
        call r8vecpc1(ict_herm, k1, p1, k1, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%fspl(1,1), &
             & iwarn, ifail)
     endif

  else   

     call ezmake_ict1(i1,ict)

     call r8vecspline(ict, k1, p1, k1, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative1_array_r8


!!
!! 2-D
!!


subroutine EZspline_derivative2_r8(spline_o, i1, i2, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: i1, i2
  real(ezspline_r8), intent(in) :: p1, p2
  real(ezspline_r8), intent(out) :: f
  integer, intent(out) :: ier

  integer ifail

  integer ict(6)
  integer ict_herm(4)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if(i1<0 .OR. i2<0) then
     ier = 12
     return
  endif
  if(max(i1,i2) > 3) then
     ier = 13
     return
  endif

  if(max(i1,i2) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict2(i1,i2,ict)

     call r8evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3),  &
          &   ict, f, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1) then
     if(i1 == 1 .AND. i2 == 0) then
        ict_herm = (/0, 1, 0, 0/) 
     else if(i1 == 0 .AND. i2 == 1) then
        ict_herm = (/0, 0, 1, 0/) 
     else if(i1 == 1 .AND. i2 == 1) then
        ict_herm = (/0, 0, 0, 1/) 
     else
        ict_herm = (/1, 0, 0, 0/) 
     endif

     if(spline_o%isLinear == 0) then
        call r8herm2ev(p1, p2, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%ilin1, spline_o%ilin2, &
             &   spline_o%fspl(1,1,1), spline_o%n1, &
             &   ict_herm, f, ifail)
     else
        call r8pc2ev(p1, p2, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%ilin1, spline_o%ilin2, &
             &   spline_o%fspl(1,1,1), spline_o%n1, &
             &   ict_herm, f, ifail)
     endif

  else

     call ezmake_ict2(i1,i2,ict)

     call r8evbicub(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, f, ifail)

  endif
  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative2_r8



subroutine EZspline_derivative2_array_r8(spline_o, i1, i2, &
     & k1, k2, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: i1, i2, k1, k2
  real(ezspline_r8), intent(in) :: p1(k1), p2(k2)
  real(ezspline_r8), intent(out) :: f(k1,k2)
  integer, intent(out) :: ier

  integer ifail

  integer ict(6)
  integer ict_herm(4)

  real(ezspline_r8), dimension(:), allocatable :: p1_cloud, p2_cloud
  integer k12
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1<0 .OR. i2<0) then
     ier = 12
     return
  endif
  if(max(i1,i2) > 3) then
     ier = 13
     return
  endif

  if(max(i1,i2) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  k12 = k1*k2
  allocate(p1_cloud(k12), p2_cloud(k12), stat=ifail)
  if(ifail /= 0) then
     ier = 32
     return
  endif

  p1_cloud = reshape( &
       & source=spread(source=p1, dim=2, ncopies=k2), &
       & shape=(/k12/))
  p2_cloud = reshape( &
       & source=spread(source=p2, dim=1, ncopies=k1), &
       & shape=(/k12/))

  if(spline_o%isHybrid == 1) then

     call ezmake_ict2(i1,i2,ict)

     call r8vecintrp2d(ict, k12, p1_cloud, p2_cloud, k12, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0) then
        ict_herm = (/0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1) then
        ict_herm = (/0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1) then
        ict_herm = (/0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8vecherm2(ict_herm, k12, p1_cloud, p2_cloud, k12, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     else
        call r8vecpc2(ict_herm, k12, p1_cloud, p2_cloud, k12, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict2(i1,i2,ict)  ! higher derivatives

     call r8vecbicub(ict, k12, p1_cloud, p2_cloud, k12, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif


  if(ifail /= 0) ier = 96

  deallocate(p1_cloud, p2_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 33
     return
  endif


end subroutine EZspline_derivative2_array_r8


subroutine EZspline_derivative2_cloud_r8(spline_o, i1, i2, &
     & k, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: i1, i2, k
  real(ezspline_r8), intent(in) :: p1(k), p2(k)
  real(ezspline_r8), intent(out) :: f(k)
  integer, intent(out) :: ier

  integer ifail

  integer ict(6)
  integer ict_herm(4)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1<0 .OR. i2<0) then
     ier = 12
     return
  endif
  if(max(i1,i2) > 3) then
     ier = 13
     return
  endif

  if(max(i1,i2) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict2(i1,i2,ict)

     call r8vecintrp2d(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0) then
        ict_herm = (/0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1) then
        ict_herm = (/0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1) then
        ict_herm = (/0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8vecherm2(ict_herm, k, p1, p2, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     else
        call r8vecpc2(ict_herm, k, p1, p2, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict2(i1,i2,ict)

     call r8vecbicub(ict, k, p1, p2, &
          & k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative2_cloud_r8

!!!
!!! 3-D
!!!


subroutine EZspline_derivative3_r8(spline_o, i1, i2, i3, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: i1, i2, i3
  real(ezspline_r8), intent(in) :: p1, p2, p3
  real(ezspline_r8), intent(out) :: f
  integer, intent(out) :: ier

  integer ifail

  integer ict(10)
  integer ict_herm(8)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1<0 .OR. i2<0 .OR. i3<0) then
     ier = 12
     return
  endif
  if(max(i1,max(i2,i3)) > 3) then
     ier = 13
     return
  endif

  if(maxval((/i1,i2,i3/)) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict3(i1,i2,i3,ict)

     call r8evintrp3d(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), size(spline_o%fspl,4), &
          &   ict, f, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0 .AND. i3 == 0) then
        ict_herm = (/0, 1, 0, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 1, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 1, 0, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 0, 0, 1, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0, 0, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8herm3ev(p1, p2, p3, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%x3(1), spline_o%n3, &
             &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
             &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             &   ict_herm, f, ifail)
     else
        call r8pc3ev(p1, p2, p3, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%x3(1), spline_o%n3, &
             &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
             &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             &   ict_herm, f, ifail)
     endif

  else

     call ezmake_ict3(i1,i2,i3,ict)

     call r8evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, f, ifail)

  endif
  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative3_r8


subroutine EZspline_derivative3_array_r8(spline_o, i1, i2, i3, &
     & k1, k2, k3, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: i1, i2, i3, k1, k2, k3
  real(ezspline_r8), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(ezspline_r8), intent(out) :: f(k1,k2,k3)
  integer, intent(out) :: ier

  integer ifail

  integer ict(10)
  integer ict_herm(8)

  real(ezspline_r8), dimension(:), allocatable :: p1_cloud, p2_cloud, p3_cloud
  integer k123
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1<0 .OR. i2<0 .OR. i3<0) then
     ier = 12
     return
  endif
  if(max(i1,max(i2,i3)) > 3) then
     ier = 13
     return
  endif

  if(maxval((/i1,i2,i3/)) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  k123 = k1*k2*k3
  allocate(p1_cloud(k123), p2_cloud(k123), p3_cloud(k123), stat=ifail)
  if(ifail /= 0) then
     ier = 32
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

  if(spline_o%isHybrid == 1) then

     call ezmake_ict3(i1,i2,i3,ict)

     call r8vecintrp3d(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0 .AND. i3 == 0) then
        ict_herm = (/0, 1, 0, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 1, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 1, 0, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 0, 0, 1, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0, 0, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8vecherm3(ict_herm, k123, p1_cloud, p2_cloud, p3_cloud, k123, &
             & f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     else
        call r8vecpc3(ict_herm, k123, p1_cloud, p2_cloud, p3_cloud, k123, &
             & f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict3(i1,i2,i3,ict)

     call r8vectricub(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif


  if(ifail /= 0) ier = 96

  deallocate(p1_cloud, p2_cloud, p3_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 33
     return
  endif


end subroutine EZspline_derivative3_array_r8

subroutine EZspline_derivative3_cloud_r8(spline_o, i1, i2, i3, &
     & k, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: i1, i2, i3, k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r8), intent(out) :: f(k)
  integer, intent(out) :: ier

  integer ifail

  integer ict(10)
  integer ict_herm(8)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r8
  if(i1<0 .OR. i2<0 .OR. i3<0) then
     ier = 12
     return
  endif
  if(max(i1,max(i2,i3)) > 3) then
     ier = 13
     return
  endif

  if(maxval((/i1,i2,i3/)) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict3(i1,i2,i3,ict)

     call r8vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0 .AND. i3 == 0) then
        ict_herm = (/0, 1, 0, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 1, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 1, 0, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 0, 0, 1, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0, 0, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call r8vecherm3(ict_herm, k, p1, p2, p3, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     else
        call r8vecpc3(ict_herm, k, p1, p2, p3, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict3(i1,i2,i3,ict)

     call r8vectricub(ict, k, p1, p2, p3, &
          & k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative3_cloud_r8
!/////
! R4 !
!/////
!
! 1-D
!

subroutine EZspline_derivative1_r4(spline_o, i1, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: i1
  real(ezspline_r4), intent(in) :: p1
  real(ezspline_r4), intent(out) :: f
  integer, intent(out) :: ier

  integer ifail

  integer ict_herm(2)
  integer ict(3)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1 < 0) then
     ier = 12
     return
  endif
  if(i1 > 3) then
     ier = 13
     return
  endif

  if(i1 >1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHermite==1 .or. spline_o%isLinear==1) then

     if(i1 == 1) then
        ict_herm = (/0, 1/)
     else
        ict_herm = (/1, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call herm1ev(p1, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%ilin1, &
             &   spline_o%fspl(1,1), &
             &   ict_herm, f, ifail)
     else
        call pc1ev(p1, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%ilin1, &
             &   spline_o%fspl(1,1), &
             &   ict_herm, f, ifail)
     endif

  else   

     call ezmake_ict1(i1,ict)

     call evspline(p1, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%ilin1, &
          &   spline_o%fspl(1,1), &
          &   ict, f, ifail)

  endif
  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative1_r4

subroutine EZspline_derivative1_array_r4(spline_o, i1, k1, p1, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: i1
  integer, intent(in) :: k1
  real(ezspline_r4), intent(in) :: p1(k1)
  real(ezspline_r4), intent(out) :: f(k1)
  integer, intent(out) :: ier

  integer ifail

  integer ict(3)
  integer ict_herm(2)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1 < 0) then
     ier = 12
     return
  endif
  if(i1 > 3) then
     ier = 13
     return
  endif

  if(i1 >1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHermite==1 .or. spline_o%isLinear==1) then

     if(i1 == 1) then
        ict_herm = (/0, 1/)
     else
        ict_herm = (/1, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call vecherm1(ict_herm, k1, p1, k1, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%fspl(1,1), &
             & iwarn, ifail)
     else
        call vecpc1(ict_herm, k1, p1, k1, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%fspl(1,1), &
             & iwarn, ifail)
     endif

  else   

     call ezmake_ict1(i1,ict)

     call vecspline(ict, k1, p1, k1, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%fspl(1,1), &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative1_array_r4


!!
!! 2-D
!!


subroutine EZspline_derivative2_r4(spline_o, i1, i2, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: i1, i2
  real(ezspline_r4), intent(in) :: p1, p2
  real(ezspline_r4), intent(out) :: f
  integer, intent(out) :: ier

  integer ifail

  integer ict(6)
  integer ict_herm(4)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  if(i1<0 .OR. i2<0) then
     ier = 12
     return
  endif
  if(max(i1,i2) > 3) then
     ier = 13
     return
  endif

  if(max(i1,i2) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict2(i1,i2,ict)

     call evintrp2d(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3),  &
          &   ict, f, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1) then
     if(i1 == 1 .AND. i2 == 0) then
        ict_herm = (/0, 1, 0, 0/) 
     else if(i1 == 0 .AND. i2 == 1) then
        ict_herm = (/0, 0, 1, 0/) 
     else if(i1 == 1 .AND. i2 == 1) then
        ict_herm = (/0, 0, 0, 1/) 
     else
        ict_herm = (/1, 0, 0, 0/) 
     endif

     if(spline_o%isLinear == 0) then
        call herm2ev(p1, p2, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%ilin1, spline_o%ilin2, &
             &   spline_o%fspl(1,1,1), spline_o%n1, &
             &   ict_herm, f, ifail)
     else
        call pc2ev(p1, p2, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%ilin1, spline_o%ilin2, &
             &   spline_o%fspl(1,1,1), spline_o%n1, &
             &   ict_herm, f, ifail)
     endif

  else

     call ezmake_ict2(i1,i2,ict)

     call evbicub(p1, p2, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%ilin1, spline_o%ilin2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   ict, f, ifail)

  endif
  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative2_r4



subroutine EZspline_derivative2_array_r4(spline_o, i1, i2, &
     & k1, k2, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: i1, i2, k1, k2
  real(ezspline_r4), intent(in) :: p1(k1), p2(k2)
  real(ezspline_r4), intent(out) :: f(k1,k2)
  integer, intent(out) :: ier

  integer ifail

  integer ict(6)
  integer ict_herm(4)

  real(ezspline_r4), dimension(:), allocatable :: p1_cloud, p2_cloud
  integer k12
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1<0 .OR. i2<0) then
     ier = 12
     return
  endif
  if(max(i1,i2) > 3) then
     ier = 13
     return
  endif

  if(max(i1,i2) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  k12 = k1*k2
  allocate(p1_cloud(k12), p2_cloud(k12), stat=ifail)
  if(ifail /= 0) then
     ier = 32
     return
  endif

  p1_cloud = reshape( &
       & source=spread(source=p1, dim=2, ncopies=k2), &
       & shape=(/k12/))
  p2_cloud = reshape( &
       & source=spread(source=p2, dim=1, ncopies=k1), &
       & shape=(/k12/))

  if(spline_o%isHybrid == 1) then

     call ezmake_ict2(i1,i2,ict)

     call vecintrp2d(ict, k12, p1_cloud, p2_cloud, k12, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0) then
        ict_herm = (/0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1) then
        ict_herm = (/0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1) then
        ict_herm = (/0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call vecherm2(ict_herm, k12, p1_cloud, p2_cloud, k12, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     else
        call vecpc2(ict_herm, k12, p1_cloud, p2_cloud, k12, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     endif


  else

     call ezmake_ict2(i1,i2,ict)

     call vecbicub(ict, k12, p1_cloud, p2_cloud, k12, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif


  if(ifail /= 0) ier = 96

  deallocate(p1_cloud, p2_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 33
     return
  endif


end subroutine EZspline_derivative2_array_r4


subroutine EZspline_derivative2_cloud_r4(spline_o, i1, i2, &
     & k, p1, p2, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: i1, i2, k
  real(ezspline_r4), intent(in) :: p1(k), p2(k)
  real(ezspline_r4), intent(out) :: f(k)
  integer, intent(out) :: ier

  integer ifail

  integer ict(6)
  integer ict_herm(4)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1<0 .OR. i2<0) then
     ier = 12
     return
  endif
  if(max(i1,i2) > 3) then
     ier = 13
     return
  endif

  if(max(i1,i2) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict2(i1,i2,ict)

     call vecintrp2d(ict, k, p1, p2, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3),  &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0) then
        ict_herm = (/0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1) then
        ict_herm = (/0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1) then
        ict_herm = (/0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call vecherm2(ict_herm, k, p1, p2, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     else
        call vecpc2(ict_herm, k, p1, p2, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%fspl(1,1,1), spline_o%n1, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict2(i1,i2,ict)

     call vecbicub(ict, k, p1, p2, &
          & k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative2_cloud_r4

!!!
!!! 3-D
!!!


subroutine EZspline_derivative3_r4(spline_o, i1, i2, i3, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: i1, i2, i3
  real(ezspline_r4), intent(in) :: p1, p2, p3
  real(ezspline_r4), intent(out) :: f
  integer, intent(out) :: ier

  integer ifail

  integer ict(10)
  integer ict_herm(8)

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1<0 .OR. i2<0 .OR. i3<0) then
     ier = 12
     return
  endif
  if(max(i1,max(i2,i3)) > 3) then
     ier = 13
     return
  endif

  if(maxval((/i1,i2,i3/)) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict3(i1,i2,i3,ict)

     call evintrp3d(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   size(spline_o%fspl,1), size(spline_o%fspl,2), &
          &   size(spline_o%fspl,3), size(spline_o%fspl,4), &
          &   ict, f, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0 .AND. i3 == 0) then
        ict_herm = (/0, 1, 0, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 1, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 1, 0, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 0, 0, 1, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0, 0, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call herm3ev(p1, p2, p3, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%x3(1), spline_o%n3, &
             &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
             &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             &   ict_herm, f, ifail)
     else
        call pc3ev(p1, p2, p3, &
             &   spline_o%x1(1), spline_o%n1, &
             &   spline_o%x2(1), spline_o%n2, &
             &   spline_o%x3(1), spline_o%n3, &
             &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
             &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             &   ict_herm, f, ifail)
     endif

  else

     call ezmake_ict3(i1,i2,i3,ict)

     call evtricub(p1, p2, p3, &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   ict, f, ifail)


  endif
  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative3_r4


subroutine EZspline_derivative3_array_r4(spline_o, i1, i2, i3, &
     & k1, k2, k3, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: i1, i2, i3, k1, k2, k3
  real(ezspline_r4), intent(in) :: p1(k1), p2(k2), p3(k3)
  real(ezspline_r4), intent(out) :: f(k1,k2,k3)
  integer, intent(out) :: ier

  integer ifail

  integer ict(10)
  integer ict_herm(8)

  real(ezspline_r4), dimension(:), allocatable :: p1_cloud, p2_cloud, p3_cloud
  integer k123
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1<0 .OR. i2<0 .OR. i3<0) then
     ier = 12
     return
  endif
  if(max(i1,max(i2,i3)) > 3) then
     ier = 13
     return
  endif

  if(maxval((/i1,i2,i3/)) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  k123 = k1*k2*k3
  allocate(p1_cloud(k123), p2_cloud(k123), p3_cloud(k123), stat=ifail)
  if(ifail /= 0) then
     ier = 32
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

  if(spline_o%isHybrid == 1) then

     call ezmake_ict3(i1,i2,i3,ict)

     call vecintrp3d(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0 .AND. i3 == 0) then
        ict_herm = (/0, 1, 0, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 1, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 1, 0, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 0, 0, 1, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0, 0, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call vecherm3(ict_herm, k123, p1_cloud, p2_cloud, p3_cloud, k123, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     else
        call vecpc3(ict_herm, k123, p1_cloud, p2_cloud, p3_cloud, k123, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict3(i1,i2,i3,ict)

     call vectricub(ict, k123, p1_cloud, p2_cloud, p3_cloud, k123, f, &
          & spline_o%n1,spline_o%x1pkg(1,1), &
          & spline_o%n2,spline_o%x2pkg(1,1), &
          & spline_o%n3,spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif


  if(ifail /= 0) ier = 96

  deallocate(p1_cloud, p2_cloud, p3_cloud, stat=ifail)
  if(ifail /= 0) then
     ier = 33
     return
  endif


end subroutine EZspline_derivative3_array_r4

subroutine EZspline_derivative3_cloud_r4(spline_o, i1, i2, i3, &
     & k, p1, p2, p3, f, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: i1, i2, i3, k
  real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)
  real(ezspline_r4), intent(out) :: f(k)
  integer, intent(out) :: ier

  integer ifail

  integer ict(10)
  integer ict_herm(8)
  integer:: iwarn=0

  ier = 0
  ifail=0
  if( .not.EZspline_allocated(spline_o) .or. spline_o%isReady /= 1) then
     ier = 94
     return
  endif

  f = 0.0_ezspline_r4
  if(i1<0 .OR. i2<0 .OR. i3<0) then
     ier = 12
     return
  endif
  if(max(i1,max(i2,i3)) > 3) then
     ier = 13
     return
  endif

  if(maxval((/i1,i2,i3/)) > 1 .AND. (spline_o%isHermite==1 .or. spline_o%isLinear==1)) then
     ier = 24; if(spline_o%isLinear==1) ier = 46
     return
  endif

  if(spline_o%isHybrid == 1) then

     call ezmake_ict3(i1,i2,i3,ict)

     call vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, spline_o%fspl(1,1,1,1), &
          & size(spline_o%fspl,1), size(spline_o%fspl,2), &
          & size(spline_o%fspl,3), size(spline_o%fspl,4), &
          & iwarn, ifail)

  else if(spline_o%isHermite==1 .or. spline_o%isLinear==1)then
     if(i1 == 1 .AND. i2 == 0 .AND. i3 == 0) then
        ict_herm = (/0, 1, 0, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 1, 0, 0, 0, 0, 0/)
     else if(i1 == 0 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 1, 0, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 0) then
        ict_herm = (/0, 0, 0, 0, 1, 0, 0, 0/)
     else if(i1 == 1 .AND. i2 == 0 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 1, 0, 0/)
     else if(i1 == 0 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 1, 0/)
     else if(i1 == 1 .AND. i2 == 1 .AND. i3 == 1) then
        ict_herm = (/0, 0, 0, 0, 0, 0, 0, 1/)
     else
        ict_herm = (/1, 0, 0, 0, 0, 0, 0, 0/)
     endif

     if(spline_o%isLinear == 0) then
        call vecherm3(ict_herm, k, p1, p2, p3, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     else
        call vecpc3(ict_herm, k, p1, p2, p3, &
             & k, f, &
             & spline_o%n1, spline_o%x1pkg(1,1), &
             & spline_o%n2, spline_o%x2pkg(1,1), &
             & spline_o%n3, spline_o%x3pkg(1,1), &
             & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
             & iwarn, ifail)
     endif

  else

     call ezmake_ict3(i1,i2,i3,ict)

     call vectricub(ict, k, p1, p2, p3, &
          & k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn, ifail)

  endif

  if(ifail /= 0) ier = 96

end subroutine EZspline_derivative3_cloud_r4
