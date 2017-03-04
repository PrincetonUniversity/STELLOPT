!/////
! R8 !
!/////
subroutine EZspline_init1_r8(spline_o, n1, BCS1, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: n1
  integer, intent(in) :: BCS1(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 2=wrong BCS1 code
  ! 3=wrong BCS2 code
  ! 99=something strange happened in EZspline_init
  integer, intent(out) :: ier
  integer i, iok
 
  call EZspline_preInit(spline_o)
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif

  spline_o%n1 = n1
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0 ! default is spline interpolation
  spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation

  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(2,n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  do i = 1, 2
 
     spline_o%bcval1min = 0.0_ezspline_r8
     spline_o%bcval1max = 0.0_ezspline_r8
     select case(BCS1(i))
     case (-1)
        spline_o%ibctype1(i) = -1
     case (0)
        spline_o%ibctype1(i) = 0
     case (1)
        spline_o%ibctype1(i) = 1
     case (2)
        spline_o%ibctype1(i) = 2
     case default
        ier = 2
        spline_o%ibctype1(i) = 0
     end select
 
  enddo
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
 
  !
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
  if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi_r8
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZspline_init1_r8
 
 
 
 
subroutine EZspline_init2_r8(spline_o, n1, n2, BCS1, BCS2, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: n1, n2
  integer, intent(in) :: BCS1(2), BCS2(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 2=wrong BCS1 code
  ! 3=wrong BCS2 code
  ! 99=something strange happened in EZspline_init
  integer, intent(out) :: ier
  integer i, iok
 
  call EZspline_preInit(spline_o)
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif

  spline_o%n1 = n1
  spline_o%n2 = n2
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0 ! default is spline interpolation
  spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(4,n1,n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  do i = 1, 2
 
     spline_o%bcval1min(1:n2) = 0.0_ezspline_r8
     spline_o%bcval1max(1:n2) = 0.0_ezspline_r8
     select case(BCS1(i))
     case (-1)
        spline_o%ibctype1(i) = -1
     case (0)
        spline_o%ibctype1(i) = 0
     case (1)
        spline_o%ibctype1(i) = 1
     case (2)
        spline_o%ibctype1(i) = 2
     case default
        ier = 2
        spline_o%ibctype1(i) = 0
     end select
 
     spline_o%bcval2min(1:n1) = 0.0_ezspline_r8
     spline_o%bcval2max(1:n1) = 0.0_ezspline_r8
     select case(BCS2(i))
     case (-1)
        spline_o%ibctype2(i) = -1
     case (0)
        spline_o%ibctype2(i) = 0
     case (1)
        spline_o%ibctype2(i) = 1
     case (2)
        spline_o%ibctype2(i) = 2
     case default
        ier = 3
        spline_o%ibctype2(i) = 0
     end select
 
  enddo
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
 
  !
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
  if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi_r8
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)
  spline_o%x2min = 0.0_ezspline_r8
  spline_o%x2max = 1.0_ezspline_r8
  if(BCS2(2)==-1) spline_o%x2max = ezspline_twopi_r8
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n2-1, ezspline_r8), i=1,spline_o%n2 ) /)
 
   spline_o%isReady = 0
 
  return
end subroutine EZspline_init2_r8
 
 
 
 
subroutine EZspline_init3_r8(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 2=wrong BCS1 code
  ! 3=wrong BCS2 code
  ! 4=wrong BCS3 code
  ! 99=something strange happened in EZspline_init
  integer, intent(out) :: ier
  integer i, iok
 
  call EZspline_preInit(spline_o)
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif

  spline_o%n1 = n1
  spline_o%n2 = n2
  spline_o%n3 = n3
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
  spline_o%klookup3 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0 ! default is spline interpolation
  spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3(n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(8,n1,n2,n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(n2, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(n2, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(n1, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(n1, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3min(n1, n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3max(n1, n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3pkg(n3, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  do i = 1, 2
 
     spline_o%bcval1min(1:n2, 1:n3) = 0.0_ezspline_r8
     spline_o%bcval1max(1:n2, 1:n3) = 0.0_ezspline_r8
     select case(BCS1(i))
     case (-1)
        spline_o%ibctype1(i) = -1
     case (0)
        spline_o%ibctype1(i) = 0
     case (1)
        spline_o%ibctype1(i) = 1
     case (2)
        spline_o%ibctype1(i) = 2
     case default
        ier = 2
        spline_o%ibctype1(i) = 0
     end select
 
     spline_o%bcval2min(1:n1, 1:n3) = 0.0_ezspline_r8
     spline_o%bcval2max(1:n1, 1:n3) = 0.0_ezspline_r8
     select case(BCS2(i))
     case (-1)
        spline_o%ibctype2(i) = -1
     case (0)
        spline_o%ibctype2(i) = 0
     case (1)
        spline_o%ibctype2(i) = 1
     case (2)
        spline_o%ibctype2(i) = 2
     case default
        ier = 3
        spline_o%ibctype2(i) = 0
     end select
 
     spline_o%bcval3min(1:n1, 1:n2) = 0.0_ezspline_r8
     spline_o%bcval3max(1:n1, 1:n2) = 0.0_ezspline_r8
     select case(BCS3(i))
     case (-1)
        spline_o%ibctype3(i) = -1
     case (0)
        spline_o%ibctype3(i) = 0
     case (1)
        spline_o%ibctype3(i) = 1
     case (2)
        spline_o%ibctype3(i) = 2
     case default
        ier = 4
        spline_o%ibctype3(i) = 0
     end select
 
  enddo
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     spline_o%ibctype3(1)=-1
     spline_o%ibctype3(2)=-1
  endif
 
 
  !
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
  if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi_r8
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)
  spline_o%x2min = 0.0_ezspline_r8
  spline_o%x2max = 1.0_ezspline_r8
  if(BCS2(2)==-1) spline_o%x2max = ezspline_twopi_r8
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n2-1, ezspline_r8), i=1,spline_o%n2 ) /)
 
  spline_o%x3min = 0.0_ezspline_r8
  spline_o%x3max = 1.0_ezspline_r8
  if(BCS3(2)==-1) spline_o%x3max = ezspline_twopi_r8
  spline_o%x3 = spline_o%x3min + (spline_o%x3max - spline_o%x3min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n3-1, ezspline_r8), i=1,spline_o%n3 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZspline_init3_r8
 
!/////
! R4 !
!/////
subroutine EZspline_init1_r4(spline_o, n1, BCS1, ier)
  use EZspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: n1
  integer, intent(in) :: BCS1(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 2=wrong BCS1 code
  ! 3=wrong BCS2 code
  ! 99=something strange happened in EZspline_init
  integer, intent(out) :: ier
  integer i, iok
 
  call EZspline_preInit(spline_o)
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif
 
  spline_o%n1 = n1
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0 ! default is spline interpolation
  spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(2,n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  do i = 1, 2
 
     spline_o%bcval1min = 0.0_ezspline_r4
     spline_o%bcval1max = 0.0_ezspline_r4
     select case(BCS1(i))
     case (-1)
        spline_o%ibctype1(i) = -1
     case (0)
        spline_o%ibctype1(i) = 0
     case (1)
        spline_o%ibctype1(i) = 1
     case (2)
        spline_o%ibctype1(i) = 2
     case default
        ier = 2
        spline_o%ibctype1(i) = 0
     end select
 
  enddo
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
 
  !
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi_r4
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZspline_init1_r4
 
 
 
 
subroutine EZspline_init2_r4(spline_o, n1, n2, BCS1, BCS2, ier)
  use EZspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: n1, n2
  integer, intent(in) :: BCS1(2), BCS2(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 2=wrong BCS1 code
  ! 3=wrong BCS2 code
  ! 99=something strange happened in EZspline_init
  integer, intent(out) :: ier
  integer i, iok
 
  call EZspline_preInit(spline_o)
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif
 
  spline_o%n1 = n1
  spline_o%n2 = n2
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0 ! default is spline interpolation
  spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(4,n1,n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  do i = 1, 2
 
     spline_o%bcval1min(1:n2) = 0.0_ezspline_r4
     spline_o%bcval1max(1:n2) = 0.0_ezspline_r4
     select case(BCS1(i))
     case (-1)
        spline_o%ibctype1(i) = -1
     case (0)
        spline_o%ibctype1(i) = 0
     case (1)
        spline_o%ibctype1(i) = 1
     case (2)
        spline_o%ibctype1(i) = 2
     case default
        ier = 2
        spline_o%ibctype1(i) = 0
     end select
 
     spline_o%bcval2min(1:n1) = 0.0_ezspline_r4
     spline_o%bcval2max(1:n1) = 0.0_ezspline_r4
     select case(BCS2(i))
     case (-1)
        spline_o%ibctype2(i) = -1
     case (0)
        spline_o%ibctype2(i) = 0
     case (1)
        spline_o%ibctype2(i) = 1
     case (2)
        spline_o%ibctype2(i) = 2
     case default
        ier = 3
        spline_o%ibctype2(i) = 0
     end select
 
  enddo
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
 
  !
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi_r4
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)
  spline_o%x2min = 0.0_ezspline_r4
  spline_o%x2max = 1.0_ezspline_r4
  if(BCS2(2)==-1) spline_o%x2max = ezspline_twopi_r4
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n2-1, ezspline_r4), i=1,spline_o%n2 ) /)
 
   spline_o%isReady = 0
 
  return
end subroutine EZspline_init2_r4
 
 
 
 
subroutine EZspline_init3_r4(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
  use EZspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 2=wrong BCS1 code
  ! 3=wrong BCS2 code
  ! 4=wrong BCS3 code
  ! 99=something strange happened in EZspline_init
  integer, intent(out) :: ier
  integer i, iok
 
  call EZspline_preInit(spline_o)
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif
 
  spline_o%n1 = n1
  spline_o%n2 = n2
  spline_o%n3 = n3
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
  spline_o%klookup3 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0 ! default is spline interpolation
  spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3(n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(8,n1,n2,n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(n2, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(n2, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(n1, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(n1, n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3min(n1, n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3max(n1, n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3pkg(n3, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  do i = 1, 2
 
     spline_o%bcval1min(1:n2, 1:n3) = 0.0_ezspline_r4
     spline_o%bcval1max(1:n2, 1:n3) = 0.0_ezspline_r4
     select case(BCS1(i))
     case (-1)
        spline_o%ibctype1(i) = -1
     case (0)
        spline_o%ibctype1(i) = 0
     case (1)
        spline_o%ibctype1(i) = 1
     case (2)
        spline_o%ibctype1(i) = 2
     case default
        ier = 2
        spline_o%ibctype1(i) = 0
     end select
 
     spline_o%bcval2min(1:n1, 1:n3) = 0.0_ezspline_r4
     spline_o%bcval2max(1:n1, 1:n3) = 0.0_ezspline_r4
     select case(BCS2(i))
     case (-1)
        spline_o%ibctype2(i) = -1
     case (0)
        spline_o%ibctype2(i) = 0
     case (1)
        spline_o%ibctype2(i) = 1
     case (2)
        spline_o%ibctype2(i) = 2
     case default
        ier = 3
        spline_o%ibctype2(i) = 0
     end select
 
     spline_o%bcval3min(1:n1, 1:n2) = 0.0_ezspline_r4
     spline_o%bcval3max(1:n1, 1:n2) = 0.0_ezspline_r4
     select case(BCS3(i))
     case (-1)
        spline_o%ibctype3(i) = -1
     case (0)
        spline_o%ibctype3(i) = 0
     case (1)
        spline_o%ibctype3(i) = 1
     case (2)
        spline_o%ibctype3(i) = 2
     case default
        ier = 4
        spline_o%ibctype3(i) = 0
     end select
 
  enddo
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     spline_o%ibctype3(1)=-1
     spline_o%ibctype3(2)=-1
  endif
 
 
  !
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi_r4
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)
  spline_o%x2min = 0.0_ezspline_r4
  spline_o%x2max = 1.0_ezspline_r4
  if(BCS2(2)==-1) spline_o%x2max = ezspline_twopi_r4
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n2-1, ezspline_r4), i=1,spline_o%n2 ) /)
 
  spline_o%x3min = 0.0_ezspline_r4
  spline_o%x3max = 1.0_ezspline_r4
  if(BCS3(2)==-1) spline_o%x3max = ezspline_twopi_r4
  spline_o%x3 = spline_o%x3min + (spline_o%x3max - spline_o%x3min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n3-1, ezspline_r4), i=1,spline_o%n3 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZspline_init3_r4
 
