!/////
! R8 !
!/////
subroutine EZlinear_init1_r8(spline_o, n1, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  integer, intent(in) :: n1
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZlinear_init
  integer, intent(out) :: ier
  integer i, iok
 
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif

  spline_o%n1 = n1
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection

  spline_o%isHermite = 0 ! no Hermite interpolation; this is EZlinear_init
  spline_o%isLinear  = 1 ! piecewise linear interpolation
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(1,n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZlinear_init1_r8
 
 
 
subroutine EZlinear_init2_r8(spline_o, n1, n2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: n1, n2
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZlinear_init
  integer, intent(out) :: ier
  integer i, iok
 
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
 
  spline_o%isHermite = 0 ! no Hermite interpolation; this is EZlinear_init
  spline_o%isLinear  = 1 ! piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(1,n1,n2), stat=iok)
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
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)
  spline_o%x2min = 0.0_ezspline_r8
  spline_o%x2max = 1.0_ezspline_r8
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n2-1, ezspline_r8), i=1,spline_o%n2 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZlinear_init2_r8
 
 
 
 
subroutine EZlinear_init3_r8(spline_o, n1, n2, n3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: n1, n2, n3
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZlinear_init
  integer, intent(out) :: ier
  integer i, iok
 
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
 
  spline_o%isHermite = 0 ! no Hermite interpolation; this is EZlinear_init
  spline_o%isLinear  = 1 ! piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3(n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(1,n1,n2,n3), stat=iok)
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
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)

  spline_o%x2min = 0.0_ezspline_r8
  spline_o%x2max = 1.0_ezspline_r8
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n2-1, ezspline_r8), i=1,spline_o%n2 ) /)
 
  spline_o%x3min = 0.0_ezspline_r8
  spline_o%x3max = 1.0_ezspline_r8
  spline_o%x3 = spline_o%x3min + (spline_o%x3max - spline_o%x3min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n3-1, ezspline_r8), i=1,spline_o%n3 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZlinear_init3_r8
 
!/////
! R4 !
!/////
subroutine EZlinear_init1_r4(spline_o, n1, ier)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  integer, intent(in) :: n1
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZlinear_init
  integer, intent(out) :: ier
  integer i, iok
 
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif
 
  spline_o%n1 = n1
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
 
  spline_o%isHermite = 0 ! no Hermite interpolation; this is EZlinear_init
  spline_o%isLinear  = 1 ! piecewise linear interpolation
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(1,n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZlinear_init1_r4
 
 
 
 
subroutine EZlinear_init2_r4(spline_o, n1, n2, ier)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: n1, n2
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZlinear_init
  integer, intent(out) :: ier
  integer i, iok
 
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
 
  spline_o%isHermite = 0 ! no Hermite interpolation; this is EZlinear_init
  spline_o%isLinear  = 1 ! piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(1,n1,n2), stat=iok)
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
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)

  spline_o%x2min = 0.0_ezspline_r4
  spline_o%x2max = 1.0_ezspline_r4
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n2-1, ezspline_r4), i=1,spline_o%n2 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZlinear_init2_r4
 
 
 
 
subroutine EZlinear_init3_r4(spline_o, n1, n2, n3, ier)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: n1, n2, n3
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZlinear_init
  integer, intent(out) :: ier
  integer i, iok
 
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
  
  spline_o%isHermite = 0 ! no Hermite interpolation; this is EZlinear_init
  spline_o%isLinear  = 1 ! piecewise linear interpolation
  spline_o%isHybrid  = 0
  spline_o%hspline = 0
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3(n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(1,n1,n2,n3), stat=iok)
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
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)

  spline_o%x2min = 0.0_ezspline_r4
  spline_o%x2max = 1.0_ezspline_r4
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n2-1, ezspline_r4), i=1,spline_o%n2 ) /)
 
  spline_o%x3min = 0.0_ezspline_r4
  spline_o%x3max = 1.0_ezspline_r4
  spline_o%x3 = spline_o%x3min + (spline_o%x3max - spline_o%x3min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n3-1, ezspline_r4), i=1,spline_o%n3 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZlinear_init3_r4
 
