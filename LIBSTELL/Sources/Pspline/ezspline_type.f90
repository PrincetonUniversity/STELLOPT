module EZspline_type
 
  integer, parameter  :: ezspline_r8 = selected_real_kind(12,100)
  integer, parameter  :: ezspline_r4 = selected_real_kind(6,37)
  real(ezspline_r8), parameter :: ezspline_twopi_r8 = 6.2831853071795865_ezspline_r8
  real(ezspline_r4), parameter :: ezspline_twopi_r4 = 6.2831853071795865_ezspline_r4
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EZspline data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  type EZspline3_r8
     !
     ! 3-d Spline/Akima Hermite/Piecewise Linear interpolations
     !
     ! Grid
     !
     real(ezspline_r8), dimension(:), allocatable :: x1, x2, x3
     !
     ! The boundary condition values (for slope and 2nd derivative).
     ! Can be optionally set by the user. Not used for periodic and
     ! not a knot boundary conditions.
     !
     real(ezspline_r8), dimension(:,:), allocatable :: bcval1min, bcval1max
     real(ezspline_r8), dimension(:,:), allocatable :: bcval2min, bcval2max
     real(ezspline_r8), dimension(:,:), allocatable :: bcval3min, bcval3max
     !
     ! Select between spline (0) and Akima spline (1); default=0 (spline)
     !
     integer :: isHermite  ! set after EZspline_init call...
     !
     ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
     ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
     !
     integer :: isLinear
     !
     ! set =0 by init routines other than EZhybrid_init which sets it =1:
     integer :: isHybrid
     !
     ! the following is set by EZhybrid_init; other EZ*_init routines clear:
     integer :: hspline(3)  ! interpolation code along each dimension
     !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
     !
     ! Grid sizes (set during EZ*_init call).
     !
     integer :: n1, n2, n3
     !
     ! Grid zone lookup method
     !
     integer :: klookup1,klookup2,klookup3
     !
     ! Type of boundary conditions (set during EZspline_init call) on left
     ! and right hand side. Possible values are:
     !
     ! -1 periodic
     ! 0 not a knot
     ! +1 1st derivative imposed
     ! +2 2nd derivative imposed
     !
     ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
     ! and not a knot boundary conditions on right-hand side. The values of
     ! the derivatives a set via  bcval1min. (See above.)
     !
     integer ibctype1(2), ibctype2(2), ibctype3(2)
     !
     ! Grid lengths. DO NOT SET.
     !
     real(ezspline_r8) :: x1min, x1max, x2min, x2max, x3min, x3max
     !
     ! Compact cubic coefficient arrays. DO NOT SET.
     !
     real(ezspline_r8), dimension(:,:,:,:), allocatable :: fspl
     !
     ! Control/Other. DO NOT SET.
     !
     integer :: isReady

     integer :: ilin1, ilin2, ilin3
     real(ezspline_r8), dimension(:,:), allocatable :: x1pkg, x2pkg, x3pkg
     !
     integer :: nguard
  end type EZspline3_r8

  type EZspline2_r8
     !
     ! 2-d Spline/Akima Hermite/Piecewise Linear interpolation
     !
     ! Grid
     !
     real(ezspline_r8), dimension(:), allocatable :: x1, x2
     !
     ! The boundary condition values (for slope and 2nd derivative).
     ! Can be optionally set by the user. Not used for periodic and
     ! not a knot boundary conditions.
     !
     real(ezspline_r8), dimension(:), allocatable :: bcval1min, bcval1max
     real(ezspline_r8), dimension(:), allocatable :: bcval2min, bcval2max
     !
     ! Select between spline (0) and Akima spline (1); default=0 (spline)
     !
     integer :: isHermite  ! set after EZspline_init call...
     !
     ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
     ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
     !
     integer :: isLinear
     !
     ! set =0 by init routines other than EZhybrid_init which sets it =1:
     integer :: isHybrid
     !
     ! the following is set by EZhybrid_init; other EZ*_init routines clear:
     integer :: hspline(2)  ! interpolation code along each dimension
     !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
     !
     ! Grid sizes (set during EZ*_init call).
     !
     integer :: n1, n2
     !
     ! Grid zone lookup method
     !
     integer :: klookup1,klookup2
     !
     ! Type of boundary conditions (set during EZspline_init call) on left
     ! and right hand side. Possible values are:
     !
     ! -1 periodic
     ! 0 not a knot
     ! +1 1st derivative imposed
     ! +2 2nd derivative imposed
     !
     ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
     ! and not a knot boundary conditions on right-hand side. The values of
     ! the derivatives are set via  bcval1min. (See above)
     !
     integer ibctype1(2), ibctype2(2)
     !
     ! Grid lengths. DO NOT SET.
     !
     real(ezspline_r8) :: x1min, x1max, x2min, x2max
     !
     ! Compact cubic coefficient arrays. DO NOT SET.
     !
     real(ezspline_r8), dimension(:,:,:), allocatable :: fspl
     !
     ! Control/Other. DO NOT SET.
     !
     integer :: isReady

     integer :: ilin1, ilin2
     real(ezspline_r8), dimension(:,:), allocatable :: x1pkg, x2pkg
     !
     integer :: nguard
  end type EZspline2_r8
 
  type EZspline1_r8
     !
     ! 1-d Spline/Akima Hermite/Piecewise Linear interpolation
     !
     ! Grid
     !
     real(ezspline_r8), dimension(:), allocatable :: x1
     !
     ! The boundary condition values (for slope and 2nd derivative).
     ! Can be optionally set by the user. Not used for periodic and
     ! not a knot boundary conditions.
     !
     real(ezspline_r8) :: bcval1min, bcval1max
     !
     ! Select between spline (0) and Akima spline (1); default=0 (spline)
     !
     integer :: isHermite  ! set after EZspline_init call...
     !
     ! set =0 for Spline or Akima; =1 for piecewise linear: this is set
     ! by EZspline_init or EZlinear_init; DO NOT SET DIRECTLY:
     !
     integer :: isLinear
     !
     ! Grid sizes (set during EZ*_init call).
     !
     integer :: n1
     !
     ! Grid zone lookup method
     !
     integer :: klookup1
     !
     ! Type of boundary conditions (set during EZspline_init call) on left
     ! and right hand side. Possible values are:
     !
     ! -1 periodic
     ! 0 not a knot
     ! +1 1st derivative imposed
     ! +2 2nd derivative imposed
     !
     ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
     ! and not a knot boundary conditions on right-hand side. The values of
     ! the derivatives are set via  bcval1min. (See above)
     !
     integer ibctype1(2)
     !
     ! Grid lengths. DO NOT SET.
     !
     real(ezspline_r8) :: x1min, x1max
     !
     ! Compact cubic coefficient arrays. DO NOT SET.
     !
     real(ezspline_r8), dimension(:,:), allocatable :: fspl
     !
     ! Control/Other. DO NOT SET.
     !
     integer :: isReady

     integer :: ilin1
     real(ezspline_r8), dimension(:,:), allocatable :: x1pkg
     !
     integer :: nguard
  end type EZspline1_r8

  type EZspline3_r4
     !
     ! 3-d Spline/Akima Hermite/Piecewise Linear interpolation
     !
     ! Grid
     !
     real(ezspline_r4), dimension(:), allocatable :: x1, x2, x3
     !
     ! The boundary condition values (for slope and 2nd derivative).
     ! Can be optionally set by the user. Not used for periodic and
     ! not a knot boundary conditions.
     !
     real(ezspline_r4), dimension(:,:), allocatable :: bcval1min, bcval1max
     real(ezspline_r4), dimension(:,:), allocatable :: bcval2min, bcval2max
     real(ezspline_r4), dimension(:,:), allocatable :: bcval3min, bcval3max
     !
     ! Select between spline (0) and Akima spline (1); default=0 (spline)
     !
     integer :: isHermite  ! set after EZspline_init call...
     !
     ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
     ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
     !
     integer :: isLinear
     !
     ! set =0 by init routines other than EZhybrid_init which sets it =1:
     integer :: isHybrid
     !
     ! the following is set by EZhybrid_init; other EZ*_init routines clear:
     integer :: hspline(3)  ! interpolation code along each dimension
     !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
     !
     ! Grid sizes (set during EZ*_init call).
     !
     integer :: n1, n2, n3
     !
     ! Grid zone lookup method
     !
     integer :: klookup1,klookup2,klookup3
     !
     ! Type of boundary conditions (set during EZspline_init call) on left
     ! and right hand side. Possible values are:
     !
     ! -1 periodic
     ! 0 not a knot
     ! +1 1st derivative imposed
     ! +2 2nd derivative imposed
     !
     ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
     ! and not a knot boundary conditions on right-hand side. The values of
     ! the derivatives a set via  bcval1min. (See above.)
     !
     integer ibctype1(2), ibctype2(2), ibctype3(2)
     !
     ! Grid lengths. DO NOT SET.
     !
     real(ezspline_r4) :: x1min, x1max, x2min, x2max, x3min, x3max
     !
     ! Compact cubic coefficient arrays. DO NOT SET.
     !
     real(ezspline_r4), dimension(:,:,:,:), allocatable :: fspl
     !
     ! Control/Other. DO NOT SET.
     !
     integer :: isReady

     integer :: ilin1, ilin2, ilin3
     real(ezspline_r4), dimension(:,:), allocatable :: x1pkg, x2pkg, x3pkg
     !
     integer :: nguard
  end type EZspline3_r4

  type EZspline2_r4
     !
     ! 2-d Spline/Akima Hermite/Piecewise Linear interpolation
     !
     ! Grid
     !
     real(ezspline_r4), dimension(:), allocatable :: x1, x2
     !
     ! The boundary condition values (for slope and 2nd derivative).
     ! Can be optionally set by the user. Not used for periodic and
     ! not a knot boundary conditions.
     !
     real(ezspline_r4), dimension(:), allocatable :: bcval1min, bcval1max
     real(ezspline_r4), dimension(:), allocatable :: bcval2min, bcval2max
     !
     ! Select between spline (0) and Akima spline (1); default=0 (spline)
     !
     integer :: isHermite  ! set after EZspline_init call...
     !
     ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
     ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
     !
     integer :: isLinear
     !
     ! set =0 by init routines other than EZhybrid_init which sets it =1:
     integer :: isHybrid
     !
     ! the following is set by EZhybrid_init; other EZ*_init routines clear:
     integer :: hspline(2)  ! interpolation code along each dimension
     !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
     !
     ! Grid sizes (set during EZ*_init call).
     !
     integer :: n1, n2
     !
     ! Grid zone lookup method
     !
     integer :: klookup1,klookup2
     !
     ! Type of boundary conditions (set during EZspline_init call) on left
     ! and right hand side. Possible values are:
     !
     ! -1 periodic
     ! 0 not a knot
     ! +1 1st derivative imposed
     ! +2 2nd derivative imposed
     !
     ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
     ! and not a knot boundary conditions on right-hand side. The values of
     ! the derivatives are set via  bcval1min. (See above)
     !
     integer ibctype1(2), ibctype2(2)
     !
     ! Grid lengths. DO NOT SET.
     !
     real(ezspline_r4) :: x1min, x1max, x2min, x2max
     !
     ! Compact cubic coefficient arrays. DO NOT SET.
     !
     real(ezspline_r4), dimension(:,:,:), allocatable :: fspl
     !
     ! Control/Other. DO NOT SET.
     !
     integer :: isReady

     integer :: ilin1, ilin2
     real(ezspline_r4), dimension(:,:), allocatable :: x1pkg, x2pkg
     !
     integer :: nguard
  end type EZspline2_r4
 
  type EZspline1_r4
     !
     ! 1-d Spline/Akima Hermite/Piecewise Linear interpolation
     !
     ! Grid
     !
     real(ezspline_r4), dimension(:), allocatable :: x1
     !
     ! The boundary condition values (for slope and 2nd derivative).
     ! Can be optionally set by the user. Not used for periodic and
     ! not a knot boundary conditions.
     !
     real(ezspline_r4) :: bcval1min, bcval1max
     !
     ! Select between spline (0) and Akima spline (1); default=0 (spline)
     !
     integer :: isHermite  ! set after EZspline_init call...
     !
     ! set =0 for Spline or Akima; =1 for piecewise linear: this is set
     ! by EZspline_init or EZlinear_init; DO NOT SET DIRECTLY:
     !
     integer :: isLinear
     !
     ! Grid sizes (set during EZ*_init call).
     !
     integer :: n1
     !
     ! Grid zone lookup method
     !
     integer :: klookup1
     !
     ! Type of boundary conditions (set during EZspline_init call) on left
     ! and right hand side. Possible values are:
     !
     ! -1 periodic
     ! 0 not a knot
     ! +1 1st derivative imposed
     ! +2 2nd derivative imposed
     !
     ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
     ! and not a knot boundary conditions on right-hand side. The values of
     ! the derivatives are set via  bcval1min. (See above)
     !
     integer ibctype1(2)
     !
     ! Grid lengths. DO NOT SET.
     !
     real(ezspline_r4) :: x1min, x1max
     !
     ! Compact cubic coefficient arrays. DO NOT SET.
     !
     real(ezspline_r4), dimension(:,:), allocatable :: fspl
     !
     ! Control/Other. DO NOT SET.
     !
     integer :: isReady

     integer :: ilin1
     real(ezspline_r4), dimension(:,:), allocatable :: x1pkg
     !
     integer :: nguard
  end type EZspline1_r4

!=========================================================================
! End type
end module EZspline_type
