module ezspline

  interface EZspline_init
     !
     ! Initialize and allocate memory. BCS1,2,3 determine the type of boundary
     ! conditions on either end of the x1,2,3 grids: eg (/0, 0/) for not-a-knot
     ! on the left and (/-1, -1/) periodic. Other BCs such as imposed slope or
     ! second derivative can also be applied by setting '1' or '2' respectively
     ! on either side. For instance (/1, 2/) for 1st derivative on the left and
     ! 2nd derivative on the right. The value of the 1st/2nd derivative must be
     ! set explicitely set through the bcval1,2,3min and bcval1,2,3max arrays.
     !
     subroutine EZspline_init3_r8(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
       use ezspline_obj
       type(EZspline3_r8) :: spline_o
       integer, intent(in) :: n1, n2, n3
       integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
       integer, intent(out) :: ier
     end subroutine EZspline_init3_r8

     subroutine EZspline_init2_r8(spline_o, n1, n2, BCS1, BCS2, ier)
       use ezspline_obj
       type(EZspline2_r8) :: spline_o
       integer, intent(in) :: n1, n2
       integer, intent(in) :: BCS1(2), BCS2(2)
       integer, intent(out) :: ier
     end subroutine EZspline_init2_r8
 
     subroutine EZspline_init1_r8(spline_o, n1, BCS1, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer, intent(in) :: n1
       integer, intent(in) :: BCS1(2)
       integer, intent(out) :: ier
     end subroutine EZspline_init1_r8
 
     subroutine EZspline_init3_r4(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: n1, n2, n3
       integer, intent(in) :: BCS1(2), BCS2(2), BCS3(2)
       integer, intent(out) :: ier
     end subroutine EZspline_init3_r4
 
     subroutine EZspline_init2_r4(spline_o, n1, n2, BCS1, BCS2, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: n1, n2
       integer, intent(in) :: BCS1(2), BCS2(2)
       integer, intent(out) :: ier
     end subroutine EZspline_init2_r4
 
     subroutine EZspline_init1_r4(spline_o, n1, BCS1, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer, intent(in) :: n1
       integer, intent(in) :: BCS1(2)
       integer, intent(out) :: ier
     end subroutine EZspline_init1_r4

  end interface


  interface EZlinear_init
     !
     ! Initialize and allocate memory for piecewise LINEAR interpolation
     ! object.  This is C0.  For C2 cubic spline or C1 Akima Hermite spline,
     ! please see EZspline_init.
     !
     ! No boundary conditions are needed for piecewise linear interpolation
     !
     subroutine EZlinear_init3_r8(spline_o, n1, n2, n3, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: n1, n2, n3
       integer, intent(out) :: ier
     end subroutine EZlinear_init3_r8
 
     subroutine EZlinear_init2_r8(spline_o, n1, n2, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: n1, n2
       integer, intent(out) :: ier
     end subroutine EZlinear_init2_r8
 
     subroutine EZlinear_init1_r8(spline_o, n1, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer, intent(in) :: n1
       integer, intent(out) :: ier
     end subroutine EZlinear_init1_r8
 
     subroutine EZlinear_init3_r4(spline_o, n1, n2, n3, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: n1, n2, n3
       integer, intent(out) :: ier
     end subroutine EZlinear_init3_r4
 
     subroutine EZlinear_init2_r4(spline_o, n1, n2, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: n1, n2
       integer, intent(out) :: ier
     end subroutine EZlinear_init2_r4
 
     subroutine EZlinear_init1_r4(spline_o, n1, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer, intent(in) :: n1
       integer, intent(out) :: ier
     end subroutine EZlinear_init1_r4
 
  end interface


  interface EZhybrid_init
     !
     ! Initialize and allocate memory for hybrid interpolation object:
     ! dimensionality > 1 only.  Interpolation method is specified separately
     ! for each dimension.  At present, Akima Hermite and Spline interpolation
     ! cannot be mixed.
     !
     ! Boundary condition arguments are optional.  They are appropriate only
     ! for the end points of dimensions for which Hermite or Spline cubic
     ! interpolation is used.
     !
     ! hspline(...) specifies the interpolation method for each dimension,
     ! according to the code: -1 for step function, 0 for piecewise linear,
     ! 1 for Akima Hermite, 2 for cubic Spline.
     !
     subroutine EZhybrid_init3_r8(spline_o, n1, n2, n3, hspline, ier, &
          BCS1, BCS2, BCS3)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: n1, n2, n3
       integer, intent(in) :: hspline(3)
       integer, intent(out) :: ier
       integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2), BCS3(2)
     end subroutine EZhybrid_init3_r8
 
     subroutine EZhybrid_init2_r8(spline_o, n1, n2, hspline, ier, &
          BCS1, BCS2)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: n1, n2
       integer, intent(in) :: hspline(2)
       integer, intent(out) :: ier
       integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2)
     end subroutine EZhybrid_init2_r8

     subroutine EZhybrid_init3_r4(spline_o, n1, n2, n3, hspline, ier, &
          BCS1, BCS2, BCS3)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: n1, n2, n3
       integer, intent(in) :: hspline(3)
       integer, intent(out) :: ier
       integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2), BCS3(2)
     end subroutine EZhybrid_init3_r4
 
     subroutine EZhybrid_init2_r4(spline_o, n1, n2, hspline, ier, &
          BCS1, BCS2)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: n1, n2
       integer, intent(in) :: hspline(2)
       integer, intent(out) :: ier
       integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2)
     end subroutine EZhybrid_init2_r4

  end interface
 
  interface EZspline_free
     !
     ! Reset and free the memory. This method must be called to avoid
     ! memory leaks after all interpolations have been computed.
     !
     subroutine EZspline_free3_r8(spline_o, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_free3_r8
 
     subroutine EZspline_free2_r8(spline_o, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_free2_r8
 
     subroutine EZspline_free1_r8(spline_o, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_free1_r8
 
     subroutine EZspline_free3_r4(spline_o, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_free3_r4
 
     subroutine EZspline_free2_r4(spline_o, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_free2_r4
 
     subroutine EZspline_free1_r4(spline_o, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_free1_r4
  end interface
 
  interface EZspline_setup
     !
     ! Compute the cubic spline coefficients. Note: the grid and the
     ! boundary conditions should be properly set prior to this call.
     !
     ! NEW optional argument: exact_dim=TRUE to requre f dimensions to
     ! match higher dimensions of spline_o%fspl exactly; default or FALSE
     ! means f dimensions can match or exceed dimensions of spline_o%fspl.
     !
     ! array arguments are now declared with f90 style dimensioning; the
     ! module interface must be used (if not feasible see ezspline_setupx.f90).
     ! 
     subroutine EZspline_setup3_r8(spline_o, f, ier, exact_dim)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       real(ezspline_r8), dimension(:,:,:), intent(in) :: f
       integer, intent(out) :: ier
       logical, intent(in), optional :: exact_dim
     end subroutine EZspline_setup3_r8
 
     subroutine EZspline_setup2_r8(spline_o, f, ier, exact_dim)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       real(ezspline_r8), dimension(:,:), intent(in) :: f
       integer, intent(out) :: ier
       logical, intent(in), optional :: exact_dim
     end subroutine EZspline_setup2_r8
 
     subroutine EZspline_setup1_r8(spline_o, f, ier, exact_dim)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       real(ezspline_r8), dimension(:), intent(in) :: f
       integer, intent(out) :: ier
       logical, intent(in), optional :: exact_dim
     end subroutine EZspline_setup1_r8
 
     subroutine EZspline_setup3_r4(spline_o, f, ier, exact_dim)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       real(ezspline_r4), dimension(:,:,:), intent(in) :: f
       integer, intent(out) :: ier
       logical, intent(in), optional :: exact_dim
     end subroutine EZspline_setup3_r4
 
     subroutine EZspline_setup2_r4(spline_o, f, ier, exact_dim)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       real(ezspline_r4), dimension(:,:), intent(in) :: f
       integer, intent(out) :: ier
       logical, intent(in), optional :: exact_dim
     end subroutine EZspline_setup2_r4
 
     subroutine EZspline_setup1_r4(spline_o, f, ier, exact_dim)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       real(ezspline_r4), dimension(:), intent(in) :: f
       integer, intent(out) :: ier
       logical, intent(in), optional :: exact_dim
     end subroutine EZspline_setup1_r4
 
  end interface
 
 
 
  interface EZspline_interp
     !
     ! Interpolation at grid point(s) p1, [p2, [p3]]. Result is returned in
     ! f. Interpolation can be sought at a single point (p1, [p2, [p3]] are
     ! scalars), on an unordered list of points (p1, [p2, [p3]] have dimension
     ! k), or on a structured grid (p1, [p2, [p3]] have dimension k1, [k2, [k3]]
     ! respectively).
     !
     subroutine EZspline_interp3_r8(spline_o, p1, p2, p3, f, ier)
       ! single point evaluation
       use ezspline_obj
       type(EZspline3_r8) spline_o
       real(ezspline_r8) :: p1, p2, p3
       real(ezspline_r8) f
       integer, intent(out) :: ier
     end subroutine EZspline_interp3_r8
 
     subroutine EZspline_interp3_array_r8(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer :: k1, k2, k3
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2), p3(k3)
       real(ezspline_r8), intent(out):: f(k1,k2,*)
       integer, intent(out) :: ier
     end subroutine EZspline_interp3_array_r8
 
     subroutine EZspline_interp3_cloud_r8(spline_o, k, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: k
       real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
       real(ezspline_r8), intent(out):: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_interp3_cloud_r8
 
     subroutine EZspline_interp2_r8(spline_o, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       real(ezspline_r8) :: p1, p2
       real(ezspline_r8) f
       integer, intent(out) :: ier
     end subroutine EZspline_interp2_r8
 
     subroutine EZspline_interp2_array_r8(spline_o, k1, k2, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer :: k1, k2
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2)
       real(ezspline_r8), intent(out):: f(k1,*)
       integer, intent(out) :: ier
     end subroutine EZspline_interp2_array_r8
 
     subroutine EZspline_interp2_cloud_r8(spline_o, k, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: k
       real(ezspline_r8), intent(in) :: p1(k), p2(k)
       real(ezspline_r8), intent(out):: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_interp2_cloud_r8
 
     subroutine EZspline_interp1_r8(spline_o, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       real(ezspline_r8) :: p1
       real(ezspline_r8) f
       integer, intent(out) :: ier
     end subroutine EZspline_interp1_r8
 
     subroutine EZspline_interp1_array_r8(spline_o, k1, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer :: k1
       real(ezspline_r8), intent(in) :: p1(k1)
       real(ezspline_r8), intent(out):: f(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_interp1_array_r8
 
     subroutine EZspline_interp3_r4(spline_o, p1, p2, p3, f, ier)
       ! single point evaluation
       use ezspline_obj
       type(EZspline3_r4) spline_o
       real(ezspline_r4) :: p1, p2, p3
       real(ezspline_r4) f
       integer, intent(out) :: ier
     end subroutine EZspline_interp3_r4
 
     subroutine EZspline_interp3_array_r4(spline_o, k1, k2, k3, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer :: k1, k2, k3
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2), p3(k3)
       real(ezspline_r4), intent(out):: f(k1,k2,*)
       integer, intent(out) :: ier
     end subroutine EZspline_interp3_array_r4
 
     subroutine EZspline_interp3_cloud_r4(spline_o, k, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: k
       real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)
       real(ezspline_r4), intent(out):: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_interp3_cloud_r4
 
     subroutine EZspline_interp2_r4(spline_o, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       real(ezspline_r4) :: p1, p2
       real(ezspline_r4) f
       integer, intent(out) :: ier
     end subroutine EZspline_interp2_r4
 
     subroutine EZspline_interp2_array_r4(spline_o, k1, k2, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer :: k1, k2
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2)
       real(ezspline_r4), intent(out):: f(k1,*)
       integer, intent(out) :: ier
     end subroutine EZspline_interp2_array_r4
 
     subroutine EZspline_interp2_cloud_r4(spline_o, k, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: k
       real(ezspline_r4), intent(in) :: p1(k), p2(k)
       real(ezspline_r4), intent(out):: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_interp2_cloud_r4
 
     subroutine EZspline_interp1_r4(spline_o, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       real(ezspline_r4) :: p1
       real(ezspline_r4) f
       integer, intent(out) :: ier
     end subroutine EZspline_interp1_r4
 
     subroutine EZspline_interp1_array_r4(spline_o, k1, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer :: k1
       real(ezspline_r4), intent(in) :: p1(k1)
       real(ezspline_r4), intent(out):: f(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_interp1_array_r4
 
  end interface
 
  interface EZspline_derivative
     !
     ! Evaluate the spline/Akima Hermite derivative of order
     ! d^{i1} d^{i2} d^{i3} f / d x1^{i1} d x2^{i2} d x2^{i2}
     ! at p1, [p2, [p3]]. The sum of i1+[i2+[i3]] should be <=2 for spline, or
     ! <=1 for Akima Hermite or Piecewise Linear. 
     ! Result is return in f. The evaluation can
     ! be sought at a single point (p1, [p2, [p3]] are scalars), on
     ! an unordered list of points (p1, [p2, [p3]] have dimension k), or
     ! on a structured grid (p1, [p2, [p3]] have dimension k1, [k2, [k3]]
     ! respectively).
     !
     !
     subroutine EZspline_derivative3_r8(spline_o, i1, i2, i3, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: i1, i2, i3
       real(ezspline_r8), intent(in) :: p1, p2, p3
       real(ezspline_r8), intent(out) :: f
       integer, intent(out) :: ier
     end subroutine EZspline_derivative3_r8
 
     subroutine EZspline_derivative3_array_r8(spline_o, i1, i2, i3, &
          & k1, k2, k3, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: i1, i2, i3, k1, k2, k3
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2), p3(k3)
       real(ezspline_r8), intent(out) :: f(k1, k2, *)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative3_array_r8
 
     subroutine EZspline_derivative3_cloud_r8(spline_o, i1, i2, i3, &
          & k, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: i1, i2, i3, k
       real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
       real(ezspline_r8), intent(out) :: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative3_cloud_r8
 
     subroutine EZspline_derivative2_r8(spline_o, i1, i2, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: i1, i2
       real(ezspline_r8), intent(in) :: p1, p2
       real(ezspline_r8), intent(out) :: f
       integer, intent(out) :: ier
     end subroutine EZspline_derivative2_r8
 
     subroutine EZspline_derivative2_array_r8(spline_o, i1, i2, k1, k2, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: i1, i2, k1, k2
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2)
       real(ezspline_r8), intent(out) :: f(k1,k2)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative2_array_r8
 
     subroutine EZspline_derivative2_cloud_r8(spline_o, i1, i2, k, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: i1, i2, k
       real(ezspline_r8), intent(in) :: p1(k), p2(k)
       real(ezspline_r8), intent(out) :: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative2_cloud_r8
 
     subroutine EZspline_derivative1_r8(spline_o, i1, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer, intent(in) :: i1
       real(ezspline_r8), intent(in) :: p1
       real(ezspline_r8), intent(out) :: f
       integer, intent(out) :: ier
     end subroutine EZspline_derivative1_r8
 
     subroutine EZspline_derivative1_array_r8(spline_o, i1, k1, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer, intent(in) :: i1, k1
       real(ezspline_r8), intent(in) :: p1(k1)
       real(ezspline_r8), intent(out) :: f(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative1_array_r8
 
     subroutine EZspline_derivative3_r4(spline_o, i1, i2, i3, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: i1, i2, i3
       real(ezspline_r4), intent(in) :: p1, p2, p3
       real(ezspline_r4), intent(out) :: f
       integer, intent(out) :: ier
     end subroutine EZspline_derivative3_r4
 
     subroutine EZspline_derivative3_array_r4(spline_o, i1, i2, i3, &
          & k1, k2, k3, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: i1, i2, i3, k1, k2, k3
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2), p3(k3)
       real(ezspline_r4), intent(out) :: f(k1, k2, *)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative3_array_r4
 
     subroutine EZspline_derivative3_cloud_r4(spline_o, i1, i2, i3, &
          & k, p1, p2, p3, f, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: i1, i2, i3, k
       real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)
       real(ezspline_r4), intent(out) :: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative3_cloud_r4
 
     subroutine EZspline_derivative2_r4(spline_o, i1, i2, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: i1, i2
       real(ezspline_r4), intent(in) :: p1, p2
       real(ezspline_r4), intent(out) :: f
       integer, intent(out) :: ier
     end subroutine EZspline_derivative2_r4
 
     subroutine EZspline_derivative2_array_r4(spline_o, i1, i2, k1, k2, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: i1, i2, k1, k2
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2)
       real(ezspline_r4), intent(out) :: f(k1,k2)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative2_array_r4
 
     subroutine EZspline_derivative2_cloud_r4(spline_o, i1, i2, k, p1, p2, f, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: i1, i2, k
       real(ezspline_r4), intent(in) :: p1(k), p2(k)
       real(ezspline_r4), intent(out) :: f(k)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative2_cloud_r4
 
     subroutine EZspline_derivative1_r4(spline_o, i1, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer, intent(in) :: i1
       real(ezspline_r4), intent(in) :: p1
       real(ezspline_r4), intent(out) :: f
       integer, intent(out) :: ier
     end subroutine EZspline_derivative1_r4
 
     subroutine EZspline_derivative1_array_r4(spline_o, i1, k1, p1, f, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer, intent(in) :: i1, k1
       real(ezspline_r4), intent(in) :: p1(k1)
       real(ezspline_r4), intent(out) :: f(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_derivative1_array_r4
 
  end interface
 
  interface EZspline_gradient
     !
     ! Return the gradient in df. When the dimensionality is 1 then
     ! df is df/dx. In more than one dimension, df has rank >=1 with
     ! the last index yielding df/dx, df/dy ... etc. Subsequent indices
     ! reflect the node positions when cloud or array evaluation is
     ! sought, as in df(k1, k2, k3, 1) for df/dx at x(k1), y(k2), z(k3).
     !
     subroutine EZspline_gradient3_r8(spline_o, p1, p2, p3, df, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       real(ezspline_r8), intent(in) :: p1, p2, p3
       real(ezspline_r8), intent(out) :: df(3)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient3_r8
 
     subroutine EZspline_gradient3_array_r8(spline_o, k1, k2, k3, &
          & p1, p2, p3, df, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: k1, k2, k3
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2), p3(k3)
       real(ezspline_r8), intent(out) :: df(k1,k2,k3,3)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient3_array_r8
 
     subroutine EZspline_gradient3_cloud_r8(spline_o, k, &
          & p1, p2, p3, df, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: k
       real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
       real(ezspline_r8), intent(out) :: df(k,3)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient3_cloud_r8
 
     subroutine EZspline_gradient2_r8(spline_o, p1, p2, df, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       real(ezspline_r8), intent(in) :: p1, p2
       real(ezspline_r8), intent(out) :: df(2)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient2_r8
 
     subroutine EZspline_gradient2_array_r8(spline_o, k1, k2, &
          & p1, p2, df, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in)  :: k1, k2
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2)
       real(ezspline_r8), intent(out) :: df(k1, k2, 2)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient2_array_r8
 
     subroutine EZspline_gradient2_cloud_r8(spline_o, k, &
          & p1, p2, df, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in)  :: k
       real(ezspline_r8), intent(in) :: p1(k), p2(k)
       real(ezspline_r8), intent(out) :: df(k, 2)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient2_cloud_r8
 
     subroutine EZspline_gradient1_r8(spline_o, p1, df, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       real(ezspline_r8), intent(in) :: p1
       real(ezspline_r8), intent(out) :: df
       integer, intent(out) :: ier
     end subroutine EZspline_gradient1_r8
 
     subroutine EZspline_gradient1_array_r8(spline_o, k1, p1, df, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       integer, intent(in) :: k1
       real(ezspline_r8), intent(in) :: p1(k1)
       real(ezspline_r8), intent(out) :: df(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient1_array_r8
 
     subroutine EZspline_gradient3_r4(spline_o, p1, p2, p3, df, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       real(ezspline_r4), intent(in) :: p1, p2, p3
       real(ezspline_r4), intent(out) :: df(3)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient3_r4
 
     subroutine EZspline_gradient3_array_r4(spline_o, k1, k2, k3, &
          & p1, p2, p3, df, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: k1, k2, k3
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2), p3(k3)
       real(ezspline_r4), intent(out) :: df(k1,k2,k3,3)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient3_array_r4
 
     subroutine EZspline_gradient3_cloud_r4(spline_o, k, &
          & p1, p2, p3, df, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: k
       real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)
       real(ezspline_r4), intent(out) :: df(k,3)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient3_cloud_r4
 
     subroutine EZspline_gradient2_r4(spline_o, p1, p2, df, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       real(ezspline_r4), intent(in) :: p1, p2
       real(ezspline_r4), intent(out) :: df(2)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient2_r4
 
     subroutine EZspline_gradient2_array_r4(spline_o, k1, k2, &
          & p1, p2, df, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in)  :: k1, k2
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2)
       real(ezspline_r4), intent(out) :: df(k1, k2, 2)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient2_array_r4
 
     subroutine EZspline_gradient2_cloud_r4(spline_o, k, &
          & p1, p2, df, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in)  :: k
       real(ezspline_r4), intent(in) :: p1(k), p2(k)
       real(ezspline_r4), intent(out) :: df(k, 2)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient2_cloud_r4
 
     subroutine EZspline_gradient1_r4(spline_o, p1, df, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       real(ezspline_r4), intent(in) :: p1
       real(ezspline_r4), intent(out) :: df
       integer, intent(out) :: ier
     end subroutine EZspline_gradient1_r4
 
     subroutine EZspline_gradient1_array_r4(spline_o, k1, p1, df, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       integer, intent(in) :: k1
       real(ezspline_r4), intent(in) :: p1(k1)
       real(ezspline_r4), intent(out) :: df(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_gradient1_array_r4
 
  end interface
 
  interface EZspline_isInDomain
     !
     ! Return error code if position (p1, [p2, [p3]]) is outside domain.
     ! The evaluation can be sought at a single point (p1, [p2, [p3]]
     ! are scalars), on an unordered list of points (p1, [p2, [p3]] have
     ! dimension k), or on a structured grid (p1, [p2, [p3]] have dimension
     ! k1, [k2, [k3]] respectively).
     !
     subroutine EZspline_isInDomain3_r8(spline_o, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r8) :: spline_o
       real(ezspline_r8), intent(in) :: p1, p2, p3
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain3_r8
 
     subroutine EZspline_isInDomain3_array_r8(spline_o, k1, k2, k3, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r8) :: spline_o
       integer, intent(in) :: k1, k2, k3
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2), p3(k3)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain3_array_r8
 
     subroutine EZspline_isInDomain3_cloud_r8(spline_o, k, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r8) :: spline_o
       integer, intent(in) :: k
       real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain3_cloud_r8
 
     subroutine EZspline_isInDomain2_r8(spline_o, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r8) :: spline_o
       real(ezspline_r8), intent(in) :: p1, p2
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain2_r8
 
     subroutine EZspline_isInDomain2_array_r8(spline_o, k1, k2, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r8) :: spline_o
       integer, intent(in) :: k1, k2
       real(ezspline_r8), intent(in) :: p1(k1), p2(k2)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain2_array_r8
 
     subroutine EZspline_isInDomain2_cloud_r8(spline_o, k, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r8) :: spline_o
       integer, intent(in) :: k
       real(ezspline_r8), intent(in) :: p1(k), p2(k)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain2_cloud_r8
 
     subroutine EZspline_isInDomain1_r8(spline_o, p1, ier)
       use ezspline_obj
       type(EZspline1_r8) :: spline_o
       real(ezspline_r8), intent(in) :: p1
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain1_r8
 
     subroutine EZspline_isInDomain1_array_r8(spline_o, k1, p1, ier)
       use ezspline_obj
       type(EZspline1_r8) :: spline_o
       integer, intent(in) :: k1
       real(ezspline_r8), intent(in) :: p1(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain1_array_r8
 
     subroutine EZspline_isInDomain3_r4(spline_o, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r4) :: spline_o
       real(ezspline_r4), intent(in) :: p1, p2, p3
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain3_r4
 
     subroutine EZspline_isInDomain3_array_r4(spline_o, k1, k2, k3, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r4) :: spline_o
       integer, intent(in) :: k1, k2, k3
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2), p3(k3)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain3_array_r4
 
     subroutine EZspline_isInDomain3_cloud_r4(spline_o, k, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r4) :: spline_o
       integer, intent(in) :: k
       real(ezspline_r4), intent(in) :: p1(k), p2(k), p3(k)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain3_cloud_r4
 
     subroutine EZspline_isInDomain2_r4(spline_o, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r4) :: spline_o
       real(ezspline_r4), intent(in) :: p1, p2
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain2_r4
 
     subroutine EZspline_isInDomain2_array_r4(spline_o, k1, k2, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r4) :: spline_o
       integer, intent(in) :: k1, k2
       real(ezspline_r4), intent(in) :: p1(k1), p2(k2)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain2_array_r4
 
     subroutine EZspline_isInDomain2_cloud_r4(spline_o, k, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r4) :: spline_o
       integer, intent(in) :: k
       real(ezspline_r4), intent(in) :: p1(k), p2(k)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain2_cloud_r4
 
     subroutine EZspline_isInDomain1_r4(spline_o, p1, ier)
       use ezspline_obj
       type(EZspline1_r4) :: spline_o
       real(ezspline_r4), intent(in) :: p1
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain1_r4
 
     subroutine EZspline_isInDomain1_array_r4(spline_o, k1, p1, ier)
       use ezspline_obj
       type(EZspline1_r4) :: spline_o
       integer, intent(in) :: k1
       real(ezspline_r4), intent(in) :: p1(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_isInDomain1_array_r4
 
  end interface
 
  interface EZspline_isGridRegular
     !
     ! Return error code if grid is not regular (not strictly  increasing).
     ! The evaluation can be sought at a single point (p1, [p2, [p3]]
     ! are scalars), on an unordered list of points (p1, [p2, [p3]] have
     ! dimension k), or on a structured grid (p1, [p2, [p3]] have dimension
     ! k1, [k2, [k3]] respectively).
     !
     subroutine EZspline_isGridRegular3_r8(spline_o, ier)
       use ezspline_obj
       type(EZspline3_r8) :: spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_isGridRegular3_r8
 
     subroutine EZspline_isGridRegular2_r8(spline_o, ier)
       use ezspline_obj
       type(EZspline2_r8) :: spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_isGridRegular2_r8
 
     subroutine EZspline_isGridRegular1_r8(spline_o, ier)
       use ezspline_obj
       type(EZspline1_r8) :: spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_isGridRegular1_r8
 
     subroutine EZspline_isGridRegular3_r4(spline_o, ier)
       use ezspline_obj
       type(EZspline3_r4) :: spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_isGridRegular3_r4
 
     subroutine EZspline_isGridRegular2_r4(spline_o, ier)
       use ezspline_obj
       type(EZspline2_r4) :: spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_isGridRegular2_r4
 
     subroutine EZspline_isGridRegular1_r4(spline_o, ier)
       use ezspline_obj
       type(EZspline1_r4) :: spline_o
       integer, intent(out) :: ier
     end subroutine EZspline_isGridRegular1_r4
 
  end interface
 
  interface EZspline_save
     !
     ! Save spline/Akima Hermite/Linear object in netcdf file 'filename'. Use
     ! EZspline_load to load spline/Akima Hermite/Linear object from netcdf 
     ! file.
     !
     ! Mod DMC March 2006 -- optionally, by giving each spline a name, 
     ! multiple splines can be saved in a single file.  Also, by specifying
     ! "fullsave=.TRUE." the spline coefficients can be saved as well in the
     ! file, so that they do not have to be recomputed at ezspline_load time.
     !
     ! When creating a single file with multiple spline objects, it is the
     ! user's responsibilitly to make sure that a different name is used for
     ! each spline that is to be saved.  Names can consist of upper or lower
     ! case letters, numerals, and "_", but must not start with a numeral.
     ! Imbedded blanks are not allowed, and the length of the name must be
     ! no more than 20 characters long-- to allow the names of the spline
     ! object elements to be appended.
     !
     subroutine EZspline_save3_r8(spline_o, filename, ier, &
          spl_name,fullsave)
       use ezspline_obj
       use ezcdf
       type(EZspline3_r8) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
       logical, intent(in), optional :: fullsave

     end subroutine EZspline_save3_r8
 
     subroutine EZspline_save2_r8(spline_o, filename, ier, &
          spl_name,fullsave)
       use ezspline_obj
       use ezcdf
       type(EZspline2_r8) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
       logical, intent(in), optional :: fullsave

     end subroutine EZspline_save2_r8
 
     subroutine EZspline_save1_r8(spline_o, filename, ier, &
          spl_name,fullsave)
       use ezspline_obj
       use ezcdf
       type(EZspline1_r8) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
       logical, intent(in), optional :: fullsave

     end subroutine EZspline_save1_r8
 
     subroutine EZspline_save3_r4(spline_o, filename, ier, &
          spl_name,fullsave)
       use ezspline_obj
       use ezcdf
       type(EZspline3_r4) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
       logical, intent(in), optional :: fullsave

     end subroutine EZspline_save3_r4
 
     subroutine EZspline_save2_r4(spline_o, filename, ier, &
          spl_name,fullsave)
       use ezspline_obj
       use ezcdf
       type(EZspline2_r4) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
       logical, intent(in), optional :: fullsave

     end subroutine EZspline_save2_r4
 
     subroutine EZspline_save1_r4(spline_o, filename, ier, &
          spl_name,fullsave)
       use ezspline_obj
       use ezcdf
       type(EZspline1_r4) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
       logical, intent(in), optional :: fullsave

     end subroutine EZspline_save1_r4
 
  end interface
 
  interface EZspline_load
     !
     ! Load spline/Akima Hermite object from netcdf file 'filename'. Use
     ! EZspline_save to save spline/Akima Hermite/Linear object in netcdf file.
     !
     ! MOD DMC March 2006-- a single NetCDF file can now contain multiple
     ! *named* spline objects.  If accessing such an object, the name must
     ! be supplied via the optional argument "spl_name".  If this is omitted,
     ! the default is to read the contents of the file which is presumed to
     ! consist of but a single spline object.

     ! Each call still opens and closes the file; if users find this 
     ! inefficient, an improvement to the control interface may be built.

     subroutine EZspline_load3_r8(spline_o, filename, ier, spl_name)
       use ezspline_obj
       use ezcdf
       type(EZspline3_r8) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
     end subroutine EZspline_load3_r8
 
     subroutine EZspline_load2_r8(spline_o, filename, ier, spl_name)
       use ezspline_obj
       use ezcdf
       type(EZspline2_r8) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
     end subroutine EZspline_load2_r8
 
     subroutine EZspline_load1_r8(spline_o, filename, ier, spl_name)
       use ezspline_obj
       use ezcdf
       type(EZspline1_r8) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
     end subroutine EZspline_load1_r8
 
     subroutine EZspline_load3_r4(spline_o, filename, ier, spl_name)
       use ezspline_obj
       use ezcdf
       type(EZspline3_r4) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
     end subroutine EZspline_load3_r4
 
     subroutine EZspline_load2_r4(spline_o, filename, ier, spl_name)
       use ezspline_obj
       use ezcdf
       type(EZspline2_r4) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
     end subroutine EZspline_load2_r4
 
     subroutine EZspline_load1_r4(spline_o, filename, ier, spl_name)
       use ezspline_obj
       use ezcdf
       type(EZspline1_r4) :: spline_o
       character*(*) :: filename
       integer, intent(out) :: ier

       character*(*), intent(in),optional :: spl_name
     end subroutine EZspline_load1_r4
  end interface
 
  interface EZspline_modulo
     !
     ! Map argument to (x1,2,3min, x1,2,3max) interval. This is useful to avoid
     ! an out-of-grid error when periodic boundary conditions are applied.
     ! The evaluation can be sought at a single point (p1, [p2, [p3]]
     ! are scalars), on an unordered list of points (p1, [p2, [p3]] have
     ! dimension k), or on a structured grid (p1, [p2, [p3]] have dimension
     ! k1, [k2, [k3]] respectively).
     !
     subroutine EZspline_modulo3_r8(spline_o, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       real(ezspline_r8) :: p1, p2, p3
       integer, intent(out) :: ier
     end subroutine EZspline_modulo3_r8
 
     subroutine EZspline_modulo_array3_r8(spline_o, k1, k2, k3, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: k1, k2, k3
       real(ezspline_r8) :: p1(k1), p2(k2), p3(k3)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_array3_r8
 
     subroutine EZspline_modulo_cloud3_r8(spline_o, k, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: k
       real(ezspline_r8) :: p1(k), p2(k), p3(k)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_cloud3_r8
 
     subroutine EZspline_modulo2_r8(spline_o, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       real(ezspline_r8) :: p1, p2
       integer, intent(out) :: ier
     end subroutine EZspline_modulo2_r8
 
     subroutine EZspline_modulo_array2_r8(spline_o, k1, k2, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: k1, k2
       real(ezspline_r8) :: p1(k1), p2(k2)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_array2_r8
 
     subroutine EZspline_modulo_cloud2_r8(spline_o, k, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r8) spline_o
       integer, intent(in) :: k
       real(ezspline_r8) :: p1(k), p2(k)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_cloud2_r8
 
     subroutine EZspline_modulo1_r8(spline_o, p1, ier)
       use ezspline_obj
       type(EZspline1_r8) spline_o
       real(ezspline_r8) :: p1
       integer, intent(out) :: ier
     end subroutine EZspline_modulo1_r8
 
     subroutine EZspline_modulo_array1_r8(spline_o, k1, p1, ier)
       use ezspline_obj
       type(EZspline3_r8) spline_o
       integer, intent(in) :: k1
       real(ezspline_r8) :: p1(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_array1_r8
 
     subroutine EZspline_modulo3_r4(spline_o, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       real(ezspline_r4) :: p1, p2, p3
       integer, intent(out) :: ier
     end subroutine EZspline_modulo3_r4
 
     subroutine EZspline_modulo_array3_r4(spline_o, k1, k2, k3, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: k1, k2, k3
       real(ezspline_r4) :: p1(k1), p2(k2), p3(k3)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_array3_r4
 
     subroutine EZspline_modulo_cloud3_r4(spline_o, k, p1, p2, p3, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: k
       real(ezspline_r4) :: p1(k), p2(k), p3(k)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_cloud3_r4
 
     subroutine EZspline_modulo2_r4(spline_o, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       real(ezspline_r4) :: p1, p2
       integer, intent(out) :: ier
     end subroutine EZspline_modulo2_r4
 
     subroutine EZspline_modulo_array2_r4(spline_o, k1, k2, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: k1, k2
       real(ezspline_r4) :: p1(k1), p2(k2)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_array2_r4
 
     subroutine EZspline_modulo_cloud2_r4(spline_o, k, p1, p2, ier)
       use ezspline_obj
       type(EZspline2_r4) spline_o
       integer, intent(in) :: k
       real(ezspline_r4) :: p1(k), p2(k)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_cloud2_r4
 
     subroutine EZspline_modulo1_r4(spline_o, p1, ier)
       use ezspline_obj
       type(EZspline1_r4) spline_o
       real(ezspline_r4) :: p1
       integer, intent(out) :: ier
     end subroutine EZspline_modulo1_r4
 
     subroutine EZspline_modulo_array1_r4(spline_o, k1, p1, ier)
       use ezspline_obj
       type(EZspline3_r4) spline_o
       integer, intent(in) :: k1
       real(ezspline_r4) :: p1(k1)
       integer, intent(out) :: ier
     end subroutine EZspline_modulo_array1_r4
 
  end interface
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EZspline object-less methods (first argument is NOT an EZspline1,2,3 type).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  interface EZspline_2NetCDF
     !
     ! Save data in netCDF file 'filename'. To save an EZspline1,2,3 object use
     ! EZspline_save method.
     !
     subroutine EZspline_2NetCDF_array3_r8(n1, n2, n3, x1, x2, x3, f, filename, ier)
       use ezspline_obj
       use ezcdf
       implicit none
       integer, intent(in) :: n1, n2, n3
       real(ezspline_r8), intent(in) :: x1(:), x2(:), x3(:), f(:, :, :)
       character*(*), intent(in) :: filename
       integer, intent(out) :: ier
     end subroutine EZspline_2NetCDF_array3_r8
 
     subroutine EZspline_2NetCDF_array2_r8(n1, n2, x1, x2, f, filename, ier)
       use ezspline_obj
       use ezcdf
       implicit none
       integer, intent(in) :: n1, n2
       real(ezspline_r8), intent(in) ::  x1(:), x2(:), f(:,:)
       character*(*), intent(in) :: filename
       integer, intent(out) :: ier
     end subroutine EZspline_2NetCDF_array2_r8
 
     subroutine EZspline_2NetCDF1_r8(n1, x1, f, filename, ier)
       use ezspline_obj
       use ezcdf
       implicit none
       integer, intent(in) :: n1
       real(ezspline_r8), intent(in) :: x1(:), f(:)
       character*(*), intent(in) :: filename
       integer, intent(out) :: ier
     end subroutine EZspline_2NetCDF1_r8
 
   subroutine EZspline_2NetCDF_cloud3_r8(n, x1, x2, x3, f, filename, ier)
     use ezspline_obj
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r8), intent(in) :: x1(:), x2(:), x3(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
   end subroutine EZspline_2NetCDF_cloud3_r8
 
   subroutine EZspline_2NetCDF_cloud2_r8(n, x1, x2, f, filename, ier)
     use ezspline_obj
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r8), intent(in) :: x1(:), x2(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
   end subroutine EZspline_2NetCDF_cloud2_r8
 
     subroutine EZspline_2NetCDF_array3_r4(n1, n2, n3, x1, x2, x3, f, filename, ier)
       use ezspline_obj
       use ezcdf
       implicit none
       integer, intent(in) :: n1, n2, n3
       real(ezspline_r4), intent(in) :: x1(:), x2(:), x3(:), f(:, :, :)
       character*(*), intent(in) :: filename
       integer, intent(out) :: ier
     end subroutine EZspline_2NetCDF_array3_r4
 
     subroutine EZspline_2NetCDF_array2_r4(n1, n2, x1, x2, f, filename, ier)
       use ezspline_obj
       use ezcdf
       implicit none
       integer, intent(in) :: n1, n2
       real(ezspline_r4), intent(in) ::  x1(:), x2(:), f(:,:)
       character*(*), intent(in) :: filename
       integer, intent(out) :: ier
     end subroutine EZspline_2NetCDF_array2_r4
 
     subroutine EZspline_2NetCDF1_r4(n1, x1, f, filename, ier)
       use ezspline_obj
       use ezcdf
       implicit none
       integer, intent(in) :: n1
       real(ezspline_r4), intent(in) :: x1(:), f(:)
       character*(*), intent(in) :: filename
       integer, intent(out) :: ier
     end subroutine EZspline_2NetCDF1_r4
 
   subroutine EZspline_2NetCDF_cloud3_r4(n, x1, x2, x3, f, filename, ier)
     use ezspline_obj
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r4), intent(in) :: x1(:), x2(:), x3(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
   end subroutine EZspline_2NetCDF_cloud3_r4
 
   subroutine EZspline_2NetCDF_cloud2_r4(n, x1, x2, f, filename, ier)
     use ezspline_obj
     use ezcdf
     implicit none
     integer, intent(in) :: n
     real(ezspline_r4), intent(in) :: x1(:), x2(:), f(:)
     character*(*), intent(in) :: filename
     integer, intent(out) :: ier
   end subroutine EZspline_2NetCDF_cloud2_r4
 
  end interface
 
 
contains
 
 
  subroutine EZspline_error(ier)
    !
    ! Error handling routine. Maps error ier code to a meaningful message.
    ! Note: does not abort nor stop if ier/=0.
    !
    implicit none
    integer, intent(in) :: ier
 
    if(ier == 0) return
    print*,'**EZspline** ERROR/WARNING #', ier,' occurred'
 
    select case(ier)
    case(1)
       print*,'**EZspline** allocation error'
    case(2)
       print*,'**EZspline** wrong BCS1 code'
    case(3)
       print*,'**EZspline** wrong BCS2 code'
    case(4)
       print*,'**EZspline** wrong BCS3 code'
    case(5)
       print*,'**EZspline** Que??'
    case(6)
       print*,'**EZspline** out of interval p1 < min(x1)'
    case(7)
       print*,'**EZspline** out of interval p1 > max(x1)'
    case(8)
       print*,'**EZspline** out of interval p2 < min(x2)'
    case(9)
       print*,'**EZspline** out of interval p2 > max(x2)'
    case(10)
       print*,'**EZspline** out of interval p3 < min(x3)'
    case(11)
       print*,'**EZspline** out of interval p3 > max(x3)'
    case(12)
       print*,'**EZspline** negative derivative order'
    case(13)
       print*,'**EZspline** derivative order too high'
    case(14)
       print*,'**EZspline** x1 grid is not strictly increasing'
    case(15)
       print*,'**EZspline** x2 grid is not strictly increasing'
    case(16)
       print*,'**EZspline** x3 grid is not strictly increasing'
    case(17)
       print*,'**EZspline** could not save spline object in file '
    case(18)
       print*,'**EZspline** memory allocation failure in coefficient setup'
 
    case(20)
       print*,'**EZspline** attempt to load spline object with wrong rank.'
    case(21)
       print*,'**EZspline** could not load spline object from file '
    case(22)
       print*,'**EZspline** loaded spline object from file but failed at coefficient set-up'
    case(23)
       print*,'**EZspline** failed to free spline object'
    case(24)
       print*,'**EZspline** 2nd order derivative not supported for Akima-Hermite (isHermite=1)'
    case(25)
       print*,'**EZspline** not supported for Akima-Hermite (isHermite=1)'
    case(26)
       print*,'**EZspline** memory allocation error in EZspline_interp'
    case(27)
       print*,'**EZspline** an error ocurred in genxpkg'
    case(28)
       print*,'**EZspline** memory allocation failure in ezspline_interp'
    case(29)
       print*,'**EZspline** memory deallocation failure in ezspline_interp'
    case(30)
       print*,'**EZspline** memory allocation error in EZspline_gradient'
    case(31)
       print*,'**EZspline** memory deallocation error in EZspline_gradient'
    case(32)
       print*,'**EZspline** memory allocation error in EZspline_derivative'
    case(33)
       print*,'**EZspline** memory deallocation error in EZspline_derivative'
    case(34)
       print*,'**EZspline** could not open netCDF file in EZspline_2netcdf'
    case(35)
       print*,'**EZspline** could not write into netCDF file in EZspline_2netcdf'
    case(36)
       print*,'**EZspline** could not read from netCDF file in EZspline_2netcdf'
    case(37)
       print*,'**EZspline** could not close netCDF file in EZspline_2netcdf'
    case(38)
       print*,'**EZspline** could not define variable (cdfDefVar) in EZspline_2netcdf'
    case(39)
       print*,'**EZspline** could not open netCDF file in EZspline_save'
    case(40)
       print*,'**EZspline** could not write into netCDF file in EZspline_save'
    case(41)
       print*,'**EZspline** could not close netCDF file in EZspline_save'
    case(42)
       print*,'**EZspline** could not define variable (cdfDefVar) in EZspline_save'
    case(43)
       print*,'**EZspline** could not open netCDF file in EZspline_load'
    case(44)
       print*,'**EZspline** could not read from netCDF file in EZspline_load'
    case(45)
       print*,'**EZspline** could not close netCDF file in EZspline_load'
    case(46)
       print*,'**EZspline** 2nd order derivative not supported for Piecewise Linear Interpolation (isLinear=1)'
    case(47)
       print*,'**EZspline** not supported for Piecewise Linear Interpolation (isLinear=1)'

    case(50)
       print*,'**EZspline** ezspline_save (optional) spline name is blank.'
    case(51)
       print*,'**EZspline** ezspline_save (optional) spline name too long (max 20 characters).'
    case(52)
       print*,'**EZspline** ezspline_save (optional) spline name contains'
       print*,'             imbedded blanks or other illegal characters.'
    case(53)
       print*,'**EZspline** attempt to write named spline object to NetCDF'
       print*,'             file with change of dimensionality or data type.'
       
    case(54)
       print*,'**EZspline** hybrid interpolation specification not in range -1:2'
       print*,'             error in EZhybrid_init.'
    case(55)
       print*,'**EZspline** hybrid interpolation cannot mix Hermite and Spline interpolation.'
       print*,'             hspline(i)=1 and hspline(j)=2 in EZhybrid_init.'
    case(56)
       print*,'**EZspline** non-default boundary condition unavailable: zonal or piecewise linear dimension.'
       print*,'             in EZhybrid_init.'

    case(57)
       print*,'**EZspline** dimension of "f" smaller than corresponding "fspl"'
       print*,'             dimension in "spline_o".'
    case(58)
       print*,'**EZspline** dimension of "f" larger than corresponding "fspl"'
       print*,'             dimension in "spline_o".'

    case(90)
       print*,'**EZspline** an error occurred after attempting to evaluate the'
       print*,'             Hermite polynomials'
    case(91)
       print*,'**EZspline** an error occurred after attempting to set up the'
       print*,'             Hermite polynomial coefficients'
    case(92)
       print*,'**EZspline** warning in EZspline_load. Looks like saved object '
       print*,'             was not properly set-up (isReady=0).'
    case(93)
       print*,'**EZspline** warning in EZspline_save. Looks like saved object '
       print*,'             was not properly set-up (isReady=0).'
    case(94)
       print*,'**EZspline** an error occurred in EZspline_interp. Did you forget'
       print*,'             to set up the cubic spline coefficients by calling'
       print*,'             call EZspline_setup(spline_o, f, ier)'
       print*,'             ?'
    case(95)
       print*,'**EZspline** some error occurred in EZspline_gradient'
    case(96)
       print*,'**EZspline** some error occurred in EZspline_derivative'
    case(97)
       print*,'**EZspline** some error occurred in EZspline_interp apparently'
       print*,'             related to a PSPLINE routine. Check if argument is '
       print*,'             outside interpolation domain by calling'
       print*,'             call EZspline_isInDomain(spline_o, [[k1, k2, k3,] .OR. k,] p1, p2, p3, ier ,ier)'
       print*,'             call EZspline_error(ier)'
    case(98)
       print*,'**EZspline** error occurred in EZspline_setup'
       print*,'  if no other explanation-- ezspline_init call never made.'
    case(99)
       print*,'**EZspline** some error occurred in EZspline_init, EZhybrid_init,  or EZlinear_init'
    case(100)
       print*,'**EZSPLINE** EZspline_init, EZhybrid_init,  or EZlinear_init -- object already allocated.'
    case(101)
       print*,'**EZSPLINE** object was never allocated.'
    case default
       print*,'**EZspline** '
    end select
 
    return
  end subroutine EZspline_error

end module EZspline
