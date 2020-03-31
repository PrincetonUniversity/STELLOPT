!
!-----------------------------------------------------------------------------------------
!
SUBROUTINE linear_interp(x,x_tables,f_tables,fx,ierr)
! 2 points linear interpolation formula
  implicit none
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  real*8, intent(in) :: x ! point where function needed
  real*8, intent(in) ::x_tables(2),f_tables(2) ! 2 points from tables F_i(x_i), i=1..2
  real*8, intent(out) :: fx !function value at x f(x)
  integer, intent(out) :: ierr
  real*8 hx

  ierr=0
  hx=(x_tables(2)-x_tables(1))
  if(hx.eq.0.0_R8) then
     ierr=1
     return
  endif
  fx = f_tables(1)+(x-x_tables(1))/hx*(f_tables(2)-f_tables(1))
  return
  
end SUBROUTINE linear_interp
!
!-----------------------------------------------------------------------------------------
!
SUBROUTINE linear_interp_2d(x,x_tables,f_tables,fx,ierr)
! 4 points 2D linear interpolation formula
  implicit none
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  real*8, intent(in) :: x(2) ! point where function needed
  real*8, intent(in) ::x_tables(4),f_tables(2,2) ! 2 points from tables F_i(x_i), i=1..2
  real*8, intent(out) :: fx !function value at x f(x)
  integer, intent(out) :: ierr
  real*8 hx1,hx2,pe,qn,pe1,qn1

  ierr=0
  hx1=(x_tables(2)-x_tables(1))
  hx2=(x_tables(4)-x_tables(3))
  if((hx1.eq.0.0_R8).or.(hx2.eq.0.0_R8)) then
     ierr=1
     return
  endif
  !pe=(e-e0)/(e1-e0)
  pe=(x(1)-x_tables(1))/hx1
  !qn=(n-n0)/(n1-n0)
  qn=(x(2)-x_tables(3))/hx2
  !1-pe=(e1-e)/(e1-e0)
  pe1=1.0_R8 - pe
  !1-qn=(n1-n)/(n1-n0)
  qn1=1.0_R8 - qn
  !four points linear extrapolation
  fx = pe1*qn1*f_tables(1,1)+pe*qn1*f_tables(2,1)+qn*pe1*f_tables(1,2)+&
       &pe*qn*f_tables(2,2)
  return
  
end SUBROUTINE linear_interp_2d
