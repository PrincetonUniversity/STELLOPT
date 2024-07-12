!-----------------------------------------------------------------------
!     Subroutine:    drkhvg
!     Authors:       Y. Suzuki
!     Date:          12/28/2002
!     Description:   Sixth order 8-stage Runge-Kutta-Huta forumla.
!
! $Revision: 118 $
! $Id: drkhvg.f90 119 2011-12-13 01:11:13Z suzuki $
!          
!    This subroutine is written by Y. Suzuki 
!        at Graduate School of Energy Science (Kyoto Univ)
!         2002/12/28
!
!    Based program is wrtten by K. Hamamatsu  
!        at Faculty of Science (Hiroshima Univ.)
!          1980/12/18
!
!  Runge-Kutta-Huta Formulas  ( Sixth order 8-stage )
!
!  " Improved Sixth-order Runge-kutta formulas and Approximate
!    Continuous Solution of Ordinary Differential Equation "
!  by D. Sarafyan:  J. Math. Anal. Appl. 40, 436-455 (1972)
!  
!  Ported by S. Lazerson (lazerson@pppl.gov) 04/25/12
!
!-----------------------------------------------------------------------
subroutine drkhvg ( x0, y0, ln, h, lnx, fun, & !(in)
  &                 y,                       & !(inout)
  &                 icount, iout             & !(out)
  &               )

!-----------------------------------------------------------------------
!     Libraries (none)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Parameters
!
!  y0(n)    real*8       initial value of y  ( start point )
!  n        integer*4    dimension of defferential equation
!  h        real*8       step size
!  f(n,nl)  real*8       work aray for calculation ( nl >= 9 )
!  y(n)     real*8       end point  ( if j=1, y(n)=y0(n) )
!  nx       integer*4    final value of x is h*(nx-1) ( end point )
!  icount   integer*4    end step
!  iout     integer*4    flag of calculation  
!          
!-----------------------------------------------------------------------
  implicit none
  integer, parameter :: ikind =  4,    &
    &                   rkind =  8,    &
    &                   dp    =  rkind
  integer(kind=ikind), intent(in) :: ln,  &
    &                                lnx
  integer(kind=ikind), intent(out) :: icount, &
    &                                 iout
  real(kind=rkind), intent(in) :: x0,    &
    &                             h,     &
    &                             y0(ln)
  real(kind=rkind), intent(out) :: y(ln,lnx)
!-----------------------------------------------------------------------
!     External Functions
!          fun           RHS of ODE
!-----------------------------------------------------------------------
  external fun
!-----------------------------------------------------------------------
!     Local Variables
!          i            Dummy Index
!-----------------------------------------------------------------------
  integer(kind=ikind) :: i
  real(kind=rkind) :: x,   &
    &                 x1,  &
    &                 c01, &
    &                 c02, &
    &                 c03, &
    &                 c04, &
    &                 c05, &
    &                 c06, &
    &                 c07, &
    &                 c08, &
    &                 c09, &
    &                 c10, &
    &                 c11, &
    &                 c12, &
    &                 c13, &
    &                 c14, &
    &                 c15, &
    &                 c16, &
    &                 c17, &
    &                 c18, &
    &                 c19, &
    &                 c20, &
    &                 c21, &
    &                 c22, &
    &                 c23, &
    &                 c24, &
    &                 c25, &
    &                 c26, &
    &                 c27
  real(kind=rkind), allocatable :: f(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------


  if(ln > 0)then
     allocate(f(ln,9))
  end if

  c01 =  h / 9.0_rkind
  c02 =  h * 0.4166666666666667e-01_rkind
  c03 =  h * 0.125_rkind
  c04 =  h / 6.0_rkind
  c05 =  h * 0.5_rkind
  c06 =  h * 0.6666666666666667_rkind
  c07 =  h / 3.0_rkind
  c08 =  h * 0.375_rkind
  c09 =  h * 0.1333333333333333e+01_rkind
  c10 =  h * 0.3333333333333333e+01_rkind
  c11 =  h * 0.7e+01_rkind
  c12 =  h * 0.9666666666666667e+01_rkind
  c13 =  h * 0.1533333333333333e+02_rkind
  c14 =  h * 0.6111111111111111_rkind
  c15 =  h * 0.1166666666666667e+01_rkind
  c16 =  h * 0.1375e+01_rkind
  c17 =  h * 0.8333333333333333_rkind
  c18 =  h * 0.4390243902439024_rkind
  c19 =  h * 0.8780487804878049_rkind
  c20 =  h * 0.1304878048780488e+01_rkind
  c21 =  h * 0.2097560975609756e+01_rkind
  c22 =  h * 0.2963414634146341e+01_rkind
  c23 =  h * 0.4317073170731707e+01_rkind
  c24 =  h * 0.3214285714285714e-01_rkind
  c25 =  h * 0.4880952380952381e-01_rkind
  c26 =  h * 0.2571428571428571_rkind
  c27 =  h * 0.3238095238095238_rkind

  iout   =  0
  icount =  1
  y      =  0.0_rkind
  y(:,1) =  y0(:)
  f(:,1) =  y0(:)


  do i=2,lnx
!
!... 1-stage
!
     x  =  x0 + real(i-2, rkind) * h
     call fun(x, f(:,1), f(:,3), iout)
     if(iout == 1) return
     f(:,2) =  c01 * f(:,3) + f(:,1)
!
!... 2-stage
!
     x1 =  x + c01
     call fun(x1, f(:,2), f(:,4), iout)
     if(iout == 1) return
     f(:,2) =  c02 * f(:,3) + c03 * f(:,4) + f(:,1)
!
!... 3-stage
!
     x1 =  x + c04
     call fun(x1, f(:,2), f(:,5), iout)
     if(iout == 1) return
     f(:,2) =  c04 * f(:,3) - c05 * f(:,4) + c06 * f(:,5) + f(:,1)
!
!... 4-stage
!
     x1 =  x + c07
     call fun(x1, f(:,2), f(:,6), iout)
     if(iout == 1) return
     f(:,2) =  c03 * f(:,3) + c08 * f(:,6) + f(:,1)
!
!... 5-stage
!
     x1 =  x + c05
     call fun(x1, f(:,2), f(:,7), iout)
     if(iout == 1) return
     f(:,2) = -c09 * f(:,3) + c10 * f(:,7) - c11 * f(:,4) - c12 * f(:,6) &
       &    +  c13 * f(:,5) +       f(:,1)
!
!... 6-stage
!
     x1 =  x + c06
     call fun(x1, f(:,2), f(:,8), iout)
     if(iout == 1) return
     f(:,2) = -c01 * f(:,3) + c03 * f(:,8) + c14 * f(:,7) - c15 * f(:,5) &
       &    +  c16 * f(:,4) +       f(:,1)
!
!... 7-stage
!
     x1 =  x + c17
     call fun(x1, f(:,2), f(:,9), iout)
     if(iout == 1) return
     f(:,2) = -c18 * f(:,8) + c19 * f(:,9) +  c20 * f(:,3) - c21 * f(:,7) &
       &    -  c22 * f(:,4) + c23 * f(:,6) +        f(:,1)
!
!... 8-stage
!
     x1 =  x + h
     call fun(x1, f(:,2), f(:,4), iout)
     if(iout == 1) return
     f(:,1) = (f(:,6) + f(:,8)) * c24 + (f(:,3) + f(:,4)) * c25          &
       &    + (f(:,5) + f(:,9)) * c26 +           f(:,7)  * c27 + f(:,1)
     icount =  i
     y(:,i) =  f(:,1)
  end do

  deallocate(f)


  return
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
end subroutine drkhvg
