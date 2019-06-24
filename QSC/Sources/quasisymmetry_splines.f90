! The subroutines in this module were adapted from subroutine in the stellarator NEO code from Graz.

module quasisymmetry_splines

  use quasisymmetry_variables

  implicit none

  real(dp), parameter :: one = 1.0d+0, zero = 0.0d+0

  type :: periodic_spline
     integer :: n
     real(dp) :: period
     real(dp), dimension (:), pointer :: s, bi, ci, di
  end type periodic_spline

contains

!====================================================

  subroutine new_periodic_spline (n_minus_1, x, y, period, spl)
    ! x is not used for anything!
    implicit none
    integer, intent (in) :: n_minus_1
    real(dp), dimension (n_minus_1), intent (in) :: x, y
    real(dp), intent (in) :: period
    type (periodic_spline), intent (out) :: spl
    integer :: n

    n = n_minus_1 + 1
    spl%n = n
    spl%period = period
    allocate (spl%s(n), spl%bi(n), spl%ci(n), spl%di(n))
    spl%s(1:n_minus_1) = y
    spl%s(n) = y(1)
    call splper(n,period/(spl%n-1),spl%s,spl%bi,spl%ci,spl%di)
  end subroutine new_periodic_spline

!====================================================

  subroutine delete_periodic_spline (spl)
    implicit none
    type (periodic_spline), intent (in out) :: spl
    spl%n = 0
    spl%period = 0.0_dp
    deallocate (spl%s, spl%bi, spl%ci, spl%di)
    nullify (spl%s)
    nullify (spl%bi)
    nullify (spl%ci)
    nullify (spl%di)
  end subroutine delete_periodic_spline

!====================================================

  function periodic_splint (x, spl)
    implicit none
    real(dp), intent (in) :: x
    type (periodic_spline), intent (in) :: spl
    real(dp) :: periodic_splint

    ! Here comes code taken from poi2d.f90
    real(dp) :: xmin, xmax, hx
    real(dp) :: dx, dxx, dxmax, x1
    integer :: ix

    xmin = 0
    xmax = spl%period
    hx = (xmax-xmin)/(spl%n-1)

    dxx = x-xmin
    dxmax = xmax - xmin
    IF(dxx .LT. ZERO) THEN
       dxx = dxx+(1+INT(abs(dxx/dxmax)))*dxmax
    ELSE IF(dxx .GT. dxmax) THEN
       dxx = dxx-(INT(abs(dxx/dxmax)))*dxmax
    END IF
    x1 = dxx/hx
    ix = INT(x1)
    dx = hx*(x1-ix)
    ix = ix+1

    ! Now comes code taken from eva2d.f90
    periodic_splint = spl%s(ix) + dx*( spl%bi(ix) + dx*( spl%ci(ix) + dx * spl%di(ix)))

  end function periodic_splint

!====================================================

SUBROUTINE splper(n,h,y,bi,ci,di)

! Makes a cubic spline of periodic function y(x)
!
! Input:  n                   number of values in y
!         h                   step size in x (aequidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  IMPLICIT NONE

  INTEGER,                     INTENT(in)  :: n
  REAL(kind=dp),               INTENT(in)  :: h
  REAL(kind=dp), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=dp), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=dp)                            :: psi, ss
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: bmx, yl
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: amx1, amx2, amx3
  INTEGER                                  :: nmx, n1, n2, i, i1

  ALLOCATE ( bmx(n), yl(n), amx1(n), amx2(n), amx3(n) )

  bmx(1) = 1.e30_dp

  nmx=n-1
  n1=nmx-1
  n2=nmx-2
  psi=3.0_dp/h/h

  CALL spfper(n,amx1,amx2,amx3)

  bmx(nmx) = (y(nmx+1)-2*y(nmx)+y(nmx-1))*psi
  bmx(1)   =(y(2)-y(1)-y(nmx+1)+y(nmx))*psi
  DO i = 3,nmx
     bmx(i-1) = (y(i)-2*y(i-1)+y(i-2))*psi
  END DO
  yl(1) = bmx(1)/amx1(1)
  DO i = 2,n1
     i1 = i-1
     yl(i) = (bmx(i)-yl(i1)*amx2(i1))/amx1(i)
  END DO
  ss = 0
  DO i = 1,n1
     ss = ss+yl(i)*amx3(i)
  END DO
  yl(nmx) = (bmx(nmx)-ss)/amx1(nmx)
  bmx(nmx) = yl(nmx)/amx1(nmx)
  bmx(n1) = (yl(n1)-amx2(n1)*bmx(nmx))/amx1(n1)
  DO i = n2,1,-1
     bmx(i) = (yl(i)-amx3(i)*bmx(nmx)-amx2(i)*bmx(i+1))/amx1(i)
  END DO
  DO i = 1,nmx
     ci(i) = bmx(i)
  END DO

  DO i = 1,n1
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2*ci(i))/3
     di(i) = (ci(i+1)-ci(i))/h/3
  END DO
  bi(nmx) = (y(n)-y(n-1))/h-h*(ci(1)+2*ci(nmx))/3
  di(nmx) = (ci(1)-ci(nmx))/h/3
!
! Fix of problems at upper periodicity boundary
!
  bi(n) = bi(1)
  ci(n) = ci(1)
  di(n) = di(1)

  DEALLOCATE ( bmx, yl, amx1, amx2, amx3 )

  RETURN
END SUBROUTINE splper


!=====================================================

SUBROUTINE spfper(np1,amx1,amx2,amx3)

! Helper routine for splfi

  IMPLICIT NONE

  INTEGER,                       INTENT(in)  :: np1
  REAL(kind=dp), DIMENSION(np1), INTENT(out) :: amx1, amx2, amx3
  REAL(kind=dp)                              :: beta, ss
  INTEGER                                    :: n, n1, i, i1

  n = np1-1

  n1 = n-1
  amx1(1) = 2
  amx2(1) = 0.5_dp
  amx3(1) = 0.5_dp
  amx1(2) = sqrt(15._dp)/2
  amx2(2) = ONE/amx1(2)
  amx3(2) = -.25_dp/amx1(2)
  beta = 3.75_dp
  DO i = 3,n1
     i1 = i-1
     beta = 4-ONE/beta
     amx1(i) = sqrt(beta)
     amx2(i) = ONE/amx1(i)
     amx3(i) = -amx3(i1)/amx1(i)/amx1(i1)
  END DO
  amx3(n1) = amx3(n1)+ONE/amx1(n1)
  amx2(n1) = amx3(n1)
  ss = 0
  DO i = 1,n1
     ss = ss+amx3(i)*amx3(i)
  END DO
  amx1(n) = sqrt(4-ss)

  RETURN
END SUBROUTINE spfper

end module quasisymmetry_splines
