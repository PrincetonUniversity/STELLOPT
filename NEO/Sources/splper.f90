!====================================================
SUBROUTINE splper(n,h,y,bi,ci,di)

! Makes a cubic spline of periodic function y(x)
!
! Input:  n                   number of values in y
!         h                   step size in x (aequidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  USE neo_precision

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
