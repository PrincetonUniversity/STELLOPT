!=====================================================
SUBROUTINE splreg(n,h,y,bi,ci,di)

! Makes a cubic spline of function y(x)
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

  REAL(kind=dp)                            :: ak1, ak2, am1, am2, c, e, c1
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: al, bt
  INTEGER                                  :: k, n2, i, i5

  ALLOCATE ( al(n), bt(n) )

  ak1 = 0
  ak2 = 0
  am1 = 0
  am2 = 0
  k = n-1
  al(1) = ak1
  bt(1) = am1
  n2 = n-2
  c = -4*h
  DO i = 1,n2
     e = -3*((y(i+2)-y(i+1))-(y(i+1)-y(i)))/h
     c1 = c-al(i)*h
     al(i+1) = h/c1
     bt(i+1) = (h*bt(i)+e)/c1
  END DO
  ci(n) = (am2+ak2*bt(k))/(1-al(k)*ak2)
  DO i = 1,k
     i5 = n-i
     ci(i5) = al(i5)*ci(i5+1)+bt(i5)
  END DO
  n2 = n-1
  DO i = 1,n2
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2*ci(i))/3
     di(i) = (ci(i+1)-ci(i))/h/3
  END DO
  DEALLOCATE ( al, bt )

  RETURN
END SUBROUTINE splreg
