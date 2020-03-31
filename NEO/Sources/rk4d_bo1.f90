
SUBROUTINE rk4d_bo1(x,y,h)
! Differential Equation Solver
  USE neo_precision
  USE sizey_bo
  USE neo_rhsbo
!
  IMPLICIT NONE
!
  REAL(kind=dp),               INTENT(inout) ::  x
  REAL(kind=dp),               INTENT(in)    ::  h
  REAL(kind=dp)                              ::  hh,h6,xh
  REAL(kind=dp), DIMENSION(:), INTENT(inout) ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE   ::  dydx,yt,dyt,dym
  INTEGER                                    ::  i
!
  ALLOCATE(dydx(ndim))
  ALLOCATE(yt(ndim))
  ALLOCATE(dyt(ndim))
  ALLOCATE(dym(ndim))
!
  hh = h / 2
  h6 = h/6
  xh = x+hh
  CALL rhs_bo(x,y,dydx)
  DO i=1,ndim
     yt(i)=y(i)+hh*dydx(i)
  END DO
  CALL rhs_bo(xh,yt,dyt)
  DO i=1,ndim
     yt(i)=y(i)+hh*dyt(i)
  END DO
  CALL rhs_bo(xh,yt,dym)
  DO i=1,ndim
     yt(i)=y(i)+h*dym(i)
     dym(i)=dyt(i)+dym(i)
  END DO
  x=x+h
  CALL rhs_bo(x,yt,dyt)
  DO i=1,ndim
     y(i)=y(i)+h6*(dydx(i)+dyt(i)+2*dym(i))
  END DO

  DEALLOCATE (dydx, yt, dyt, dym)

  RETURN
END SUBROUTINE rk4d_bo1
