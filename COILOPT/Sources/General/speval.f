      SUBROUTINE speval(n,xk,c,t,f,ifail)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: n, ifail
      REAL(rprec), DIMENSION(0:*) :: xk, c
      REAL(rprec), DIMENSION(*) :: f
      REAL(rprec) :: t
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, l, ip1
      REAL(rprec) :: a10, a11, a12, a20, a21
      REAL(rprec) :: b02, b12, b03, b13, b23, b04, b14, b24, b34
      REAL(rprec) :: d10, d11, d12, d20, d21, d30
      REAL(rprec) :: tm0, tm1, tm2, tp1, tp2, tp3

!     evaluate the cubic spline with knots xk and b-spline
!     coefficients c at the point t. return function value
!     in f(1) and derivatives in f(2)-f(4).

      ifail=1
      DO j=1,4
         f(j)=0
      END DO
      IF((t.lt.xk(4)) .or. (t.gt.xk(n+1))) RETURN
      i=4
      ip1=n+1
   20 l=(i+ip1)/2
      IF(ip1-i.le.1) GO TO 40
      IF(t.lt.xk(l)) GO TO 30
      i=l
      GO TO 20
   30 ip1=l
      GO TO 20
   40 tm2=t-xk(i-2)
      tm1=t-xk(i-1)
      tm0=t-xk(i)
      tp1=xk(i+1)-t
      tp2=xk(i+2)-t
      tp3=xk(i+3)-t
      d10=tp1+tm0
      d11=tp1+tm1
      d20=tp2+tm0
      d12=tp1+tm2
      d21=tp2+tm1
      d30=tp3+tm0
      b12=tp1/d10
      b02=tm0/d10
      b23=tp1*b12/d11
      b13=tm1*b12/d11+tp2*b02/d20
      b03=tm0*b02/d20
      b34=tp1*b23/d12
      b24=tm2*b23/d12+tp2*b13/d21
      b14=tm1*b13/d21+tp3*b03/d30
      b04=tm0*b03/d30
      f(1)=c(i-3)*b34+c(i-2)*b24+c(i-1)*b14+c(i)*b04
      a12=(c(i-2)-c(i-3))/d12
      a11=(c(i-1)-c(i-2))/d21
      a10=(c(i)-c(i-1))/d30
      f(2)=3*(a12*b23+a11*b13+a10*b03)
      a21=(a11-a12)/d11
      a20=(a10-a11)/d20
      f(3)=6*(a21*b12+a20*b02)
      f(4)=6*(a20-a21)/d10
      ifail=0

      END SUBROUTINE speval
