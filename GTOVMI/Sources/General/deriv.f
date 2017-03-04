      FUNCTION deriv(xx,yy,n)

!df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
! Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.

       USE precision
       IMPLICIT NONE
       INTEGER :: n,n2
       REAL(rprec) :: xx(n),yy(n)
       REAL(rprec),DIMENSION(n) :: x, y, dydx, x01, x02, x12
       REAL(rprec) :: deriv(n)
       INTERFACE shiftx
        FUNCTION shiftx(x,n,m)
         IMPLICIT NONE 
         INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
         INTEGER :: m, n
         REAL(rprec) :: x(n)
         REAL(rprec),DIMENSION(n) :: x1
         REAL(rprec) :: shiftx(n)
        END FUNCTION shiftx
       END INTERFACE shiftx
       x=xx; y=yy; dydx=0
       x12 = x - shiftx(x,n,-1) 
       x01 = shiftx(x,n,1) -x
       x02 = shiftx(x,n,1) - shiftx(x,n,-1)
       dydx = shiftx(y,n,1) * (x12 / (x01*x02)) + !	Middle points
     .	y * (1./x12 - 1./x01) - 
     .	shiftx(y,n,-1) * (x01 / (x02 * x12))
       dydx(1) = y(1) * (x01(2)+x02(2))/(x01(2)*x02(2)) - !	First point
     .  y(2) * x02(2)/(x01(2)*x12(2)) + 
     .  y((3)) * x01(2)/(x02(2)*x12(2)) 
       n2 = n- 1
       dydx(n) = -y(n-2) * x12(n2)/(x01(n2)*x02(n2)) + !	Last point
     .  y(n-1) * x02(n2)/(x01(n2)*x12(n2)) - 
     .  y(n) * (x02(n2)+x12(n2)) / (x02(n2)*x12(n2))
       deriv=dydx
       RETURN
      END FUNCTION deriv
