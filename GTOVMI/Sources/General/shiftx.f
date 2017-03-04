      FUNCTION shiftx(x,n,m)
! shift elememts of x(n) bu m ; USEful for calculating dx
       USE precision
       IMPLICIT NONE 
       INTEGER :: m, n
       REAL(rprec) :: x(n)
       REAL(rprec),DIMENSION(n) :: x1
       REAL(rprec) :: shiftx(n)
       x1=0
       IF (m .eq. 0) THEN
        x1=x
       ELSEIF (m .gt. 0) THEN
        x1(m+1:n)=x(1:n-m) 
        x1(1:m)=x(n-m+1:n)
       ELSE
        x1(1:n+m)=x(1-m:n) 
        x1(n+m+1:n)=x(1:-m)
       ENDIF
       shiftx=x1
       RETURN
      END FUNCTION shiftx


