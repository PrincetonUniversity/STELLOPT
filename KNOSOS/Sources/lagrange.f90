!Performs interpolations and extrapolations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE LAGRANGE_C(xa,ya,n,x,y,order) 

!-----------------------------------------------------------------------------------------------
!Given xa(1:n) and the complex array ya(1:n) of length n, which tabulate a function 
!(with the xai’s in order),and given a value of x, returns a lagrange-interpolated y
!-----------------------------------------------------------------------------------------------  

  IMPLICIT NONE
  !Input
  INTEGER n,order
  REAL*8 x,xa(n)
  COMPLEX*16 ya(n)
  !Output
  COMPLEX*16 y
  !Others
  LOGICAL REVERSE
  INTEGER k,khi,klo,j 
  REAL*8 xl,tempxa(n)
  COMPLEX*16 tempya(n)

  !If xa is not in growing order, reverse
  REVERSE=.FALSE.
  IF(xa(2).LT.xa(1)) THEN
     REVERSE=.TRUE.
     DO k=1,n
        tempxa(k)=xa(n-k+1)
        tempya(k)=ya(n-k+1)
     END DO
     xa=tempxa
     ya=tempya
  END IF

  klo=1 
  khi=n
1 IF (khi-klo.GT.1) then  
     k=(khi+klo)/2
     IF(xa(k).GT.x) then 
        khi=k
     ELSE 
        klo=k
     END IF
     GOTO 1
  END IF

  IF (klo.GT.(n-order)) klo=n-order
  
  y=0.d0
  DO k = 0, order
     xl=1.d0
     DO j = 0, order
        IF (j.NE.k) xl=xl*(x-xa(j+klo))/(xa(klo+k)-xa(j+klo))
     END DO
     y=y+xl*ya(klo+k)
  END DO
  
  IF(REVERSE) THEN
     DO k=1,n
        tempxa(k)=xa(n-k+1)
        tempya(k)=ya(n-k+1)
     END DO
     xa=tempxa
     ya=tempya
  END IF


END SUBROUTINE LAGRANGE_C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE LAGRANGE(xa,ya,n,x,y,order) 

!-----------------------------------------------------------------------------------------------
!Given xa(1:n) and the real array ya(1:n) of length n, which tabulate a function 
!(with the xai’s in order),and given a value of x, returns a lagrange-interpolated y
!-----------------------------------------------------------------------------------------------
  
  IMPLICIT NONE
  !Input
  INTEGER n,order
  REAL*8 x,xa(n),ya(n)
  !Output
  REAL*8 y
  !Others
  LOGICAL REVERSE
  INTEGER k,khi,klo,j 
  REAL*8 xl,tempxa(n),tempya(n)

  !If xa is not in growing order, reverse
  REVERSE=.FALSE.
  IF(xa(2).LT.xa(1)) THEN
     REVERSE=.TRUE.
     DO k=1,n
        tempxa(k)=xa(n-k+1)
        tempya(k)=ya(n-k+1)
     END DO
     xa=tempxa
     ya=tempya
  END IF

  klo=1 
  khi=n
1 IF (khi-klo.GT.1) then  
     k=(khi+klo)/2
     IF(xa(k).GT.x) then 
        khi=k
     ELSE 
        klo=k
     END IF
     GOTO 1
  END IF

  IF (klo.GT.(n-order)) klo=n-order
  
  y=0.d0
  DO k = 0, order
     xl=1.d0
     DO j = 0, order
        IF (j.NE.k) xl=xl*(x-xa(j+klo))/(xa(klo+k)-xa(j+klo))
     END DO
     y=y+xl*ya(klo+k)
  END DO
  
  IF(REVERSE) THEN
     DO k=1,n
        tempxa(k)=xa(n-k+1)
        tempya(k)=ya(n-k+1)
     END DO
     xa=tempxa
     ya=tempya
  END IF


END SUBROUTINE LAGRANGE
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE BILAGRANGE(x1a,x2a,ya,m,n,x1,x2,y,order) 

!-----------------------------------------------------------------------------------------------
!Given x1a(1:n) and x2a(1:m) and the real array ya(1:n,1:m), which tabulate a function 
!(with the xai’s in order),and given a value of x1, x2, returns a lagrange-interpolated y
!-----------------------------------------------------------------------------------------------

  IMPLICIT NONE
  !Input
  INTEGER m,n,order
  REAL*8 x1,x2,x1a(m),x2a(n),ya(m,n) 
  !Output
  REAL*8 y
  !Others
  INTEGER j,k
  INTEGER, PARAMETER :: NN=100 ! Maximum expected value of n and m.
  REAL*8 ytmp(NN),yytmp(NN)

  !Interpolates in one direction and creates a 1D grid
  DO j=1,m
     DO k=1,n 
        ytmp(k)=ya(j,k)
     END DO
     CALL LAGRANGE(x2a,ytmp,n,x2,yytmp(j),order) 
  END DO
  !Interpolates in the other direction
  CALL LAGRANGE(x1a,yytmp,m,x1,y,order)

  
END SUBROUTINE BILAGRANGE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE TRILAGRANGE(x1a,x2a,x3a,ya,m,n,o,x1,x2,x3,y,order) 

!-----------------------------------------------------------------------------------------------
!Given x1a(1:n), x2a(1:m) and x3a(1:o) and the real array ya(1:n,1:m,1:o), which tabulate a 
!function  (with the xai’s in order),and given a value of x1, x2, x3, returns a lagrange-
!interpolated y
!-----------------------------------------------------------------------------------------------

  IMPLICIT NONE
  !Input
  INTEGER m,n,o,order
  REAL*8 x1,x2,x3,x1a(m),x2a(n),x3a(o),ya(m,n,o) 
  !Output
  REAL*8 y
  !Others
  INTEGER j,k,l
  REAL*8 ytmp(o),yytmp(m,n)

  !Interpolates in one direction and creates a 2D grid
  DO j=1,m
     DO k=1,n 
        DO l=1,o
           ytmp(l)=ya(j,k,l)
         END DO
         CALL LAGRANGE(x3a,ytmp(1:o),o,x3,yytmp(j,k),order)
     END DO
  END DO
  !Interpolates in the other two directions
  CALL BILAGRANGE(x1a,x2a,yytmp(1:m,1:n),m,n,x1,x2,y,order)

END SUBROUTINE TRILAGRANGE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
