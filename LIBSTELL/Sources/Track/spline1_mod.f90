MODULE SPLINE1_MOD
!-------------------------------------------------------------------------------
!SPLINE1-SPLINE interpolation in 1d
!
!SPLINE1_MOD is an F90 module of cubic interpolating spline routines in 1d
!
!References:
!
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!  W.A.Houlberg, F90 free form 8/2004
!
!Contains PUBLIC routines:
!
!  SPLINE1_FIT         -get the coefficients
!  SPLINE1_EVAL        -evaluate the spline
!  SPLINE1_INTERP      -interpolate from one grid to another
!  SPLINE1_INTEG       -integrate the spline fit
!
!Comments:
!
!  Spline interpolation routines are C2 (f, f', and f'' are continuous)
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  SPLINE1_SEARCH         !find the indices that bracket an abscissa value
                         !  called from SPLINE1_EVAL

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE SPLINE1_FIT(n,x,f, &
                       K_BC1,K_BCN)
!-------------------------------------------------------------------------------
!SPLINE1_FIT gets the coefficients for a 1d cubic interpolating spline
!
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!  D.McCune, W.A.Houlberg 3/2004
!  W.A.Houlberg, F90 free form 8/2004
!
!Comments:
!  For x(i).le.x.le.x(i+1)
!    s(x) = f(1,i) + f(2,i)*(x-x(i)) + f(3,i)*(x-x(i))**2/2!
!          +f(4,i)*(x-x(i))**3/3!
!  The cubic spline is twice differentiable (C2)
!  The BCs default to not-a-knot conditions, K_BC1=0 and K_BCN=0
!  Two errors fixed in 3/2004 version:
!    Logic for not-a-knot on f(3,1) and f(3,n) after back substitution
!    Treatment of asymmetrix coefficient matrix for s''=0 BC
!  For s''=0 the off-diagonal elements nearest the corners are not symmetric,
!    i.e., el(1,2) /= el(2,1) and el(n-1,n) /= el(n,n-1).  The correct values:
!    el(1,2)=0, el(2,1)=x(2)-x(1), el(n,n-1)=0, el(n-1,n)=x(n)-x(n-1) are fixed
!    by adding variables to hold the el21 and elnn1 values.
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n                      !number of data points or knots [>=2]

REAL(KIND=rspec), INTENT(IN) :: &
  x(:)                   !abscissas of the knots in increasing order [arb]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  f(:,:)                 !ordinates of the knots [arb]
                         !f(2,1)=input value of s'(x1) for K_BC1=1 [arb]
                         !f(2,n)=input value of s'(xn) for K_BCN=1 [arb]
                         !f(3,1)=input value of s''(x1) for K_BC1=2 [arb]
                         !f(3,n)=input value of s''(xn) for K_BCN=2 [arb]
                         !f(2,i)=output s'(x(i)) [arb]
                         !f(3,i)=output s''(x(i)) [arb]
                         !f(4,i)=output s'''(x(i)) [arb]

!Declaration of optional input variables
INTEGER, INTENT(IN), OPTIONAL :: &
  K_BC1,               & !option for BC at x(1) [-]
                         !=-1 periodic, ignore K_BCN
                         !=0 not-a-knot (default)
                         !=1 s'(x1) = input value of f(2,1)
                         !=2 s''(x1) = input value of f(3,1)
                         !=3 s'(x1) = 0.0
                         !=4 s''(x1) = 0.0
                         !=5 match first derivative to first 2 points
                         !=6 match second derivative to first 3 points
                         !=7 match third derivative to first 4 points
                         !=else use not-a-knot
  K_BCN                  !option for boundary condition at x(n) [-]
                         !=0 not-a-knot (default)
                         !=1 s'(x1) = input value of f(2,1)
                         !=2 s''(x1) = input value of f(3,1)
                         !=3 s'(x1) = 0.0
                         !=4 s''(x1) = 0.0
                         !=5 match first derivative to first 2 points
                         !=6 match second derivative to first 3 points
                         !=7 match third derivative to first 4 points
                         !=else use knot-a-knot

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,ib,imax,imin,kbc1,kbcn

REAL(KIND=rspec) :: &
  a1,an,b1,bn,el21,elnn1,hn,q,t,wk(1:n)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Set left BC option, and leftmost node
!Default to not-a-knot
kbc1=0

!Otherwise use requested value
IF(PRESENT(K_BC1)) THEN

  IF(K_BC1 >= -1 .AND. &    
     K_BC1 <= 7) kbc1=K_BC1

ENDIF

!Include first node for all but not-a-knot
imin=1

!Not-a-knot condition removes first node
IF(kbc1 == 0) imin=2

!Set left BC values
!Default for not-a-knot
a1=0
b1=0

IF(kbc1 == 1) THEN

  !First derivative specified
  a1=f(2,1)

ELSEIF(kbc1 == 2) THEN

  !Second derivative specified
  b1=f(3,1)

ELSEIF(kbc1 == 5) THEN

  !Match first derivative to first two points
  a1=(f(1,2)-f(1,1))/(x(2)-x(1))

ELSEIF(kbc1 == 6) THEN

  !Match second derivative to first three points
  b1=2*((f(1,3)-f(1,2))/(x(3)-x(2))-(f(1,2)-f(1,1))/(x(2)-x(1)))/(x(3)-x(1))

ENDIF

!Set right BC option, and rightmost node
!Default to not-a-knot
kbcn=0

!Otherwise use requested value
IF(PRESENT(K_BCN)) THEN

  IF(K_BCN >= -1 .AND. &
     K_BCN <= 7) kbcn=K_BCN

ENDIF

!Include last node for all but not-a-knot
imax=n

!Not-a-knot condition removes last node
IF(kbcn == 0) imax=n-1

!Set right BC values
!Default for not-a-knot
an=0
bn=0

IF(kbcn == 1) THEN

  !First derivative specified
  an=f(2,n)

ELSEIF(kbcn == 2) THEN

  !Second derivative specified
  bn=f(3,n)

ELSEIF(kbcn == 5) THEN

  !Match first derivative to last two points
  an=(f(1,n)-f(1,n-1))/(x(n)-x(n-1))

ELSEIF(kbcn == 6) THEN

  !Match second derivative to first three points
  bn=2*((f(1,n)-f(1,n-1))/(x(n)-x(n-1)) &   
     -(f(1,n-1)-f(1,n-2))/(x(n-1)-x(n-2)))/(x(n)-x(n-2))

ENDIF

!Clear derivatives, f(2:4,1:n), and work array
f(2:4,1:n)=0
wk(1:n)=0

!-------------------------------------------------------------------------------
!Evaluate coefficients
!-------------------------------------------------------------------------------
!Two and three nodes are special cases
IF(n == 2) THEN

  !Coefficients for n=2
  f(2,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
  f(3,1)=0
  f(4,1)=0
  f(2,2)=f(2,1)
  f(3,2)=0
  f(4,2)=0

!Set up tridiagonal system for A*y=B where y(i) are the second
!  derivatives at the knots
!  f(2,i) are the diagonal elements of A
!  f(4,i) are the off-diagonal elements of A
!  f(3,i) are the B elements/3, and will become c/3 upon solution

ELSEIF(n > 2) THEN

  f(4,1)=x(2)-x(1)
  f(3,2)=(f(1,2)-f(1,1))/f(4,1)

  DO i=2,n-1 !Over nodes

    f(4,i)=x(i+1)-x(i)
    f(2,i)=2*(f(4,i-1)+f(4,i))
    f(3,i+1)=(f(1,i+1)-f(1,i))/f(4,i)
    f(3,i)=f(3,i+1)-f(3,i)

  ENDDO !Over nodes

  !Save elements for non-symmetric matrix cases
  el21=f(4,1)
  elnn1=f(4,n-1)

!Apply left BC
  IF(kbc1 == -1) THEN

    !Periodic
    f(2,1)=2*(f(4,1)+f(4,n-1))
    f(3,1)=(f(1,2)-f(1,1))/f(4,1)-(f(1,n)-f(1,n-1))/f(4,n-1)
    wk(1)=f(4,n-1)
    wk(2:n-3)=0
    wk(n-2)=f(4,n-2)
    wk(n-1)=f(4,n-1)

  ELSEIF((kbc1 == 1) .OR. &
         (kbc1 == 3) .OR. & 
         (kbc1 == 5)) THEN

    !First derivative conditions
    f(2,1)=2*f(4,1)
    f(3,1)=(f(1,2)-f(1,1))/f(4,1)-a1

  ELSEIF((kbc1 == 2) .OR. (kbc1 == 4) .OR. (kbc1 == 6)) THEN

    !Second derivative conditions
    f(2,1)=2*f(4,1)
    f(3,1)=f(4,1)*b1/3
    f(4,1)=0

  ELSEIF(kbc1 == 7) THEN

    !Third derivative condition
    f(2,1)=-f(4,1)
    f(3,1)=f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
    f(3,1)=f(3,1)*f(4,1)**2/(x(4)-x(1))

  ELSE

    !Not-a-knot condition
    f(2,2)=f(4,1)+2*f(4,2)
    f(3,2)=f(3,2)*f(4,2)/(f(4,1)+f(4,2))

  ENDIF

  !Apply right BC
  IF((kbcn == 1) .OR. &
     (kbcn == 3) .OR. & 
     (kbcn == 5)) THEN

    !First derivative conditions
    f(2,n)=2*f(4,n-1)
    f(3,n)=-(f(1,n)-f(1,n-1))/f(4,n-1)+an

  ELSEIF((kbcn == 2) .OR. &
         (kbcn == 4) .OR. &
         (kbcn == 6)) THEN

    !Second derivative conditions
    f(2,n)=2*f(4,n-1)
    f(3,n)=f(4,n-1)*bn/3
    elnn1=0

  ELSEIF(kbcn == 7) THEN

    !Third derivative condition
    f(2,n)=-f(4,n-1)
    f(3,n)=f(3,n-1)/(x(n)-x(n-2))-f(3,n-2)/(x(n-1)-x(n-3))
    f(3,n)=-f(3,n)*f(4,n-1)**2/(x(n)-x(n-3))

  ELSEIF(kbc1 /= -1) THEN

    !Not-a-knot condition
    f(2,n-1)=2*f(4,n-2)+f(4,n-1)
    f(3,n-1)=f(3,n-1)*f(4,n-2)/(f(4,n-1)+f(4,n-2))

  ENDIF

  !Limit solution for only three points in domain
  IF(n == 3) THEN

    f(3,1)=0
    f(3,n)=0

  ENDIF

  !Solve system of equations for second derivatives at the knots
  IF(kbc1 == -1) THEN

    !Periodic BC - requires special treatment at ends
    !Forward elimination
    DO i=2,n-2 !Over nodes in forward elimination

      t=f(4,i-1)/f(2,i-1)
      f(2,i)=f(2,i)-t*f(4,i-1)
      f(3,i)=f(3,i)-t*f(3,i-1)
      wk(i)=wk(i)-t*wk(i-1)
      q=wk(n-1)/f(2,i-1)
      wk(n-1)=-q*f(4,i-1)
      f(2,n-1)=f(2,n-1)-q*wk(i-1)
      f(3,n-1)=f(3,n-1)-q*f(3,i-1)

    ENDDO !Over nodes in forward elimination

    !Correct the n-1 element
    wk(n-1)=wk(n-1)+f(4,n-2)

    !Complete the forward elimination
    !wk(n-1) and wk(n-2) are the off-diag elements of the lower corner
    t=wk(n-1)/f(2,n-2)
    f(2,n-1)=f(2,n-1)-t*wk(n-2)
    f(3,n-1)=f(3,n-1)-t*f(3,n-2)

    !Back substitution
    f(3,n-1)=f(3,n-1)/f(2,n-1)
    f(3,n-2)=(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)

    DO ib=3,n-1 !Over nodes in back substitution

      i=n-ib
      f(3,i)=(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)

    ENDDO !Over nodes in back substitution

    f(3,n)=f(3,1)

  ELSE

    !Non-periodic BC
    !Forward elimination
    !For not-a-knot and s''=0 BCs the off-diagonal end elements are not equal
    DO i=imin+1,imax !Over nodes in forward elimination

      IF((i == n-1) .AND. &
         (imax == n-1)) THEN

        t=(f(4,i-1)-f(4,i))/f(2,i-1)

      ELSE

        IF(i == 2) THEN

          t=el21/f(2,i-1)

        ELSEIF(i == n) THEN

          t=elnn1/f(2,i-1)

        ELSE

          t=f(4,i-1)/f(2,i-1)

        ENDIF

      ENDIF

      IF((i == imin+1) .AND. &
         (imin == 2)) THEN

        f(2,i)=f(2,i)-t*(f(4,i-1)-f(4,i-2))

      ELSE

        f(2,i)=f(2,i)-t*f(4,i-1)

      ENDIF

      f(3,i)=f(3,i)-t*f(3,i-1)

    ENDDO !Over nodes in forward elimination

    !Back substitution
    f(3,imax)=f(3,imax)/f(2,imax)

    DO ib=1,imax-imin !Over nodes in back substitution

      i=imax-ib

      IF((i == 2) .AND. &
         (imin == 2)) THEN

        f(3,i)=(f(3,i)-(f(4,i)-f(4,i-1))*f(3,i+1))/f(2,i)

      ELSE

        f(3,i)=(f(3,i)-f(4,i)*f(3,i+1))/f(2,i)

      ENDIF

    ENDDO !Over nodes in back substitution

    !Reset d array to step size
    f(4,1)=x(2)-x(1)
    f(4,n-1)=x(n)-x(n-1)

    !Set f(3,1) for not-a-knot; 5 changed to 7 in 3/2004
    IF((kbc1 <= 0) .OR. &
       (kbc1 > 7)) THEN

      f(3,1)=(f(3,2)*(f(4,1)+f(4,2))-f(3,3)*f(4,1))/f(4,2)

    ENDIF

    !Set f(3,n) for not-a-knot; 5 changed to 7 in 3/2004
    IF((kbcn <= 0) .OR. &
       (kbcn > 7)) THEN

      f(3,n)=f(3,n-1)+(f(3,n-1)-f(3,n-2))*f(4,n-1)/f(4,n-2)

    ENDIF

  ENDIF

  !f(3,i) is now the sigma(i) of the text and f(4,i) is the step size
  !Compute polynomial coefficients
  DO i=1,n-1 !Over nodes

    f(2,i)=(f(1,i+1)-f(1,i))/f(4,i)-f(4,i)*(f(3,i+1)+2*f(3,i))
    f(4,i)=(f(3,i+1)-f(3,i))/f(4,i)
    f(3,i)=6*f(3,i)
    f(4,i)=6*f(4,i)

  ENDDO !Over nodes

  IF(kbc1 == -1) THEN

    !Periodic BC
    f(2,n)=f(2,1)
    f(3,n)=f(3,1)
    f(4,n)=f(4,1)

  ELSE

    !All other BCs
    hn=x(n)-x(n-1)
    f(2,n)=f(2,n-1)+hn*(f(3,n-1)+hn*f(4,n-1)/2)
    f(3,n)=f(3,n-1)+hn*f(4,n-1)
    f(4,n)=f(4,n-1)

    IF((kbcn == 1) .OR. &
       (kbcn == 3) .OR. &
       (kbcn == 5)) THEN

      !First derivative BC
      f(2,n)=an

    ELSEIF((kbcn == 2) .OR. &
           (kbcn == 4) .OR. &
           (kbcn == 6)) THEN

      !Second derivative BC
      f(3,n)=bn

    ENDIF

  ENDIF

ENDIF

END SUBROUTINE SPLINE1_FIT

SUBROUTINE SPLINE1_EVAL(k_vopt,n,u,x,f, &
                        i, &
                        value)
!-------------------------------------------------------------------------------
!SPLINE1_EVAL evaluates the cubic spline function and its derivatives
!
!References:
!  W.A.Houlberg, P.I. Strand, D.McCune 8/2001
!  W.A.Houlberg, F90 free form 8/2004
!
!Comments:
!  s=f(1,i)+f(2,i)*dx+f(3,i)*dx**2/2!+f(4,i)*dx**3/3!
!  s'=f(2,i)+f(3,i)*dx+f(4,i)*dx**2/2!
!  s''=f(3,i)+f(4,i)*dx
!    where dx=u-x(i) and x(i).lt.u.lt.x(i+1)
!  If u <= x(1) then i=1 is used
!  If u >= x(n) then i=n is used
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_vopt(:),           & !k_vopt(1)=calculate the function [0=off]
                         !k_vopt(2)=calculate the first derivative [0=off]
                         !k_vopt(3)=calculate the second derivative [0=off]
  n                      !number of data points [-]

REAL(KIND=rspec), INTENT(IN) :: &
  u,                   & !abscissa at which the spline is to be evaluated [arb]
  f(:,:),              & !array containing the data ordinates [arb]
  x(:)                   !array containing the data abcissas [arb]

!Declaration of input/output variables
INTEGER, INTENT(INOUT) :: &
  i                      !in=guess for target lower bound if 1<i<n [-]
                         !out=target lower bound [-]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  value(:)               !ordering is a subset of the sequence under k_vopt [arb]
                         !value(1)=function or lowest order derivative requested [arb]
                         !value(2)=next order derivative requested [arb]
                         !value(3)=second derivative if all k_vopt are non-zero [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  j

REAL(KIND=rspec) :: &
  dx

!-------------------------------------------------------------------------------
!Find target cell
!-------------------------------------------------------------------------------
CALL SPLINE1_SEARCH(x,n,u, &
                    i)

!-------------------------------------------------------------------------------
!Set values outside of range to endpoints
!-------------------------------------------------------------------------------
IF(i <= 0) THEN

  i=1
  dx=0

ELSEIF(i >= n) THEN

  i=n
  dx=0

ELSE

 dx=u-x(i)

ENDIF

!-------------------------------------------------------------------------------
!Evaluate spline
!-------------------------------------------------------------------------------
j=0

!Value
IF(k_vopt(1) /= 0) THEN

  j=j+1
  value(j)=f(1,i)+dx*(f(2,i)+dx/2*(f(3,i)+f(4,i)*dx/3))

ENDIF

!First derivative
IF(k_vopt(2) /= 0) THEN

  j=j+1
  value(j)=f(2,i)+dx*(f(3,i)+f(4,i)*dx/2)

ENDIF

!Second derivative
IF(k_vopt(3) /= 0) THEN

  j=j+1
  value(j)=f(3,i)+dx*f(4,i)

ENDIF

END SUBROUTINE SPLINE1_EVAL

SUBROUTINE SPLINE1_INTERP(k_vopt,n0,x0,f,n1,x1,value, &   
                          iflag,message, &
                          K_BC1,K_BCN)
!-------------------------------------------------------------------------------
!SPLINE1_INTERP performs a spline interpolation from the x0 mesh to the x1 mesh
!
!References:
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!  W.A.Houlberg, F90 free form 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_vopt(:),           & !option for output values
                         !k_vopt(1)=calculate the function [0=off]
                         !k_vopt(2)=calculate the first derivative [0=off]
                         !k_vopt(3)=calculate the second derivative [0=off]
  n0,                  & !number of source abscissas [-]
  n1                     !number of target abscissas [-]

REAL(KIND=rspec), INTENT(IN) :: &
  x0(:),               & !source abscissas (in increasing order) [arb]
  x1(:)                  !target abscissas (in increasing order) [arb]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &
  f(:,:)                 !f(1,n0)=input source values [arb]
                         !f(2,1)=input value of s'(x1) for K_BC1=1 [arb]
                         !f(2,n0)=input value of s'(xn) for K_BCN=1 [arb]
                         !f(3,1)=input value of s''(x1) for K_BC1=2 [arb]
                         !f(3,n0)=input value of s''(xn) for K_BCN=2 [arb]
                         !f(4,n)=output arrays of n0 spline coefficients
                         !f(2,i)=output s'(x0(i))/1!
                         !f(3,i)=output s''(x0(i))/2!
                         !f(4,i)=output s'''(x0(i))/3!

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]
     
INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  value(:,:)             !ordering is a subset of the sequence under k_vopt [arb]
                         !value(1,j)=function or lowest order derivative requested [arb]
                         !value(2,j)=next order derivative requested [arb]
                         !value(3,j)=second derivative if all k_vopt are non-zero [arb]


!Declaration of optional input variables
INTEGER, INTENT(IN), OPTIONAL :: &
  K_BC1,               & !option for BC at x0(1) [-]
                         !=-1 periodic, ignore K_BCN
                         !=0 not-a-knot (default)
                         !=1 s'(x1) = input value of f(2,1)
                         !=2 s''(x1) = input value of f(3,1)
                         !=3 s'(x1) = 0.0
                         !=4 s''(x1) = 0.0
                         !=5 match first derivative to first 2 points
                         !=6 match second derivative to first 3 points
                         !=7 match third derivative to first 4 points
                         !=else use not-a-knot
  K_BCN                  !option for boundary condition at x0(n0) [-]
                         !=0 not-a-knot (default)
                         !=1 s'(x1) = input value of f(2,1)
                         !=2 s''(x1) = input value of f(3,1)
                         !=3 s'(x1) = 0.0
                         !=4 s''(x1) = 0.0
                         !=5 match first derivative to first 2 points
                         !=6 match second derivative to first 3 points
                         !=7 match third derivative to first 4 points
                         !=else use knot-a-knot

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
IF(n1 <= 0) THEN

  iflag=-1
  message='SPLINE1_INTERP/WARNING(1):no points in output array'
  GOTO 9999

ENDIF

!Make sure there are at least two nodes in the input grid
IF(n0 < 2) THEN

  iflag=1
  message='SPLINE1_INTERP/ERROR(2):<2 points in source array'
  GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Get spline coefficients
!-------------------------------------------------------------------------------
IF(PRESENT(K_BC1) .AND. &
   PRESENT(K_BCN)) THEN

  CALL SPLINE1_FIT(n0,x0,f, &
                   K_BC1=K_BC1, &
                   K_BCN=K_BCN)

  ELSEIF(PRESENT(K_BC1)) THEN

    CALL SPLINE1_FIT(n0,x0,f, &
                     K_BC1=K_BC1)

  ELSEIF(PRESENT(K_BCN)) THEN

    CALL SPLINE1_FIT(n0,x0,f, &
                     K_BCN=K_BCN)

  ELSE

    CALL SPLINE1_FIT(n0,x0,f)

ENDIF

!-------------------------------------------------------------------------------
!Interpolate onto x1 mesh
!-------------------------------------------------------------------------------
i=1

DO j=1,n1 !Over nodes

  CALL SPLINE1_EVAL(k_vopt,n0,x1(j),x0,f,i,value(1:3,j))

ENDDO !Over nodes

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
 9999 CONTINUE

END SUBROUTINE SPLINE1_INTERP

SUBROUTINE SPLINE1_INTEG(k_order,n0,x0,f,n1,x1, &
                         value,iflag,message)
!-------------------------------------------------------------------------------
!SPLINE1_INTEG evaluates the integral f(x)*x**k_order, where f(x) is a cubic
!  spline function and k_order is 0 or 1
!
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall (1977) 76
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!  W.A.Houlberg, F90 free form 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_order,             & !exponent in integral (s(x)*x**k_order) [-]
                         !=0 or 1 only
  n0,                  & !number of source nodes [-]
  n1                     !number of output nodes (may be 1) [-]

REAL(KIND=rspec), INTENT(IN) :: &
  x0(:),               & !source abcissas [arb]
  x1(:),               & !output abcissas at which the integral is wanted [arb]
  f(:,:)                 !f(1,n0)=source ordinates [arb]
                         !f(2-4,n0)=spline coefficients computed by SPLINE1_FIT [arb]


!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  value(:)               !integral of f(x)*x**k_order from x0(1) to x1(i) [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,j,ido,jdo

REAL(KIND=rspec) :: &
  add,dx,sum,xnew,xold

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
IF(n1 <= 0) GOTO 9999

value(1:n1)=0

!Integral value
sum=0

!Source index and abscissa
j=1
xold=x0(j)

!-------------------------------------------------------------------------------
!Integrate over target abscissas
!-------------------------------------------------------------------------------
DO i=1,n1 !Over target abscissas

  !Find dx to next source or target abscissa
  ido=0

  DO WHILE(ido == 0) !Over source abscissas

    IF(x1(i) < x0(j+1)) THEN

      !Hit target abscissa
      xnew=x1(i)
      ido=1
      jdo=0

    ELSEIF(x1(i) == x0(j+1)) THEN

      !Hit both source and target abscissas
      xnew=x1(i)
      ido=1
      jdo=1

    ELSEIF(x1(i) > x0(j+1)) THEN

      !Hit source abscissa
      xnew=x0(j+1)
      ido=0
      jdo=1

    ENDIF

    !Integrate over dx
    dx=xnew-xold

    IF(k_order == 1) THEN

      !Integral x s(x)dx
      add= dx*(xold*f(1,j) &
          +dx*((xold*f(2,j)+f(1,j))/2 &
          +dx*((xold*f(3,j)/2+f(2,j))/3 &
          +dx*((xold*f(4,j)/3+f(3,j))/8 &
          +dx*f(4,j)/30))))

    ELSE

      !Integral s(x)dx
      add=dx*(f(1,j)+dx*(f(2,j)+dx*(f(3,j)+dx*f(4,j)/4)/3)/2)

    ENDIF

    !Add increment and update endpoint
    sum=sum+add
    xold=xnew

    !Check whether to increment source index
    IF(jdo == 1) j=j+1

    !Check whether target index is in range
    IF(j > n0) THEN

      iflag=1
      message='LINEAR1_INTEG/ERROR:target node out of range'
      GOTO 9999

    ENDIF

    !Set integral value
    value(i)=sum

  ENDDO !Over source abscissas

ENDDO !Over target abscissas
   
!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE SPLINE1_INTEG

SUBROUTINE SPLINE1_SEARCH(x,n,xl, &
                          jlo)
!-------------------------------------------------------------------------------
!SPLINE1_SEARCH is a correlated table search routine to find the indices of the
!  array x that bound xl
!
!References:
!  W.A.Houlberg, P.I.Strand, D.McCune 8/2001
!
!Comments:
!  This is similar to the Numerical Recipes routine HUNT
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n                      !number of abscissas [-]

REAL(KIND=rspec), INTENT(IN) :: &
  xl,                  & !target value [arb]
  x(:)                   !monotonically increasing array of abscissas [arb]

!Declaration of input/output variables
INTEGER, INTENT(INOUT) :: &
  jlo                    !input starting lower index [-]
                         !<1     binary search
                         !=1,n-1 use value
                         !>n-1   binary search
                         !output starting lower index [-]
                         !=0     xl < x(1) 
                         !=1     x(1) <= xl <= x(2)
                         !=2,n-1 x(jlo) < xl <= x(jlo+1)
                         !=n     x(jlo) > x(n)

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  inc,jhi,jmid

!-------------------------------------------------------------------------------
!Check lower end of array, first two points
!-------------------------------------------------------------------------------
IF(xl < x(1)) THEN

  !Target is below node 1
  jlo=0

ELSEIF(xl <= x(2)) THEN

  !Target is between nodes 1 and 2 inclusive
  jlo=1
   
!-------------------------------------------------------------------------------
!Check middle range
!-------------------------------------------------------------------------------
ELSEIF(xl <= x(n)) THEN

  !Target is between nodes 2 and n
  IF(jlo < 1 .OR. &
     jlo > (n-1)) THEN

    !jlo from previous call is unusable
    jlo=2
    jhi=n

  !Bracket target value
  ELSE

    !Start with jlo from previous call
    inc=1

    IF(xl > x(jlo)) THEN

      !Search up
      jhi=jlo+1

      DO WHILE(xl > x(jhi))

        inc=2*inc
        jlo=jhi
        jhi=MIN(jlo+inc,n)

      ENDDO

    ELSE

      !Search down
      jhi=jlo
      jlo=jlo-1

      DO WHILE(xl <= x(jlo))

        inc=inc+inc
        jhi=jlo
        jlo=MAX(jlo-inc,1)

      ENDDO

    ENDIF

  ENDIF

  !Bisection
  DO WHILE(jhi-jlo > 1)

    jmid=(jhi+jlo)/2

    IF(xl > x(jmid)) THEN

      jlo=jmid

    ELSE

      jhi=jmid

    ENDIF

  ENDDO

!-------------------------------------------------------------------------------
!Target is above node n
!-------------------------------------------------------------------------------
ELSE

  jlo=n

ENDIF

END SUBROUTINE SPLINE1_SEARCH

END MODULE SPLINE1_MOD
