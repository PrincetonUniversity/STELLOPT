MODULE LINEAR1_MOD
!-------------------------------------------------------------------------------
!LINEAR1-LINEAR interpolation in 1d
!
!LINEAR1_MOD is an F90 module of linear interpolating routines in 1d
!
!References:
!
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!
!  LINEAR1_INTERP      -interpolate from one grid to another
!  LINEAR1_INTEG       -integrate the linear fit
!
!Comments:
!
!  Linear interpolation routines are C0 (only f is continuous)
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
! Public procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE LINEAR1_INTERP(n0,x0,y0,n1,x1, &
                          y1,iflag,message)
!-------------------------------------------------------------------------------
!W_LIN_INTERP performs a linear interpolation from the x0 mesh to the x1 mesh
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  If the target mesh is outside the domain of the input mesh, the end data
!    value is returned
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n0,                  & !number of source abscissas [-]
  n1                     !number of target abscissas [-]

REAL(KIND=rspec), INTENT(IN) :: &
  x0(:),               & !source abscissas (in increasing order) [arb]
  x1(:),               & !target abscissas [arb]
  y0(:)                  !source values [arb]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &
  y1(:)                  !target values [arb]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,il

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
IF(n1 <= 0) THEN

  iflag=-1
  message='LINEAR1_INTERP/WARNING(1):no points in output array'
  GOTO 9999

ENDIF

!Make sure there are at least two nodes in the input grid
IF(n0 < 2) THEN

  iflag=1
  message='LINEAR1_INTERP/ERROR(2):<2 points in source array'
  GOTO 9999

ENDIF

!Set starting index of source grid
il=1

!-------------------------------------------------------------------------------
!Interpolate from x0 to x1
!-------------------------------------------------------------------------------
DO i=1,n1 !Over index of target grid

  10 IF(x1(i) < x0(1)) THEN

    !Target is below data range, use innermost data value
    y1(i)=y0(1)
    iflag=-1
    message='LINEAR1_INTERP(3)/WARNING:x<x(1), use end point'

  ELSEIF(x1(i) == x0(1)) THEN

    !Target and source nodes coincide
    y1(i)=y0(1)

  ELSEIF(x1(i) > x0(il+1)) THEN

    !Beyond next source node
    !Step x0 grid forward and loop
    IF(il < n0-1) THEN

      il=il+1
      GOTO 10

    ELSE

      !Target is above data range, set to last value
      y1(i)=y0(n0)
      iflag=-1
      message='LINEAR1_INTERP(4)/WARNING:x>x(n0), use end point'

    ENDIF

  ELSE

    !Between the proper set of source nodes, interpolate
    y1(i)=y0(il)+(y0(il+1)-y0(il))*(x1(i)-x0(il))/(x0(il+1)-x0(il))

  ENDIF

ENDDO !Over index of target grid

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE
      
END SUBROUTINE LINEAR1_INTERP

SUBROUTINE LINEAR1_INTEG(k_order,n0,x0,f,n1,x1, &
                         value,iflag,message)
!-------------------------------------------------------------------------------
!LINEAR1_INTEG evaluates the integral f(x)*x**k_order, where f(x) is a linear
!  function and k_order is 0 or 1
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  k_order,             & !exponent of weighting factor (s(x)*x**k_order) [-]
                         !=0 or 1 only
  n0,                  & !number of source abcissas [-]
  n1                     !number of output abcissas (may be 1) [-]

REAL(KIND=rspec), INTENT(IN) :: &
  f(:),                & !source ordinates [arb]
  x0(:),               & !source abcissas [arb]
  x1(:)                  !output abcissas [arb]

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
  add,dx,f2,sum,xnew,xold

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!If no target values are requested, return
IF(n1 <= 0) THEN

  iflag=-1
  message='LINEAR1_INTEG(1)/WARNING:no target points'
  GOTO 9999

ENDIF

!If first target point is below range, return
IF(x1(1) < 0.99999_rspec*x0(1)) THEN

  iflag=1
  message='LINEAR1_INTEG(2)/WARNING:target below range'
  GOTO 9999

ENDIF

!If first target point is above range, return
IF(x1(n1) > 1.00001_rspec*x0(n0)) THEN

  iflag=1
  message='LINEAR1_INTEG(3)/WARNING:target above range'
  GOTO 9999

ENDIF

value(1:n1)=0

!Integral value
sum=0

!Source node, value and derivative
j=1
xold=x0(j)
f2=(f(2)-f(1))/(x0(2)-x0(1))

!-------------------------------------------------------------------------------
!Integrate over target nodes
!-------------------------------------------------------------------------------
DO i=1,n1 !Over target nodes

!Find dx to next source or target nodes
  ido=0

  DO WHILE(ido == 0) !Over source nodes

    IF(x1(i) < x0(j+1)) THEN

      !Hit target nodes
      xnew=x1(i)
      ido=1
      jdo=0

    ELSEIF(x1(i) == x0(j+1)) THEN

      !Hit both source and target nodes
      xnew=x1(i)
      ido=1
      jdo=1

    ELSEIF(x1(i) > x0(j+1)) THEN

      !Hit source nodes
      xnew=x0(j+1)
      ido=0
      jdo=1

    ENDIF

!Integrate over dx
    dx=xnew-xold

    IF(k_order == 1) THEN

      !Integral x s(x)dx
      add=dx*(xold*f(j)+dx*((xold*f2+f(j))/2))

    ELSE

      !Integral s(x)dx
      add=dx*(f(j)+dx*f2/2)

    ENDIF

!Add increment and update endpoint
    sum=sum+add
    xold=xnew

    IF(jdo == 1) THEN

      !Increment source node and derivative
      j=j+1

      IF(j == n0) THEN

        j=n0-1
        f2=0

      ELSE

        f2=(f(j+1)-f(j))/(x0(j+1)-x0(j))

      ENDIF

    ENDIF

!Set integral value
    value(i)=sum

  ENDDO !Over source nodes

ENDDO !Over target nodes
   
!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE LINEAR1_INTEG

END MODULE LINEAR1_MOD
