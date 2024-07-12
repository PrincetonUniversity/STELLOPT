      SUBROUTINE r8v_spline(k_bc1,k_bcn,n,x,f,wk)
!***********************************************************************
!V_SPLINE evaluates the coefficients for a 1d cubic interpolating spline
!References:
!  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!    Computations, Prentice-Hall, 1977, p.76
!  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
!    1996, p.251
!  W.A.Houlberg, D.McCune 3/2000
!Input:
!  k_bc1-option for BC at x1 = x(1)
!       =-1 periodic, ignore k_bcn
!       =0 not-a-knot
!       =1 s'(x1) = input value of f(2,1)
!       =2 s''(x1) = input value of f(3,1)
!       =3 s'(x1) = 0.0
!       =4 s''(x1) = 0.0
!       =5 match first derivative to first 2 points
!       =6 match second derivative to first 3 points
!       =7 match third derivative to first 4 points
!       =else use not-a-knot
!  k_bcn-option for boundary condition at xn = x(n)
!       =0 not-a-knot
!       =1 s'(xn) = input value of f(2,n)
!       =2 s''(xn) = input value of f(3,n)
!       =3 s'(xn) = 0.0
!       =4 s''(xn) = 0.0
!       =5 match first derivative to last 2 points
!       =6 match second derivative to lasst 3 points
!       =7 match third derivative to last 4 points
!       =else use knot-a-knot
!  n-number of data points or knots-(n.ge.2)
!  x(n)-abscissas of the knots in strictly increasing order
!  f(1,i)-ordinates of the knots
!  f(2,1)-input value of s'(x1) for k_bc1=1
!  f(2,n)-input value of s'(xn) for k_bcn=1
!  f(3,1)-input value of s''(x1) for k_bc1=2
!  f(3,n)-input value of s''(xn) for k_bcn=2
!  wk(n)-scratch work area for periodic BC
!Output:
!  f(2,i)=s'(x(i))
!  f(3,i)=s''(x(i))
!  f(4,i)=s'''(x(i))
!Comments:
!  s(x)=f(1,i)+f(2,i)*(x-x(i))+f(3,i)*(x-x(i))**2/2!
!       +f(4,i)*(x-x(i))**3/3! for x(i).le.x.le.x(i+1)
!  W_SPLINE can be used to evaluate the spline and its derivatives
!  The cubic spline is twice differentiable (C2)
!
!  modifications -- dmc 18 Feb 2010:
!    Deal with n.lt.4 -- the general tridiagonal spline method
!    does not have the right formulation for n.eq.3 "not a knot" or periodic
!    boundary conditions, nor for n.eq.2 with any boundary conditions.
!
!    Apply boundary conditions even for n=2, when the "spline" is really
!    just a single cubic polynomial.
!    In this case, several boundary condition (BC) options are mapped to
!      BC option 5, 1st divided difference.  If 5 is used on both sides
!      of an n=2 "spline" you get a linear piece which is what the old
!      code always gave, regardless of BC option settings.  The following
!      BC controls are mapped to 5:
!        periodic (-1)
!        not a knot (0) (for n=2 no grid point exists for knot location).
!        option (5) is preserved
!        options 6 and 7 -- mapped to (5); higher divided differences
!          need n>2; in fact 7 needs n>3; for n=3 option 6 is substituted.
!
!    The boundary condition screening is done at the start of the code;
!    passed controls k_bc1 and k_bcn are mapped to i_bc1 and i_bcn.
!
!    ALSO: for n=3, "not a knot" from both left and right needs special
!      interpretation, since the 2 boundary conditions overlap.  The
!      chosen interpretation is a parabolic fit to the 3 given data points.
!      and so f''' = 0 and f'' = constant.  If "not a knot" is used on
!      one side only, the solution is a single cubic piece and special
!      code is also needed.
!    ALSO: for n=3, "periodic" boundary condition needs special code; this
!      is added.
!
!  bugfixes -- dmc 24 Feb 2004:
!    (a) fixed logic for not-a-knot:
!          !    Set f(3,1) for not-a-knot
!                    IF(k_bc1.le.0.or.k_bc1.gt.7) THEN ...
!        instead of
!          !    Set f(3,1) for not-a-knot
!                    IF(k_bc1.le.0.or.k_bc1.gt.5) THEN ...
!        and similarly for logic after cmt
!          !    Set f(3,n) for not-a-knot
!        as required since k_bc*=6 and k_bc*=7 are NOT not-a-knot BCs.
!
!    (b) the BCs to fix 2nd derivative at end points did not work if that
!        2nd derivative were non-zero.  The reason is that in those cases
!        the off-diagonal matrix elements nearest the corners are not
!        symmetric; i.e. elem(1,2).ne.elem(2,1) and
!        elem(n-1,n).ne.elem(n,n-1) where I use "elem" to refer to
!        the tridiagonal matrix elements.  The correct values for the
!        elements is:   elem(1,2)=0, elem(2,1)=x(2)-x(1)
!                       elem(n,n-1)=0, elem(n-1,n)=x(n)-x(n-1)
!        the old code in effect had these as all zeroes.  Since this
!        meant the wrong set of linear equations was solved, the
!        resulting spline had a discontinuity in its 1st derivative
!        at x(2) and x(n-1).  Fixed by introducing elem21 and elemnn1
!        to represent the non-symmetric lower-diagonal values.  Since
!        elem21 & elemnn1 are both on the lower diagonals, logic to
!        use them occurs in the non-periodic forward elimination loop
!        only.  DMC 24 Feb 2004.
!***********************************************************************
!Declaration of input variables
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER        k_bc1,                   k_bcn,
     &               n
      REAL*8           x(*),                    wk(*),
     &               f(4,*)
!Declaration in local variables
      INTEGER        i,                       ib,
     &               imax,                    imin
      REAL*8           a1,                      an,
     &               b1,                      bn,
     &               q,                       t,
     &               hn
      REAL*8           elem21,                  elemnn1    ! (dmc)
 
      integer :: i_bc1,i_bcn  ! screened BC controls
 
      integer :: iord1,iord2  ! used for n=2 only
      REAL*8 :: h,f0,fh         ! used for n=2,3 only
      REAL*8 :: h1,h2,h3,dels   ! used for n=3 special cases
      REAL*8 :: f1,f2,f3,aa,bb  ! used for n=3
 
      integer :: i3knots      ! for n=3, number of not-a-knot BCs
      integer :: i3perio      ! for n=3, periodic BC
 
!------------------------------------------------------------
!  screen the BC options (DMC Feb. 2010...)
 
      i_bc1=k_bc1
      i_bcn=k_bcn
 
      if((i_bc1.lt.-1).or.(i_bc1.gt.7)) i_bc1=0  ! outside [-1:7] -> not-a-knot
      if((i_bcn.lt.0).or.(i_bcn.gt.7)) i_bcn=0   ! outside [0:7] -> not-a-knot
 
      if(i_bc1.eq.-1) i_bcn=-1  ! periodic BC
 
      i3knots=0
      i3perio=0
      if(n.eq.3) then
         i_bc1=min(6,i_bc1)
         i_bcn=min(6,i_bcn)
         if(i_bc1.eq.0) i3knots = i3knots + 1
         if(i_bcn.eq.0) i3knots = i3knots + 1
         if(i_bc1.eq.-1) i3perio = 1
      endif
 
      if(n.eq.2) then
         if(i_bc1.eq.-1) then
            i_bc1=5
            i_bcn=5
         endif
         if((i_bc1.eq.0).or.(i_bc1.gt.5)) i_bc1=5
         if((i_bcn.eq.0).or.(i_bcn.gt.5)) i_bcn=5
 
         if((i_bc1.eq.1).or.(i_bc1.eq.3).or.(i_bc1.eq.5)) then
            iord1=1  ! 1st derivative match on LHS
         else
            iord1=2  ! 2nd derivative match on LHS
         endif
 
         if((i_bcn.eq.1).or.(i_bcn.eq.3).or.(i_bcn.eq.5)) then
            iord2=1  ! 1st derivative match on RHS
         else
            iord2=2  ! 2nd derivative match on RHS
         endif
      endif
 
!Set default range
      imin=1
      imax=n
!Set first and second BC values
      a1=0.0_r8
      b1=0.0_r8
      an=0.0_r8
      bn=0.0_r8
      IF(i_bc1.eq.1) THEN
        a1=f(2,1)
      ELSEIF(i_bc1.eq.2) THEN
        b1=f(3,1)
      ELSEIF(i_bc1.eq.5) THEN
        a1=(f(1,2)-f(1,1))/(x(2)-x(1))
      ELSEIF(i_bc1.eq.6) THEN
        b1=2.0_r8*((f(1,3)-f(1,2))/(x(3)-x(2))
     &         -(f(1,2)-f(1,1))/(x(2)-x(1)))/(x(3)-x(1))
      ENDIF
      IF(i_bcn.eq.1) THEN
        an=f(2,n)
      ELSEIF(i_bcn.eq.2) THEN
        bn=f(3,n)
      ELSEIF(i_bcn.eq.5) THEN
        an=(f(1,n)-f(1,n-1))/(x(n)-x(n-1))
      ELSEIF(i_bcn.eq.6) THEN
        bn=2.0_r8*((f(1,n)-f(1,n-1))/(x(n)-x(n-1))
     &         -(f(1,n-1)-f(1,n-2))/(x(n-1)-x(n-2)))/(x(n)-x(n-2))
      ENDIF
!Clear f(2:4,n)
      f(2,n)=0.0_r8
      f(3,n)=0.0_r8
      f(4,n)=0.0_r8
      IF(n.eq.2) THEN
        if((i_bc1.eq.5).and.(i_bcn.eq.5)) then
!Coefficients for n=2 (this was the original code)
           f(2,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
           f(3,1)=0.0_r8
           f(4,1)=0.0_r8
           f(2,2)=f(2,1)
           f(3,2)=0.0_r8
           f(4,2)=0.0_r8
        else if((iord1.eq.1).and.(iord2.eq.1)) then
           ! LHS: match a1 for 1st deriv; RHS: match an for 1st deriv.
           f(2,1)=a1
           f(2,2)=an
           h = (x(2)-x(1))
           f0 = f(1,1)
           fh = f(1,2)
 
           ! setting xx = x-x(1),
           ! f = c1*xx**3 + c2*xx**2 + a1*xx + f0
           !   -->  c1*h**3   + c2*h**2 = fh - f0 - a1*h
           !   and  3*c1*h**2 + 2*c2*h  = an - a1
           ! f' = 3*c1*xx*2 + 2*c2*xx + a1
           ! f'' = 6*c1*xx + 2*c2
           ! f''' = 6*c1
 
           ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2
 
           f(3,1) = (3*(fh-f0)/(h*h) - (2*a1 + an)/h)*2       ! 2*c2
           f(4,1) = (-2*(fh-f0)/(h*h*h) + (a1 + an)/(h*h))*6  ! 6*c1
 
           f(4,2) = f(4,1)
           f(3,2) = f(4,1)*h + f(3,1)
 
        else if((iord1.eq.1).and.(iord2.eq.2)) then
           ! LHS: match a1 for 1st deriv; RHS: match bn for 2nd deriv.
           f(2,1)=a1
           f(3,2)=bn
           h = (x(2)-x(1))
           f0 = f(1,1)
           fh = f(1,2)
 
           ! setting xx = x-x(1),
           ! f = c1*xx**3 + c2*xx**2 + a1*xx + f0
           !   -->  c1*h**3   + c2*h**2 = fh - f0 - a1*h
           !   and  6*c1*h    + 2*c2    = bn
           ! f' = 3*c1*xx*2 + 2*c2*xx + a1
           ! f'' = 6*c1*xx + 2*c2
           ! f''' = 6*c1
 
           ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2
 
           f(3,1) = (-bn/4 + 3*(fh-f0)/(2*h*h) - 3*a1/(2*h))*2       ! 2*c2
           f(4,1) = (bn/(4*h) - (fh-f0)/(2*h*h*h) + a1/(2*h*h))*6    ! 6*c1
 
           f(4,2) = f(4,1)
           f(2,2) = f(4,1)*h*h/2 + f(3,1)*h + a1
        else if((iord1.eq.2).and.(iord2.eq.1)) then
           ! LHS: match b1 for 2nd deriv; RHS: match an for 1st deriv.
           f(3,1)=b1
           f(2,2)=an
           h = (x(2)-x(1))
           f0 = f(1,1)
           fh = f(1,2)
 
           ! setting xx = x-x(1),
           ! f = c1*xx**3 + (b1/2)*xx**2 + c3*xx + f0
           !   -->  c1*h**3   + c3*h = fh - f0 - b1*h**2/2
           !   and  3*c1*h**2 + c3   = an - b1*h
           ! f' = 3*c1*xx*2 + b1*xx + c3
           ! f'' = 6*c1*xx + b1
           ! f''' = 6*c1
 
           ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2
 
           f(2,1) = 3*(fh-f0)/(2*h) - b1*h/4 - an/2                  ! c3
           f(4,1) = (an/(2*h*h) - (fh-f0)/(2*h*h*h) - b1/(4*h))*6    ! 6*c1
 
           f(4,2) = f(4,1)
           f(3,2) = f(4,1)*h + f(3,1)
        else if((iord1.eq.2).and.(iord2.eq.2)) then
           ! LHS: match b1 for 2nd deriv; RHS: match bn for 2nd deriv.
           f(3,1)=b1
           f(3,2)=bn
           h = (x(2)-x(1))
           f0 = f(1,1)
           fh = f(1,2)
 
           ! setting xx = x-x(1),
           ! f = c1*xx**3 + (b1/2)*xx**2 + c3*xx + f0
           !   -->  c1*h**3   + c3*h = fh - f0 - b1*h**2/2
           !   and  6*c1*h           = bn - b1
           ! f' = 3*c1*xx*2 + b1*xx + c3
           ! f'' = 6*c1*xx + b1
           ! f''' = 6*c1
 
           ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2
 
           f(2,1) = (fh-f0)/h -b1*h/3 -bn*h/6                 ! c3
           f(4,1) = (bn-b1)/h                                 ! 6*c1
 
           f(4,2) = f(4,1)
           f(2,2) = f(4,1)*h*h/2 + b1*h + f(2,1)
        endif
 
      ELSE IF(i3perio.eq.1) then
!Special case: nx=3 periodic spline
         h1=x(2)-x(1)
         h2=x(3)-x(2)
         h=h1+h2
 
         dels=(f(1,3)-f(1,2))/h2 - (f(1,2)-f(1,1))/h1
 
         f(2,1)= (f(1,2)-f(1,1))/h1 + (h1*dels)/h
         f(3,1)= -6*dels/h
         f(4,1)= 12*dels/(h1*h)
 
         f(2,2)= (f(1,3)-f(1,2))/h2 - (h2*dels)/h
         f(3,2)= 6*dels/h
         f(4,2)= -12*dels/(h2*h)
 
         f(2,3)=f(2,1)
         f(3,3)=f(3,1)
         f(4,3)=f(4,1)
 
 
      ELSE IF(i3knots.eq.2) then
!Special case: nx=3, not-a-knot on both sides
         h1=x(2)-x(1)
         h2=x(3)-x(2)
         h=h1+h2
                                ! this is just a quadratic fit through 3 pts
         f1=f(1,1)-f(1,2)
         f2=f(1,3)-f(1,2)
 
!  quadratic around origin at (x(2),f(1,2))
!          aa*h1**2 - bb*h1 = f1
!          aa*h2**2 + bb*h2 = f2
 
         aa = (f2*h1 + f1*h2)/(h1*h2*h)
         bb = (f2*h1*h1 - f1*h2*h2)/(h1*h2*h)
 
         f(4,1:3)=0.0_r8  ! f''' = 0 (quadratic polynomial)
         f(3,1:3)=2*aa  ! f'' = const
 
         f(2,1)=bb-2*aa*h1
         f(2,2)=bb
         f(2,3)=bb+2*aa*h2
 
      ELSE IF(i3knots.eq.1) then
!Special cases: nx=3, not-a-knot on single side
         if((i_bc1.eq.1).or.(i_bc1.eq.3).or.(i_bc1.eq.5)) then
                                ! f' LHS condition; not-a-knot RHS
!  a1 = f' LHS BC
            h2=x(2)-x(1)
            h3=x(3)-x(1)
 
            f2=f(1,2)-f(1,1)
            f3=f(1,3)-f(1,1)
 
!  find cubic aa*xx**3 + bb*xx**2 + a1*xx
!    satisfying aa*h2**3 + bb*h2**2 + a1*h2 = f2
!           and aa*h3**3 + bb*h3**2 + a1*h3 = f3
 
            aa=a1/(h2*h3) + f3/(h3*h3*(h3-h2)) - f2/(h2*h2*(h3-h2))
            bb=-a1*(h3*h3-h2*h2)/(h2*h3*(h3-h2))
     >           + f2*h3/(h2*h2*(h3-h2)) - f3*h2/(h3*h3*(h3-h2))
 
            f(2,1)=a1
            f(3,1)=2*bb
            f(4,1)=6*aa
 
            f(2,2)=3*aa*h2*h2 + 2*bb*h2 + a1
            f(3,2)=6*aa*h2 + 2*bb
            f(4,2)=6*aa
 
            f(2,3)=3*aa*h3*h3 + 2*bb*h3 + a1
            f(3,3)=6*aa*h3 + 2*bb
            f(4,3)=6*aa
 
         else if((i_bc1.eq.2).or.(i_bc1.eq.4).or.(i_bc1.eq.6)) then
                                ! f'' LHS condition; not-a-knot RHS
!  b1 = f'' LHS BC
            h2=x(2)-x(1)
            h3=x(3)-x(1)
 
            f2=f(1,2)-f(1,1)
            f3=f(1,3)-f(1,1)
 
!  find cubic aa*xx**3 + (b1/2)*xx**2 + bb*xx
!    satisfying aa*h2**3 + bb*h2 = f2 -(b1/2)*h2**2
!           and aa*h3**3 + bb*h3 = f3 -(b1/2)*h3**2
 
            aa= -(b1/2)*(h3-h2)/(h3*h3-h2*h2)
     >           -f2/(h2*(h3*h3-h2*h2)) + f3/(h3*(h3*h3-h2*h2))
            bb= -(b1/2)*h2*h3*(h3-h2)/(h3*h3-h2*h2)
     >           +f2*h3*h3/(h2*(h3*h3-h2*h2))
     >           -f3*h2*h2/(h3*(h3*h3-h2*h2))
 
            f(2,1)=bb
            f(3,1)=b1
            f(4,1)=6*aa
 
            f(2,2)=3*aa*h2*h2 + b1*h2 + bb
            f(3,2)=6*aa*h2 + b1
            f(4,2)=6*aa
 
            f(2,3)=3*aa*h3*h3 + b1*h3 + bb
            f(3,3)=6*aa*h3 + b1
            f(4,3)=6*aa
 
         else if((i_bcn.eq.1).or.(i_bcn.eq.3).or.(i_bcn.eq.5)) then
                                ! f' RHS condition; not-a-knot LHS
!  an = f' RHS BC
            h2=x(2)-x(3)
            h3=x(1)-x(3)
 
            f2=f(1,2)-f(1,3)
            f3=f(1,1)-f(1,3)
 
!  find cubic aa*xx**3 + bb*xx**2 + an*xx
!    satisfying aa*h2**3 + bb*h2**2 + an*h2 = f2
!           and aa*h3**3 + bb*h3**2 + an*h3 = f3
 
            aa=an/(h2*h3) + f3/(h3*h3*(h3-h2)) - f2/(h2*h2*(h3-h2))
            bb=-an*(h3*h3-h2*h2)/(h2*h3*(h3-h2))
     >           + f2*h3/(h2*h2*(h3-h2)) - f3*h2/(h3*h3*(h3-h2))
 
            f(2,3)=an
            f(3,3)=2*bb
            f(4,3)=6*aa
 
            f(2,2)=3*aa*h2*h2 + 2*bb*h2 + an
            f(3,2)=6*aa*h2 + 2*bb
            f(4,2)=6*aa
 
            f(2,1)=3*aa*h3*h3 + 2*bb*h3 + an
            f(3,1)=6*aa*h3 + 2*bb
            f(4,1)=6*aa
 
         else if((i_bcn.eq.2).or.(i_bcn.eq.4).or.(i_bcn.eq.6)) then
                                ! f'' RHS condition; not-a-knot LHS
!  bn = f'' RHS BC
            h2=x(2)-x(3)
            h3=x(1)-x(3)
 
            f2=f(1,2)-f(1,3)
            f3=f(1,1)-f(1,3)
 
!  find cubic aa*xx**3 + (bn/2)*xx**2 + bb*xx
!    satisfying aa*h2**3 + bb*h2 = f2 -(bn/2)*h2**2
!           and aa*h3**3 + bb*h3 = f3 -(bn/2)*h3**2
 
            aa= -(bn/2)*(h3-h2)/(h3*h3-h2*h2)
     >           -f2/(h2*(h3*h3-h2*h2)) + f3/(h3*(h3*h3-h2*h2))
            bb= -(bn/2)*h2*h3*(h3-h2)/(h3*h3-h2*h2)
     >           +f2*h3*h3/(h2*(h3*h3-h2*h2))
     >           -f3*h2*h2/(h3*(h3*h3-h2*h2))
 
            f(2,3)=bb
            f(3,3)=bn
            f(4,3)=6*aa
 
            f(2,2)=3*aa*h2*h2 + bn*h2 + bb
            f(3,2)=6*aa*h2 + bn
            f(4,2)=6*aa
 
            f(2,1)=3*aa*h3*h3 + bn*h3 + bb
            f(3,1)=6*aa*h3 + bn
            f(4,1)=6*aa
 
         endif
      ELSE IF(n.gt.2) THEN
!Set up tridiagonal system for A*y=B where y(i) are the second
!  derivatives at the knots
!  f(2,i) are the diagonal elements of A
!  f(4,i) are the off-diagonal elements of A
!  f(3,i) are the B elements/3, and will become c/3 upon solution
        f(4,1)=x(2)-x(1)
        f(3,2)=(f(1,2)-f(1,1))/f(4,1)
        DO i=2,n-1
          f(4,i)=x(i+1)-x(i)
          f(2,i)=2.0_r8*(f(4,i-1)+f(4,i))
          f(3,i+1)=(f(1,i+1)-f(1,i))/f(4,i)
          f(3,i)=f(3,i+1)-f(3,i)
        ENDDO
!
!  (dmc): save now:
!
        elem21=f(4,1)
        elemnn1=f(4,n-1)
!
!  BC's
!    Left
        IF(i_bc1.eq.-1) THEN
          f(2,1)=2.0_r8*(f(4,1)+f(4,n-1))
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-(f(1,n)-f(1,n-1))/f(4,n-1)
          wk(1)=f(4,n-1)
          DO i=2,n-3
            wk(i)=0.0_r8
          ENDDO
          wk(n-2)=f(4,n-2)
          wk(n-1)=f(4,n-1)
        ELSEIF(i_bc1.eq.1.or.i_bc1.eq.3.or.i_bc1.eq.5) THEN
          f(2,1)=2.0_r8*f(4,1)
          f(3,1)=(f(1,2)-f(1,1))/f(4,1)-a1
        ELSEIF(i_bc1.eq.2.or.i_bc1.eq.4.or.i_bc1.eq.6) THEN
          f(2,1)=2.0_r8*f(4,1)
          f(3,1)=f(4,1)*b1/3.0_r8
          f(4,1)=0.0_r8  ! upper diagonal only (dmc: cf elem21)
        ELSEIF(i_bc1.eq.7) THEN
          f(2,1)=-f(4,1)
          f(3,1)=f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
          f(3,1)=f(3,1)*f(4,1)**2/(x(4)-x(1))
        ELSE                             ! not a knot:
          imin=2
          f(2,2)=f(4,1)+2.0_r8*f(4,2)
          f(3,2)=f(3,2)*f(4,2)/(f(4,1)+f(4,2))
        ENDIF
!    Right
        IF(i_bcn.eq.1.or.i_bcn.eq.3.or.i_bcn.eq.5) THEN
          f(2,n)=2.0_r8*f(4,n-1)
          f(3,n)=-(f(1,n)-f(1,n-1))/f(4,n-1)+an
        ELSEIF(i_bcn.eq.2.or.i_bcn.eq.4.or.i_bcn.eq.6) THEN
          f(2,n)=2.0_r8*f(4,n-1)
          f(3,n)=f(4,n-1)*bn/3.0_r8
!xxx          f(4,n-1)=0.0  ! dmc: preserve f(4,n-1) for back subst.
          elemnn1=0.0_r8  !  lower diaganol only (dmc)
        ELSEIF(i_bcn.eq.7) THEN
          f(2,n)=-f(4,n-1)
          f(3,n)=f(3,n-1)/(x(n)-x(n-2))-f(3,n-2)/(x(n-1)-x(n-3))
          f(3,n)=-f(3,n)*f(4,n-1)**2/(x(n)-x(n-3))
        ELSEIF(i_bc1.ne.-1) THEN         ! not a knot:
          imax=n-1
          f(2,n-1)=2.0_r8*f(4,n-2)+f(4,n-1)
          f(3,n-1)=f(3,n-1)*f(4,n-2)/(f(4,n-1)+f(4,n-2))
        ENDIF
        IF(i_bc1.eq.-1) THEN
!Solve system of equations for second derivatives at the knots
!  Periodic BC
!    Forward elimination
          DO i=2,n-2
            t=f(4,i-1)/f(2,i-1)
            f(2,i)=f(2,i)-t*f(4,i-1)
            f(3,i)=f(3,i)-t*f(3,i-1)
            wk(i)=wk(i)-t*wk(i-1)
            q=wk(n-1)/f(2,i-1)
            wk(n-1)=-q*f(4,i-1)
            f(2,n-1)=f(2,n-1)-q*wk(i-1)
            f(3,n-1)=f(3,n-1)-q*f(3,i-1)
          ENDDO
!    Correct the n-1 element
          wk(n-1)=wk(n-1)+f(4,n-2)
!    Complete the forward elimination
!    wk(n-1) and wk(n-2) are the off-diag elements of the lower corner
          t=wk(n-1)/f(2,n-2)
          f(2,n-1)=f(2,n-1)-t*wk(n-2)
          f(3,n-1)=f(3,n-1)-t*f(3,n-2)
!    Back substitution
          f(3,n-1)=f(3,n-1)/f(2,n-1)
          f(3,n-2)=(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)
          DO ib=3,n-1
            i=n-ib
            f(3,i)=(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)
          ENDDO
          f(3,n)=f(3,1)
        ELSE
!  Non-periodic BC
!    Forward elimination
!    For Not-A-Knot BC the off-diagonal end elements are not equal
          DO i=imin+1,imax
            IF((i.eq.n-1).and.(imax.eq.n-1)) THEN
              t=(f(4,i-1)-f(4,i))/f(2,i-1)
            ELSE
              if(i.eq.2) then
                 t=elem21/f(2,i-1)
              else if(i.eq.n) then
                 t=elemnn1/f(2,i-1)
              else
                 t=f(4,i-1)/f(2,i-1)
              endif
            ENDIF
            IF((i.eq.imin+1).and.(imin.eq.2)) THEN
              f(2,i)=f(2,i)-t*(f(4,i-1)-f(4,i-2))
            ELSE
              f(2,i)=f(2,i)-t*f(4,i-1)
            ENDIF
            f(3,i)=f(3,i)-t*f(3,i-1)
          ENDDO
!    Back substitution
          f(3,imax)=f(3,imax)/f(2,imax)
          DO ib=1,imax-imin
            i=imax-ib
            IF((i.eq.2).and.(imin.eq.2)) THEN
              f(3,i)=(f(3,i)-(f(4,i)-f(4,i-1))*f(3,i+1))/f(2,i)
            ELSE
              f(3,i)=(f(3,i)-f(4,i)*f(3,i+1))/f(2,i)
            ENDIF
          ENDDO
!    Reset d array to step size
          f(4,1)=x(2)-x(1)
          f(4,n-1)=x(n)-x(n-1)
!    Set f(3,1) for not-a-knot
          IF(i_bc1.le.0.or.i_bc1.gt.7) THEN
            f(3,1)=(f(3,2)*(f(4,1)+f(4,2))-f(3,3)*f(4,1))/f(4,2)
          ENDIF
!    Set f(3,n) for not-a-knot
          IF(i_bcn.le.0.or.i_bcn.gt.7) THEN
            f(3,n)=f(3,n-1)+(f(3,n-1)-f(3,n-2))*f(4,n-1)/f(4,n-2)
          ENDIF
        ENDIF
!f(3,i) is now the sigma(i) of the text and f(4,i) is the step size
!Compute polynomial coefficients
        DO i=1,n-1
          f(2,i)=
     >        (f(1,i+1)-f(1,i))/f(4,i)-f(4,i)*(f(3,i+1)+2.0_r8*f(3,i))
          f(4,i)=(f(3,i+1)-f(3,i))/f(4,i)
          f(3,i)=6.0_r8*f(3,i)
          f(4,i)=6.0_r8*f(4,i)
        ENDDO
        IF(i_bc1.eq.-1) THEN
          f(2,n)=f(2,1)
          f(3,n)=f(3,1)
          f(4,n)=f(4,1)
        ELSE
           hn=x(n)-x(n-1)
           f(2,n)=f(2,n-1)+hn*(f(3,n-1)+0.5_r8*hn*f(4,n-1))
           f(3,n)=f(3,n-1)+hn*f(4,n-1)
           f(4,n)=f(4,n-1)
           IF(i_bcn.eq.1.or.i_bcn.eq.3.or.i_bcn.eq.5) THEN
              f(2,n)=an
           ELSE IF(i_bcn.eq.2.or.i_bcn.eq.4.or.i_bcn.eq.6) THEN
              f(3,n)=bn
           ENDIF
        ENDIF
      ENDIF
      RETURN
      END
