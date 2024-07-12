c  cspline -- dmc 15 Feb 1999
c
c  a standard interface to the 1d spline setup routine
c    modified dmc 3 Mar 2000 -- to use Wayne Houlberg's v_spline code.
c    new BC options added.
c
      subroutine cspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   wk,iwk,ilinx,ier)
c
      real x(nx)                        ! x axis (in)
      real fspl(4,nx)                   ! spline data (in/out)
      integer ibcxmin                   ! x(1) BC flag (in, see comments)
      real bcxmin                       ! x(1) BC data (in, see comments)
      integer ibcxmax                   ! x(nx) BC flag (in, see comments)
      real bcxmax                       ! x(nx) BC data (in, see comments)
      real wk(iwk)                      ! workspace of size at least nx
      integer ilinx                     ! even spacing flag (out)
      integer ier                       ! output, =0 means OK
c
c  ** note wk(...) array is not used unless ibcxmin=-1 (periodic spline
c  evaluation)
c
c  this routine computes spline coefficients for a 1d spline --
c  evaluation of the spline can be done by cspeval.for subroutines
c  or directly by inline code.
c
c  the input x axis x(1...nx) must be strictly ascending, i.e.
c  x(i+1).gt.x(i) is required for i=1 to nx-1.  This is checked and
c  ier=1 is set and the routine exits if the test is not satisfied.
c
c  on output, ilinx=1 is set if, to a reasonably close tolerance,
c  all grid spacings x(i+1)-x(i) are equal.  This allows a speedier
c  grid lookup algorithm on evaluation of the spline.  If on output
c  ilinx=2, this means the spline x axis is not evenly spaced.
c
c  the input data for the spline are given in f[j] = fspl(1,j).  The
c  output data are the spline coefficients fspl(2,j),fspl(3,j), and
c  fspl(4,j), j=1 to nx.  The result is a spline s(x) satisfying the
c  boundary conditions and with the properties
c
c     s(x(j)) = fspl(1,j)
c     s'(x) is continuous even at the grid points x(j)
c     s''(x) is continuous even at the grid points x(j)
c
c  the formula for evaluation of s(x) is:
c
c     let dx = x-x(i), where x(i).le.x.le.x(i+1).  Then,
c     s(x)=fspl(1,i) + dx*(fspl(2,i) +dx*(fspl(3,i) + dx*fspl(4,i)))
c
c  ==>boundary conditions.  Complete specification of a 1d spline
c  requires specification of boundary conditions at x(1) and x(nx).
c
c  this routine provides 4 options:
c
c -1 ***** PERIODIC BC
c  ibcxmin=-1  --  periodic boundary condition.  This means the
c    boundary conditions s'(x(1))=s'(x(nx)) and s''(x(1))=s''(x(nx))
c    are imposed.  Note that s(x(1))=s(x(nx)) (i.e. fspl(1,1)=fspl(1,nx))
c    is not required -- that is determined by the fspl array input data.
c    The periodic boundary condition is to be preferred for periodic
c    data.  When splining periodic data f(x) with period P, the relation
c    x(nx)=x(1)+n*P, n = the number of periods (usually 1), should hold.
c    (ibcxmax, bcxmin, bcxmax are ignored).
c
c  if a periodic boundary condition is set, this covers both boundaries.
c  for the other types of boundary conditions, the type of condition
c  chosen for the x(1) boundary need not be the same as the type chosen
c  for the x(nx) boundary.
c
c  0 ***** NOT A KNOT BC
c  ibcxmin=0 | ibcxmax=0 -- this specifies a "not a knot" boundary
c    condition -- see cubsplb.for.  This is a common way for inferring
c    a "good" spline boundary condition automatically from data in the
c    vicinity of the boundary.  (bcxmin | bcxmax are ignored).
c
c  1 ***** BC:  SPECIFIED SLOPE
c  ibcxmin=1 | ibcxmax=1 -- boundary condition is to have s'(x(1)) |
c    s'(x(nx)) match the passed value (bcxmin | bcxmax).
c
c  2 ***** BC:  SPECIFIED 2nd DERIVATIVE
c  ibcxmin=2 | ibcxmax=2 -- boundary condition is to have s''(x(1)) |
c    s''(x(nx)) match the passed value (bcxmin | bcxmax).
c
c  3 ***** BC:  SPECIFIED SLOPE = 0.0
c  ibcxmin=3 | ibcxmax=3 -- boundary condition is to have s'(x(1)) |
c    s'(x(nx)) equal to ZERO.
c
c  4 ***** BC:  SPECIFIED 2nd DERIVATIVE = 0.0
c  ibcxmin=4 | ibcxmax=4 -- boundary condition is to have s''(x(1)) |
c    s''(x(nx)) equal to ZERO.
c
c  5 ***** BC:  1st DIVIDED DIFFERENCE
c  ibcxmin=5 | ibcxmax=5 -- boundary condition is to have s'(x(1)) |
c    s'(x(nx)) equal to the slope from the 1st|last 2 points
c
c  6 ***** BC:  2nd DIVIDED DIFFERENCE
c  ibcxmin=6 | ibcxmax=6 -- boundary condition is to have s''(x(1)) |
c    s''(x(nx)) equal to the 2nd derivative from the 1st|last 3 points
c
c  7 ***** BC:  3rd DIVIDED DIFFERENCE
c  ibcxmin=7 | ibcxmax=7 -- boundary condition is to have s'''(x(1)) |
c    s'''(x(nx)) equal to the 3rd derivative from the 1st|last 4 points
c
c---------------------------------------------------------------------
      data half/0.5/
      data sixth/0.166666666666666667/
c
c  error checks
c
      ier = 0
      if(nx.lt.2) then
         write(6,'('' ?cspline:  at least 2 x points required.'')')
         ier=1
      endif
      call ibc_ck(ibcxmin,'cspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'cspline','xmax',0,7,ier)
c
c  x axis check
c
      call splinck(x,nx,ilinx,1.0e-3,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?cspline:  x axis not strict ascending'')')
      endif
c
      if(ibcxmin.eq.-1) then
         inum=nx
         if(iwk.lt.inum) then
            write(6,1009) inum,iwk,nx
 1009       format(
     >      ' ?cspline:  workspace too small.  need:  ',i6,' got:  ',i6/
     >      '  (need = nx, nx=',i6)
            ier=3
         endif
      endif
c
      if(ier.ne.0) return
c
c  OK -- evaluate spline
c
      if(ibcxmin.eq.1) then
         fspl(2,1)=bcxmin
      else if(ibcxmin.eq.2) then
         fspl(3,1)=bcxmin
      endif
c
      if(ibcxmax.eq.1) then
         fspl(2,nx)=bcxmax
      else if(ibcxmax.eq.2) then
         fspl(3,nx)=bcxmax
      endif
c
      call v_spline(ibcxmin,ibcxmax,nx,x,fspl,wk)
c
      do i=1,nx
         fspl(3,i)=half*fspl(3,i)
         fspl(4,i)=sixth*fspl(4,i)
      enddo
c
      return
      end
