      subroutine mkspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ilinx,ier)
C
C  make a 2-coefficient 1d spline
C
C  only 2 coefficients, the data and its 2nd derivative, are needed to
C  fully specify a spline.  See e.g. Numerical Recipies in Fortran-77
C  (2nd edition) chapter 3, section on cubic splines.
C
C  input:
      integer nx                        ! no. of data points
      real x(nx)                        ! x axis data, strict ascending order
C
C  input/output:
      real fspl(2,nx)                   ! f(1,*):  data in; f(2,*):  coeffs out
C
C     f(1,j) = f(x(j))  on input (unchanged on output)
C     f(2,j) = f''(x(j)) (of interpolating spline) (on output).
C
C  ...boundary conditions...
C
C  input:
C
      integer ibcxmin                   ! b.c. flag @ x=xmin=x(1)
      real bcxmin                       ! b.c. data @xmin
C
      integer ibcxmax                   ! b.c. flag @ x=xmax=x(nx)
      real bcxcmax                      ! b.c. data @xmax
C
C  ibcxmin=-1 -- periodic boundary condition
C                (bcxmin,ibcxmax,bcxmax are ignored)
C
C                the output spline s satisfies
C                s'(x(1))=s'(x(nx)) ..and.. s''(x(1))=s''(x(nx))
C
C  if non-periodic boundary conditions are used, then the xmin and xmax
C  boundary conditions can be specified independently:
C
C  ibcxmin (ibcxmax) = 0 -- this specifies a "not a knot" boundary
C                condition, see "cubsplb.for".  This is a common way
C                for inferring a "good" spline boundary condition
C                automatically from data in the vicinity of the
C                boundary.  ... bcxmin (bcxmax) are ignored.
C
C  ibcxmin (ibcxmax) = 1 -- boundary condition is to have s'(x(1))
C                ( s'(x(nx)) ) match the passed value bcxmin (bcxmax).
C
C  ibcxmin (ibcxmax) = 2 -- boundary condition is to have s''(x(1))
C                ( s''(x(nx)) ) match the passed value bcxmin (bcxmax).
C
C  ibcxmin (ibcxmax) = 3 -- boundary condition is to have s'(x(1))=0.0
C                ( s'(x(nx))=0.0 )
C
C  ibcxmin (ibcxmax) = 4 -- boundary condition is to have s''(x(1))=0.0
C                ( s''(x(nx))=0.0 )
C
C  ibcxmin (ibcxmax) = 5 -- boundary condition is to have s'(x(1))
C                ( s'(x(nx)) ) match the 1st divided difference
C                e.g. at x(1):  d(1)/h(1), where
C                           d(j)=f(1,j+1)-f(1,j)
C                           h(j)=x(j+1)-x(j)
C
C  ibcxmin (ibcxmax) = 6 -- BC is to have s''(x(1)) ( s''(x(nx)) )
C                match the 2nd divided difference
C                e.g. at x(1):
C                     e(1) = [d(2)/h(2) - d(1)/h(1)]/(0.5*(h(1)+h(2)))
C
C  ibcxmin (ibcxmax) = 7 -- BC is to have s'''(x(1)) ( s'''(x(nx)) )
C                match the 3rd divided difference
C                e.g. at x(1): [e(2)-e(1)]/(0.33333*(h(1)+h(2)+h(3)))
C
C  output:
C
      integer ilinx                     ! =1: hint, x axis is ~evenly spaced
C
C  let dx[avg] = (x(nx)-x(1))/(nx-1)
C  let dx[j] = x(j+1)-x(j), for all j satisfying 1.le.j.lt.nx
C
C  if for all such j, abs(dx[j]-dx[avg]).le.(1.0e-3*dx[avg]) then
C  ilinx=1 is returned, indicating the data is (at least nearly)
C  evenly spaced.  Even spacing is useful, for speed of zone lookup,
C  when evaluating a spline.
C
C  if the even spacing condition is not satisfied, ilinx=2 is returned.
C
      integer ier                       ! exit code, 0=OK
C
C  an error code is returned if the x axis is not strict ascending,
C  or if nx.lt.4, or if an invalid boundary condition specification was
C  input.
C
C------------------------------------
C
C  this routine calls traditional 4 coefficient spline software, and
C  translates the result to 2 coefficient form.
C
C  this could be done more efficiently but we decided out of conservatism
C  to use the traditional software.
C
C------------------------------------
C  workspaces -- f90 dynamically allocated
C
      real, dimension(:,:), allocatable :: fspl4 ! traditional 4-spline
      real, dimension(:), allocatable :: wk ! cspline workspace
C
C------------------------------------
C
      allocate(fspl4(4,nx),wk(nx))
C
C  make the traditional call
C
      do i=1,nx
         fspl4(1,i)=fspl(1,i)
         fspl(2,i)=0.0                  ! for now
      enddo
C
      inwk=nx
C
C  boundary conditions imposed by cspline...
C
      call cspline(x,nx,fspl4,ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   wk,inwk,ilinx,ier)
C
      if(ier.eq.0) then
C
C  copy the output -- careful of end point.
C
         do i=1,nx-1
            fspl(2,i)=2.0*fspl4(3,i)
         enddo
         fspl(2,nx)=2.0*fspl4(3,nx-1) +
     >        (x(nx)-x(nx-1))*6.0*fspl4(4,nx-1)
      endif
C
      deallocate(fspl4,wk)
C
      return
      end
