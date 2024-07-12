c  bpsplinb -- dmc 30 May 1996
c
c  set up coefficients for bicubic spline with following BC's:
c  * LHS and RHS BC under user control (see comments)
c  * derivatives periodic in second coordinate (use pspline.for)
c
      subroutine bpsplinb(x,inx,th,inth,fspl,inf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,
     >                    wk,nwk,ilinx,ilinth,ier)
c
      implicit NONE
      integer inx,inth,inf3,nwk
      real x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
      integer ibcxmin,ibcxmax
      real bcxmin(inth),bcxmax(inth)
c
c  input:
c    x(1...inx) -- abscissae, first dimension of data
c   th(1...inth) -- abscissae, second (periodic) dimension of data
c   fspl(1,1,1..inx,1..inth) -- function values
c   inf3 -- fspl dimensioning, inf3.ge.inx required.
c
c  boundary conditions input:
c   ibcxmin -- indicator for boundary condition at x(1):
c    bcxmin(...) -- boundary condition data
c     =0 -- use "not a knot", bcxmin(...) ignored
c     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
c     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
c     =5 -- match 1st derivative to 1st divided difference
c     =6 -- match 2nd derivative to 2nd divided difference
c     =7 -- match 3rd derivative to 3rd divided difference
c           (for more detailed definition of BCs 5-7, see the
c           comments of subroutine mkspline)
c   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
c
c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpolation as with ibcxmin, bcxmin)
c
c  output:
c   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
c   ...fspl(1,1,*,*) is not replaced.
c
      integer ilinx,ilinth,ier
c
c   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
c   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
c
c   ier -- completion code, 0 for normal
c
c  workspace:
c   wk -- must be at least 5*max(inx,inth) large
c  nwk -- size of workspace
c
c---------------------------------
c  compute bicubic spline of 2d function, given values at the grid
c  grid crossing points, f(1,1,i,j)=f(x(i),th(j)).
c
c  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                       and th btw th(j) and th(j+1), dt=th-th(j),
c
c      spline =
c        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
c   +dt*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
c   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
c   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))
c
c      where d2=dt**2 and d3=dt**3.
c
c---------------------------------
c
      integer ierx,ierth
      real zdum(1)
c
c---------------------------------
c
      ier=0
      if(nwk.lt.5*max(inx,inth)) then
         write(6,'('' ?bpsplinb:  workspace too small.'')')
         ier=1
      endif
      if(inx.lt.2) then
         write(6,'('' ?bpsplinb:  at least 2 x points required.'')')
         ier=1
      endif
      if(inth.lt.2) then
         write(6,'('' ?bpsplinb:  need at least 2 theta points.'')')
         ier=1
      endif
c
      call ibc_ck(ibcxmin,'bcspline','xmin',0,7,ier)
      call ibc_ck(ibcxmax,'bcspline','xmax',0,7,ier)
c
c  check ilinx & x vector
c
      call splinck(x,inx,ilinx,1.0e-3,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?bpsplinb:  x axis not strict ascending'')')
      endif
c
c  check ilinth & th vector
c
      call splinck(th,inth,ilinth,1.0e-3,ierth)
      if(ierth.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?bpsplinb:  th axis not strict ascending'')')
      endif
c
      if(ier.ne.0) return
c
c------------------------------------
      zdum=0.0
c
      call bcspline(x,inx,th,inth,fspl,inf3,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   -1,zdum,-1,zdum,
     >   wk,nwk,ilinx,ilinth,ier)
c
      end
