c  bpspline -- dmc 30 May 1996
c
c  set up coefficients for bicubic spline with following BC's:
c  * LHS and RHS BCs -- d3f/dx3 from 3rd divided difference, 1st coordinate
c  * derivatives periodic in second coordinate
c
c  see similar routine, bpsplinb, for more control over x boundary
c  conditions...
c
      subroutine bpspline(x,inx,th,inth,fspl,inf3,
     >                    wk,nwk,ilinx,ilinth,ier)
c
      implicit NONE
      integer inx,inth,inf3,nwk
      real x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
      integer ilinx,ilinth
      integer ier
c
c  input:
c    x(1...inx) -- abscissae, first dimension of data
c   th(1...inth) -- abscissae, second (periodic) dimension of data
c   fspl(1,1,1..inx,1..inth) -- function values
c   inf3 -- fspl dimensioning, inf3.ge.inx required.
c
c  output:
c   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
c   ...fspl(1,1,*,*) is not replaced.
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
      integer ierx,ierth
      real zdum(1)
c---------------------------------
c
      ier=0
      if(nwk.lt.5*max(inx,inth)) then
         write(6,'('' ?bpspline:  workspace too small.'')')
         ier=1
      endif
      if(inx.lt.2) then
         write(6,'('' ?bpspline:  at least 2 x points required.'')')
         ier=1
      endif
      if(inth.lt.2) then
         write(6,'('' ?bpspline:  need at least 2 theta points.'')')
         ier=1
      endif
c
c  check ilinx & x vector
c
      call splinck(x,inx,ilinx,1.0e-3,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?bpspline:  x axis not strict ascending'')')
      endif
c
c  check ilinth & th vector
c
      call splinck(th,inth,ilinth,1.0e-3,ierth)
      if(ierth.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?bpspline:  th axis not strict ascending'')')
      endif
c
      if(ier.ne.0) return
c
c------------------------------------
c
      zdum=0.0
c
c  not-a-knot BCs for x, periodic for theta
c
      call bcspline(x,inx,th,inth,fspl,inf3,
     >   7,zdum,7,zdum,
     >   -1,zdum,-1,zdum,
     >   wk,nwk,ilinx,ilinth,ier)
c
      return
      end
