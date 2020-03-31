c  tpspline -- dmc 20 Jan 1999
c
c  set up coefficients for bicubic spline with following BC's:
c  * 1st dimension:  3rd divided difference
c  * derivatives periodic in second coordinate
c  * derivatives periodic in third coordinate
c
c  for more control over boundary conditions, use tpspline or tcspline.
c
c  for evaluation of interpolant, see tcspeval.for
c
      subroutine r8tpspline(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >                    wk,nwk,ilinx,ilinth,ilinph,ier)
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer inx,inth,inph,inf4,inf5,nwk
      REAL*8 x(inx),th(inth),ph(inph)
      REAL*8 fspl(4,4,4,inf4,inf5,inph),wk(nwk)
      integer ilinx,ilinth,ilinph,ier
c
c  input:
c    x(1...inx) -- abscissae, first dimension of data
c   th(1...inth) -- abscissae, second (periodic) dimension of data
c   ph(1...inph) -- abscissae, third (periodic) dimension of data
c   fspl(1,1,1,1..inx,1..inth,1..inph) -- function values
c   inf4 -- fspl dimensioning, inf4.ge.inx required.
c   inf5 -- fspl dimensioning, inf5.ge.inth required.
c
c  output:
c   fspl(*,*,*,1..inx,1..inth,1..inph) -- bicubic spline coeffs (4x4)
c   ...fspl(1,1,1,*,*,*) is not replaced.
c
c   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
c   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
c   ilinph-- =1 on output if ph(inph) evenly spaced (tol=1e-3)
c
c   ier -- completion code, 0 for normal
c
c  workspace:
c   wk -- must be at least 5*max(inx,inth,inph) large
c  nwk -- size of workspace
c
c---------------------------------
c  compute tricubic spline of 3d function, given values at the
c  grid crossing points, f(1,1,1,i,j,k)=f(x(i),th(j),ph(k)).
c
c  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                       and th btw th(j) and th(j+1), dt=th-th(j),
c                       and ph btw ph(k) and ph(k+1), dp=ph-ph(k),
c
c   spline =
c       f(1,1,1,i,j,k)+dx*f(2,1,1,i,j,k)+dx2*f(3,1,1,i,j,k)+dx3*f(4,1,1,i,j,k)
c  +dt*(f(1,2,1,i,j,k)+dx*f(2,2,1,i,j,k)+dx2*f(3,2,1,i,j,k)+dx3*f(4,2,1,i,j,k))
c +dt2*(f(1,3,1,i,j,k)+dx*f(2,3,1,i,j,k)+dx2*f(3,3,1,i,j,k)+dx3*f(4,3,1,i,j,k))
c +dt3*(f(1,4,1,i,j,k)+dx*f(2,4,1,i,j,k)+dx2*f(3,4,1,i,j,k)+dx3*f(4,4,1,i,j,k))
c        +dp*(
c       f(1,1,2,i,j,k)+dx*f(2,1,2,i,j,k)+dx2*f(3,1,2,i,j,k)+dx3*f(4,1,2,i,j,k)
c  +dt*(f(1,2,2,i,j,k)+dx*f(2,2,2,i,j,k)+dx2*f(3,2,2,i,j,k)+dx3*f(4,2,2,i,j,k))
c +dt2*(f(1,3,2,i,j,k)+dx*f(2,3,2,i,j,k)+dx2*f(3,3,2,i,j,k)+dx3*f(4,3,2,i,j,k))
c +dt3*(f(1,4,2,i,j,k)+dx*f(2,4,2,i,j,k)+dx2*f(3,4,2,i,j,k)+dx3*f(4,4,2,i,j,k)))
c        +dp2*(
c       f(1,1,3,i,j,k)+dx*f(2,1,3,i,j,k)+dx2*f(3,1,3,i,j,k)+dx3*f(4,1,3,i,j,k)
c  +dt*(f(1,2,3,i,j,k)+dx*f(2,2,3,i,j,k)+dx2*f(3,2,3,i,j,k)+dx3*f(4,2,3,i,j,k))
c +dt2*(f(1,3,3,i,j,k)+dx*f(2,3,3,i,j,k)+dx2*f(3,3,3,i,j,k)+dx3*f(4,3,3,i,j,k))
c +dt3*(f(1,4,3,i,j,k)+dx*f(2,4,3,i,j,k)+dx2*f(3,4,3,i,j,k)+dx3*f(4,4,3,i,j,k)))
c        +dp3*(
c       f(1,1,4,i,j,k)+dx*f(2,1,4,i,j,k)+dx2*f(3,1,4,i,j,k)+dx3*f(4,1,4,i,j,k)
c  +dt*(f(1,2,4,i,j,k)+dx*f(2,2,4,i,j,k)+dx2*f(3,2,4,i,j,k)+dx3*f(4,2,4,i,j,k))
c +dt2*(f(1,3,4,i,j,k)+dx*f(2,3,4,i,j,k)+dx2*f(3,3,4,i,j,k)+dx3*f(4,3,4,i,j,k))
c +dt3*(f(1,4,4,i,j,k)+dx*f(2,4,4,i,j,k)+dx2*f(3,4,4,i,j,k)+dx3*f(4,4,4,i,j,k)))
c
c      where dx2=dx**2 and dx3=dx**3.
c      where dt2=dt**2 and dt3=dt**3.
c      where dp2=dp**2 and dp3=dp**3.
c
c---------------------------------
      integer ierx,ierth,ierph
      REAL*8 zdumth(inth),zdumx(inx)
c---------------------------------
c
      ier=0
      if(nwk.lt.5*max(inx,inth,inph)) then
         write(6,'('' ?tpspline:  workspace too small.'')')
         ier=1
      endif
      if(inx.lt.2) then
         write(6,'('' ?tpspline:  at least 2 x points required.'')')
         ier=1
      endif
      if(inth.lt.2) then
         write(6,'('' ?tpspline:  need at least 2 theta points.'')')
         ier=1
      endif
      if(inph.lt.2) then
         write(6,'('' ?tpspline:  need at least 2 phi points.'')')
         ier=1
      endif
c
c  check ilinx & x vector
c
      call r8splinck(x,inx,ilinx,1.0E-3_r8,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?tpspline:  x axis not strict ascending'')')
      endif
c
c  check ilinth & th vector
c
      call r8splinck(th,inth,ilinth,1.0E-3_r8,ierth)
      if(ierth.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?tpspline:  theta axis not strict ascending'')')
      endif
c
c  check ilinth & th vector
c
      call r8splinck(ph,inph,ilinph,1.0E-3_r8,ierph)
      if(ierph.ne.0) ier=4
c
      if(ier.eq.4) then
         write(6,'('' ?tpspline:  phi axis not strict ascending'')')
      endif
c
      if(ier.ne.0) return
c
c------------------------------------
c
      call r8tcspline(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >   7,zdumth,7,zdumth,inth,
     >   -1,zdumx,-1,zdumx,inx,
     >   -1,zdumx,-1,zdumx,inx,
     >   wk,nwk,ilinx,ilinth,ilinph,ier)
c
      return
      end
