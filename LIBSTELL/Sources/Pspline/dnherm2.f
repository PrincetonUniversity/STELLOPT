      subroutine dnherm2(x,nx,y,ny,fherm,nf2,ilinx,iliny,ier)
C
C  create a data set for Hermite interpolation, based on simple
C  numerical differentiation using the given grid.  The grid should
C  be "approximately" evenly spaced.
C
C  2d routine
C
C  input:
C
C  nf2.gt.nx expected!
C
      integer nx,ny,nf2                 ! array dimensions
      real x(nx)                        ! x coordinate array
      real y(ny)                        ! y coordinate array
      real fherm(0:3,nf2,ny)            ! data/Hermite array
C
C  fherm(0,i,j) = function value f at x(i),y(j)  **on input**
C
C  fherm(1,i,j) = derivative df/dx at x(i),y(j)  **on output**
C  fherm(2,i,j) = derivative df/dy at x(i),y(j)  **on output**
C  fherm(3,i,j) = derivative d2f/dxdy at x(i),y(j)  **on output**
C
C  addl output:
C    ilinx=1 if x axis is evenly spaced
C    iliny=1 if y axis is evenly spaced
C    ier=0 if no error:
C      x, y must both be strict ascending
C      nf2.ge.nx is required.
C
C----------------------------
c
      ier=0
c
      call splinck(x,nx,ilinx,1.0e-3,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?dnherm2:  x axis not strict ascending'')')
      endif
c
      call splinck(y,ny,iliny,1.0e-3,iery)
      if(iery.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?dnherm2:  y axis not strict ascending'')')
      endif
c
      if(nf2.lt.nx) then
         ier=4
         write(6,*) '?dnherm2:  fherm array dimension too small.'
      endif
C
      if(ier.ne.0) return
C
      do iy=1,ny
c
         iyp=min(ny,iy+1)
         iym=max(1,iy-1)
c
         do ix=1,nx
c
            ixp=min(nx,ix+1)
            ixm=max(1,ix-1)
c
c  x diffs in vicinity
c
            zd=(fherm(0,ixp,iy)-fherm(0,ixm,iy))/(x(ixp)-x(ixm))
c
            fherm(1,ix,iy)=zd
c
c  y diffs in vicinity
c
            zd=(fherm(0,ix,iyp)-fherm(0,ix,iym))/(y(iyp)-y(iym))
c
            fherm(2,ix,iy)=zd
c
c  xy cross derivs (except at corners, iedge=2)
c
            fherm(3,ix,iy)=(fherm(0,ixp,iyp)-fherm(0,ixm,iyp)
     >         -fherm(0,ixp,iym)+fherm(0,ixm,iym))/
     >         ((x(ixp)-x(ixm))*(y(iyp)-y(iym)))
c
         enddo
      enddo
C
      return
      end
