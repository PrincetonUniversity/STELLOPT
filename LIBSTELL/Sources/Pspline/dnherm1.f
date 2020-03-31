      subroutine dnherm1(x,nx,fherm,ilinx,ier)
C
C  create a data set for Hermite interpolation, based on simple
C  numerical differentiation using the given grid.
C
C  1d routine
C
C  input:
C
      integer nx                        ! array dimensions
      real x(nx)                        ! x coordinate array
      real fherm(0:1,nx)                ! data/Hermite array
C
C  fherm(0,i) = function value f at x(i)       **on input**
C
C  fherm(1,i) = derivative df/dx at x(i)       **on output**
C
C addl output:
C  ilinx=1 if x is "evenly spaced" ier=0 if no errors
C
C  ** x must be strict ascending **
C
C----------------------------
C
      ztol=1.0e-3
      call splinck(x,nx,ilinx,ztol,ier)
      if(ier.ne.0) then
         write(6,*) '?dnherm1:  x axis not strict ascending.'
         return
      endif
C
      do ix=1,nx
c
c  x div. diffs in vicinity
c
         ixp=min(nx,ix+1)
         ixm=max(1,ix-1)
         zd=(fherm(0,ixp)-fherm(0,ixm))/(x(ixp)-x(ixm))
c
         fherm(1,ix)=zd
      enddo
C
      return
      end
