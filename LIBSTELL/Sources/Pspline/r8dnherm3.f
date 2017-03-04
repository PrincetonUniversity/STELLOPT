      subroutine r8dnherm3(x,nx,y,ny,z,nz,fherm,nf2,nf3,
     >   ilinx,iliny,ilinz,ier)
C
C  create a data set for Hermite interpolation, based on simple
C  numerical differentiation using the given grid.  The grid should
C  be "approximately" evenly spaced.
C
C  3d routine
C
C  input:
C
C  nf2.get.nx and nf3.ge.ny  expected!
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iliny,ilinz,ier,ilinx,ierx,iery,ierz,iz,izp,izm,iy
      INTEGER iyp,iym,ix,ixp,ixm
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zd
!============
      integer nx,ny,nz,nf2,nf3          ! array dimensions
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 y(ny)                        ! y coordinate array
      REAL*8 z(nz)                        ! z coordinate array
      REAL*8 fherm(0:7,nf2,nf3,nz)        ! data/Hermite array
C
C  fherm(0,i,j,k) = function value f at x(i),y(j),z(k)  **on input**
C
C  fherm(1,i,j,k) = derivative df/dx at x(i),y(j),z(k)  **on output**
C  fherm(2,i,j,k) = derivative df/dy at x(i),y(j),z(k)  **on output**
C  fherm(3,i,j,k) = derivative df/dz at x(i),y(j),z(k)  **on output**
C  fherm(4,i,j,k) = derivative d2f/dxdy at x(i),y(j),z(k)  **on output**
C  fherm(5,i,j,k) = derivative d2f/dxdz at x(i),y(j),z(k)  **on output**
C  fherm(6,i,j,k) = derivative d2f/dydz at x(i),y(j),z(k)  **on output**
C  fherm(7,i,j,k) = derivative d3f/dxdydz at x(i),y(j),z(k)  **on output**
C
C  additional outputs:
C
C    ilinx = 1  if x is evenly spaced (approximately)
C    iliny = 1  if y is evenly spaced (approximately)
C    ilinz = 1  if z is evenly spaced (approximately)
C
C    ier = 0 if there are no errors
C
C  note possible errors:  x y or z NOT strict ascending
C  nf2.lt.nx   .or.   nf3.lt.ny
C
C----------------------------
C
C
C  error check
C
      ier=0
c
      call r8splinck(x,nx,ilinx,1.0E-3_r8,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?dnherm3:  x axis not strict ascending'')')
      endif
c
      call r8splinck(y,ny,iliny,1.0E-3_r8,iery)
      if(iery.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?dnherm3:  y axis not strict ascending'')')
      endif
c
      call r8splinck(z,nz,ilinz,1.0E-3_r8,ierz)
      if(ierz.ne.0) ier=4
c
      if(ier.eq.4) then
         write(6,'('' ?dnherm3:  z axis not strict ascending'')')
      endif
c
      if(nf2.lt.nx) then
         ier=5
         write(6,*) '?dnherm3:  fherm (x) array dimension too small.'
      endif
c
      if(nf3.lt.ny) then
         ier=6
         write(6,*) '?dnherm3:  fherm (y) array dimension too small.'
      endif
C
      if(ier.ne.0) return
C
C  error check OK
C
      do iz=1,nz
c
         izp=min(nz,iz+1)
         izm=max(1,iz-1)
c
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
               zd=(fherm(0,ixp,iy,iz)-fherm(0,ixm,iy,iz))/
     >            (x(ixp)-x(ixm))
c
               fherm(1,ix,iy,iz)=zd
c
c  y diffs in vicinity
c
               zd=(fherm(0,ix,iyp,iz)-fherm(0,ix,iym,iz))/
     >            (y(iyp)-y(iym))
c
               fherm(2,ix,iy,iz)=zd
c
c  z diffs in vicinity
c
               zd=(fherm(0,ix,iy,izp)-fherm(0,ix,iy,izm))/
     >            (z(izp)-z(izm))
c
               fherm(3,ix,iy,iz)=zd
c
c  xy cross derivs
c
               fherm(4,ix,iy,iz)=
     >            (fherm(0,ixp,iyp,iz)-fherm(0,ixm,iyp,iz)
     >            -fherm(0,ixp,iym,iz)+fherm(0,ixm,iym,iz))/
     >            ((x(ixp)-x(ixm))*(y(iyp)-y(iym)))
c
c  xz cross derivs
c
               fherm(5,ix,iy,iz)=
     >            (fherm(0,ixp,iy,izp)-fherm(0,ixm,iy,izp)
     >            -fherm(0,ixp,iy,izm)+fherm(0,ixm,iy,izm))/
     >            ((x(ixp)-x(ixm))*(z(izp)-z(izm)))
c
c  yz cross derivs
c
               fherm(6,ix,iy,iz)=
     >            (fherm(0,ix,iyp,izp)-fherm(0,ix,iym,izp)
     >            -fherm(0,ix,iyp,izm)+fherm(0,ix,iym,izm))/
     >            ((y(iyp)-y(iym))*(z(izp)-z(izm)))
c
c  xyz cross deriv
c
               fherm(7,ix,iy,iz)=
     >            ((fherm(0,ixp,iyp,izp)-fherm(0,ixp,iym,izp)
     >            -fherm(0,ixp,iyp,izm)+fherm(0,ixp,iym,izm))-
     >            (fherm(0,ixm,iyp,izp)-fherm(0,ixm,iym,izp)
     >            -fherm(0,ixm,iyp,izm)+fherm(0,ixm,iym,izm)))/
     >            ((x(ixp)-x(ixm))*(y(iyp)-y(iym))*(z(izp)-z(izm)))
c
            enddo
         enddo
      enddo
C
      return
      end
