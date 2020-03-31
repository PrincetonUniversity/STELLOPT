      subroutine r8mktricub(x,nx,y,ny,z,nz,f,nf2,nf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x,
     >                    ibcymin,bcymin,ibcymax,bcymax,inb1y,
     >                    ibczmin,bczmin,ibczmax,bczmax,inb1z,
     >                    ilinx,iliny,ilinz,ier)
c
c  setup a tricubic spline; store coefficients in compact form
c  (as per suggestion of L. Zakharov, PPPL, Feb. 1999)
C  8 coeffs per (x,y,z) grid point:
C          f,fxx,fyy,fzz,fxxyy,fxxzz,fyyzz,fxxyyzz
C
C  dmc -- modified Feb 2004 -- rewritten to compute coefficients
C  directly rather than by conversion from the non-compact representation
C  (to reduce cpu and memory cost)
C
C
C  input:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer nx                        ! length of x vector
      integer ny                        ! length of y vector
      integer nz                        ! length of z vector
      REAL*8 x(nx)                        ! x vector, strict ascending
      REAL*8 y(ny)                        ! y vector, strict ascending
      REAL*8 z(nz)                        ! z vector, strict ascending
c
      integer nf2                       ! 2nd dim. of f array, nf2.ge.nx
      integer nf3                       ! 3rd dim. of f array, nf3.ge.ny
c
c  input/output:
c
      REAL*8 f(8,nf2,nf3,nz)              ! data and spline coefficients
c
C  on input:  f(1,i,j,k) = f(x(i),y(j),z(k))
C  on output:  f(1,i,j,k) unchanged
C              f(2,i,j,k) = d2f/dx2(x(i),y(j),z(k))
C              f(3,i,j,k) = d2f/dy2(x(i),y(j),z(k))
C              f(4,i,j,k) = d2f/dz2(x(i),y(j),z(k))
C              f(5,i,j,k) = d4f/dx2dy2(x(i),y(j),z(k))
C              f(6,i,j,k) = d4f/dx2dz2(x(i),y(j),z(k))
C              f(7,i,j,k) = d4f/dy2dz2(x(i),y(j),z(k))
C              f(8,i,j,k) = d6f/dx2dy2dz2(x(i),y(j),z(k))
C
C  there is a rather Hermite like interpolation formula to go with
C  this-- see evtricub.for.  Also the bicubic formula is given in
C  mkbicubw.for; the tricubic formula is precisely analogous.
C
C  boundary condition data
C  inputs:
      integer inb1x                     ! 1st dim of xmin & xmax bc arrays
      integer inb1y                     ! 1st dim of ymin & ymax bc arrays
      integer inb1z                     ! 1st dim of zmin & zmax bc arrays
C
      integer ibcxmin,ibcxmax           ! BC type flag @xmin, xmax
      integer ibcymin,ibcymax           ! BC type flag @ymin, ymax
      integer ibczmin,ibczmax           ! BC type flag @zmin, zmax
C
      REAL*8 bcxmin(inb1x,nz),bcxmax(inb1x,nz) ! xmin & xmax BC data, ny x nz
      REAL*8 bcymin(inb1y,nz),bcymax(inb1y,nz) ! ymin & ymax BC data, nx x nz
      REAL*8 bczmin(inb1z,ny),bczmax(inb1z,ny) ! zmin & zmax BC data, nx x ny
c
c  where BC data is not required, dummy scalars may be passed.
C  the ibc* flags determine whether BC data isneeded.
c
c  BC data:  bcxmin & bcxmax:  BC vs. y,z @xmin,xmax
C            bcymin & bcymax:  BC vs. x,z @ymin,ymax
C            bczmin & bczmax:  BC vs. x,y @zmin,zmax
C
c   ibcxmin -- indicator for boundary condition at xmin=x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- use periodic boundary condition
c     =0 -- use "not a knot"
c     =1 -- match slope, specified at x(1),y(iy),z(iz) by bcxmin(iy,iz)
c     =2 -- match 2nd derivative, specified at x(1),y(iy),z(iz)
c           by bcxmin(iy,iz
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all y(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all y(j)
c     =5 -- df/dx BC from 1st divided difference
c     =6 -- d2f/dx2 BC from 2nd divided difference (parabolic fit)
c     =7 -- d3f/dx3 BC from 3rd divided difference (cubic fit)
c   ***NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
c
c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
c            and ibcxmax, bcxmax are ignored.
c   inb1x -- 1st dimension of bcxmin, bcxmax: if ibcxmin or ibcxmax .gt. 0
c            this must be .ge. ny.
c
c   interpretation of ibcymin,bcymin,ibcymax,bcymax,inb1y
c     is same as with ibcxmin,...
c
c   interpretation of ibczmin,bczmin,ibczmax,bczmax,inb1z
c     is same as with ibcxmin,...
c
c   the explicit bdy condition arrays are referenced only if the
c     corresponding "ibc" flag values are set to 1 or 2.
c
c  output:
      integer ilinx                     ! x vector equal spacing flag
      integer iliny                     ! y vector equal spacing flag
      integer ilinz                     ! z vector equal spacing flag
c
c   ilinx -- =1 on output if x(nx) pts are nearly evenly spaced (tol=1e-3)
c   iliny -- =1 on output if y(ny) evenly spaced (tol=1e-3)
c   ilinz -- =1 on output if z(nz) evenly spaced (tol=1e-3)
c
      integer ier                       ! exit code
c   ier -- completion code, 0 for normal
c
C-----------------------------------------------------
c  workspace **dynamic allocation**
C  f90 dynamic array
C
      REAL*8, dimension(:,:,:), allocatable :: fbicub ! bicubic subsection
      REAL*8, dimension(:,:), allocatable :: fwk ! work array
      REAL*8, dimension(:), allocatable :: bcx1,bcx2,bcy1,bcy2 ! BCs for mkbicub
c
      REAL*8, dimension(:,:,:,:), allocatable :: fcorr ! correction spline
      REAL*8, dimension(:,:), allocatable :: bcc1,bcc2 ! correction BCs
c
      integer iflg,ierx,iery,ierz
      integer ix,iy,iz
c
      REAL*8 ztol
      REAL*8 zbc1,zbc2,hz
      integer ibc1,ibc2
c
      data ztol/1.0E-3_r8/
c-----------------------------------------------------
c
      ier=0
c
      iflg=0
c
c  check z bdy condition "linearity"
c
      if(ibczmin.ne.-1) then
         if((ibczmin.eq.1).or.(ibczmin.eq.2)) then
            do iy=1,ny
               do ix=1,nx
                  if(bczmin(ix,iy).ne.0.0_r8) iflg=1
               enddo
            enddo
         endif
         if((ibczmax.eq.1).or.(ibczmax.eq.2)) then
            do iy=1,ny
               do ix=1,nx
                  if(bczmax(ix,iy).ne.0.0_r8) iflg=1
               enddo
            enddo
         endif
      endif
c
      if(nx.lt.2) then
         write(6,'('' ?mktricub:  at least 2 x points required.'')')
         ier=1
      endif
      if(ny.lt.2) then
         write(6,'('' ?mktricub:  need at least 2 y points.'')')
         ier=1
      endif
      if(nz.lt.2) then
         write(6,'('' ?mktricub:  need at least 2 z points.'')')
         ier=1
      endif
c
      if((ibcxmin.eq.1).or.(ibcxmax.eq.1).or.(ibcxmin.eq.2).or.
     >   (ibcxmax.eq.2)) then
         if(inb1x.lt.ny) then
            ier=1
            write(6,
     >'('' ?mktricub:  1st dim of bcxmin/max arrays .lt. ny'')')
         endif
      endif
c
      if((ibcymin.eq.1).or.(ibcymax.eq.1).or.(ibcymin.eq.2).or.
     >   (ibcymax.eq.2)) then
         if(inb1y.lt.nx) then
            ier=1
            write(6,
     >'('' ?mktricub:  1st dim of bcymin/max arrays .lt. nx'')')
         endif
      endif
c
      if((ibczmin.eq.1).or.(ibczmax.eq.1).or.(ibczmin.eq.2).or.
     >   (ibczmax.eq.2)) then
         if(inb1z.lt.nx) then
            ier=1
            write(6,
     >'('' ?mktricub:  1st dim of bczmin/max arrays .lt. nx'')')
         endif
      endif
c
      call ibc_ck(ibcxmin,'mktricub','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'mktricub','xmax',0,7,ier)
c
      call ibc_ck(ibcymin,'mktricub','ymin',-1,7,ier)
      if(ibcymin.ge.0) call ibc_ck(ibcymax,'mktricub','ymax',0,7,ier)
c
      call ibc_ck(ibczmin,'mktricub','zmin',-1,7,ier)
      if(ibczmax.ge.0) call ibc_ck(ibczmax,'mktricub','zmax',0,7,ier)
c
c  check ilinx & x vector
c
      call r8splinck(x,nx,ilinx,ztol,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?mktricub:  x axis not strict ascending'')')
      endif
c
c  check iliny & y vector
c
      call r8splinck(y,ny,iliny,ztol,iery)
      if(iery.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?mktricub:  y axis not strict ascending'')')
      endif
c
c  check ilinz & z vector
c
      call r8splinck(z,nz,ilinz,ztol,ierz)
      if(ierz.ne.0) ier=4
c
      if(ier.eq.4) then
         write(6,'('' ?mktricub:  z axis not strict ascending'')')
      endif
c
      if(ier.ne.0) return
c
c------------------------------------
c  1.  compute (x,y) bicubic splines using mkbicub
c
      allocate(fbicub(4,nx,ny))
      allocate(bcx1(ny),bcx2(ny),bcy1(nx),bcy2(nx))
      bcx1=0.0; bcx2=0.0; bcy1=0.0; bcy2=0.0_r8
c
      do iz=1,nz
         if(ibcxmin.ne.-1) then
            if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) then
               bcx1(1:ny)=bcxmin(1:ny,iz)
            endif
            if((ibcxmax.eq.1).or.(ibcxmax.eq.2)) then
               bcx2(1:ny)=bcxmax(1:ny,iz)
            endif
         endif
         if(ibcymin.ne.-1) then
            if((ibcymin.eq.1).or.(ibcymin.eq.2)) then
               bcy1(1:nx)=bcymin(1:nx,iz)
            endif
            if((ibcymax.eq.1).or.(ibcymax.eq.2)) then
               bcy2(1:nx)=bcymax(1:nx,iz)
            endif
         endif
c
         fbicub(1,1:nx,1:ny) = f(1,1:nx,1:ny,iz)
c
         call r8mkbicub(x,nx,y,ny,fbicub,nx,
     >      ibcxmin,bcx1,ibcxmax,bcx2,
     >      ibcymin,bcy1,ibcymax,bcy2,
     >      ilinx,iliny,ier)
         if(ier.ne.0) return
c
         f(2:3,1:nx,1:ny,iz) = fbicub(2:3,1:nx,1:ny)  ! fxx, fyy
         f(5,1:nx,1:ny,iz) = fbicub(4,1:nx,1:ny)      ! fxxyy
c
      enddo
c
      deallocate(fbicub,bcx1,bcx2,bcy1,bcy2)
c
c  2.  homogeneous spline in z direction; inhomogeneous BC imposed later
c      if necessary
c
      zbc1=0.0_r8
      zbc2=0.0_r8
      ibc1=ibczmin
      ibc2=ibczmax
      if(iflg.eq.1) then
         if((ibczmin.eq.1).or.(ibczmin.eq.2)) ibc1=0
         if((ibczmax.eq.1).or.(ibczmax.eq.2)) ibc2=0
      endif
c
      allocate(fwk(2,nz))
c
      do iy=1,ny
         do ix=1,nx
 
            fwk(1,1:nz) = f(1,ix,iy,1:nz)
            call r8mkspline(z,nz,fwk,
     >         ibc1,zbc1,ibc2,zbc2,ilinz,ier)
            if(ier.ne.0) return
            f(4,ix,iy,1:nz) = fwk(2,1:nz) ! fzz
 
            fwk(1,1:nz) = f(2,ix,iy,1:nz)
            call r8mkspline(z,nz,fwk,
     >         ibc1,zbc1,ibc2,zbc2,ilinz,ier)
            if(ier.ne.0) return
            f(6,ix,iy,1:nz) = fwk(2,1:nz) ! fxxzz
 
            fwk(1,1:nz) = f(3,ix,iy,1:nz)
            call r8mkspline(z,nz,fwk,
     >         ibc1,zbc1,ibc2,zbc2,ilinz,ier)
            if(ier.ne.0) return
            f(7,ix,iy,1:nz) = fwk(2,1:nz) ! fyyzz
 
            fwk(1,1:nz) = f(5,ix,iy,1:nz)
            call r8mkspline(z,nz,fwk,
     >         ibc1,zbc1,ibc2,zbc2,ilinz,ier)
            if(ier.ne.0) return
            f(8,ix,iy,1:nz) = fwk(2,1:nz) ! fxxyyzz
 
         enddo
      enddo
c
      deallocate(fwk)
c
      if(iflg.eq.1) then
c
c  3. inhomogeneous BC correction
c
         allocate(fwk(2,max(nx,ny,nz)))
         allocate(bcc1(nx,ny),bcc2(nx,ny))
         allocate(fcorr(4,nx,ny,nz))
c
c  correction BCs
c
         do iy=1,ny
            do ix=1,nx
               bcc1(ix,iy)=0.0_r8
               if(ibczmin.eq.1) then
                  hz=z(2)-z(1)
                  bcc1(ix,iy)=(f(1,ix,iy,2)-f(1,ix,iy,1))/hz +
     >               hz*(-2*f(4,ix,iy,1)-f(4,ix,iy,2))/6
                  bcc1(ix,iy)=bczmin(ix,iy)-bcc1(ix,iy)
               else if(ibczmin.eq.2) then
                  bcc1(ix,iy)=bczmin(ix,iy)-f(4,ix,iy,1)
               endif
            enddo
         enddo
c
         do iy=1,ny
            do ix=1,nx
               bcc2(ix,iy)=0.0_r8
               if(ibczmax.eq.1) then
                  hz=z(2)-z(1)
                  bcc2(ix,iy)=(f(1,ix,iy,2)-f(1,ix,iy,1))/hz +
     >               hz*(-2*f(4,ix,iy,1)-f(4,ix,iy,2))/6
                  bcc2(ix,iy)=bczmax(ix,iy)-bcc2(ix,iy)
               else if(ibczmax.eq.2) then
                  bcc2(ix,iy)=bczmax(ix,iy)-f(4,ix,iy,1)
               endif
            enddo
         enddo
c
         fwk(1,1:nz)=0.0_r8  ! values are all zero, only BC is set...
         do iy=1,ny
            do ix=1,nx
               call r8mkspline(z,nz,fwk,
     >            ibczmin,bcc1(ix,iy),ibczmax,bcc2(ix,iy),ilinz,ier)
               if(ier.ne.0) return
               fcorr(1,ix,iy,1:nz)=fwk(2,1:nz)  ! fzz-correction
            enddo
         enddo
c
c  higher order corrections
c
         zbc1=0.0_r8
         zbc2=0.0_r8
c
         do iz=1,nz
            do iy=1,ny
               fwk(1,1:nx)=fcorr(1,1:nx,iy,iz)
               call r8mkspline(x,nx,fwk,
     >            ibcxmin,zbc1,ibcxmax,zbc2,ilinx,ier)
               if(ier.ne.0) return
               fcorr(2,1:nx,iy,iz)=fwk(2,1:nx)  ! fxxzz-correction
            enddo
         enddo
c
         do iz=1,nz
            do ix=1,nx
               fwk(1,1:ny)=fcorr(1,ix,1:ny,iz)
               call r8mkspline(y,ny,fwk,
     >            ibcymin,zbc1,ibcymax,zbc2,iliny,ier)
               if(ier.ne.0) return
               fcorr(3,ix,1:ny,iz)=fwk(2,1:ny)  ! fyyzz-correction
 
               fwk(1,1:ny)=fcorr(2,ix,1:ny,iz)
               call r8mkspline(y,ny,fwk,
     >            ibcymin,zbc1,ibcymax,zbc2,iliny,ier)
               if(ier.ne.0) return
               fcorr(4,ix,1:ny,iz)=fwk(2,1:ny)  ! fxxyyzz-correction
            enddo
         enddo
c
c  apply correction
c
         do iz=1,nz
            do iy=1,ny
               do ix=1,nx
                  f(4,ix,iy,iz)=f(4,ix,iy,iz)+fcorr(1,ix,iy,iz)
                  f(6,ix,iy,iz)=f(6,ix,iy,iz)+fcorr(2,ix,iy,iz)
                  f(7,ix,iy,iz)=f(7,ix,iy,iz)+fcorr(3,ix,iy,iz)
                  f(8,ix,iy,iz)=f(8,ix,iy,iz)+fcorr(4,ix,iy,iz)
               enddo
            enddo
         enddo
c
         deallocate(fwk,fcorr,bcc1,bcc2)
c
      endif
c
c  that's all
c
      return
      end
 
