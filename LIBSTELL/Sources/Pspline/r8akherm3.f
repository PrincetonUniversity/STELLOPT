      subroutine r8akherm3(x,nx,y,ny,z,nz,fherm,nf2,nf3,
     >   ilinx,iliny,ilinz,ier)
C
C  create a data set for Hermite interpolation, based on Akima's method
C  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
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
      INTEGER iliny,ilinz,ier,ilinx
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
C  a default boundary condition is used, based on divided differences
C  in the edge region.  For more control of BC, use akherm3p...
C
      call r8akherm3p(x,nx,y,ny,z,nz,fherm,nf2,nf3,
     >   ilinx,iliny,ilinz,0,0,0,ier)
C
      return
      end
C----------------------------
      subroutine r8akherm3p(x,nx,y,ny,z,nz,fherm,nf2,nf3,
     >   ilinx,iliny,ilinz,ipx,ipy,ipz,ier)
C
C  create a data set for Hermite interpolation, based on Akima's method
C  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
C
C  3d routine --
C  with independently settable boundary condition options:
C     ipx|ipy|ipz  =0  -- default, boundary conditions from divided diffs
C     ipx|ipy|ipz  =1  -- periodic boundary condition
C     ipx|ipy|ipz  =2  -- user supplied df/dx|df/dy|df/dz
C
C  input:
C
C  nf2.get.nx and nf3.ge.ny  expected!
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iliny,ilinz,ier,ilinx,ierx,iery,ierz,ierbc,iz,iy,ix
      INTEGER icx,incx,ix1,icy,incy,iy1,icz,incz,iz1,izm2,izm1
      INTEGER izmm2,izmm1,izp2,izp1,izpp2,izpp1,iym2,iym1,iymm2
      INTEGER iymm1,iyp2,iyp1,iypp2,iypp1,ixm2,ixm1,ixmm2,ixmm1
      INTEGER iflagx,ixp2,ixp1,ixpp2,ixpp1,iflagy,iflagz,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cxp,cxm,cxpp,cxmm,cxtrap0,cxtrapn,cyp,cym,cypp,cymm
      REAL*8 cytrap0,cytrapn,czp,czm,czpp,czmm,cztrap0,cztrapn
      REAL*8 zansr,cxm2,cxp2,cym2,cyp2,fzp,fzm
!============
      integer nx,ny,nz,nf2,nf3          ! array dimensions
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 y(ny)                        ! y coordinate array
      REAL*8 z(nz)                        ! z coordinate array
      REAL*8 fherm(0:7,nf2,nf3,nz)        ! data/Hermite array
C
      integer ipx                       ! =0: dflt; =1: periodic df/dx
      integer ipy                       ! =0: dflt; =1: periodic df/dy
      integer ipz                       ! =0: dflt; =1: periodic df/dz
C
C  if ipx=2:  fherm(1,1,1:ny,1:nz) & fherm(1,nx,1:ny,1:nz) = INPUT df/dx BCs
C  if ipy=2:  fherm(2,1:nx,1,1:nz) & fherm(2,1:nx,ny,1:nz) = INPUT df/dy BCs
C  if ipz=2:  fherm(3,1:nx,1:ny,1) & fherm(3,1:nx,1:ny,nz) = INPUT df/dz BCs
C
C  for the rest of the grid points...
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
C---------------------
C  local...
C
      REAL*8 wx(2),wy(2),wz(2),exy(2,2),exz(2,2),eyz(2,2),fxyz(2,2,2)
C
      REAL*8, dimension(:,:,:), allocatable ::  ftmp
c
      REAL*8 xx(0:nx+1)
      REAL*8 yy(0:ny+1)
      REAL*8 zz(0:nz+1)
C
C----------------------------
C
C  error check
C
      ier=0
c
      call r8splinck(x,nx,ilinx,1.0E-3_r8,ierx)
      if(ierx.ne.0) ier=ier+1
c
      if(ierx.ne.0) then
         write(6,'('' ?akherm3:  x axis not strict ascending'')')
      endif
c
      call r8splinck(y,ny,iliny,1.0E-3_r8,iery)
      if(iery.ne.0) ier=ier+1
c
      if(iery.ne.0) then
         write(6,'('' ?akherm3:  y axis not strict ascending'')')
      endif
c
      call r8splinck(z,nz,ilinz,1.0E-3_r8,ierz)
      if(ierz.ne.0) ier=ier+1
c
      if(ierz.ne.0) then
         write(6,'('' ?akherm3:  z axis not strict ascending'')')
      endif
c
      if(nf2.lt.nx) then
         ier=ier+1
         write(6,*) '?akherm3:  fherm (x) array dimension too small.'
      endif
c
      if(nf3.lt.ny) then
         ier=ier+1
         write(6,*) '?akherm3:  fherm (y) array dimension too small.'
      endif
C
      ierbc=0
      call ibc_ck(ipx,'akherm3','X Bdy Cond',0,2,ierbc)
      ier=ier+ierbc
C
      ierbc=0
      call ibc_ck(ipy,'akherm3','Y Bdy Cond',0,2,ierbc)
      ier=ier+ierbc
C
      ierbc=0
      call ibc_ck(ipz,'akherm3','Z Bdy Cond',0,2,ierbc)
      ier=ier+ierbc
C
      if(ier.ne.0) return
C
C  error check OK
C
C---------------------------------------
C
C  get a  temporary array for f -- will extend out by 1 zone in
C  each direction so that numerical derivative evaluation can be
C  done without a lot of special case logic near the edges...
C
      allocate(ftmp(0:nx+1,0:ny+1,0:nz+1))
C
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               ftmp(ix,iy,iz)=fherm(0,ix,iy,iz)
            enddo
         enddo
      enddo
C
C  also create expanded axes grids...
C
      xx(1:nx)=x
      yy(1:ny)=y
      zz(1:nz)=z
c
      xx(0)=2*x(1)-x(2)
      xx(nx+1)=2*x(nx)-x(nx-1)
      yy(0)=2*y(1)-y(2)
      yy(ny+1)=2*y(ny)-y(ny-1)
      zz(0)=2*z(1)-z(2)
      zz(nz+1)=2*z(nz)-z(nz-1)
C
C---------------------------------------
C
C  handle boundary conditions and create rows of extrapolated points
C  in ftmp.  first do ftmp(0,1:ny), ftmp(nx+1,1:ny),
C                then ftmp(1:nx,0), ftmp(1:nx,ny+1),
c                then ... fill in the corners
c
c  also, for ipx.le.1 fill in the bdy fherm(1,*) values;
c        for ipy.le.1 fill in the bdy fherm(2,*) values
c
c  x bc's
c
      do iz=1,nz
         do iy=1,ny
c
            cxp=(ftmp(2,iy,iz)-ftmp(1,iy,iz))/(xx(2)-xx(1))
            cxm=(ftmp(nx,iy,iz)-ftmp(nx-1,iy,iz))/(xx(nx)-xx(nx-1))
c
            if(ipx.eq.1) then
c
c  periodic BC
c
               if(nx.gt.2) then
                  cxpp=(ftmp(3,iy,iz)-ftmp(2,iy,iz))/(xx(3)-xx(2))
                  cxmm=(ftmp(nx-1,iy,iz)-ftmp(nx-2,iy,iz))/
     >                 (xx(nx-1)-xx(nx-2))
c
                  call r8akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,1,iy,iz))
                  fherm(1,nx,iy,iz)=fherm(1,1,iy,iz)
c
               else
                  fherm(1,1,iy,iz)=cxp  ! =cxm, nx=2
                  fherm(1,nx,iy,iz)=fherm(1,1,iy,iz)
               endif
c
               cxtrap0=cxm
               cxtrapn=cxp
c
            else if(ipx.eq.0) then
C
C  default BC -- standard numeric extrapolation
C
               if(nx.gt.2) then
                  cxpp=(ftmp(3,iy,iz)-ftmp(2,iy,iz))/(xx(3)-xx(2))
                  fherm(1,1,iy,iz)=1.5_r8*cxp-0.5_r8*cxpp
C
                  cxmm=(ftmp(nx-1,iy,iz)-ftmp(nx-2,iy,iz))/
     >                 (xx(nx-1)-xx(nx-2))
                  fherm(1,nx,iy,iz)=1.5_r8*cxm-0.5_r8*cxmm
c
               else
                  fherm(1,1,iy,iz)=cxp  ! =cxm, nx=2
                  fherm(1,nx,iy,iz)=fherm(1,1,iy,iz)
               endif
C
C  extrapolate to slope to ghost points just past bdy...
C
               cxtrap0=2.0_r8*fherm(1,1,iy,iz)-cxp
               cxtrapn=2.0_r8*fherm(1,nx,iy,iz)-cxm
C
            else
C
C  BC supplied by user.  Also use this for extrapolation...
C
               cxtrap0=2.0_r8*fherm(1,1,iy,iz)-cxp
               cxtrapn=2.0_r8*fherm(1,nx,iy,iz)-cxm
C
            endif
C
            ftmp(0,iy,iz)=ftmp(1,iy,iz)-cxtrap0*(xx(1)-xx(0))
            ftmp(nx+1,iy,iz)=ftmp(nx,iy,iz)+cxtrapn*(xx(nx+1)-xx(nx))
C
         enddo
      enddo
c
c  y bc's
c
      do iz=1,nz
         do ix=1,nx
c
            cyp=(ftmp(ix,2,iz)-ftmp(ix,1,iz))/(yy(2)-yy(1))
            cym=(ftmp(ix,ny,iz)-ftmp(ix,ny-1,iz))/(yy(ny)-yy(ny-1))
c
            if(ipy.eq.1) then
c
c  periodic BC
c
               if(ny.gt.2) then
                  cypp=(ftmp(ix,3,iz)-ftmp(ix,2,iz))/(yy(3)-yy(2))
                  cymm=(ftmp(ix,ny-1,iz)-ftmp(ix,ny-2,iz))/
     >                 (yy(ny-1)-yy(ny-2))
c
                  call r8akherm0(cymm,cym,cyp,cypp,wy,fherm(2,ix,1,iz))
                  fherm(2,ix,ny,iz)=fherm(2,ix,1,iz)
c
               else
                  fherm(2,ix,1,iz) = cyp  ! =cym, ny=2
                  fherm(2,ix,ny,iz)=fherm(2,ix,1,iz)
               endif
               cytrap0=cym
               cytrapn=cyp
c
            else if(ipy.eq.0) then
C
C  default BC -- standard numeric extrapolation
C
               if(ny.gt.2) then
                  cypp=(ftmp(ix,3,iz)-ftmp(ix,2,iz))/(yy(3)-yy(2))
                  fherm(2,ix,1,iz)=1.5_r8*cyp-0.5_r8*cypp
C
                  cymm=(ftmp(ix,ny-1,iz)-ftmp(ix,ny-2,iz))/
     >                 (yy(ny-1)-yy(ny-2))
                  fherm(2,ix,ny,iz)=1.5_r8*cym-0.5_r8*cymm
c
               else
                  fherm(2,ix,1,iz) = cyp  ! =cym, ny=2
                  fherm(2,ix,ny,iz)=fherm(2,ix,1,iz)
               endif
C
C  extrapolate to slope to ghost points just past bdy...
C
               cytrap0=2.0_r8*fherm(2,ix,1,iz)-cyp
               cytrapn=2.0_r8*fherm(2,ix,ny,iz)-cym
C
            else
C
C  BC supplied by user.  Also use this for extrapolation...
C
               cytrap0=2.0_r8*fherm(2,ix,1,iz)-cyp
               cytrapn=2.0_r8*fherm(2,ix,ny,iz)-cym
C
            endif
C
            ftmp(ix,0,iz)=ftmp(ix,1,iz)-cytrap0*(yy(1)-yy(0))
            ftmp(ix,ny+1,iz)=ftmp(ix,ny,iz)+cytrapn*(yy(ny+1)-yy(ny))
C
         enddo
      enddo
c
c  z bc's
c
      do iy=1,ny
         do ix=1,nx
c
            czp=(ftmp(ix,iy,2)-ftmp(ix,iy,1))/(zz(2)-zz(1))
            czm=(ftmp(ix,iy,nz)-ftmp(ix,iy,nz-1))/(zz(nz)-zz(nz-1))
c
            if(ipy.eq.1) then
c
c  periodic BC
c
               if(nz.gt.2) then
                  czpp=(ftmp(ix,iy,3)-ftmp(ix,iy,2))/(zz(3)-zz(2))
                  czmm=(ftmp(ix,iy,nz-1)-ftmp(ix,iy,nz-2))/
     >                 (zz(nz-1)-zz(nz-2))
c
                  call r8akherm0(czmm,czm,czp,czpp,wz,fherm(3,ix,iy,1))
                  fherm(3,ix,iy,nz)=fherm(3,ix,iy,1)
c
               else
                  fherm(3,ix,iy,1) = czp
                  fherm(3,ix,iy,nz)=fherm(3,ix,iy,1)
               endif
c
               cztrap0=czm
               cztrapn=czp
c
            else if(ipy.eq.0) then
C
C  default BC -- standard numeric extrapolation
C
               if(nz.gt.2) then
                  czpp=(ftmp(ix,iy,3)-ftmp(ix,iy,2))/(zz(3)-zz(2))
                  fherm(3,ix,iy,1)=1.5_r8*czp-0.5_r8*czpp
C
                  czmm=(ftmp(ix,iy,nz-1)-ftmp(ix,iy,nz-2))/
     >                 (zz(nz-1)-zz(nz-2))
                  fherm(3,ix,iy,nz)=1.5_r8*czm-0.5_r8*czmm
c
               else
                  fherm(3,ix,iy,1) = czp
                  fherm(3,ix,iy,nz)=fherm(3,ix,iy,1)
               endif
C
C  extrapolate to slope to ghost points just past bdy...
C
               cztrap0=2.0_r8*fherm(3,ix,iy,1)-czp
               cztrapn=2.0_r8*fherm(3,ix,iy,nz)-czm
C
            else
C
C  BC supplied by user.  Also use this for extrapolation...
C
               cztrap0=2.0_r8*fherm(3,ix,iy,1)-czp
               cztrapn=2.0_r8*fherm(3,ix,iy,nz)-czm
C
            endif
C
            ftmp(ix,iy,0)=ftmp(ix,iy,1)-cztrap0*(zz(1)-zz(0))
            ftmp(ix,iy,nz+1)=ftmp(ix,iy,nz)+cztrapn*(zz(nz+1)-zz(nz))
C
         enddo
      enddo
c
c  and for the edges...
c
      do iz=1,nz
         do ix=0,1
            icx=ix*(nx+1)
            incx=1-2*ix
            ix1=icx+incx
            do iy=0,1
               icy=iy*(ny+1)
               incy=1-2*iy
               iy1=icy+incy
C
               ftmp(icx,icy,iz)=ftmp(icx,iy1,iz)+ftmp(ix1,icy,iz)-
     >              ftmp(ix1,iy1,iz)
C
            enddo
         enddo
      enddo
c
      do iy=1,ny
         do iz=0,1
            icz=iz*(nz+1)
            incz=1-2*iz
            iz1=icz+incz
            do ix=0,1
               icx=ix*(nx+1)
               incx=1-2*ix
               ix1=icx+incx
C
               ftmp(icx,iy,icz)=ftmp(icx,iy,iz1)+ftmp(ix1,iy,icz)-
     >              ftmp(ix1,iy,iz1)
C
            enddo
         enddo
      enddo
c
      do ix=1,nx
         do iy=0,1
            icy=iy*(ny+1)
            incy=1-2*iy
            iy1=icy+incy
            do iz=0,1
               icz=iz*(nz+1)
               incz=1-2*iz
               iz1=icz+incz
C
               ftmp(ix,icy,icz)=ftmp(ix,icy,iz1)+ftmp(ix,iy1,icz)-
     >              ftmp(ix,iy1,iz1)
C
            enddo
         enddo
      enddo
c
c  and the corners...
c
      do iz=0,1
         icz=iz*(nz+1)
         incz=1-2*iz
         iz1=icz+incz
         do iy=0,1
            icy=iy*(ny+1)
            incy=1-2*iy
            iy1=icy+incy
            do ix=0,1
               icx=ix*(nx+1)
               incx=1-2*ix
               ix1=icx+incx
C
               ftmp(icx,icy,icz)=(
     >              2*(ftmp(icx,icy,iz1)+ftmp(icx,iy1,icz)+
     >              ftmp(ix1,icy,icz)) -
     >              (ftmp(icx,iy1,iz1)+ftmp(ix1,icy,iz1)+
     >              ftmp(ix1,iy1,icz)))/3
C
            enddo
         enddo
      enddo
c
c-----------------------------------------------------
c  boundary done; now fill in the interior
c
c
      do iz=1,nz
         izm2=iz
         izm1=izm2-1
c
         izmm2=iz-1
         izmm1=izmm2-1
c
         izp2=iz+1
         izp1=izp2-1
c
         izpp2=iz+2
         izpp1=izpp2-1
c
         do iy=1,ny
            iym2=iy
            iym1=iym2-1
c
            iymm2=iy-1
            iymm1=iymm2-1
c
            iyp2=iy+1
            iyp1=iyp2-1
c
            iypp2=iy+2
            iypp1=iypp2-1
c
            do ix=1,nx
c
c  x div. diffs in vicinity
c
               ixm2=ix
               ixm1=ixm2-1
c
               ixmm2=ix-1
               ixmm1=ixmm2-1
c
               iflagx=0
               cxm=(ftmp(ixm2,iy,iz)-ftmp(ixm1,iy,iz))/
     >            (xx(ixm2)-xx(ixm1))
               if(ix.gt.1) then
                  cxmm=(ftmp(ixmm2,iy,iz)-ftmp(ixmm1,iy,iz))/
     >               (xx(ixmm2)-xx(ixmm1))
               else
                  if(ipx.eq.1) then
                     cxmm=(ftmp(nx-1,iy,iz)-ftmp(nx-2,iy,iz))/
     >                  (xx(nx-1)-xx(nx-2))
                  else
                     iflagx=1
                  endif
               endif
c
               ixp2=ix+1
               ixp1=ixp2-1
c
               ixpp2=ix+2
               ixpp1=ixpp2-1
c
               cxp=(ftmp(ixp2,iy,iz)-ftmp(ixp1,iy,iz))/
     >            (xx(ixp2)-xx(ixp1))
               if(ix.lt.nx) then
                  cxpp=(ftmp(ixpp2,iy,iz)-ftmp(ixpp1,iy,iz))/
     >               (xx(ixpp2)-xx(ixpp1))
               else
                  if(ipx.eq.1) then
                     cxpp=(ftmp(3,iy,iz)-ftmp(2,iy,iz))/
     >                  (xx(3)-xx(2))
                  else
                     cxpp=cxp+(cxm-cxmm)
                  endif
               endif
c
               if(iflagx.eq.1) then
                  cxmm=cxm+(cxp-cxpp)
               endif
c
c  Akima weightings + df/dx for interior pts
c
               call r8akherm0(cxmm,cxm,cxp,cxpp,wx,zansr)
               if((ix.gt.1).and.(ix.lt.nx)) fherm(1,ix,iy,iz)=zansr
c
c  y div. diffs in vicinity
c
               iflagy=0
               cym=(ftmp(ix,iym2,iz)-ftmp(ix,iym1,iz))/
     >            (yy(iym2)-yy(iym1))
               if(iy.gt.1) then
                  cymm=(ftmp(ix,iymm2,iz)-ftmp(ix,iymm1,iz))/
     >               (yy(iymm2)-yy(iymm1))
               else
                  if(ipy.eq.1) then
                     cymm=(ftmp(ix,ny-1,iz)-ftmp(ix,ny-2,iz))/
     >                  (yy(ny-1)-yy(ny-2))
                  else
                     iflagy=1
                  endif
               endif
c
               cyp=(ftmp(ix,iyp2,iz)-ftmp(ix,iyp1,iz))/
     >            (yy(iyp2)-yy(iyp1))
               if(iy.lt.ny) then
                  cypp=(ftmp(ix,iypp2,iz)-ftmp(ix,iypp1,iz))/
     >               (yy(iypp2)-yy(iypp1))
               else
                  if(ipy.eq.1) then
                     cypp=(ftmp(ix,3,iz)-ftmp(ix,2,iz))/(yy(3)-yy(2))
                  else
                     cypp=cyp+(cym-cymm)
                  endif
               endif
c
               if(iflagy.eq.1) then
                  cymm=cym+(cyp-cypp)
               endif
c
c  Akima weightings + df/dy for interior pts
c
               call r8akherm0(cymm,cym,cyp,cypp,wy,zansr)
               if((iy.gt.1).and.(iy.lt.ny)) fherm(2,ix,iy,iz)=zansr
c
c  z div. diffs in vicinity
c
               iflagz=0
               czm=(ftmp(ix,iy,izm2)-ftmp(ix,iy,izm1))/
     >            (zz(izm2)-zz(izm1))
               if(iz.gt.1) then
                  czmm=(ftmp(ix,iy,izmm2)-ftmp(ix,iy,izmm1))/
     >               (zz(izmm2)-zz(izmm1))
               else
                  if(ipz.eq.1) then
                     czmm=(ftmp(ix,iy,nz-1)-ftmp(ix,iy,nz-2))/
     >                  (zz(nz-1)-zz(nz-2))
                  else
                     iflagz=1
                  endif
               endif
c
               czp=(ftmp(ix,iy,izp2)-ftmp(ix,iy,izp1))/
     >            (zz(izp2)-zz(izp1))
               if(iz.lt.nz) then
                  czpp=(ftmp(ix,iy,izpp2)-ftmp(ix,iy,izpp1))/
     >               (zz(izpp2)-zz(izpp1))
               else
                  if(ipz.eq.1) then
                     czpp=(ftmp(ix,iy,3)-ftmp(ix,iy,2))/(zz(3)-zz(2))
                  else
                     czpp=czp+(czm-czmm)
                  endif
               endif
c
               if(iflagz.eq.1) then
                  czmm=czm+(czp-czpp)
               endif
c
c  Akima weightings + df/dy for interior pts
c
               call r8akherm0(czmm,czm,czp,czpp,wz,zansr)
               if((iz.gt.1).and.(iz.lt.nz)) fherm(3,ix,iy,iz)=zansr
c
c  cross derivatives (2nd order divided differences)
c
               cxm2=(ftmp(ixm2,iym1,iz)-ftmp(ixm1,iym1,iz))/
     >            (xx(ixm2)-xx(ixm1))
               exy(1,1)=(cxm-cxm2)/(yy(iym2)-yy(iym1))
c
               cxm2=(ftmp(ixm2,iyp2,iz)-ftmp(ixm1,iyp2,iz))/
     >            (xx(ixm2)-xx(ixm1))
               exy(1,2)=(cxm2-cxm)/(yy(iyp2)-yy(iyp1))
c
               cxp2=(ftmp(ixp2,iym1,iz)-ftmp(ixp1,iym1,iz))/
     >            (xx(ixp2)-xx(ixp1))
               exy(2,1)=(cxp-cxp2)/(yy(iym2)-yy(iym1))
c
               cxp2=(ftmp(ixp2,iyp2,iz)-ftmp(ixp1,iyp2,iz))/
     >            (xx(ixp2)-xx(ixp1))
               exy(2,2)=(cxp2-cxp)/(yy(iyp2)-yy(iyp1))
c
c---------
c
               cxm2=(ftmp(ixm2,iy,izm1)-ftmp(ixm1,iy,izm1))/
     >            (xx(ixm2)-xx(ixm1))
               exz(1,1)=(cxm-cxm2)/(zz(izm2)-zz(izm1))
c
               cxm2=(ftmp(ixm2,iy,izp2)-ftmp(ixm1,iy,izp2))/
     >            (xx(ixm2)-xx(ixm1))
               exz(1,2)=(cxm2-cxm)/(zz(izp2)-zz(izp1))
c
               cxp2=(ftmp(ixp2,iy,izm1)-ftmp(ixp1,iy,izm1))/
     >            (xx(ixp2)-xx(ixp1))
               exz(2,1)=(cxp-cxp2)/(zz(izm2)-zz(izm1))
c
               cxp2=(ftmp(ixp2,iy,izp2)-ftmp(ixp1,iy,izp2))/
     >            (xx(ixp2)-xx(ixp1))
               exz(2,2)=(cxp2-cxp)/(zz(izp2)-zz(izp1))
c
c---------
c
               cym2=(ftmp(ix,iym2,izm1)-ftmp(ix,iym1,izm1))/
     >            (yy(iym2)-yy(iym1))
               eyz(1,1)=(cym-cym2)/(zz(izm2)-zz(izm1))
c
               cym2=(ftmp(ix,iym2,izp2)-ftmp(ix,iym1,izp2))/
     >            (yy(iym2)-yy(iym1))
               eyz(1,2)=(cym2-cym)/(zz(izp2)-zz(izp1))
c
               cyp2=(ftmp(ix,iyp2,izm1)-ftmp(ix,iyp1,izm1))/
     >            (yy(iyp2)-yy(iyp1))
               eyz(2,1)=(cyp-cyp2)/(zz(izm2)-zz(izm1))
c
               cyp2=(ftmp(ix,iyp2,izp2)-ftmp(ix,iyp1,izp2))/
     >            (yy(iyp2)-yy(iyp1))
               eyz(2,2)=(cyp2-cyp)/(zz(izp2)-zz(izp1))
c
c  dxdydz 3rd order cross deriv.
c
               k=1
c
               fzp=(ftmp(ixm2,iym2,izm2)-ftmp(ixm2,iym1,izm2)-
     >            ftmp(ixm1,iym2,izm2)+ftmp(ixm1,iym1,izm2))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iym2)-yy(iym1)))
               fzm=(ftmp(ixm2,iym2,izm1)-ftmp(ixm2,iym1,izm1)-
     >            ftmp(ixm1,iym2,izm1)+ftmp(ixm1,iym1,izm1))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iym2)-yy(iym1)))
               fxyz(1,1,k)=(fzp-fzm)/(zz(izm2)-zz(izm1))
c
               fzp=(ftmp(ixm2,iyp2,izm2)-ftmp(ixm2,iyp1,izm2)-
     >            ftmp(ixm1,iyp2,izm2)+ftmp(ixm1,iyp1,izm2))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iyp2)-yy(iyp1)))
               fzm=(ftmp(ixm2,iyp2,izm1)-ftmp(ixm2,iyp1,izm1)-
     >            ftmp(ixm1,iyp2,izm1)+ftmp(ixm1,iyp1,izm1))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iyp2)-yy(iyp1)))
               fxyz(1,2,k)=(fzp-fzm)/(zz(izm2)-zz(izm1))
c
               fzp=(ftmp(ixp2,iym2,izm2)-ftmp(ixp2,iym1,izm2)-
     >            ftmp(ixp1,iym2,izm2)+ftmp(ixp1,iym1,izm2))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iym2)-yy(iym1)))
               fzm=(ftmp(ixp2,iym2,izm1)-ftmp(ixp2,iym1,izm1)-
     >            ftmp(ixp1,iym2,izm1)+ftmp(ixp1,iym1,izm1))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iym2)-yy(iym1)))
               fxyz(2,1,k)=(fzp-fzm)/(zz(izm2)-zz(izm1))
c
               fzp=(ftmp(ixp2,iyp2,izm2)-ftmp(ixp2,iyp1,izm2)-
     >            ftmp(ixp1,iyp2,izm2)+ftmp(ixp1,iyp1,izm2))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iyp2)-yy(iyp1)))
               fzm=(ftmp(ixp2,iyp2,izm1)-ftmp(ixp2,iyp1,izm1)-
     >            ftmp(ixp1,iyp2,izm1)+ftmp(ixp1,iyp1,izm1))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iyp2)-yy(iyp1)))
               fxyz(2,2,k)=(fzp-fzm)/(zz(izm2)-zz(izm1))
c
               k=2
c
               fzp=(ftmp(ixm2,iym2,izp2)-ftmp(ixm2,iym1,izp2)-
     >            ftmp(ixm1,iym2,izp2)+ftmp(ixm1,iym1,izp2))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iym2)-yy(iym1)))
               fzm=(ftmp(ixm2,iym2,izp1)-ftmp(ixm2,iym1,izp1)-
     >            ftmp(ixm1,iym2,izp1)+ftmp(ixm1,iym1,izp1))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iym2)-yy(iym1)))
               fxyz(1,1,k)=(fzp-fzm)/(zz(izp2)-zz(izp1))
c
               fzp=(ftmp(ixm2,iyp2,izp2)-ftmp(ixm2,iyp1,izp2)-
     >            ftmp(ixm1,iyp2,izp2)+ftmp(ixm1,iyp1,izp2))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iyp2)-yy(iyp1)))
               fzm=(ftmp(ixm2,iyp2,izp1)-ftmp(ixm2,iyp1,izp1)-
     >            ftmp(ixm1,iyp2,izp1)+ftmp(ixm1,iyp1,izp1))/
     >            ((xx(ixm2)-xx(ixm1))*(yy(iyp2)-yy(iyp1)))
               fxyz(1,2,k)=(fzp-fzm)/(zz(izp2)-zz(izp1))
c
               fzp=(ftmp(ixp2,iym2,izp2)-ftmp(ixp2,iym1,izp2)-
     >            ftmp(ixp1,iym2,izp2)+ftmp(ixp1,iym1,izp2))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iym2)-yy(iym1)))
               fzm=(ftmp(ixp2,iym2,izp1)-ftmp(ixp2,iym1,izp1)-
     >            ftmp(ixp1,iym2,izp1)+ftmp(ixp1,iym1,izp1))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iym2)-yy(iym1)))
               fxyz(2,1,k)=(fzp-fzm)/(zz(izp2)-zz(izp1))
c
               fzp=(ftmp(ixp2,iyp2,izp2)-ftmp(ixp2,iyp1,izp2)-
     >            ftmp(ixp1,iyp2,izp2)+ftmp(ixp1,iyp1,izp2))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iyp2)-yy(iyp1)))
               fzm=(ftmp(ixp2,iyp2,izp1)-ftmp(ixp2,iyp1,izp1)-
     >            ftmp(ixp1,iyp2,izp1)+ftmp(ixp1,iyp1,izp1))/
     >            ((xx(ixp2)-xx(ixp1))*(yy(iyp2)-yy(iyp1)))
               fxyz(2,2,k)=(fzp-fzm)/(zz(izp2)-zz(izp1))
c
c  the values
c
               fherm(4,ix,iy,iz)=
     >            (wx(1)*(wy(1)*exy(1,1)+wy(2)*exy(1,2))+
     >             wx(2)*(wy(1)*exy(2,1)+wy(2)*exy(2,2)))/
     >             ((wx(1)+wx(2))*(wy(1)+wy(2)))
c
               fherm(5,ix,iy,iz)=
     >            (wx(1)*(wz(1)*exz(1,1)+wz(2)*exz(1,2))+
     >             wx(2)*(wz(1)*exz(2,1)+wz(2)*exz(2,2)))/
     >             ((wx(1)+wx(2))*(wz(1)+wz(2)))
c
               fherm(6,ix,iy,iz)=
     >            (wy(1)*(wz(1)*eyz(1,1)+wz(2)*eyz(1,2))+
     >             wy(2)*(wz(1)*eyz(2,1)+wz(2)*eyz(2,2)))/
     >             ((wy(1)+wy(2))*(wz(1)+wz(2)))
c
               fherm(7,ix,iy,iz)=(wz(1)*
     >            (wx(1)*(wy(1)*fxyz(1,1,1)+wy(2)*fxyz(1,2,1))+
     >            wx(2)*(wy(1)*fxyz(2,1,1)+wy(2)*fxyz(2,2,1))) + wz(2)*
     >            (wx(1)*(wy(1)*fxyz(1,1,2)+wy(2)*fxyz(1,2,2))+
     >            wx(2)*(wy(1)*fxyz(2,1,2)+wy(2)*fxyz(2,2,2))) )/
     >            ((wx(1)+wx(2))*(wy(1)+wy(2))*(wz(1)+wz(2)))
c
            enddo
         enddo
      enddo
C
      deallocate (ftmp)
      return
      end
