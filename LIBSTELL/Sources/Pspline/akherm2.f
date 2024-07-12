      subroutine akherm2(x,nx,y,ny,fherm,nf2,ilinx,iliny,ier)
C
C  create a data set for Hermite interpolation, based on Akima's method
C  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
C
C  input:
C
      integer nx,ny,nf1                 ! array dimensions
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
C  a default boundary condition is used, based on divided differences
C  in the edge region.  For more control of BC, use akherm2p...
C
      call akherm2p(x,nx,y,ny,fherm,nf2,ilinx,iliny,0,0,ier)
C
      return
      end
C----------------------------
      subroutine akherm2p(x,nx,y,ny,fherm,nf2,ilinx,iliny,ipx,ipy,ier)
C
C  create a data set for Hermite interpolation, based on Akima's method
C  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
C
C  with independently settable boundary condition options:
C     ipx or ipy  =0  -- default, boundary conditions from divided diffs
C     ipx or ipy  =1  -- periodic boundary condition
C     ipx or ipy  =2  -- user supplied df/dx or df/dy
C  input:
C
      integer nx,ny,nf1                 ! array dimensions
      real x(nx)                        ! x coordinate array
      real y(ny)                        ! y coordinate array
      real fherm(0:3,nf2,ny)            ! data/Hermite array
C
      integer ipx                       ! =1 if df/dx periodic in x
      integer ipy                       ! =1 if df/dy periodic in y
C
C  fherm(0,1:nx,1:ny) supplied by user; this routine computes the
C  rest of the elements, but note:
C
C    if ipx=2:  fherm(1,1,1:ny) & fherm(1,nx,1:ny) = INPUT df/dx BCs
C    if ipy=2:  fherm(2,1:nx,1) & fherm(2,1:nx,ny) = INPUT df/dy BCs
C
C  on output, at all grid points (i,j) covering (1:nx,1:ny):
C  fherm1(1,i,j) -- df/dx at the grid point
C  fherm1(2,i,j) -- df/dy at the grid point
C  fherm1(3,i,j) -- d2f/dxdy at the grid point
C
C---------------------
C  local...
C
      real wx(2),wy(2),e(2,2)
c
      real, dimension(:,:), allocatable ::  ftmp
c
      real xx(0:nx+1)
      real yy(0:ny+1)
c
c---------------------
c
c  error checks
c
      ier=0
c
      call splinck(x,nx,ilinx,1.0e-3,ierx)
      if(ierx.ne.0) ier=ier+1
c
      if(ierx.ne.0) then
         write(6,'('' ?akherm2:  x axis not strict ascending'')')
      endif
c
      call splinck(y,ny,iliny,1.0e-3,iery)
      if(iery.ne.0) ier=ier+1
c
      if(iery.ne.0) then
         write(6,'('' ?akherm2:  y axis not strict ascending'')')
      endif
c
      if(nf2.lt.nx) then
         ier=ier+1
         write(6,*) '?akherm2:  fherm array dimension too small.'
      endif
C
      ierbc=0
      call ibc_ck(ipx,'akherm2','X Bdy Cond',0,2,ierbc)
      ier=ier+ierbc
C
      ierbc=0
      call ibc_ck(ipy,'akherm2','Y Bdy Cond',0,2,ierbc)
      ier=ier+ierbc
C
      if(ier.ne.0) return
C
C---------------------------------------
C
C  get a  temporary array for f -- will extend out by 1 zone in
C  each direction so that numerical derivative evaluation can be
C  done without a lot of special case logic near the edges...
C
      allocate(ftmp(0:nx+1,0:ny+1))
C
      do iy=1,ny
         do ix=1,nx
            ftmp(ix,iy)=fherm(0,ix,iy)
         enddo
      enddo
C
C  also create expanded axes grids...
C
      xx(1:nx)=x
      yy(1:ny)=y
      xx(0)=2*x(1)-x(2)
      xx(nx+1)=2*x(nx)-x(nx-1)
      yy(0)=2*y(1)-y(2)
      yy(ny+1)=2*y(ny)-y(ny-1)
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
      do iy=1,ny
c
         cxp=(ftmp(2,iy)-ftmp(1,iy))/(xx(2)-xx(1))
         cxm=(ftmp(nx,iy)-ftmp(nx-1,iy))/(xx(nx)-xx(nx-1))
c
         if(ipx.eq.1) then
c
c  periodic BC
c
            if(nx.gt.2) then
               cxpp=(ftmp(3,iy)-ftmp(2,iy))/(xx(3)-xx(2))
               cxmm=(ftmp(nx-1,iy)-ftmp(nx-2,iy))/(xx(nx-1)-xx(nx-2))
c
               call akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,1,iy))
               fherm(1,nx,iy)=fherm(1,1,iy)
            else
               fherm(1,1,iy) = cxp  ! =cxm, nx=2
               fherm(1,nx,iy) = fherm(1,1,iy)
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
               cxpp=(ftmp(3,iy)-ftmp(2,iy))/(xx(3)-xx(2))
               fherm(1,1,iy)=1.5*cxp-0.5*cxpp
C
               cxmm=(ftmp(nx-1,iy)-ftmp(nx-2,iy))/(xx(nx-1)-xx(nx-2))
               fherm(1,nx,iy)=1.5*cxm-0.5*cxmm
c
            else
               fherm(1,1,iy) = cxp  ! =cxm, nx=2
               fherm(1,nx,iy) = fherm(1,1,iy)
            endif
C
C  extrapolate to slope to ghost points just past bdy...
C
            cxtrap0=2.0*fherm(1,1,iy)-cxp
            cxtrapn=2.0*fherm(1,nx,iy)-cxm
C
         else
C
C  BC supplied by user.  Also use this for extrapolation...
C
            cxtrap0=2.0*fherm(1,1,iy)-cxp
            cxtrapn=2.0*fherm(1,nx,iy)-cxm
C
         endif
C
         ftmp(0,iy)=ftmp(1,iy)-cxtrap0*(xx(1)-xx(0))
         ftmp(nx+1,iy)=ftmp(nx,iy)+cxtrapn*(xx(nx+1)-xx(nx))
C
      enddo
c
c  y bc's
c
      do ix=1,nx
c
         cyp=(ftmp(ix,2)-ftmp(ix,1))/(yy(2)-yy(1))
         cym=(ftmp(ix,ny)-ftmp(ix,ny-1))/(yy(ny)-yy(ny-1))
c
         if(ipy.eq.1) then
c
c  periodic BC
c
            if(ny.gt.2) then
               cypp=(ftmp(ix,3)-ftmp(ix,2))/(yy(3)-yy(2))
               cymm=(ftmp(ix,ny-1)-ftmp(ix,ny-2))/(yy(ny-1)-yy(ny-2))
c
               call akherm0(cymm,cym,cyp,cypp,wy,fherm(2,ix,1))
               fherm(2,ix,ny)=fherm(2,ix,1)
c
            else
               fherm(2,ix,1) = cyp  ! =cym, ny=2
               fherm(2,ix,ny)=fherm(2,ix,1)
            endif
c
            cytrap0=cym
            cytrapn=cyp
c
         else if(ipy.eq.0) then
C
C  default BC -- standard numeric extrapolation
C
            if(ny.gt.2) then
               cypp=(ftmp(ix,3)-ftmp(ix,2))/(yy(3)-yy(2))
               fherm(2,ix,1)=1.5*cyp-0.5*cypp
C
               cymm=(ftmp(ix,ny-1)-ftmp(ix,ny-2))/(yy(ny-1)-yy(ny-2))
               fherm(2,ix,ny)=1.5*cym-0.5*cymm
c
            else
               fherm(2,ix,1) = cyp  ! =cym, ny=2
               fherm(2,ix,ny)=fherm(2,ix,1)
            endif
C
C  extrapolate to slope to ghost points just past bdy...
C
            cytrap0=2.0*fherm(2,ix,1)-cyp
            cytrapn=2.0*fherm(2,ix,ny)-cym
C
         else
C
C  BC supplied by user.  Also use this for extrapolation...
C
            cytrap0=2.0*fherm(2,ix,1)-cyp
            cytrapn=2.0*fherm(2,ix,ny)-cym
C
         endif
C
         ftmp(ix,0)=ftmp(ix,1)-cytrap0*(yy(1)-yy(0))
         ftmp(ix,ny+1)=ftmp(ix,ny)+cytrapn*(yy(ny+1)-yy(ny))
C
      enddo
C
C  and do something for the corners...
C
      do ix=0,1
         do iy=0,1
            icx=ix*(nx+1)
            icy=iy*(ny+1)
            incx=1-2*ix
            incy=1-2*iy

            ix1=icx+incx
            iy1=icy+incy

            ftmp(icx,icy)=ftmp(icx,iy1)+ftmp(ix1,icy)-ftmp(ix1,iy1)
C
         enddo
      enddo
C
C----------------------------------------------------------------
C  OK, now ready to compute all the interior coefficients and the
C  rest of the edge coefficients as well...
C
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
            cxm=(ftmp(ixm2,iy)-ftmp(ixm1,iy))/(xx(ixm2)-xx(ixm1))
            if(ix.gt.1) then
               cxmm=(ftmp(ixmm2,iy)-ftmp(ixmm1,iy))/
     >            (xx(ixmm2)-xx(ixmm1))
            else
               if(ipx.eq.1) then
                  cxmm=(ftmp(nx-1,iy)-ftmp(nx-2,iy))/
     >               (xx(nx-1)-xx(nx-2))
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
            cxp=(ftmp(ixp2,iy)-ftmp(ixp1,iy))/(xx(ixp2)-xx(ixp1))
            if(ix.lt.nx) then
               cxpp=(ftmp(ixpp2,iy)-ftmp(ixpp1,iy))/
     >            (xx(ixpp2)-xx(ixpp1))
            else
               if(ipx.eq.1) then
                  cxpp=(ftmp(3,iy)-ftmp(2,iy))/(xx(3)-xx(2))
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
            call akherm0(cxmm,cxm,cxp,cxpp,wx,zansr)
            if((ix.gt.1).and.(ix.lt.nx)) fherm(1,ix,iy)=zansr
c
c  y div. diffs in vicinity
c
            iflagy=0
            cym=(ftmp(ix,iym2)-ftmp(ix,iym1))/(yy(iym2)-yy(iym1))
            if(iy.gt.1) then
               cymm=(ftmp(ix,iymm2)-ftmp(ix,iymm1))/
     >            (yy(iymm2)-yy(iymm1))
            else
               if(ipy.eq.1) then
                  cymm=(ftmp(ix,ny-1)-ftmp(ix,ny-2))/
     >               (yy(ny-1)-yy(ny-2))
               else
                  iflagy=1
               endif
            endif
c
            cyp=(ftmp(ix,iyp2)-ftmp(ix,iyp1))/(yy(iyp2)-yy(iyp1))
            if(iy.lt.ny) then
               cypp=(ftmp(ix,iypp2)-ftmp(ix,iypp1))/
     >            (yy(iypp2)-yy(iypp1))
            else
               if(ipy.eq.1) then
                  cypp=(ftmp(ix,3)-ftmp(ix,2))/(yy(3)-yy(2))
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
            call akherm0(cymm,cym,cyp,cypp,wy,zansr)
            if((iy.gt.1).and.(iy.lt.ny)) fherm(2,ix,iy)=zansr
c
c  cross derivatives (2nd order divided differences)
c
            cxm2=(ftmp(ixm2,iym1)-ftmp(ixm1,iym1))/
     >         (xx(ixm2)-xx(ixm1))
            e(1,1)=(cxm-cxm2)/(yy(iym2)-yy(iym1))
c
            cxm2=(ftmp(ixm2,iyp2)-ftmp(ixm1,iyp2))/
     >         (xx(ixm2)-xx(ixm1))
            e(1,2)=(cxm2-cxm)/(yy(iyp2)-yy(iyp1))
c
            cxp2=(ftmp(ixp2,iym1)-ftmp(ixp1,iym1))/
     >         (xx(ixp2)-xx(ixp1))
            e(2,1)=(cxp-cxp2)/(yy(iym2)-yy(iym1))
c
            cxp2=(ftmp(ixp2,iyp2)-ftmp(ixp1,iyp2))/
     >         (xx(ixp2)-xx(ixp1))
            e(2,2)=(cxp2-cxp)/(yy(iyp2)-yy(iyp1))
c
c  the values
c
            fherm(3,ix,iy)=(wx(1)*(wy(1)*e(1,1)+wy(2)*e(1,2))+
     >         wx(2)*(wy(1)*e(2,1)+wy(2)*e(2,2)))/
     >         ((wx(1)+wx(2))*(wy(1)+wy(2)))
c
         enddo
      enddo
C
      deallocate (ftmp)
      return
      end
