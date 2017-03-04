      subroutine pc2ev(xget,yget,x,nx,y,ny,ilinx,iliny,
     >                   f,inf2,ict,fval,ier)
C
C  evaluate a piecewise bilinear interpolant on a rectilinear
C  grid -- this is C1 in both directions.
C
C  this subroutine calls two subroutines:
C     herm2xy  -- find cell containing (xget,yget)
C     pc2fcn -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
      real xget,yget                    ! target of this interpolation
      real x(nx)                        ! ordered x grid
      real y(ny)                        ! ordered y grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
C
      real f(inf2,ny)                   ! function data
C
C       f 2nd dimension inf2 must be .ge. nx
C       contents of f:
C
C  f(i,j) = f @ x(i),y(j)
C
      integer ict(4)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return d2f/dxdy (0, don't)
C
C output arguments:
C =================
C
      real fval(*)                      ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the fourth output (depends on ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,1]
C   on output fval= [f,df/dx,df/dy,d2f/dxdy]
C
C    on input ict = [1,0,0,0]
C   on output fval= [f] ... elements 2 & 3 & 4 never referenced
C
C    on input ict = [0,1,1,0]
C   on output fval= [df/dx,df/dy] ... element 3 & 4 never referenced
C
C    on input ict = [0,0,1,0]
C   on output fval= [df/dy] ... elements 2 & 3 & 4 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i,j                       ! cell indices
C
C  normalized displacement from (x(i),y(j)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C
      real xparam,yparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real hx,hy
      real hxi,hyi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C
C---------------------------------------------------------------------
C
      call herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
      if(ier.ne.0) return
c
      call pc2fcn(ict,1,1,
     >   fval,i,j,xparam,yparam,hx,hxi,hy,hyi,f,inf2,ny)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate piecewise bilinear function interpolation -- 2d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine pc2fcn(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   fin,inf2,ny)
C
      integer ict(4)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec)         ! target cells (i,j)
      real xparam(ivec),yparam(ivec)
                          ! normalized displacements from (i,j) corners
C
      real hx(ivec),hy(ivec)            ! grid spacing, and
      real hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))
C
      real fin(inf2,ny)                 ! interpolant data (cf "pc2ev")
C
      real fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine
C  pc2ev comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C
      real sum
      integer v
C
C   ...in x direction
C
      do v=1,ivec
         i=ii(v)
         j=jj(v)
c
         xp=xparam(v)
         xpi=1.0-xp
C
C   ...in y direction
C
         yp=yparam(v)
         ypi=1.0-yp
C
         iadr=0
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            sum=ypi*(xpi*fin(i,j)+xp*fin(i+1,j))
     >         + yp*(xpi*fin(i,j+1)+xp*fin(i+1,j+1))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
C
            sum=ypi*(fin(i+1,j)-fin(i,j))
     >         + yp*(fin(i+1,j+1)-fin(i,j+1))
            fval(v,iadr)=sum*hxi(v)
C
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
C
            sum=xpi*(fin(i,j+1)-fin(i,j))
     >         + xp*(fin(i+1,j+1)-fin(i+1,j))
            fval(v,iadr)=sum*hyi(v)
         endif
C
         if(ict(4).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
C
            sum=fin(i+1,j+1)-fin(i,j+1)-fin(i+1,j)+fin(i,j)
            fval(v,iadr)=sum*hxi(v)*hyi(v)
         endif
C
      enddo                             ! vector loop
C
      return
      end
C---------------------------------------------------------------------
C  evaluate piecewise bilinear function interpolation -- 2d fcn
C   --vectorized-- dmc 10 Feb 1999
C    --optimized for VARIATION along x axis ONLY--
C
      subroutine pc2fcnx(ict,ivec,ivecd,
     >   fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >   fin,inf2,ny)
C
      integer ict(4)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj               ! target cells (i,j)
      real xparam(ivec),yparam
                          ! normalized displacements from (i,j) corners
C
      real hx(ivec),hy                  ! grid spacing, and
      real hxi(ivec),hyi                ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))
C
      real fin(inf2,ny)                 ! interpolant data (cf "pc2ev")
C
      real fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine
C  pc2ev comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C
      real sum
      integer v
C
C   ...in x direction
C
      j=jj
C
C   ...in y direction
C
      yp=yparam
      ypi=1.0-yp
C
      do v=1,ivec
         i=ii(v)
c
         xp=xparam(v)
         xpi=1.0-xp
C
         iadr=0
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            sum=ypi*(xpi*fin(i,j)+xp*fin(i+1,j))
     >         + yp*(xpi*fin(i,j+1)+xp*fin(i+1,j+1))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
C
            sum=ypi*(fin(i+1,j)-fin(i,j))
     >         + yp*(fin(i+1,j+1)-fin(i,j+1))
            fval(v,iadr)=sum*hxi(v)
C
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
C
            sum=xpi*(fin(i,j+1)-fin(i,j))
     >         + xp*(fin(i+1,j+1)-fin(i+1,j))
            fval(v,iadr)=sum*hyi
         endif
C
         if(ict(4).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
C
            sum=fin(i+1,j+1)-fin(i,j+1)-fin(i+1,j)+fin(i,j)
            fval(v,iadr)=sum*hxi(v)*hyi
         endif
C
      enddo                             ! vector loop
C
      return
      end
