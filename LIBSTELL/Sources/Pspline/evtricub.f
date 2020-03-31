      subroutine evtricub(xget,yget,zget,x,nx,y,ny,z,nz,
     >                   ilinx,iliny,ilinz,
     >                   f,inf2,inf3,ict,fval,ier)
C
C  use mktricub to set up spline coefficients...
C
C  evaluate a 3d cubic Spline interpolant on a rectilinear
C  grid -- this is C2 in all directions.
C
C  this subroutine calls two subroutines:
C     herm3xyz  -- find cell containing (xget,yget,zget)
C     fvtricub  -- evaluate the spline function (w/derivatives if req.)
C
C  input arguments:
C  ================
C
      real xget,yget,zget               ! target of this interpolation
      real x(nx)                        ! ordered x grid
      real y(ny)                        ! ordered y grid
      real z(nz)                        ! ordered z grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
      integer ilinz                     ! ilinz=1 => assume z evenly spaced
C
      real f(0:7,inf2,inf3,nz)          ! function data
C
C       f 2nd dimension inf2 must be .ge. nx; 3rd dim inf3 .ge. ny
C       contents of f:
C
C  f(0,i,j,k) = f @ x(i),y(j),z(k)
C  f(1,i,j,k) = d2f/dx2 @ x(i),y(j),z(k)
C  f(2,i,j,k) = d2f/dy2 @ x(i),y(j),z(k)
C  f(3,i,j,k) = d2f/dz2 @ x(i),y(j),z(k)
C  f(4,i,j,k) = d4f/dx2dy2 @ x(i),y(j),z(k)
C  f(5,i,j,k) = d4f/dx2dz2 @ x(i),y(j),z(k)
C  f(6,i,j,k) = d4f/dy2dz2 @ x(i),y(j),z(k)
C  f(7,i,j,k) = d6f/dx2dy2dz2 @ x(i),y(j),z(k)
C
      integer ict(10)                   ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return df/dz  (0, don't)
C  ict(5)=1 -- return d2f/dx2  (0, don't)
C  ict(6)=1 -- return d2f/dy2  (0, don't)
C  ict(7)=1 -- return d2f/dz2  (0, don't)
C  ict(8)=1 -- return d2f/dxdy (0, don't)
C  ict(9)=1 -- return d2f/dxdz (0, don't)
C  ict(10)=1-- return d2f/dydz (0, don't)
c
c  (new dmc Dec 2005 -- higher derivatives available)
c    ict(1)=3 --> 3rd derivative, .le.2 diff. in any coordinate
c      ict(2:8) select: fxxy fxxz fxyy fxyz fxzz fyyz fyzz
c      ->note ict(1)=3, ict(5)=1 gives fxyz = d3f/dxdydz
c    ict(1)=-3 --> 3rd derivative, 3 in one coordinate
c      ict(2:4) select: fxxx fyyy fzzz
c    ict(1)=4 --> 3rd derivative, .le.2 diff. in any coordinate
c      ict(2:7) select: fxxyy fxxyz fxxzz fxyyz fxyzz fyyzz
c    ict(1)=-4 --> 3rd derivative, 3 in one coordinate
c      ict(2:7) select: fxxxy fxxxz fxyyy fxzzz fyyyz fyzzz
c    ict(1)=5 --> 3rd derivative, .le.2 diff. in any coordinate
c      ict(2:4) select: fxxyyz fxxyzz fxyyzz
c    ict(1)=-5 --> 3rd derivative, 3 in one coordinate
c      ict(2:10) select:  fxxxyy fxxxyz fxxxzz fxxyyy fxxzzz
c                         fxyyyz fxyzzz fyyyzz fzzzyy
c    ict(1)=6 --> 3rd derivative, .le.2 diff. in any coordinate
c      fxxyyzz
c    ict(1)=-6 --> 3rd derivative, 3 in one coordinate
c      ict(2:10) select: fxxxyyy fxxxyyz fxxxyzz fxxxyyz
c                        fxxyyyz fxxyzzz fxyyyzz fxyyzzz fyyyzzz
c    ict(1)=-7 --> 7th derivative
c      ict(2:7) select: fxxxyyyz fxxxyyzz fxxxyzzz
c                       fxxyyyzz fxxyyzzz fxyyyzzz
c    ict(1)=-8 --> 8th derivative
c      ict(2:4) select: fxxxyyyzz fxxxyyzzz fxxyyyzzz
c    ict(1)=-9 --> 9th derivative:  fxxxyyyzzz
c
C
C output arguments:
C =================
C
      real fval(*)                     ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the 4th output (depends on ict(...) spec)
C  fval(5-10) receive 5th thru 10th outputs (if required by ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,1,0,0,0,0,0,0,0]
C   on output fval= [f,df/dx,df/dy,df/dz]
C
C    on input ict = [1,0,0,0,0,0,0,0,0,0,0]
C   on output fval= [f] ... elements 2-10 never referenced
C
C    on input ict = [0,1,1,0,0,0,0,0,0,0,0]
C   on output fval= [df/dx,df/dy] ... elements 3-10 never referenced
C
C    on input ict = [0,0,0,0,1,0,0,0,0,0,0]
C   on output fval= [d2f/dx2] ... elements 2-10 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i(1),j(1),k(1)                     ! cell indices
C
C  normalized displacement from (x(i),y(j),z(k)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C    zparam=0 @z(k)  zparam=1 @z(k+1)
C
      real xparam(1),yparam(1),zparam(1)
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real hx(1),hy(1),hz(1)
      real hxi(1),hyi(1),hzi(1)
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C  0 .le. zparam .le. 1
C
C---------------------------------------------------------------------
C  use lookup routine as in Hermite interpolation
C
      i(1)=0
      j(1)=0
      k(1)=0
      call herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz,
     >     i(1),j(1),k(1),xparam(1),yparam(1),zparam(1),
     >     hx(1),hxi(1),hy(1),hyi(1),hz(1),hzi(1),ier)
      if(ier.ne.0) return
c
      call fvtricub(ict,1,1,
     >   fval,i,j,k,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   f,inf2,inf3,nz)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine fvtricub(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
C
C  use mktricub to set up spline coefficients...
C
      integer ict(10)                   ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec),kk(ivec) ! target cells (i,j,k)
      real xparam(ivec),yparam(ivec),zparam(ivec)
                          ! normalized displacements from (i,j,k) corners
C
      real hx(ivec),hy(ivec),hz(ivec)   ! grid spacing, and
      real hxi(ivec),hyi(ivec),hzi(ivec) ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
C
      real fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "evtricub")
C
      real fval(ivecd,*)               ! output returned
C
C  for detailed description of fin, ict and fval see subroutine evtricub
C  comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj,kk and parameter inputs
C     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C
      integer v
C
      real sum
      real, parameter :: sixth = 0.166666666666666667
C
C---------------
C
      z36th=sixth*sixth
      z216th=sixth*sixth*sixth
C
      iadr=0
      if(abs(ict(1)).le.2) then
C
C  0, 1st, 2nd derivatives...
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value...
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >              xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >              xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >              xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >              xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C     
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*(
     >            zpi*(
     >              cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz2*(
     >            czi*(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz2*(
     >            czi*(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy2*hz2*(
     >            czi*(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))
     >              +(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))
     >              +(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy2*(
     >            zpi*(
     >              -(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))
     >              +(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))
     >              +(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hz2*(
     >            czi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*(
     >            zpi*(
     >              cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz2*(
     >            czi*(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy2*hz2*(
     >            czi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy2*hz2*(
     >            czi*(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C     
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(0,i,j,k)  +fin(0,i,j+1,k))+
     >              xp*(-fin(0,i+1,j,k)+fin(0,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(0,i,j,k+1)  +fin(0,i,j+1,k+1))+
     >              xp*(-fin(0,i+1,j,k+1)+fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi(v)*hx2*(
     >            zpi*(
     >              cxi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              cx*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              cx*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*(
     >            zpi*(
     >              xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >              xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >              xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi(v)*hz2*(
     >            czi*(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy(v)*(
     >            zpi*(
     >              cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi(v)*hx2*hz2*(
     >            czi*(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy(v)*hz2*(
     >            czi*(
     >              xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy(v)*hz2*(
     >            czi*(
     >              cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >              xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >              xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi(v)*(
     >            -(
     >              cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi(v)*(
     >            -(
     >              xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >              xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >              xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*(
     >            czdi*(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*hzi(v)*(
     >            -(
     >              cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz(v)*(
     >            czdi*(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz(v)*(
     >            czdi*(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy2*hz(v)*(
     >            czdi*(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dx2:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz2*(
     >            czi*(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dy2:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+
     >              xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+
     >              xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz2*(
     >            czi*(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dz2:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*(
     >            zpi*(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hxi(v)*hyi(v)*(
     >            zpi*(
     >              (fin(0,i,j,k)  -fin(0,i,j+1,k))-
     >              (fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >            +zp*(
     >              (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))-
     >              (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi(v)*hx(v)*(
     >            zpi*(
     >              cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy(v)*(
     >            zpi*(
     >              -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))
     >              +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))
     >              +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hyi(v)*hz2*(
     >            czi*(
     >              (fin(3,i,j,k)  -fin(3,i,j+1,k))-
     >              (fin(3,i+1,j,k)-fin(3,i+1,j+1,k)))
     >            +cz*(
     >              (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))-
     >              (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy(v)*(
     >            zpi*(
     >              cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi(v)*hx(v)*hz2*(
     >            czi*(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy(v)*hz2*(
     >            czi*(
     >              -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy(v)*hz2*(
     >            czi*(
     >              cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(9).eq.1) then
C
C  d2f/dxdz:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hxi(v)*hzi(v)*(
     >            (
     >              (ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) -
     >              (ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            -(
     >              (ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) -
     >              (ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi(v)*(
     >            -(
     >              cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy2*hzi(v)*(
     >            (
     >              (cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) -
     >              (cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            -(
     >              (cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) -
     >              (cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hz(v)*(
     >            czdi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*hzi(v)*(
     >            -(
     >              cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz(v)*(
     >            czdi*(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy2*hz(v)*(
     >            czdi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy2*hz(v)*(
     >            czdi*(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(10).eq.1) then
C
C  d2f/dydz:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hyi(v)*hzi(v)*(
     >            (
     >              xpi*(fin(0,i,j,k)  -fin(0,i,j+1,k))+
     >              xp*(fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >            -(
     >              xpi*(fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))+
     >              xp*(fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi(v)*hx2*hzi(v)*(
     >            (
     >              cxi*(fin(1,i,j,k)  -fin(1,i,j+1,k))+
     >              cx*(fin(1,i+1,j,k)-fin(1,i+1,j+1,k)))
     >            -(
     >              cxi*(fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+
     >              cx*(fin(1,i+1,j,k+1)-fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hzi(v)*(
     >            -(
     >              xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >              xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >              xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi(v)*hz(v)*(
     >            czdi*(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy(v)*hzi(v)*(
     >            -(
     >              cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi(v)*hx2*hz(v)*(
     >            czdi*(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy(v)*hz(v)*(
     >            czdi*(
     >              xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy(v)*hz(v)*(
     >            czdi*(
     >              cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  3rd derivatives (.le.2 in each coordinate)
C
      else if(ict(1).eq.3) then
         if(ict(2).eq.1) then
C                               ! fxxy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*(
     >            zpi*(
     >              xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi(v)*(
     >            czi*(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy(v)*hz2*(
     >            czi*(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi(v)*(
     >            -(
     >              xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*(
     >            czdi*(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz(v)*(
     >            czdi*(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))
     >              +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))
     >              +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*(
     >            czi*(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz2*(
     >            czi*(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hxi(v)*hyi(v)*hzi(v)*(
     >            -(
     >              (fin(0,i,j,k)  -fin(0,i,j+1,k))-
     >              (fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >            +(
     >              (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))-
     >              (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi(v)*hx(v)*hzi(v)*(
     >            -(
     >              cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy(v)*hzi(v)*(
     >            -(
     >              -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))
     >              +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))
     >              +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hyi(v)*hz(v)*(
     >            czdi*(
     >              (fin(3,i,j,k)  -fin(3,i,j+1,k))-
     >              (fin(3,i+1,j,k)-fin(3,i+1,j+1,k)))
     >            +czd*(
     >              (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))-
     >              (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy(v)*hzi(v)*(
     >            -(
     >              cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi(v)*hx(v)*hz(v)*(
     >            czdi*(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy(v)*hz(v)*(
     >            czdi*(
     >              -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy(v)*hz(v)*(
     >            czdi*(
     >              cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)

C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >            zpi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*(
     >            zpi*(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+
     >              xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+
     >              xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi(v)*(
     >            -(
     >              cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*(
     >            czdi*(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz(v)*(
     >            czdi*(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi(v)*(
     >            zpi*(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*(
     >            zpi*(
     >              xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy(v)*(
     >            zpi*(
     >              cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  3rd derivatives (3 in each coordinate)
C
      else if(ict(1).eq.-3) then
         if(ict(2).eq.1) then
C                               ! fxxx
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))
     >              +(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))
     >              +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >            zpi*(
     >              -(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))
     >              +(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))
     >              +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*(
     >            czi*(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz2*hxi(v)*(
     >            czi*(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fyyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+
     >              xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+
     >              xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi(v)*(
     >            zpi*(
     >              cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi(v)*(
     >            czi*(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz2*hyi(v)*(
     >            czi*(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi(v)*(
     >            -(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi(v)*(
     >            -(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*hzi(v)*(
     >            -(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  4th derivatives (.le.2 in each coordinate)
C
      else if(ict(1).eq.4) then
         if(ict(2).eq.1) then
C                               ! fxxyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+
     >              xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >              xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >              xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hyi(v)*hzi(v)*(
     >            -(
     >              xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hzi(v)*(
     >            -(
     >              xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hyi(v)*(
     >            czdi*(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy(v)*hz(v)*(
     >            czdi*(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hxi(v)*hzi(v)*(
     >            -(
     >              -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))
     >              +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))
     >              +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi(v)*(
     >            -(
     >              cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hxi(v)*(
     >            czdi*(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz(v)*(
     >            czdi*(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hyi(v)*hxi(v)*(
     >            zpi*(
     >               ( +fin(3,i,j,k)  -fin(3,i,j+1,k))
     >              +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +zp*(
     >               ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))
     >              +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi(v)*(
     >            zpi*(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hxi(v)*(
     >            zpi*(
     >              -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy(v)*(
     >            zpi*(
     >              cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  4th derivatives (3 in a coordinate)
C
      else if(ict(1).eq.-4) then
         if(ict(2).eq.1) then
C                               ! fxxxy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hyi(v)*hxi(v)*(
     >            zpi*(
     >              (  fin(1,i,j,k)  -fin(1,i,j+1,k))+
     >              ( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              (  fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+
     >              ( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hxi(v)*(
     >            zpi*(
     >              -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >               (cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >               (cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi(v)*hxi(v)*(
     >            czi*(
     >              (  fin(5,i,j,k)  -fin(5,i,j+1,k))+
     >              ( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              (  fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))+
     >              ( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy(v)*hz2*hxi(v)*(
     >            czi*(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >               (cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >               (cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hzi(v)*hxi(v)*(
     >             (
     >              +(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))
     >              -(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))
     >              +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi(v)*hxi(v)*(
     >             (
     >              +(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))
     >              -(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))
     >              +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hxi(v)*(
     >            czdi*(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz(v)*hxi(v)*(
     >            czdi*(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hxi(v)*hyi(v)*(
     >            zpi*(
     >               ( fin(2,i,j,k)  -fin(2,i,j+1,k))
     >              +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +zp*(
     >               ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1))
     >              +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi(v)*(
     >            zpi*(
     >              cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cxd*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cxd*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*hyi(v)*(
     >            czi*(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +cz*(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz2*hyi(v)*(
     >            czi*(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
               sum=hxi(v)*hzi(v)*(
     >            -(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi(v)*(
     >            -(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*hzi(v)*(
     >            -(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*hzi(v)*(
     >            -(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fyyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hyi(v)*hzi(v)*(
     >            -(
     >              xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+
     >              xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+
     >              xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi(v)*hzi(v)*(
     >            -(
     >              cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hyi(v)*(
     >            czdi*(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz(v)*hyi(v)*(
     >            czdi*(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
               sum=hyi(v)*hzi(v)*(
     >            -(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi(v)*hzi(v)*(
     >            -(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hzi(v)*(
     >            -(
     >              xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy(v)*hzi(v)*(
     >            -(
     >              cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  5th derivatives (.le.2 in each coordinate)
C
      else if(ict(1).eq.5) then
         if(ict(2).eq.1) then
C                               ! fxxyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+
     >              xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*(
     >            czdi*(
     >              xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >              xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >              xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*(
     >            zpi*(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  5th derivatives (3 in a coordinate)
C
      else if(ict(1).eq.-5) then
         if(ict(2).eq.1) then
C                               ! fxxxyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))
     >              +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))
     >              +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*(
     >            czi*(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hyi(v)*hzi(v)*hxi(v)*(
     >            -(
     >              -(-fin(1,i,j,k)  +fin(1,i,j+1,k))
     >              +( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +(
     >              -(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))
     >              +( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hzi(v)*hxi(v)*(
     >            -(
     >              -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))
     >              +(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))
     >              +(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hyi(v)*hxi(v)*(
     >            czdi*(
     >              -(-fin(5,i,j,k)  +fin(5,i,j+1,k))
     >              +( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              -(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))
     >              +( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy(v)*hz(v)*hxi(v)*(
     >            czdi*(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))
     >              +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))
     >              +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >            zpi*(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxyyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)

C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              xp*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              xp*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi(v)*(
     >            czi*(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi(v)*(
     >            -(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxyyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hxi(v)*hzi(v)*hyi(v)*(
     >            -(
     >               ( fin(2,i,j,k)  -fin(2,i,j+1,k))
     >              +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +(
     >               ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1))
     >              +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi(v)*hyi(v)*(
     >            -(
     >              cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cxd*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cxd*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hxi(v)*hyi(v)*(
     >            czdi*(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +czd*(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz(v)*hyi(v)*(
     >            czdi*(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fxyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
               sum=hyi(v)*hxi(v)*hzi(v)*(
     >            -(
     >               ( +fin(3,i,j,k)  -fin(3,i,j+1,k))
     >              +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +(
     >               ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))
     >              +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi(v)*hzi(v)*(
     >            -(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hxi(v)*hzi(v)*(
     >            -(
     >              -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy(v)*hzi(v)*(
     >            -(
     >              cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(9).eq.1) then
C                               ! fyyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi(v)*(
     >            zpi*(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(10).eq.1) then
C                               ! fyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
               sum=hzi(v)*(
     >            -(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi(v)*(
     >            -(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  6th derivatives (2 in each coordinate)
C
      else if(ict(1).eq.6) then
C                               ! fxxyyzz
         iadr=iadr+1
         do v=1,ivec
            i=ii(v)
            j=jj(v)
            k=kk(v)
C
C   ...in x direction
C
            xp=xparam(v)
            xpi=1.0-xp
C
C   ...and in y direction
C
            yp=yparam(v)
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam(v)
            zpi=1.0-zp
C
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >           xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >         +zp*(
     >           xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >           xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         enddo
      endif
C
C----------------------------------
C  6th derivatives (3 in a coordinate)
C
      if(ict(1).eq.-6) then
         if(ict(2).eq.1) then
C                               ! fxxxyyy
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz(v)*hz(v)
C
               sum=hyi(v)*hxi(v)*(
     >            zpi*(
     >              ( fin(4,i,j,k)  -fin(4,i,j+1,k))
     >             +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1))
     >             +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi(v)*hxi(v)*(
     >            czi*(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hxi(v)*hzi(v)*(
     >            -(
     >              -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))
     >              +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))
     >              +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hxi(v)*(
     >            czdi*(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hxi(v)*hyi(v)*(
     >            zpi*(
     >               ( fin(5,i,j,k)  -fin(5,i,j+1,k))
     >              +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >               ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))
     >              +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hxi(v)*(
     >            zpi*(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))
     >              +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))
     >              +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxxzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy(v)*hy(v)
C
               sum=hxi(v)*hzi(v)*(
     >            -(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*hzi(v)*(
     >            -(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxyyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hzi(v)*hyi(v)*(
     >            -(
     >              xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              xp*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              xp*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hyi(v)*(
     >            czdi*(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxxyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
               sum=hyi(v)*hzi(v)*(
     >            -(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hzi(v)*(
     >            -(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fxyyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hxi(v)*hyi(v)*(
     >            zpi*(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +zp*(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi(v)*(
     >            zpi*(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(9).eq.1) then
C                               ! fxyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
               sum=hxi(v)*hzi(v)*(
     >            -(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi(v)*(
     >            -(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(10).eq.1) then
C                               ! fyyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi(v)*hzi(v)*(
     >            -(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi(v)*hzi(v)*(
     >            -(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  7th derivatives
C
      else if(abs(ict(1)).eq.7) then
         if(ict(2).eq.1) then
C                               ! fxxxyyyz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi

               czd=3.0*zp2-1.0
               czdi=-3.0*zpi2+1.0
C
               sum=hyi(v)*hxi(v)*hzi(v)*(
     >            -(
     >              ( fin(4,i,j,k)  -fin(4,i,j+1,k))
     >             +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +(
     >              ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1))
     >             +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz(v)*hyi(v)*hxi(v)*(
     >            czdi*(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +czd*(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi

               cyd=3.0*yp2-1.0
               cydi=-3.0*ypi2+1.0
C
               sum=hxi(v)*hyi(v)*hzi(v)*(
     >            -(
     >               ( fin(5,i,j,k)  -fin(5,i,j+1,k))
     >              +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >               ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))
     >              +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy(v)*hxi(v)*hzi(v)*(
     >            -(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))
     >              +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))
     >              +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxyyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hyi(v)*(
     >            zpi*(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >           +zp*(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
               sum=hzi(v)*(
     >           -(
     >              xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >              xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >           +(
     >              xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >              xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxyyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hyi(v)*hzi(v)*(
     >            -(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi(v)*hzi(v)*(
     >            -(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  8th derivatives
C
      else if(abs(ict(1)).eq.8) then
         if(ict(2).eq.1) then
C                               ! fxxxyyyzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in z direction
C
               zp=zparam(v)
               zpi=1.0-zp
C
               sum=hyi(v)*hxi(v)*(
     >            zpi*(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +zp*(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...and in y direction
C
               yp=yparam(v)
               ypi=1.0-yp
C
               sum=hxi(v)*hzi(v)*(
     >            -(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxyyyzzz
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj(v)
               k=kk(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi(v)*hzi(v)*(
     >           -(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >           +(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
               fval(v,iadr)=sum
C
            enddo
         endif
C
C----------------------------------
C  9th derivative
C
      else if(abs(ict(1)).eq.9) then
C                               ! fxxxyyyzzz
         iadr=iadr+1
         do v=1,ivec
            i=ii(v)
            j=jj(v)
            k=kk(v)
C
            sum=hyi(v)*hxi(v)*hzi(v)*(
     >            -(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         enddo
      endif
C
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C    --optimized for VARIATION along x axis ONLY--
C
      subroutine fvtricubx(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
C
C  use mktricub to set up spline coefficients...
C
      integer ict(10)                   ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj,kk            ! target cells (i,j,k)
      real xparam(ivec),yparam,zparam
                          ! normalized displacements from (i,j,k) corners
C
      real hx(ivec),hy,hz               ! grid spacing, and
      real hxi(ivec),hyi,hzi            ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
C
      real fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "evtricub")
C
      real fval(ivecd,*)               ! output returned
C
C  for detailed description of fin, ict and fval see subroutine evtricub
C  comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj,kk and parameter inputs
C     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C
      integer v
C
      real sum
      real, parameter :: sixth = 0.166666666666666667
C
C---------------
C
      z36th=sixth*sixth
      z216th=sixth*sixth*sixth
C
      iadr=0
      if(abs(ict(1)).le.2) then
C
C  0, 1st, 2nd derivatives...
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value...
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >              xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >              xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >              xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >              xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C     
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*(
     >            zpi*(
     >              cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz2*(
     >            czi*(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz2*(
     >            czi*(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy2*hz2*(
     >            czi*(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))
     >              +(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))
     >              +(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy2*(
     >            zpi*(
     >              -(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))
     >              +(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))
     >              +(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hz2*(
     >            czi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*(
     >            zpi*(
     >              cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz2*(
     >            czi*(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy2*hz2*(
     >            czi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy2*hz2*(
     >            czi*(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C     
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(0,i,j,k)  +fin(0,i,j+1,k))+
     >              xp*(-fin(0,i+1,j,k)+fin(0,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(0,i,j,k+1)  +fin(0,i,j+1,k+1))+
     >              xp*(-fin(0,i+1,j,k+1)+fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi*hx2*(
     >            zpi*(
     >              cxi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              cx*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              cx*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*(
     >            zpi*(
     >              xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >              xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >              xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi*hz2*(
     >            czi*(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy*(
     >            zpi*(
     >              cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi*hx2*hz2*(
     >            czi*(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy*hz2*(
     >            czi*(
     >              xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy*hz2*(
     >            czi*(
     >              cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >              xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >              xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi*(
     >            -(
     >              cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi*(
     >            -(
     >              xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >              xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >              xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*(
     >            czdi*(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*hzi*(
     >            -(
     >              cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz*(
     >            czdi*(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz*(
     >            czdi*(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy2*hz*(
     >            czdi*(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dx2:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz2*(
     >            czi*(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dy2:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+
     >              xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+
     >              xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz2*(
     >            czi*(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dz2:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp

            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*(
     >            zpi*(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C
C  d2f/dxdy:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hyi*(
     >            zpi*(
     >              (fin(0,i,j,k)  -fin(0,i,j+1,k))-
     >              (fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >            +zp*(
     >              (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))-
     >              (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi*hx(v)*(
     >            zpi*(
     >              cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy*(
     >            zpi*(
     >              -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))
     >              +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))
     >              +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hyi*hz2*(
     >            czi*(
     >              (fin(3,i,j,k)  -fin(3,i,j+1,k))-
     >              (fin(3,i+1,j,k)-fin(3,i+1,j+1,k)))
     >            +cz*(
     >              (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))-
     >              (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy*(
     >            zpi*(
     >              cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi*hx(v)*hz2*(
     >            czi*(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy*hz2*(
     >            czi*(
     >              -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy*hz2*(
     >            czi*(
     >              cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(9).eq.1) then
C
C  d2f/dxdz:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hzi*(
     >            (
     >              (ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) -
     >              (ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >            -(
     >              (ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) -
     >              (ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi*(
     >            -(
     >              cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy2*hzi*(
     >            (
     >              (cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) -
     >              (cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >            -(
     >              (cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) -
     >              (cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hz*(
     >            czdi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*hzi*(
     >            -(
     >              cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz*(
     >            czdi*(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy2*hz*(
     >            czdi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy2*hz*(
     >            czdi*(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(10).eq.1) then
C
C  d2f/dydz:
C
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*hzi*(
     >            (
     >              xpi*(fin(0,i,j,k)  -fin(0,i,j+1,k))+
     >              xp*(fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >            -(
     >              xpi*(fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))+
     >              xp*(fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi*hx2*hzi*(
     >            (
     >              cxi*(fin(1,i,j,k)  -fin(1,i,j+1,k))+
     >              cx*(fin(1,i+1,j,k)-fin(1,i+1,j+1,k)))
     >            -(
     >              cxi*(fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+
     >              cx*(fin(1,i+1,j,k+1)-fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hzi*(
     >            -(
     >              xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >              xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >              xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi*hz*(
     >            czdi*(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy*hzi*(
     >            -(
     >              cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi*hx2*hz*(
     >            czdi*(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy*hz*(
     >            czdi*(
     >              xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx2*hy*hz*(
     >            czdi*(
     >              cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  3rd derivatives (.le.2 in each coordinate)
C
      else if(ict(1).eq.3) then
         if(ict(2).eq.1) then
C                               ! fxxy
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*(
     >            zpi*(
     >              xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi*(
     >            czi*(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy*hz2*(
     >            czi*(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >              xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >              xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi*(
     >            -(
     >              xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >              xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >              xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*(
     >            czdi*(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz*(
     >            czdi*(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyy
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))
     >              +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))
     >              +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*(
     >            czi*(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz2*(
     >            czi*(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxyz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hyi*hzi*(
     >            -(
     >              (fin(0,i,j,k)  -fin(0,i,j+1,k))-
     >              (fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >            +(
     >              (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))-
     >              (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
               sum=sum+sixth*hyi*hx(v)*hzi*(
     >            -(
     >              cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hy*hzi*(
     >            -(
     >              -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))
     >              +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))
     >              +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hxi(v)*hyi*hz*(
     >            czdi*(
     >              (fin(3,i,j,k)  -fin(3,i,j+1,k))-
     >              (fin(3,i+1,j,k)-fin(3,i+1,j+1,k)))
     >            +czd*(
     >              (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))-
     >              (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy*hzi*(
     >            -(
     >              cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+z36th*hyi*hx(v)*hz*(
     >            czdi*(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hxi(v)*hy*hz*(
     >            czdi*(
     >              -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z216th*hx(v)*hy*hz*(
     >            czdi*(
     >              cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy

C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >            zpi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*(
     >            zpi*(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyyz
            j=jj
            k=kk
C
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+
     >              xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+
     >              xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi*(
     >            -(
     >              cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*(
     >            czdi*(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz*(
     >            czdi*(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi*(
     >            zpi*(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*(
     >            zpi*(
     >              xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy*(
     >            zpi*(
     >              cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  3rd derivatives (3 in each coordinate)
C
      else if(ict(1).eq.-3) then
         if(ict(2).eq.1) then
C                               ! fxxx
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
               j=jj
               k=kk
C
C   ...and in y direction
C
               yp=yparam
               ypi=1.0-yp
               yp2=yp*yp
               ypi2=ypi*ypi
C
               cy=yp*(yp2-1.0)
               cyi=ypi*(ypi2-1.0)
               hy2=hy*hy
C
C   ...and in z direction
C
               zp=zparam
               zpi=1.0-zp
               zp2=zp*zp
               zpi2=zpi*zpi
C
               cz=zp*(zp2-1.0)
               czi=zpi*(zpi2-1.0)
               hz2=hz*hz
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))
     >              +(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))
     >              +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >            zpi*(
     >              -(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))
     >              +(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))
     >              +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*(
     >            czi*(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz2*hxi(v)*(
     >            czi*(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fyyy
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+
     >              xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+
     >              xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi*(
     >            zpi*(
     >              cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi*(
     >            czi*(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz2*hyi*(
     >            czi*(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >              xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >              xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi*(
     >            -(
     >              cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi*(
     >            -(
     >              xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >              xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >              xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy2*hzi*(
     >            -(
     >              cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  4th derivatives (.le.2 in each coordinate)
C
      else if(ict(1).eq.4) then
         if(ict(2).eq.1) then
C                               ! fxxyy
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+
     >              xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*(
     >            czi*(
     >              xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >              xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >              xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxyz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*hzi*(
     >            -(
     >              xpi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >              xp*( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >              xp*( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hzi*(
     >            -(
     >              xpi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >              xp*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >              xp*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hyi*(
     >            czdi*(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy*hz*(
     >            czdi*(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*(
     >            zpi*(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxyyz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hzi*(
     >            -(
     >              -(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))
     >              +(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))
     >              +(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi*(
     >            -(
     >              cxdi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              cxd*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >              cxd*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hxi(v)*(
     >            czdi*(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz*(
     >            czdi*(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hyi*hxi(v)*(
     >            zpi*(
     >               ( +fin(3,i,j,k)  -fin(3,i,j+1,k))
     >              +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +zp*(
     >               ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))
     >              +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi*(
     >            zpi*(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hxi(v)*(
     >            zpi*(
     >              -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy*(
     >            zpi*(
     >              cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=(
     >            zpi*(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*(
     >            zpi*(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  4th derivatives (3 in a coordinate)
C
      else if(ict(1).eq.-4) then
         if(ict(2).eq.1) then
C                               ! fxxxy
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hyi*hxi(v)*(
     >            zpi*(
     >              (  fin(1,i,j,k)  -fin(1,i,j+1,k))+
     >              ( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +zp*(
     >              (  fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+
     >              ( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hxi(v)*(
     >            zpi*(
     >              -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >               (cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >               (cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi*hxi(v)*(
     >            czi*(
     >              (  fin(5,i,j,k)  -fin(5,i,j+1,k))+
     >              ( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +cz*(
     >              (  fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))+
     >              ( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy*hz2*hxi(v)*(
     >            czi*(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >               (cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >               (cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hzi*hxi(v)*(
     >             (
     >              +(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))
     >              -(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))
     >              +(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi*hxi(v)*(
     >             (
     >              +(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))
     >              -(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >            +(
     >              -(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))
     >              +(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hxi(v)*(
     >            czdi*(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy2*hz*hxi(v)*(
     >            czdi*(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyyy
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hyi*(
     >            zpi*(
     >               ( fin(2,i,j,k)  -fin(2,i,j+1,k))
     >              +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +zp*(
     >               ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1))
     >              +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi*(
     >            zpi*(
     >              cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cxd*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cxd*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*hyi*(
     >            czi*(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +cz*(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz2*hyi*(
     >            czi*(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hzi*(
     >            -(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi*(
     >            -(
     >              cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*hzi*(
     >            -(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >            +(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy2*hzi*(
     >            -(
     >              cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fyyyz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*hzi*(
     >            -(
     >              xpi*(-fin(2,i,j,k)  +fin(2,i,j+1,k))+
     >              xp*( -fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(2,i,j,k+1)  +fin(2,i,j+1,k+1))+
     >              xp*( -fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi*hzi*(
     >            -(
     >              cxi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cx*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +(
     >              cxi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cx*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hyi*(
     >            czdi*(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hz*hyi*(
     >            czdi*(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*hzi*(
     >            -(
     >              xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >              xp*( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >              xp*( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi*hzi*(
     >            -(
     >              cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cx*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >              cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cx*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hzi*(
     >            -(
     >              xpi*(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))+
     >              xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))+
     >              xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx2*hy*hzi*(
     >            -(
     >              cxi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  5th derivatives (.le.2 in each coordinate)
C
      else if(ict(1).eq.5) then
         if(ict(2).eq.1) then
C                               ! fxxyyz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >              xp*( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))+
     >              xp*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*(
     >            czdi*(
     >              xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >              xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >              xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*(
     >            zpi*(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*(
     >            zpi*(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  5th derivatives (3 in a coordinate)
C
      else if(ict(1).eq.-5) then
         if(ict(2).eq.1) then
C                               ! fxxxyy
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))
     >              +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))
     >              +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hxi(v)*(
     >            czi*(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +cz*(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hyi*hzi*hxi(v)*(
     >            -(
     >              -(-fin(1,i,j,k)  +fin(1,i,j+1,k))
     >              +( -fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >            +(
     >              -(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))
     >              +( -fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hzi*hxi(v)*(
     >            -(
     >              -(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))
     >              +(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))
     >              +(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hyi*hxi(v)*(
     >            czdi*(
     >              -(-fin(5,i,j,k)  +fin(5,i,j+1,k))
     >              +( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +czd*(
     >              -(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))
     >              +( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+z36th*hy*hz*hxi(v)*(
     >            czdi*(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))
     >              +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))
     >              +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*(
     >            zpi*(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxyyy
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              xp*( -fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              xp*( -fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi*(
     >            czi*(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >              xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >              xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hzi*(
     >            -(
     >              xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >              xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >              xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxyyyz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hzi*hyi*(
     >            -(
     >               ( fin(2,i,j,k)  -fin(2,i,j+1,k))
     >              +(-fin(2,i+1,j,k)+fin(2,i+1,j+1,k)))
     >            +(
     >               ( fin(2,i,j,k+1)  -fin(2,i,j+1,k+1))
     >              +(-fin(2,i+1,j,k+1)+fin(2,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi*hyi*(
     >            -(
     >              cxdi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              cxd*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              cxd*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hxi(v)*hyi*(
     >            czdi*(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +czd*(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hz*hyi*(
     >            czdi*(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +czd*(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fxyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hyi*hxi(v)*hzi*(
     >            -(
     >               ( +fin(3,i,j,k)  -fin(3,i,j+1,k))
     >              +( -fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >            +(
     >               ( +fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))
     >              +( -fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi*hzi*(
     >            -(
     >              cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              cxd*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              cxd*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hxi(v)*hzi*(
     >            -(
     >              -(cydi*fin(6,i,j,k) +cyd*fin(6,i,j+1,k))
     >              +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(6,i,j,k+1) +cyd*fin(6,i,j+1,k+1))
     >              +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
               sum=sum+z36th*hx(v)*hy*hzi*(
     >            -(
     >              cxdi*(cydi*fin(7,i,j,k) +cyd*fin(7,i,j+1,k))+
     >              cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(cydi*fin(7,i,j,k+1) +cyd*fin(7,i,j+1,k+1))+
     >              cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(9).eq.1) then
C                               ! fyyyzz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +zp*(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi*(
     >            zpi*(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(10).eq.1) then
C                               ! fyyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hzi*(
     >            -(
     >              xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >              xp*( ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(ypi*fin(6,i,j,k+1) +yp*fin(6,i,j+1,k+1))+
     >              xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hzi*(
     >            -(
     >              cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cx*( ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cx*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  6th derivatives (2 in each coordinate)
C
      else if(ict(1).eq.6) then
C                               ! fxxyyzz
         j=jj
         k=kk
C
C   ...and in y direction
C
         yp=yparam
         ypi=1.0-yp
C
C   ...and in z direction
C
         zp=zparam
         zpi=1.0-zp
C
         iadr=iadr+1
         do v=1,ivec
            i=ii(v)
C
C   ...in x direction
C
            xp=xparam(v)
            xpi=1.0-xp
C
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >           xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >         +zp*(
     >           xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >           xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         enddo
      endif
C
C----------------------------------
C  6th derivatives (3 in a coordinate)
C
      if(ict(1).eq.-6) then
         if(ict(2).eq.1) then
C                               ! fxxxyyy
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi
C
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz*hz
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hyi*hxi(v)*(
     >            zpi*(
     >              ( fin(4,i,j,k)  -fin(4,i,j+1,k))
     >             +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +zp*(
     >              ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1))
     >             +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz2*hyi*hxi(v)*(
     >            czi*(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +cz*(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*hzi*(
     >            -(
     >              -(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))
     >              +( ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(4,i,j,k+1) +yp*fin(4,i,j+1,k+1))
     >              +(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hxi(v)*(
     >            czdi*(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +czd*(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*hyi*(
     >            zpi*(
     >               ( fin(5,i,j,k)  -fin(5,i,j+1,k))
     >              +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +zp*(
     >               ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))
     >              +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hxi(v)*(
     >            zpi*(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))
     >              +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))
     >              +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxxzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi
C
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy*hy
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*hzi*(
     >            -(
     >              -(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))
     >              +(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))
     >              +(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy2*hxi(v)*hzi*(
     >            -(
     >              -(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))
     >              +(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >            +(
     >              -(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))
     >              +(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxyyyz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hzi*hyi*(
     >            -(
     >              xpi*(-fin(4,i,j,k)  +fin(4,i,j+1,k))+
     >              xp*(-fin(4,i+1,j,k) +fin(4,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(4,i,j,k+1)  +fin(4,i,j+1,k+1))+
     >              xp*(-fin(4,i+1,j,k+1) +fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hyi*(
     >            czdi*(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +czd*(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxxyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*hzi*(
     >            -(
     >              xpi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >              xp*( -fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >              xp*( -fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hzi*(
     >            -(
     >              xpi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >              xp*( cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              xpi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >              xp*( cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fxyyyzz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hyi*(
     >            zpi*(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +zp*(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi*(
     >            zpi*(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +zp*(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(9).eq.1) then
C                               ! fxyyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hzi*(
     >            -(
     >              -(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))
     >              +(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))
     >              +(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hzi*(
     >            -(
     >              cxdi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >              cxd*(ypi*fin(7,i+1,j,k) +yp*fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >              cxd*( ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(10).eq.1) then
C                               ! fyyyzzz
            j=jj
            k=kk
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0)
               cxi=xpi*(xpi2-1.0)
               hx2=hx(v)*hx(v)
C
               sum=hyi*hzi*(
     >            -(
     >              xpi*(-fin(6,i,j,k)  +fin(6,i,j+1,k))+
     >              xp*( -fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +(
     >              xpi*(-fin(6,i,j,k+1)  +fin(6,i,j+1,k+1))+
     >              xp*( -fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx2*hyi*hzi*(
     >            -(
     >              cxi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cx*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +(
     >              cxi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cx*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  7th derivatives
C
      else if(abs(ict(1)).eq.7) then
         if(ict(2).eq.1) then
C                               ! fxxxyyyz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
            zp2=zp*zp
            zpi2=zpi*zpi

            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hyi*hxi(v)*hzi*(
     >            -(
     >              ( fin(4,i,j,k)  -fin(4,i,j+1,k))
     >             +(-fin(4,i+1,j,k)+fin(4,i+1,j+1,k)))
     >            +(
     >              ( fin(4,i,j,k+1)  -fin(4,i,j+1,k+1))
     >             +(-fin(4,i+1,j,k+1)+fin(4,i+1,j+1,k+1))))
C
               sum=sum+sixth*hz*hyi*hxi(v)*(
     >            czdi*(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +czd*(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*(
     >            zpi*(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +zp*(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
            yp2=yp*yp
            ypi2=ypi*ypi

            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*hyi*hzi*(
     >            -(
     >               ( fin(5,i,j,k)  -fin(5,i,j+1,k))
     >              +(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >            +(
     >               ( fin(5,i,j,k+1)  -fin(5,i,j+1,k+1))
     >              +(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
               sum=sum+sixth*hy*hxi(v)*hzi*(
     >            -(
     >              -(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))
     >              +(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >            +(
     >              -(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))
     >              +(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxyyyzz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*(
     >            zpi*(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >           +zp*(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxyyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hzi*(
     >           -(
     >              xpi*(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))+
     >              xp*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >           +(
     >              xpi*(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))+
     >              xp*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxyyyzzz
            j=jj
            k=kk
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
               xp2=xp*xp
               xpi2=xpi*xpi

               cxd=3.0*xp2-1.0
               cxdi=-3.0*xpi2+1.0
C
               sum=hxi(v)*hyi*hzi*(
     >            -(
     >               ( fin(6,i,j,k)  -fin(6,i,j+1,k))
     >              +(-fin(6,i+1,j,k)+fin(6,i+1,j+1,k)))
     >            +(
     >               ( fin(6,i,j,k+1)  -fin(6,i,j+1,k+1))
     >              +(-fin(6,i+1,j,k+1)+fin(6,i+1,j+1,k+1))))
C
               sum=sum+sixth*hx(v)*hyi*hzi*(
     >            -(
     >              cxdi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              cxd*(-fin(7,i+1,j,k) +fin(7,i+1,j+1,k)))
     >            +(
     >              cxdi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              cxd*(-fin(7,i+1,j,k+1) +fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
C----------------------------------
C  8th derivatives
C
      else if(abs(ict(1)).eq.8) then
         if(ict(2).eq.1) then
C                               ! fxxxyyyzz
            j=jj
            k=kk
C
C   ...and in z direction
C
            zp=zparam
            zpi=1.0-zp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hyi*hxi(v)*(
     >            zpi*(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +zp*(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyzzz
            j=jj
            k=kk
C
C   ...and in y direction
C
            yp=yparam
            ypi=1.0-yp
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               sum=hxi(v)*hzi*(
     >            -(
     >              -(ypi*fin(7,i,j,k) +yp*fin(7,i,j+1,k))
     >              +(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(7,i,j,k+1) +yp*fin(7,i,j+1,k+1))
     >              +(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxyyyzzz
            j=jj
            k=kk
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
C   ...in x direction
C
               xp=xparam(v)
               xpi=1.0-xp
C
               sum=hyi*hzi*(
     >           -(
     >              xpi*(-fin(7,i,j,k)  +fin(7,i,j+1,k))+
     >              xp*( -fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >           +(
     >              xpi*(-fin(7,i,j,k+1)  +fin(7,i,j+1,k+1))+
     >              xp*( -fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
               fval(v,iadr)=sum
C
            enddo
         endif
C
C----------------------------------
C  9th derivative
C
      else if(abs(ict(1)).eq.9) then
C                               ! fxxxyyyzzz
         j=jj
         k=kk
C
         iadr=iadr+1
         do v=1,ivec
            i=ii(v)
C
            sum=hyi*hxi(v)*hzi*(
     >            -(
     >              ( fin(7,i,j,k)  -fin(7,i,j+1,k))
     >             +(-fin(7,i+1,j,k)+fin(7,i+1,j+1,k)))
     >            +(
     >              ( fin(7,i,j,k+1)  -fin(7,i,j+1,k+1))
     >             +(-fin(7,i+1,j,k+1)+fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         enddo
      endif
C
      return
      end
