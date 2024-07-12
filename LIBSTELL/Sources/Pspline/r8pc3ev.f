      subroutine r8pc3ev(xget,yget,zget,x,nx,y,ny,z,nz,
     >                   ilinx,iliny,ilinz,
     >                   f,inf2,inf3,ict,fval,ier)
C
C  evaluate a trilinear interpolant on a rectilinear grid
C  derivatives are available, but, not continuous across grid planes.
C
C  this subroutine calls two subroutines:
C     herm3xyz  -- find cell containing (xget,yget,zget)
C     pc3fcn -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ny,nz,inf2,inf3,nx
!============
      REAL*8 xget,yget,zget               ! target of this interpolation
      REAL*8 x(nx)                        ! ordered x grid
      REAL*8 y(ny)                        ! ordered y grid
      REAL*8 z(nz)                        ! ordered z grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
      integer ilinz                     ! ilinz=1 => assume z evenly spaced
C
      REAL*8 f(inf2,inf3,nz)              ! function data
C
C       f 2nd dimension inf2 must be .ge. nx; 3rd dim inf3 .ge. ny
C       contents of f:
C
C  f(i,j,k) = f @ x(i),y(j),z(k)
C
      integer ict(8)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return df/dz  (0, don't)
C  ict(5)=1 -- return d2f/dxdy  (0, don't)
C  ict(6)=1 -- return d2f/dxdz  (0, don't)
C  ict(7)=1 -- return d2f/dydz  (0, don't)
C  ict(8)=1 -- return d3f/dxdydz  (0, don't)
C
C output arguments:
C =================
C
      REAL*8 fval(*)                      ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the 4th output (depends on ict(...) spec)
C  fval(5-8) receive 5th thru 8th outputs (if required by ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,1,0,0,0,0]
C   on output fval= [f,df/dx,df/dy,df/dz]
C
C    on input ict = [1,0,0,0,0,0,0,0]
C   on output fval= [f] ... elements 2-8 never referenced
C
C    on input ict = [0,1,1,0,0,0,0,0]
C   on output fval= [df/dx,df/dy] ... elements 3-8 never referenced
C
C    on input ict = [0,0,0,0,1,0,0,0]
C   on output fval= [d2f/dxdy] ... elements 2-8 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i,j,k                     ! cell indices
C
C  normalized displacement from (x(i),y(j)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C    zparam=0 @z(k)  zparam=1 @z(k+1)
C
      REAL*8 xparam,yparam,zparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      REAL*8 hx,hy,hz
      REAL*8 hxi,hyi,hzi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C  0 .le. zparam .le. 1
C
C---------------------------------------------------------------------
C
      call r8herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz,
     >   i,j,k,xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi,ier)
      if(ier.ne.0) return
c
      call r8pc3fcn(ict,1,1,
     >   fval,i,j,k,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   f,inf2,inf3,nz)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate trilinear function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine r8pc3fcn(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inf3,nz,inf2,i,j,k,iadr
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xp,xpi,yp,ypi,zp,zpi
!============
      integer ict(8)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec),kk(ivec) ! target cells (i,j,k)
      REAL*8 xparam(ivec),yparam(ivec),zparam(ivec)
                          ! normalized displacements from (i,j,k) corners
C
      REAL*8 hx(ivec),hy(ivec),hz(ivec)   ! grid spacing, and
      REAL*8 hxi(ivec),hyi(ivec),hzi(ivec) ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
C
      REAL*8 fin(inf2,inf3,nz)            ! interpolant data (cf "pc3ev")
C
      REAL*8 fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine pc3ev
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
      REAL*8 sum
      integer v
C
      do v=1,ivec
         i=ii(v)
         j=jj(v)
         k=kk(v)
C
C   ...in x direction
C
         xp=xparam(v)
         xpi=1.0_r8-xp
C
C   ...in y direction
C
         yp=yparam(v)
         ypi=1.0_r8-yp
C
C   ...in z direction
C
         zp=zparam(v)
         zpi=1.0_r8-zp
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
            sum=zpi*(
     >           xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+
     >            xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >          +zp*(
     >           xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+
     >            xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
            sum=zpi*(
     >              -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))
     >              +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >          +zp*(
     >              -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))
     >              +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
            sum=zpi*(
     >           xpi*(-fin(i,j,k)  +fin(i,j+1,k))+
     >            xp*(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >          +zp*(
     >           xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+
     >            xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hyi(v)
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            iadr=iadr+1
            sum=   -(
     >           xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+
     >            xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >             +(
     >           xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+
     >            xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hzi(v)
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
            sum=zpi*(
     >              -(-fin(i,j,k)  +fin(i,j+1,k))
     >              +(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >          +zp*(
     >              -(-fin(i,j,k+1)  +fin(i,j+1,k+1))
     >              +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)*hyi(v)
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dxdz:
C
            iadr=iadr+1
            sum=  -(
     >              -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))
     >              +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))
     >              +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)*hzi(v)
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dydz:
C
            iadr=iadr+1
            sum=  -(
     >           xpi*(-fin(i,j,k)  +fin(i,j+1,k))+
     >            xp*(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >            +(
     >           xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+
     >            xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hyi(v)*hzi(v)
         endif
C
         if(ict(8).eq.1) then
C
C  d3f/dxdydz:
C
            iadr=iadr+1
            sum=  -(
     >              -(-fin(i,j,k)  +fin(i,j+1,k))
     >              +(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >            +(
     >              -(-fin(i,j,k+1)  +fin(i,j+1,k+1))
     >              +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)*hyi(v)*hzi(v)
         endif
C
      enddo                             ! vector loop
C
      return
      end
C---------------------------------------------------------------------
C  evaluate trilinear function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C    --optimized for VARIATION along x axis ONLY--
C
      subroutine r8pc3fcnx(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inf3,nz,inf2,j,k,i,iadr
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 yp,ypi,zp,zpi,xp,xpi
!============
      integer ict(8)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj,kk ! target cells (i,j,k)
      REAL*8 xparam(ivec),yparam,zparam
                          ! normalized displacements from (i,j,k) corners
C
      REAL*8 hx(ivec),hy,hz               ! grid spacing, and
      REAL*8 hxi(ivec),hyi,hzi            ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
C
      REAL*8 fin(inf2,inf3,nz)            ! interpolant data (cf "pc3ev")
C
      REAL*8 fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine pc3ev
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
      REAL*8 sum
      integer v
C
      j=jj
      k=kk
C
C   ...in y direction
C
      yp=yparam
      ypi=1.0_r8-yp
C
C   ...in z direction
C
      zp=zparam
      zpi=1.0_r8-zp
 
      do v=1,ivec
         i=ii(v)
C
C   ...in x direction
C
         xp=xparam(v)
         xpi=1.0_r8-xp
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
            sum=zpi*(
     >           xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+
     >            xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >          +zp*(
     >           xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+
     >            xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
            sum=zpi*(
     >              -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))
     >              +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >          +zp*(
     >              -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))
     >              +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
            sum=zpi*(
     >           xpi*(-fin(i,j,k)  +fin(i,j+1,k))+
     >            xp*(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >          +zp*(
     >           xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+
     >            xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hyi
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            iadr=iadr+1
            sum=   -(
     >           xpi*(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))+
     >            xp*(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >             +(
     >           xpi*(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))+
     >            xp*(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hzi
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
            sum=zpi*(
     >              -(-fin(i,j,k)  +fin(i,j+1,k))
     >              +(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >          +zp*(
     >              -(-fin(i,j,k+1)  +fin(i,j+1,k+1))
     >              +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)*hyi
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dxdz:
C
            iadr=iadr+1
            sum=  -(
     >              -(ypi*fin(i,j,k)  +yp*fin(i,j+1,k))
     >              +(ypi*fin(i+1,j,k)+yp*fin(i+1,j+1,k)))
     >            +(
     >              -(ypi*fin(i,j,k+1)  +yp*fin(i,j+1,k+1))
     >              +(ypi*fin(i+1,j,k+1)+yp*fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)*hzi
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dydz:
C
            iadr=iadr+1
            sum=  -(
     >           xpi*(-fin(i,j,k)  +fin(i,j+1,k))+
     >            xp*(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >            +(
     >           xpi*(-fin(i,j,k+1)  +fin(i,j+1,k+1))+
     >            xp*(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hyi*hzi
         endif
C
         if(ict(8).eq.1) then
C
C  d3f/dxdydz:
C
            iadr=iadr+1
            sum=  -(
     >              -(-fin(i,j,k)  +fin(i,j+1,k))
     >              +(-fin(i+1,j,k)+fin(i+1,j+1,k)))
     >            +(
     >              -(-fin(i,j,k+1)  +fin(i,j+1,k+1))
     >              +(-fin(i+1,j,k+1)+fin(i+1,j+1,k+1)))
C
            fval(v,iadr)=sum*hxi(v)*hyi*hzi
         endif
C
      enddo                             ! vector loop
C
      return
      end
