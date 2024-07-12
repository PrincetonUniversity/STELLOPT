      subroutine herm3ev(xget,yget,zget,x,nx,y,ny,z,nz,
     >                   ilinx,iliny,ilinz,
     >                   f,inf2,inf3,ict,fval,ier)
C
C  evaluate a 3d cubic Hermite interpolant on a rectilinear
C  grid -- this is C1 in all directions.
C
C  this subroutine calls two subroutines:
C     herm3xyz  -- find cell containing (xget,yget,zget)
C     herm3fcn -- evaluate interpolant function and (optionally) derivatives
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
C  f(1,i,j,k) = df/dx @ x(i),y(j),z(k)
C  f(2,i,j,k) = df/dy @ x(i),y(j),z(k)
C  f(3,i,j,k) = df/dz @ x(i),y(j),z(k)
C  f(4,i,j,k) = d2f/dxdy @ x(i),y(j),z(k)
C  f(5,i,j,k) = d2f/dxdz @ x(i),y(j),z(k)
C  f(6,i,j,k) = d2f/dydz @ x(i),y(j),z(k)
C  f(7,i,j,k) = d3f/dxdydz @ x(i),y(j),z(k)
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
      real fval(*)                      ! output data
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
      real xparam,yparam,zparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real hx,hy,hz
      real hxi,hyi,hzi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C  0 .le. zparam .le. 1
C
C---------------------------------------------------------------------
C
      call herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz,
     >   i,j,k,xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi,ier)
      if(ier.ne.0) return
c
      call herm3fcn(ict,1,1,
     >   fval,i,j,k,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   f,inf2,inf3,nz)
C
      return
      end
C---------------------------------------------------------------------
c  herm3xyz -- look up x-y-z zone
c
c  this is the "first part" of herm3ev, see comments, above.
c
      subroutine herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,
     >   ilinx,iliny,ilinz,
     >   i,j,k,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,ier)
c
c  input of herm3xyz
c  ================
c
      integer nx,ny,nz                  ! coord. grid dimensions
c
      real xget,yget,zget               ! target point
      real x(nx),y(ny),z(nz)            ! indep. coords. (ascending order)
c
      integer ilinx                     ! =1:  x evenly spaced
      integer iliny                     ! =1:  y evenly spaced
      integer ilinz                     ! =1:  z evenly spaced
c
c  output of herm3xyz
c  =================
      integer i,j,k                     ! index to cell containing target pt
c          on exit:  1.le.i.le.nx-1   1.le.j.le.ny-1  1.le.k.le.nz-1
c
c  normalized position w/in (i,j) cell (btw 0 and 1):
c
      real xparam                       ! (xget-x(i))/(x(i+1)-x(i))
      real yparam                       ! (yget-y(j))/(y(j+1)-y(j))
      real zparam                       ! (zget-z(k))/(z(k+1)-z(k))
c
c  grid spacing
c
      real hx                           ! hx = x(i+1)-x(i)
      real hy                           ! hy = y(j+1)-y(j)
      real hz                           ! hz = z(k+1)-z(k)
c
c  inverse grid spacing:
c
      real hxi                          ! 1/hx = 1/(x(i+1)-x(i))
      real hyi                          ! 1/hy = 1/(y(j+1)-y(j))
      real hzi                          ! 1/hz = 1/(z(k+1)-z(k))
c
      integer ier                       ! return ier.ne.0 on error
c
c------------------------------------
c
      ier=0
c
c  range check
c
      zxget=xget
      zyget=yget
      zzget=zget
      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0e-7*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
            write(6,1001) xget,x(1),x(nx)
 1001       format(' ?herm3ev:  xget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((xget.lt.x(1)-0.5*zxtol).or.
     >         (xget.gt.x(nx)+0.5*zxtol))
     >      write(6,1011) xget,x(1),x(nx)
 1011       format(' %herm3ev:  xget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(xget.lt.x(1)) then
               zxget=x(1)
            else
               zxget=x(nx)
            endif
         endif
      endif
      if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
         zytol=4.0e-7*max(abs(y(1)),abs(y(ny)))
         if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
            ier=1
            write(6,1002) yget,y(1),y(ny)
 1002       format(' ?herm3ev:  yget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((yget.lt.y(1)-0.5*zytol).or.
     >         (yget.gt.y(ny)+0.5*zytol))
     >      write(6,1012) yget,y(1),y(ny)
 1012       format(' %herm3ev:  yget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(yget.lt.y(1)) then
               zyget=y(1)
            else
               zyget=y(ny)
            endif
         endif
      endif
      if((zget.lt.z(1)).or.(zget.gt.z(nz))) then
         zztol=4.0e-7*max(abs(z(1)),abs(z(nz)))
         if((zget.lt.z(1)-zztol).or.(zget.gt.z(nz)+zztol)) then
            ier=1
            write(6,1003) zget,z(1),z(nz)
 1003       format(' ?herm3ev:  zget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((zget.lt.z(1)-0.5*zztol).or.
     >         (zget.gt.z(nz)+0.5*zztol))
     >      write(6,1013) zget,z(1),z(nz)
 1013       format(' %herm3ev:  zget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(zget.lt.z(1)) then
               zzget=z(1)
            else
               zzget=z(nz)
            endif
         endif
      endif
      if(ier.ne.0) return
c
c  now find interval in which target point lies..
c
      nxm=nx-1
      nym=ny-1
      nzm=nz-1
c
      if(ilinx.eq.1) then
         ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
         i=min(nxm, ii)
         if(zxget.lt.x(i)) then
            i=i-1
         else if(zxget.gt.x(i+1)) then
            i=i+1
         endif
      else
         if((1.le.i).and.(i.lt.nxm)) then
            if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
               continue  ! already have the zone
            else
               call zonfind(x,nx,zxget,i)
            endif
         else
            i=nx/2
            call zonfind(x,nx,zxget,i)
         endif
      endif
c
      if(iliny.eq.1) then
         jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
         j=min(nym, jj)
         if(zyget.lt.y(j)) then
            j=j-1
         else if(zyget.gt.y(j+1)) then
            j=j+1
         endif
      else
         if((1.le.j).and.(j.lt.nym)) then
            if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
               continue  ! already have the zone
            else
               call zonfind(y,ny,zyget,j)
            endif
         else
            j=ny/2
            call zonfind(y,ny,zyget,j)
         endif
      endif
c
      if(ilinz.eq.1) then
         kk=1+nzm*(zzget-z(1))/(z(nz)-z(1))
         k=min(nzm, kk)
         if(zzget.lt.z(k)) then
            k=k-1
         else if(zzget.gt.z(k+1)) then
            k=k+1
         endif
      else
         if((1.le.k).and.(k.lt.nzm)) then
            if((z(k).le.zzget).and.(zzget.le.z(k+1))) then
               continue  ! already have the zone
            else
               call zonfind(z,nz,zzget,k)
            endif
         else
            k=nz/2
            call zonfind(z,nz,zzget,k)
         endif
      endif
c
      hx=(x(i+1)-x(i))
      hy=(y(j+1)-y(j))
      hz=(z(k+1)-z(k))
c
      hxi=1.0/hx
      hyi=1.0/hy
      hzi=1.0/hz
c
      xparam=(zxget-x(i))*hxi
      yparam=(zyget-y(j))*hyi
      zparam=(zzget-z(k))*hzi
c
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine herm3fcn(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
C
      integer ict(8)                    ! requested output control
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
      real fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "herm3ev")
C
      real fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine herm3ev
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
C  Hermite cubic basis functions
C  -->for function value matching
C     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
C   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
C
C   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.0*x+3.0)
C   abar(x)=1-a(x)
C   a'(x)=-abar'(x)          = 6.0*x*(1.0-x)
C
C  -->for derivative matching
C     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
C   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
C
C     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
C     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
C
      real sum
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
         xpi=1.0-xp
         xp2=xp*xp
         xpi2=xpi*xpi
         ax=xp2*(3.0-2.0*xp)
         axbar=1.0-ax
         bx=-xp2*xpi
         bxbar=xpi2*xp
C
C   ...in y direction
C
         yp=yparam(v)
         ypi=1.0-yp
         yp2=yp*yp
         ypi2=ypi*ypi
         ay=yp2*(3.0-2.0*yp)
         aybar=1.0-ay
         by=-yp2*ypi
         bybar=ypi2*yp
C
C   ...in z direction
C
         zp=zparam(v)
         zpi=1.0-zp
         zp2=zp*zp
         zpi2=zpi*zpi
         az=zp2*(3.0-2.0*zp)
         azbar=1.0-az
         bz=-zp2*zpi
         bzbar=zpi2*zp
C
         iadr=0
C
C  derivatives:
C
         axp=6.0*xp*xpi
         axbarp=-axp
         bxp=xp*(3.0*xp-2.0)
         bxbarp=xpi*(3.0*xpi-2.0)
C
         ayp=6.0*yp*ypi
         aybarp=-ayp
         byp=yp*(3.0*yp-2.0)
         bybarp=ypi*(3.0*ypi-2.0)
C
         azp=6.0*zp*zpi
         azbarp=-azp
         bzp=zp*(3.0*zp-2.0)
         bzbarp=zpi*(3.0*zpi-2.0)
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            sum=azbar*(
     >         axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
C
            sum=sum+hx(v)*(
     >         azbar*(
     >         bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*(
     >         azbar*(
     >         axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz(v)*(
     >         bzbar*(
     >         axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hy(v)*(
     >         azbar*(
     >         bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hz(v)*(
     >         bzbar*(
     >         bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*hz(v)*(
     >         bzbar*(
     >         axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hy(v)*hz(v)*(
     >         bzbar*(
     >         bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C     df/dx:
C
            iadr=iadr+1
C
            sum=hxi(v)*(
     >          azbar*(
     >         axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         azbar*(
     >         bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy(v)*(
     >         azbar*(
     >         axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hz(v)*(
     >         bzbar*(
     >         axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*(
     >         azbar*(
     >         bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz(v)*(
     >         bzbar*(
     >         bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy(v)*hz(v)*(
     >         bzbar*(
     >         axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*hz(v)*(
     >         bzbar*(
     >         bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
C
            sum=hyi(v)*(
     >          azbar*(
     >         axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*hx(v)*(
     >         azbar*(
     >         bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         azbar*(
     >         axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*hz(v)*(
     >         bzbar*(
     >         axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*(
     >         azbar*(
     >         bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hyi(v)*hz(v)*(
     >         bzbar*(
     >         bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz(v)*(
     >         bzbar*(
     >         axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hz(v)*(
     >         bzbar*(
     >         bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            iadr=iadr+1
C
            sum=hzi(v)*(
     >          azbarp*(
     >         axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi(v)*hx(v)*(
     >         azbarp*(
     >         bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi(v)*hy(v)*(
     >         azbarp*(
     >         axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi(v)*hx(v)*hy(v)*(
     >         azbarp*(
     >         bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*(
     >         bzbarp*(
     >         bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*(
     >         bzbarp*(
     >         axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hy(v)*(
     >         bzbarp*(
     >         bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
C
            sum=hxi(v)*hyi(v)*(
     >          azbar*(
     >         axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*(
     >         azbar*(
     >         bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*(
     >         azbar*(
     >         axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hyi(v)*hz(v)*(
     >         bzbar*(
     >         axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         azbar*(
     >         bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*hz(v)*(
     >         bzbar*(
     >         bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hz(v)*(
     >         bzbar*(
     >         axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz(v)*(
     >         bzbar*(
     >         bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dxdz:
C
            iadr=iadr+1
C
            sum=hxi(v)*hzi(v)*(
     >          azbarp*(
     >         axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi(v)*(
     >         azbarp*(
     >         bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy(v)*hzi(v)*(
     >         azbarp*(
     >         axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*(
     >         bzbarp*(
     >         axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*hzi(v)*(
     >         azbarp*(
     >         bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy(v)*(
     >         bzbarp*(
     >         axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy(v)*(
     >         bzbarp*(
     >         bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dydz:
C
            iadr=iadr+1
C
            sum=hyi(v)*hzi(v)*(
     >          azbarp*(
     >         axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*hzi(v)*hx(v)*(
     >         azbarp*(
     >         bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi(v)*(
     >         azbarp*(
     >         axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*(
     >         bzbarp*(
     >         axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hzi(v)*(
     >         azbarp*(
     >         bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hyi(v)*(
     >         bzbarp*(
     >         bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*(
     >         bzbarp*(
     >         bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(8).eq.1) then
C
C  d3f/dxdydz:
C
            iadr=iadr+1
C
            sum=hxi(v)*hyi(v)*hzi(v)*(
     >          azbarp*(
     >         axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*hzi(v)*(
     >         azbarp*(
     >         bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hzi(v)*(
     >         azbarp*(
     >         axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hyi(v)*(
     >         bzbarp*(
     >         axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi(v)*(
     >         azbarp*(
     >         bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi(v)*(
     >         bzbarp*(
     >         bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*(
     >         bzbarp*(
     >         axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
      enddo                             ! vector loop
C
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C    --optimized for VARIATION along x axis ONLY--
C
      subroutine herm3fcnx(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
C
      integer ict(8)                    ! requested output control
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
      real fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "herm3ev")
C
      real fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine herm3ev
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
C  Hermite cubic basis functions
C  -->for function value matching
C     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
C   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
C
C   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.0*x+3.0)
C   abar(x)=1-a(x)
C   a'(x)=-abar'(x)          = 6.0*x*(1.0-x)
C
C  -->for derivative matching
C     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
C   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
C
C     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
C     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
C
      real sum
      integer v
C
      j=jj
      k=kk
C
C   ...in y direction
C
      yp=yparam
      ypi=1.0-yp
      yp2=yp*yp
      ypi2=ypi*ypi
      ay=yp2*(3.0-2.0*yp)
      aybar=1.0-ay
      by=-yp2*ypi
      bybar=ypi2*yp
C
C   ...in z direction
C
      zp=zparam
      zpi=1.0-zp
      zp2=zp*zp
      zpi2=zpi*zpi
      az=zp2*(3.0-2.0*zp)
      azbar=1.0-az
      bz=-zp2*zpi
      bzbar=zpi2*zp
C
C  derivatives
C
      ayp=6.0*yp*ypi
      aybarp=-ayp
      byp=yp*(3.0*yp-2.0)
      bybarp=ypi*(3.0*ypi-2.0)
C
      azp=6.0*zp*zpi
      azbarp=-azp
      bzp=zp*(3.0*zp-2.0)
      bzbarp=zpi*(3.0*zpi-2.0)

      do v=1,ivec
         i=ii(v)
C
C   ...in x direction
C
         xp=xparam(v)
         xpi=1.0-xp
         xp2=xp*xp
         xpi2=xpi*xpi
         ax=xp2*(3.0-2.0*xp)
         axbar=1.0-ax
         bx=-xp2*xpi
         bxbar=xpi2*xp
C
         iadr=0
C
C  derivatives:
C
         axp=6.0*xp*xpi
         axbarp=-axp
         bxp=xp*(3.0*xp-2.0)
         bxbarp=xpi*(3.0*xpi-2.0)
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            sum=azbar*(
     >         axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
C
            sum=sum+hx(v)*(
     >         azbar*(
     >         bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*(
     >         azbar*(
     >         axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz*(
     >         bzbar*(
     >         axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hy*(
     >         azbar*(
     >         bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hz*(
     >         bzbar*(
     >         bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*hz*(
     >         bzbar*(
     >         axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hy*hz*(
     >         bzbar*(
     >         bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C     df/dx:
C
            iadr=iadr+1
C
            sum=hxi(v)*(
     >          azbar*(
     >         axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         azbar*(
     >         bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy*(
     >         azbar*(
     >         axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hz*(
     >         bzbar*(
     >         axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*(
     >         azbar*(
     >         bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz*(
     >         bzbar*(
     >         bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy*hz*(
     >         bzbar*(
     >         axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*hz*(
     >         bzbar*(
     >         bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
C
            sum=hyi*(
     >          azbar*(
     >         axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*hx(v)*(
     >         azbar*(
     >         bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         azbar*(
     >         axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*hz*(
     >         bzbar*(
     >         axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*(
     >         azbar*(
     >         bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hyi*hz*(
     >         bzbar*(
     >         bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz*(
     >         bzbar*(
     >         axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hz*(
     >         bzbar*(
     >         bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            iadr=iadr+1
C
            sum=hzi*(
     >          azbarp*(
     >         axbar*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            ax*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbar*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            ax*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi*hx(v)*(
     >         azbarp*(
     >         bxbar*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bx*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bx*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi*hy*(
     >         azbarp*(
     >         axbar*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            ax*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbar*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            ax*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         axbar*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            ax*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            ax*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi*hx(v)*hy*(
     >         azbarp*(
     >         bxbar*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bx*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bx*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*(
     >         bzbarp*(
     >         bxbar*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bx*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bx*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*(
     >         bzbarp*(
     >         axbar*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            ax*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            ax*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hy*(
     >         bzbarp*(
     >         bxbar*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bx*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bx*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
C
            sum=hxi(v)*hyi*(
     >          azbar*(
     >         axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  az*(
     >         axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*(
     >         azbar*(
     >         bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*(
     >         azbar*(
     >         axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + az*(
     >         axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hyi*hz*(
     >         bzbar*(
     >         axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         azbar*(
     >         bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + az*(
     >         bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*hz*(
     >         bzbar*(
     >         bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hz*(
     >         bzbar*(
     >         axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bz*(
     >         axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hz*(
     >         bzbar*(
     >         bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bz*(
     >         bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dxdz:
C
            iadr=iadr+1
C
            sum=hxi(v)*hzi*(
     >          azbarp*(
     >         axbarp*(aybar*fin(0,i,j,k)  +ay*fin(0,i,j+1,k))+
     >            axp*(aybar*fin(0,i+1,j,k)+ay*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbarp*(aybar*fin(0,i,j,k+1)  +ay*fin(0,i,j+1,k+1))+
     >            axp*(aybar*fin(0,i+1,j,k+1)+ay*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi*(
     >         azbarp*(
     >         bxbarp*(aybar*fin(1,i,j,k)  +ay*fin(1,i,j+1,k))+
     >            bxp*(aybar*fin(1,i+1,j,k)+ay*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(aybar*fin(1,i,j,k+1)  +ay*fin(1,i,j+1,k+1))+
     >            bxp*(aybar*fin(1,i+1,j,k+1)+ay*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy*hzi*(
     >         azbarp*(
     >         axbarp*(bybar*fin(2,i,j,k)  +by*fin(2,i,j+1,k))+
     >            axp*(bybar*fin(2,i+1,j,k)+by*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbarp*(bybar*fin(2,i,j,k+1)  +by*fin(2,i,j+1,k+1))+
     >            axp*(bybar*fin(2,i+1,j,k+1)+by*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*(
     >         bzbarp*(
     >         axbarp*(aybar*fin(3,i,j,k)  +ay*fin(3,i,j+1,k))+
     >            axp*(aybar*fin(3,i+1,j,k)+ay*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(aybar*fin(3,i,j,k+1)  +ay*fin(3,i,j+1,k+1))+
     >            axp*(aybar*fin(3,i+1,j,k+1)+ay*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*hzi*(
     >         azbarp*(
     >         bxbarp*(bybar*fin(4,i,j,k)  +by*fin(4,i,j+1,k))+
     >            bxp*(bybar*fin(4,i+1,j,k)+by*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(bybar*fin(4,i,j,k+1)  +by*fin(4,i,j+1,k+1))+
     >            bxp*(bybar*fin(4,i+1,j,k+1)+by*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         bxbarp*(aybar*fin(5,i,j,k)  +ay*fin(5,i,j+1,k))+
     >            bxp*(aybar*fin(5,i+1,j,k)+ay*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(aybar*fin(5,i,j,k+1)  +ay*fin(5,i,j+1,k+1))+
     >            bxp*(aybar*fin(5,i+1,j,k+1)+ay*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hy*(
     >         bzbarp*(
     >         axbarp*(bybar*fin(6,i,j,k)  +by*fin(6,i,j+1,k))+
     >            axp*(bybar*fin(6,i+1,j,k)+by*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(bybar*fin(6,i,j,k+1)  +by*fin(6,i,j+1,k+1))+
     >            axp*(bybar*fin(6,i+1,j,k+1)+by*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hy*(
     >         bzbarp*(
     >         bxbarp*(bybar*fin(7,i,j,k)  +by*fin(7,i,j+1,k))+
     >            bxp*(bybar*fin(7,i+1,j,k)+by*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(bybar*fin(7,i,j,k+1)  +by*fin(7,i,j+1,k+1))+
     >            bxp*(bybar*fin(7,i+1,j,k+1)+by*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dydz:
C
            iadr=iadr+1
C
            sum=hyi*hzi*(
     >          azbarp*(
     >         axbar*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            ax*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbar*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            ax*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*hzi*hx(v)*(
     >         azbarp*(
     >         bxbar*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bx*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bx*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi*(
     >         azbarp*(
     >         axbar*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            ax*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbar*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            ax*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*(
     >         bzbarp*(
     >         axbar*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            ax*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            ax*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hzi*(
     >         azbarp*(
     >         bxbar*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bx*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbar*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bx*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*hyi*(
     >         bzbarp*(
     >         bxbar*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bx*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bx*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         axbar*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            ax*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbar*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            ax*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hx(v)*(
     >         bzbarp*(
     >         bxbar*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bx*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbar*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bx*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
         if(ict(8).eq.1) then
C
C  d3f/dxdydz:
C
            iadr=iadr+1
C
            sum=hxi(v)*hyi*hzi*(
     >          azbarp*(
     >         axbarp*(aybarp*fin(0,i,j,k)  +ayp*fin(0,i,j+1,k))+
     >            axp*(aybarp*fin(0,i+1,j,k)+ayp*fin(0,i+1,j+1,k)))
     >          +  azp*(
     >         axbarp*(aybarp*fin(0,i,j,k+1)  +ayp*fin(0,i,j+1,k+1))+
     >            axp*(aybarp*fin(0,i+1,j,k+1)+ayp*fin(0,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*hzi*(
     >         azbarp*(
     >         bxbarp*(aybarp*fin(1,i,j,k)  +ayp*fin(1,i,j+1,k))+
     >            bxp*(aybarp*fin(1,i+1,j,k)+ayp*fin(1,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(aybarp*fin(1,i,j,k+1)  +ayp*fin(1,i,j+1,k+1))+
     >            bxp*(aybarp*fin(1,i+1,j,k+1)+ayp*fin(1,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hzi*(
     >         azbarp*(
     >         axbarp*(bybarp*fin(2,i,j,k)  +byp*fin(2,i,j+1,k))+
     >            axp*(bybarp*fin(2,i+1,j,k)+byp*fin(2,i+1,j+1,k)))
     >          + azp*(
     >         axbarp*(bybarp*fin(2,i,j,k+1)  +byp*fin(2,i,j+1,k+1))+
     >            axp*(bybarp*fin(2,i+1,j,k+1)+byp*fin(2,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*hyi*(
     >         bzbarp*(
     >         axbarp*(aybarp*fin(3,i,j,k)  +ayp*fin(3,i,j+1,k))+
     >            axp*(aybarp*fin(3,i+1,j,k)+ayp*fin(3,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(aybarp*fin(3,i,j,k+1)  +ayp*fin(3,i,j+1,k+1))+
     >            axp*(aybarp*fin(3,i+1,j,k+1)+ayp*fin(3,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hzi*(
     >         azbarp*(
     >         bxbarp*(bybarp*fin(4,i,j,k)  +byp*fin(4,i,j+1,k))+
     >            bxp*(bybarp*fin(4,i+1,j,k)+byp*fin(4,i+1,j+1,k)))
     >          + azp*(
     >         bxbarp*(bybarp*fin(4,i,j,k+1)  +byp*fin(4,i,j+1,k+1))+
     >            bxp*(bybarp*fin(4,i+1,j,k+1)+byp*fin(4,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hyi*(
     >         bzbarp*(
     >         bxbarp*(aybarp*fin(5,i,j,k)  +ayp*fin(5,i,j+1,k))+
     >            bxp*(aybarp*fin(5,i+1,j,k)+ayp*fin(5,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(aybarp*fin(5,i,j,k+1)  +ayp*fin(5,i,j+1,k+1))+
     >            bxp*(aybarp*fin(5,i+1,j,k+1)+ayp*fin(5,i+1,j+1,k+1)))
     >         )
C
            sum=sum+hxi(v)*(
     >         bzbarp*(
     >         axbarp*(bybarp*fin(6,i,j,k)  +byp*fin(6,i,j+1,k))+
     >            axp*(bybarp*fin(6,i+1,j,k)+byp*fin(6,i+1,j+1,k)))
     >          + bzp*(
     >         axbarp*(bybarp*fin(6,i,j,k+1)  +byp*fin(6,i,j+1,k+1))+
     >            axp*(bybarp*fin(6,i+1,j,k+1)+byp*fin(6,i+1,j+1,k+1)))
     >         )
C
            sum=sum+(
     >         bzbarp*(
     >         bxbarp*(bybarp*fin(7,i,j,k)  +byp*fin(7,i,j+1,k))+
     >            bxp*(bybarp*fin(7,i+1,j,k)+byp*fin(7,i+1,j+1,k)))
     >          + bzp*(
     >         bxbarp*(bybarp*fin(7,i,j,k+1)  +byp*fin(7,i,j+1,k+1))+
     >            bxp*(bybarp*fin(7,i+1,j,k+1)+byp*fin(7,i+1,j+1,k+1)))
     >         )
C
            fval(v,iadr)=sum
         endif
C
      enddo                             ! vector loop
C
      return
      end
