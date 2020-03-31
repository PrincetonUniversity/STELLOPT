c     tcspeval -- eval tricubic spline function and/or derivatives
c
      subroutine tcspeval(xget,yget,zget,iselect,fval,
     >                    x,nx,y,ny,z,nz,ilinx,iliny,ilinz,f,inf4,inf5,
     >                    ier)
c
      integer iselect(10)
      integer ilinx,iliny,ilinz,nx,ny,nz,inf4,inf5,ier
c
      real xget,yget,zget
      real fval(*)
      real x(nx),y(ny),z(nz),f(4,4,4,inf4,inf5,nz)
c
c  modification -- dmc 11 Jan 1999 -- remove SAVE stmts; break routine
C    into these parts:
C
C    tcspevxyz -- find grid cell of target pt.
C    tcspevfn -- evaluate function using output of tcspevxyz
C
C    in cases where multiple functions are defined on the same grid,
C    time can be saved by using tcspevxyz once and then tcspevfn
C    multiple times.
c
c  input:
c     (xget,yget,zget)   location where interpolated value is desired
c                   x(1).le.xget.le.x(nx) expected
c                   y(1).le.yget.le.y(ny) expected
c                   z(1).le.zget.le.z(nz) expected
c
c     iselect       select desired output
c
c                     iselect(1)=1 -- want function value (f) itself
c                     iselect(2)=1 -- want  df/dx
c                     iselect(3)=1 -- want  df/dy
c                     iselect(4)=1 -- want  df/dz
c                     iselect(5)=1 -- want  d2f/dx2
c                     iselect(6)=1 -- want  d2f/dy2
c                     iselect(7)=1 -- want  d2f/dz2
c                     iselect(8)=1 -- want  d2f/dxdy
c                     iselect(9)=1 -- want  d2f/dxdz
c                     iselect(10)=1 -- want  d2f/dydz
c
c
c              example:  iselect(1)=iselect(2)=iselect(3)=iselect(4)=1
c                            f, df/dx, df/dy, and df/dz all evaluated
c                        iselect(5)=iselect(6)=iselect(7)=0
c                        iselect(8)=iselect(9)=iselect(10)=0
c                            2nd derivatives not evaluated.
c
c  (new dmc Dec 2005 -- higher derivatives available)
c    iselect(1)=3 --> 3rd derivative, .le.2 diff. in any coordinate
c      iselect(2:8) select: fxxy fxxz fxyy fxyz fxzz fyyz fyzz
c      ->note iselect(1)=3, iselect(5)=1 gives fxyz = d3f/dxdydz
c    iselect(1)=-3 --> 3rd derivative, 3 in one coordinate
c      iselect(2:4) select: fxxx fyyy fzzz
c    iselect(1)=4 --> 3rd derivative, .le.2 diff. in any coordinate
c      iselect(2:7) select: fxxyy fxxyz fxxzz fxyyz fxyzz fyyzz
c    iselect(1)=-4 --> 3rd derivative, 3 in one coordinate
c      iselect(2:7) select: fxxxy fxxxz fxyyy fxzzz fyyyz fyzzz
c    iselect(1)=5 --> 3rd derivative, .le.2 diff. in any coordinate
c      iselect(2:4) select: fxxyyz fxxyzz fxyyzz
c    iselect(1)=-5 --> 3rd derivative, 3 in one coordinate
c      iselect(2:10) select:  fxxxyy fxxxyz fxxxzz fxxyyy fxxzzz
c                             fxyyyz fxyzzz fyyyzz fzzzyy
c    iselect(1)=6 --> 3rd derivative, .le.2 diff. in any coordinate
c      fxxyyzz
c    iselect(1)=-6 --> 3rd derivative, 3 in one coordinate
c      iselect(2:10) select: fxxxyyy fxxxyyz fxxxyzz fxxxyyz
c                            fxxyyyz fxxyzzz fxyyyzz fxyyzzz fyyyzzz
c    iselect(1)=-7 --> 7th derivative
c      iselect(2:7) select: fxxxyyyz fxxxyyzz fxxxyzzz
c                           fxxyyyzz fxxyyzzz fxyyyzzz
c    iselect(1)=-8 --> 8th derivative
c      iselect(2:4) select: fxxxyyyzz fxxxyyzzz fxxyyyzzz
c    iselect(1)=-9 --> 9th derivative:  fxxxyyyzzz
c
c-------
c
c     x(1...nx)     independent coordinate x, strict ascending
c     y(1...ny)     independent coordinate y, strict ascending
c     z(1...nz)     independent coordinate y, strict ascending
c
c     ilinx  --  =1: flag that x is linearly spaced
c
c                   see fval (output) description.
c
c     x(1...nx)     independent coordinate x, strict ascending
c     y(1...ny)     independent coordinate y, strict ascending
c     z(1...nz)     independent coordinate y, strict ascending
c
c     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
c     iliny  --  =1: flag that y is linearly spaced (avoid search for speed)
c     ilinz  --  =1: flat that z is linearly spaced (avoid search for speed)
c
c  **CAUTION** actual even spacing of x, y, z is NOT CHECKED HERE!
c
c
c     f             the function values (at grid points) and spline coefs
c
c  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                             and y btw y(j) and y(j+1), dy=y-y(j),
c                             and z btw z(k) and z(k+1), dz=z-z(k)
c
c  do m=1,4
c   p(m) =
c    f(1,1,m,i,j,k)+dx*f(2,1,m,i,j,k)+dx**2*f(3,1,m,i,j,k)+dx**3*f(4,1,m,i,j,k)
c   +dy*(
c   f(1,2,m,i,j,k)+dx*f(2,2,m,i,j,k)+dx**2*f(3,2,m,i,j,k)+dx**3*f(4,2,m,i,j,k))
c   +dy**2*(
c   f(1,3,m,i,j,k)+dx*f(2,3,m,i,j,k)+dx**2*f(3,3,m,i,j,k)+dx**3*f(4,3,m,i,j,k))
c   +dy**3*(
c   f(1,4,m,i,j,k)+dx*f(2,4,m,i,j,k)+dx**2*f(3,4,m,i,j,k)+dx**3*f(4,4,m,i,j,k))
c  enddo
c  answer = p(1)+dz*p(2)+dz**2*p(3)+dz**3*p(4)
c
c      where d2=dy**2 and d3=dy**3.
c
c  nb dmc Feb 1999 -- p loops unrolled, by hand, to aid vector compilers
c
c  output:
c      up to 10 elements of fval, ordered as follows:
c        fval(1)=function value or lowest order derivative requested
c        fval(2)=next order derivative
c             etc
c        the ordering is a subset of the sequence given under the "iselect"
c        description; the first M elements of fval are used, where M = the
c        number of non-zero elements of iselect.
c
c      ier = 0 -- successful completion; = 1 -- an error occurred.
c
c-------------------------------------------------------------------
c  local
c
      integer :: i(1),j(1),k(1)
c
      real dx(1),dy(1),dz(1)
c
c--------------------------
c
      i(1)=0
      j(1)=0
      k(1)=0
c
      call tcspevxyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz,
     >   i(1),j(1),k(1),dx(1),dy(1),dz(1),ier)
      if(ier.ne.0) return
c
      call tcspevfn(iselect,1,1,fval,i,j,k,dx,dy,dz,f,inf4,inf5,nz)
c
      return
      end
c
c-------------------------------------------------------------------------
c  tcspevxyz -- look up x-y zone
c
c  this is the "first part" of tcspeval, see comments, above.
c
      subroutine tcspevxyz(xget,yget,zget,x,nx,y,ny,z,nz,
     >   ilinx,iliny,ilinz,
     >   i,j,k,dx,dy,dz,ier)
c
      integer nx,ny,nz                  ! array dimensions
c
      real xget,yget,zget               ! target point
      real x(nx),y(ny),z(nz)            ! indep. coords.
c
      integer ilinx                     ! =1:  assume x evenly spaced
      integer iliny                     ! =1:  assume y evenly spaced
      integer ilinz                     ! =1:  assume z evenly spaced
c
c  output of tcspevxyz
c
      integer i,j,k                     ! index to cell containing target pt
      real dx,dy,dz                     ! displacement of target pt w/in cell
                                        ! dx=x-x(i)  dy=y-y(j)  dz=z-z(k)
C
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
c
      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0e-7*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
            write(6,1001) xget,x(1),x(nx)
 1001       format(' ?tcspeval:  xget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((xget.lt.x(1)-0.5*zxtol).or.
     >         (xget.gt.x(nx)+0.5*zxtol))
     >      write(6,1011) xget,x(1),x(nx)
 1011       format(' %tcspeval:  xget=',1pe15.8,' beyond range ',
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
 1002       format(' ?tcspeval:  yget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((yget.lt.y(1)-0.5*zytol).or.
     >         (yget.gt.y(ny)+0.5*zytol))
     >      write(6,1012) yget,y(1),y(ny)
 1012       format(' %tcspeval:  yget=',1pe15.8,' beyond range ',
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
 1003       format(' ?tcspeval:  zget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((zget.lt.z(1)-0.5*zztol).or.
     >         (zget.gt.z(nz)+0.5*zztol))
     >      write(6,1013) zget,z(1),z(nz)
 1013       format(' %tcspeval:  zget=',1pe15.8,' beyond range ',
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
      dx=zxget-x(i)
      dy=zyget-y(j)
      dz=zzget-z(k)
c
      return
      end
c------------------------------------------------------------------------
c  tcspevfn -- OK now evaluate the tricubic spline
c     evaluate at a vector set of target locations as specified by
c     input vectors (iv,jv,kv), (dxv,dyv,dzv)
c
      subroutine tcspevfn(ict,ivec,ivd,fval,iv,jv,kv,dxv,dyv,dzv,
     >   f,inf4,inf5,nz)
c
c  input:
      integer ict(10)                   ! selector:
c        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
c        ict(2)=1 for df/dx   ""
c        ict(3)=1 for df/dy   ""
c        ict(4)=1 for df/dz   ""
c        ict(5)=1 for d2f/dx2
c        ict(6)=1 for d2f/dy2
c        ict(7)=1 for d2f/dz2
c        ict(8)=1 for d2f/dxdy
c        ict(9)=1 for d2f/dxdz
c        ict(10)=1 for d2f/dydz
c
      integer ivec,ivd                  ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c    ivd -- 1st dimension of fval, .ge. ivec
c
c output:
      real fval(ivd,*)                     ! output array
c
c  for vector elements v,  (iv(v),jv(v),kv(v),dxv(v),dyv(v),dzv(v))
c
c  fval(v,1) = first item requested by ict(...),
c  fval(v,2) = 2nd item requested,  ...etc...
c
c  input:
      integer iv(ivec),jv(ivec),kv(ivec) ! grid cell indices
      real dxv(ivec),dyv(ivec),dzv(ivec) ! displacements w/in cell
c
      integer inf4                      ! 4th dimension of f, .le.nx
      integer inf5                      ! 5th dimension of f, .le.ny
      real f(4,4,4,inf4,inf5,nz)        ! tricubic fcn spline coeffs array
c
c  usage example:
c
c  1.  for each element (xx(v),yy(v),zz(v)) in a vector of (x,y,z)
c    triples, find the x,y,z zone indices and displacements with
c    to the "lower left corner" of the zone; store these in vectors
c    iv,jv,kv and dxv,dyv,dzv
c
c  2.  set ict(1)=0, ict(2)=1, ict(3)=1, ict(4)=1 & the rest zero --
c      to get only the 1st derivatives.
c
c  3.  ivec is the length of the vector; ivd is the 1st dimension
c      of the array fval to receive the output
c
c      real fval(ivd,10)
c      real xv(ivd),yv(ivd),zv(ivd)
c      integer iv(ivd),jv(ivd),kv(ivd)
c      real dxv(ivd),dyv(ivd),dzv(ivd)
c      integer ict(10)
c
c      real fspline(4,4,4,nx,ny,nz)  ! spline coeffs
c      data ict/0,1,1,1,0,0,0,0,0,0/    ! this call:  want 1st
c                               ! derivatives only ... these will
c                               ! be output to
c                               ! fval(*,1) fval(*,2) fval(*,3)
c      ...
c      do iv=1,ivec
c        ...                    ! find indices and displacements
c      enddo
c      call tcspevfn(ict,ivec,ivd,fval,iv,jv,kv,dxv,dyv,dzv,
c     >              fspline,nx,ny,nz)
c
c-------------------------------------------------------------------
c
c  local --
c
ccc      real p(4)  use p1,p2,p3,p4 now
c
      integer v                         ! vector index
c
c  OK can now do evaluations
c
      iaval=0  ! fval addressing
c
      if(abs(ict(1)).le.2) then
         if((ict(1).gt.0).or.(ict(1).eq.-1)) then
c  evaluate f
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1=
     >            f(1,1,1,i,j,k)+dy*(f(1,2,1,i,j,k)+dy*(f(1,3,1,i,j,k)+
     >            dy*f(1,4,1,i,j,k)))
     >        +dx*(f(2,1,1,i,j,k)+dy*(f(2,2,1,i,j,k)+dy*(f(2,3,1,i,j,k)+
     >            dy*f(2,4,1,i,j,k)))
     >        +dx*(f(3,1,1,i,j,k)+dy*(f(3,2,1,i,j,k)+dy*(f(3,3,1,i,j,k)+
     >            dy*f(3,4,1,i,j,k)))
     >        +dx*(f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+
     >            dy*f(4,4,1,i,j,k)))
     >            )))
               p2=
     >            f(1,1,2,i,j,k)+dy*(f(1,2,2,i,j,k)+dy*(f(1,3,2,i,j,k)+
     >            dy*f(1,4,2,i,j,k)))
     >        +dx*(f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+
     >            dy*f(2,4,2,i,j,k)))
     >        +dx*(f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+
     >            dy*f(3,4,2,i,j,k)))
     >        +dx*(f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k)))
     >            )))
               p3=
     >            f(1,1,3,i,j,k)+dy*(f(1,2,3,i,j,k)+dy*(f(1,3,3,i,j,k)+
     >            dy*f(1,4,3,i,j,k)))
     >        +dx*(f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+
     >            dy*f(2,4,3,i,j,k)))
     >        +dx*(f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k)))
     >        +dx*(f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k)))
     >            )))
               p4=
     >            f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+
     >            dy*f(1,4,4,i,j,k)))
     >        +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            )))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if((ict(2).gt.0).and.(ict(1).ne.-1)) then
c  evaluate df/dx
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1=
     >            f(2,1,1,i,j,k)+dy*(f(2,2,1,i,j,k)+dy*(f(2,3,1,i,j,k)+
     >            dy*f(2,4,1,i,j,k)))
     >         +2.*dx*
     >            (f(3,1,1,i,j,k)+dy*(f(3,2,1,i,j,k)+dy*(f(3,3,1,i,j,k)+
     >            dy*f(3,4,1,i,j,k)))
     >         +1.5*dx*
     >            (f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+
     >            dy*f(4,4,1,i,j,k)))
     >            ))
               p2=
     >            f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+
     >            dy*f(2,4,2,i,j,k)))
     >         +2.*dx*
     >            (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+
     >            dy*f(3,4,2,i,j,k)))
     >         +1.5*dx*
     >            (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k)))
     >            ))
               p3=
     >            f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+
     >            dy*f(2,4,3,i,j,k)))
     >         +2.*dx*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k)))
     >         +1.5*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k)))
     >            ))
               p4=
     >            f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >         +2.*dx*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >         +1.5*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            ))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if((ict(3).gt.0).and.(ict(1).ne.-1)) then
c  evaluate df/dy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1=
     >     f(1,2,1,i,j,k)+
     >         dy*(2.0*f(1,3,1,i,j,k)+dy*3.0*f(1,4,1,i,j,k))
     >+dx*(f(2,2,1,i,j,k)+
     >         dy*(2.0*f(2,3,1,i,j,k)+dy*3.0*f(2,4,1,i,j,k))
     >+dx*(f(3,2,1,i,j,k)+
     >         dy*(2.0*f(3,3,1,i,j,k)+dy*3.0*f(3,4,1,i,j,k))
     >+dx*(f(4,2,1,i,j,k)+
     >         dy*(2.0*f(4,3,1,i,j,k)+dy*3.0*f(4,4,1,i,j,k))
     >         )))
               p2=
     >     f(1,2,2,i,j,k)+
     >         dy*(2.0*f(1,3,2,i,j,k)+dy*3.0*f(1,4,2,i,j,k))
     >+dx*(f(2,2,2,i,j,k)+
     >         dy*(2.0*f(2,3,2,i,j,k)+dy*3.0*f(2,4,2,i,j,k))
     >+dx*(f(3,2,2,i,j,k)+
     >         dy*(2.0*f(3,3,2,i,j,k)+dy*3.0*f(3,4,2,i,j,k))
     >+dx*(f(4,2,2,i,j,k)+
     >         dy*(2.0*f(4,3,2,i,j,k)+dy*3.0*f(4,4,2,i,j,k))
     >         )))
               p3=
     >     f(1,2,3,i,j,k)+
     >         dy*(2.0*f(1,3,3,i,j,k)+dy*3.0*f(1,4,3,i,j,k))
     >+dx*(f(2,2,3,i,j,k)+
     >         dy*(2.0*f(2,3,3,i,j,k)+dy*3.0*f(2,4,3,i,j,k))
     >+dx*(f(3,2,3,i,j,k)+
     >         dy*(2.0*f(3,3,3,i,j,k)+dy*3.0*f(3,4,3,i,j,k))
     >+dx*(f(4,2,3,i,j,k)+
     >         dy*(2.0*f(4,3,3,i,j,k)+dy*3.0*f(4,4,3,i,j,k))
     >         )))
               p4=
     >     f(1,2,4,i,j,k)+
     >         dy*(2.0*f(1,3,4,i,j,k)+dy*3.0*f(1,4,4,i,j,k))
     >+dx*(f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >+dx*(f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >+dx*(f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >         )))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if((ict(4).gt.0).and.(ict(1).ne.-1)) then
c  evaluate df/dz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p2=
     >            f(1,1,2,i,j,k)+dy*(f(1,2,2,i,j,k)+dy*(f(1,3,2,i,j,k)+
     >            dy*f(1,4,2,i,j,k)))
     >        +dx*(f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+
     >            dy*f(2,4,2,i,j,k)))
     >        +dx*(f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+
     >            dy*f(3,4,2,i,j,k)))
     >        +dx*(f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k)))
     >         )))
               p3=
     >            f(1,1,3,i,j,k)+dy*(f(1,2,3,i,j,k)+dy*(f(1,3,3,i,j,k)+
     >            dy*f(1,4,3,i,j,k)))
     >        +dx*(f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+
     >            dy*f(2,4,3,i,j,k)))
     >        +dx*(f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k)))
     >        +dx*(f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k)))
     >         )))
               p4=
     >            f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+
     >            dy*f(1,4,4,i,j,k)))
     >        +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >         )))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
c
         if((ict(5).gt.0).or.(ict(1).eq.-1)) then
c  evaluate d2f/dx2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1= 2.0*
     >            (f(3,1,1,i,j,k)+dy*(f(3,2,1,i,j,k)+dy*(f(3,3,1,i,j,k)+
     >            dy*f(3,4,1,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+
     >            dy*f(4,4,1,i,j,k))))
               p2= 2.0*
     >            (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+
     >            dy*f(3,4,2,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k))))
               p3= 2.0*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 2.0*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if((ict(6).gt.0).or.(ict(1).eq.-1)) then
c  evaluate d2f/dy2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1=
     >            2.0*f(1,3,1,i,j,k)+6.0*dy*f(1,4,1,i,j,k)
     >            +dx*(2.0*f(2,3,1,i,j,k)+6.0*dy*f(2,4,1,i,j,k)
     >            +dx*(2.0*f(3,3,1,i,j,k)+6.0*dy*f(3,4,1,i,j,k)
     >            +dx*(2.0*f(4,3,1,i,j,k)+6.0*dy*f(4,4,1,i,j,k))))
               p2=
     >            2.0*f(1,3,2,i,j,k)+6.0*dy*f(1,4,2,i,j,k)
     >            +dx*(2.0*f(2,3,2,i,j,k)+6.0*dy*f(2,4,2,i,j,k)
     >            +dx*(2.0*f(3,3,2,i,j,k)+6.0*dy*f(3,4,2,i,j,k)
     >            +dx*(2.0*f(4,3,2,i,j,k)+6.0*dy*f(4,4,2,i,j,k))))
               p3=
     >            2.0*f(1,3,3,i,j,k)+6.0*dy*f(1,4,3,i,j,k)
     >            +dx*(2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >            +dx*(2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >            +dx*(2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k))))
               p4=
     >            2.0*f(1,3,4,i,j,k)+6.0*dy*f(1,4,4,i,j,k)
     >            +dx*(2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >            +dx*(2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >            +dx*(2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if((ict(7).gt.0).or.(ict(1).eq.-1)) then
c  evaluate df2/dz2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p3=
     >            f(1,1,3,i,j,k)+dy*(f(1,2,3,i,j,k)+dy*(f(1,3,3,i,j,k)+
     >            dy*f(1,4,3,i,j,k)))
     >        +dx*(f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+
     >            dy*f(2,4,3,i,j,k)))
     >        +dx*(f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k)))
     >        +dx*(f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k)))
     >            )))
               p4=
     >            f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+
     >            dy*f(1,4,4,i,j,k)))
     >        +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            )))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
c
         if((ict(8).gt.0).and.(ict(1).ne.-1)) then
c  evaluate d2f/dxdy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1=
     >     f(2,2,1,i,j,k)+
     >         dy*(2.0*f(2,3,1,i,j,k)+dy*3.0*f(2,4,1,i,j,k))
     >            +2.*dx*(
     >     f(3,2,1,i,j,k)+
     >         dy*(2.0*f(3,3,1,i,j,k)+dy*3.0*f(3,4,1,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,1,i,j,k)+
     >         dy*(2.0*f(4,3,1,i,j,k)+dy*3.0*f(4,4,1,i,j,k))
     >            ))
               p2=
     >     f(2,2,2,i,j,k)+
     >         dy*(2.0*f(2,3,2,i,j,k)+dy*3.0*f(2,4,2,i,j,k))
     >            +2.*dx*(
     >     f(3,2,2,i,j,k)+
     >         dy*(2.0*f(3,3,2,i,j,k)+dy*3.0*f(3,4,2,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,2,i,j,k)+
     >         dy*(2.0*f(4,3,2,i,j,k)+dy*3.0*f(4,4,2,i,j,k))
     >            ))
               p3=
     >     f(2,2,3,i,j,k)+
     >         dy*(2.0*f(2,3,3,i,j,k)+dy*3.0*f(2,4,3,i,j,k))
     >            +2.*dx*(
     >     f(3,2,3,i,j,k)+
     >         dy*(2.0*f(3,3,3,i,j,k)+dy*3.0*f(3,4,3,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,3,i,j,k)+
     >         dy*(2.0*f(4,3,3,i,j,k)+dy*3.0*f(4,4,3,i,j,k))
     >            ))
               p4=
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +2.*dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            ))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if((ict(9).gt.0).and.(ict(1).ne.-1)) then
c  evaluate d2f/dxdz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p2=
     >            f(2,1,2,i,j,k)+dy*(f(2,2,2,i,j,k)+dy*(f(2,3,2,i,j,k)+
     >            dy*f(2,4,2,i,j,k)))
     >        +2.*dx*
     >            (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+
     >            dy*f(3,4,2,i,j,k)))
     >        +1.5*dx*
     >            (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k)))
     >            ))
               p3=
     >            f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+
     >            dy*f(2,4,3,i,j,k)))
     >        +2.*dx*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k)))
     >        +1.5*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k)))
     >            ))
               p4=
     >            f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +2.*dx*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +1.5*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            ))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
c
         if((ict(10).gt.0).and.(ict(1).ne.-1)) then
c  evaluate d2f/dydz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p2=
     >     f(1,2,2,i,j,k)+
     >         dy*(2.0*f(1,3,2,i,j,k)+dy*3.0*f(1,4,2,i,j,k))
     >            +dx*(
     >     f(2,2,2,i,j,k)+
     >         dy*(2.0*f(2,3,2,i,j,k)+dy*3.0*f(2,4,2,i,j,k))
     >            +dx*(
     >     f(3,2,2,i,j,k)+
     >         dy*(2.0*f(3,3,2,i,j,k)+dy*3.0*f(3,4,2,i,j,k))
     >            +dx*(
     >     f(4,2,2,i,j,k)+
     >         dy*(2.0*f(4,3,2,i,j,k)+dy*3.0*f(4,4,2,i,j,k))
     >            )))
               p3=
     >     f(1,2,3,i,j,k)+
     >         dy*(2.0*f(1,3,3,i,j,k)+dy*3.0*f(1,4,3,i,j,k))
     >            +dx*(
     >     f(2,2,3,i,j,k)+
     >         dy*(2.0*f(2,3,3,i,j,k)+dy*3.0*f(2,4,3,i,j,k))
     >            +dx*(
     >     f(3,2,3,i,j,k)+
     >         dy*(2.0*f(3,3,3,i,j,k)+dy*3.0*f(3,4,3,i,j,k))
     >            +dx*(
     >     f(4,2,3,i,j,k)+
     >         dy*(2.0*f(4,3,3,i,j,k)+dy*3.0*f(4,4,3,i,j,k))
     >            )))
               p4=
     >     f(1,2,4,i,j,k)+
     >         dy*(2.0*f(1,3,4,i,j,k)+dy*3.0*f(1,4,4,i,j,k))
     >            +dx*(
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            )))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
c
         if(ict(1).eq.-1) then
c  evaluate d4f/dx2dy2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1=
     >         4.0*f(3,3,1,i,j,k)+12.0*dy*f(3,4,1,i,j,k)
     >         +dx*(12.0*f(4,3,1,i,j,k)+36.0*dy*f(4,4,1,i,j,k))
               p2=
     >         4.0*f(3,3,2,i,j,k)+12.0*dy*f(3,4,2,i,j,k)
     >         +dx*(12.0*f(4,3,2,i,j,k)+36.0*dy*f(4,4,2,i,j,k))
               p3=
     >         4.0*f(3,3,3,i,j,k)+12.0*dy*f(3,4,3,i,j,k)
     >         +dx*(12.0*f(4,3,3,i,j,k)+36.0*dy*f(4,4,3,i,j,k))
               p4=
     >         4.0*f(3,3,4,i,j,k)+12.0*dy*f(3,4,4,i,j,k)
     >         +dx*(12.0*f(4,3,4,i,j,k)+36.0*dy*f(4,4,4,i,j,k))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
c
         if(ict(1).eq.-1) then
c  evaluate d4f/dx2dz2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p3= 2.0*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 2.0*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
c
         if(ict(1).eq.-1) then
c  evaluate d4f/dy2dz2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p3=
     >            2.0*f(1,3,3,i,j,k)+6.0*dy*f(1,4,3,i,j,k)
     >            +dx*(2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >            +dx*(2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >            +dx*(2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k))))
               p4=
     >            2.0*f(1,3,4,i,j,k)+6.0*dy*f(1,4,4,i,j,k)
     >            +dx*(2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >            +dx*(2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >            +dx*(2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
c
         if(ict(1).eq.-1) then
c  evaluate d6f/dx2dy2dz2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p3=
     >         4.0*f(3,3,3,i,j,k)+12.0*dy*f(3,4,3,i,j,k)
     >         +dx*(12.0*f(4,3,3,i,j,k)+36.0*dy*f(4,4,3,i,j,k))
               p4=
     >         4.0*f(3,3,4,i,j,k)+12.0*dy*f(3,4,4,i,j,k)
     >         +dx*(12.0*f(4,3,4,i,j,k)+36.0*dy*f(4,4,4,i,j,k))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
c
      endif
C
C----------------------------------
C  3rd derivatives (.le.2 in each coordinate)
C
      if(ict(1).eq.3) then
         if(ict(2).eq.1) then
C                               ! fxxy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
               p1= 2.0*
     >            (f(3,2,1,i,j,k)+2.0*dy*(f(3,3,1,i,j,k)+
     >            1.5*dy*f(3,4,1,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,1,i,j,k)+2.0*dy*(f(4,3,1,i,j,k)+
     >            1.5*dy*f(4,4,1,i,j,k)))
               p2= 2.0*
     >            (f(3,2,2,i,j,k)+2.0*dy*(f(3,3,2,i,j,k)+
     >            1.5*dy*f(3,4,2,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,2,i,j,k)+2.0*dy*(f(4,3,2,i,j,k)+
     >            1.5*dy*f(4,4,2,i,j,k)))
               p3= 2.0*
     >            (f(3,2,3,i,j,k)+2.0*dy*(f(3,3,3,i,j,k)+
     >            1.5*dy*f(3,4,3,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,3,i,j,k)+2.0*dy*(f(4,3,3,i,j,k)+
     >            1.5*dy*f(4,4,3,i,j,k)))
               p4= 2.0*
     >            (f(3,2,4,i,j,k)+2.0*dy*(f(3,3,4,i,j,k)+
     >            1.5*dy*f(3,4,4,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2= 2.0*
     >            (f(3,1,2,i,j,k)+dy*(f(3,2,2,i,j,k)+dy*(f(3,3,2,i,j,k)+
     >            dy*f(3,4,2,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k))))
               p3= 2.0*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 2.0*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=
     >            2.0*f(2,3,1,i,j,k)+6.0*dy*f(2,4,1,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,1,i,j,k)+6.0*dy*f(3,4,1,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,1,i,j,k)+6.0*dy*f(4,4,1,i,j,k)))
               p2=
     >            2.0*f(2,3,2,i,j,k)+6.0*dy*f(2,4,2,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,2,i,j,k)+6.0*dy*f(3,4,2,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,2,i,j,k)+6.0*dy*f(4,4,2,i,j,k)))
               p3=
     >            2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k)))
               p4=
     >            2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k)))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=
     >     f(2,2,2,i,j,k)+
     >         dy*(2.0*f(2,3,2,i,j,k)+dy*3.0*f(2,4,2,i,j,k))
     >            +2.*dx*(
     >     f(3,2,2,i,j,k)+
     >         dy*(2.0*f(3,3,2,i,j,k)+dy*3.0*f(3,4,2,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,2,i,j,k)+
     >         dy*(2.0*f(4,3,2,i,j,k)+dy*3.0*f(4,4,2,i,j,k))
     >            ))
               p3=
     >     f(2,2,3,i,j,k)+
     >         dy*(2.0*f(2,3,3,i,j,k)+dy*3.0*f(2,4,3,i,j,k))
     >            +2.*dx*(
     >     f(3,2,3,i,j,k)+
     >         dy*(2.0*f(3,3,3,i,j,k)+dy*3.0*f(3,4,3,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,3,i,j,k)+
     >         dy*(2.0*f(4,3,3,i,j,k)+dy*3.0*f(4,4,3,i,j,k))
     >            ))
               p4=
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +2.*dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            ))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >            f(2,1,3,i,j,k)+dy*(f(2,2,3,i,j,k)+dy*(f(2,3,3,i,j,k)+
     >            dy*f(2,4,3,i,j,k)))
     >        +2.*dx*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k)))
     >        +1.5*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k)))
     >            ))
               p4=
     >            f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +2.*dx*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +1.5*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            ))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=
     >            2.0*f(1,3,2,i,j,k)+6.0*dy*f(1,4,2,i,j,k)
     >            +dx*(2.0*f(2,3,2,i,j,k)+6.0*dy*f(2,4,2,i,j,k)
     >            +dx*(2.0*f(3,3,2,i,j,k)+6.0*dy*f(3,4,2,i,j,k)
     >            +dx*(2.0*f(4,3,2,i,j,k)+6.0*dy*f(4,4,2,i,j,k))))
               p3=
     >            2.0*f(1,3,3,i,j,k)+6.0*dy*f(1,4,3,i,j,k)
     >            +dx*(2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >            +dx*(2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >            +dx*(2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k))))
               p4=
     >            2.0*f(1,3,4,i,j,k)+6.0*dy*f(1,4,4,i,j,k)
     >            +dx*(2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >            +dx*(2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >            +dx*(2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >     f(1,2,3,i,j,k)+
     >         dy*(2.0*f(1,3,3,i,j,k)+dy*3.0*f(1,4,3,i,j,k))
     >            +dx*(
     >     f(2,2,3,i,j,k)+
     >         dy*(2.0*f(2,3,3,i,j,k)+dy*3.0*f(2,4,3,i,j,k))
     >            +dx*(
     >     f(3,2,3,i,j,k)+
     >         dy*(2.0*f(3,3,3,i,j,k)+dy*3.0*f(3,4,3,i,j,k))
     >            +dx*(
     >     f(4,2,3,i,j,k)+
     >         dy*(2.0*f(4,3,3,i,j,k)+dy*3.0*f(4,4,3,i,j,k))
     >            )))
               p4=
     >     f(1,2,4,i,j,k)+
     >         dy*(2.0*f(1,3,4,i,j,k)+dy*3.0*f(1,4,4,i,j,k))
     >            +dx*(
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            )))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  3rd derivatives (3 in each coordinate)
C
      if(ict(1).eq.-3) then
         if(ict(2).eq.1) then
C                               ! fxxx
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1= 6.0*
     >            (f(4,1,1,i,j,k)+dy*(f(4,2,1,i,j,k)+dy*(f(4,3,1,i,j,k)+
     >            dy*f(4,4,1,i,j,k))))
               p2= 6.0*
     >            (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k))))
               p3= 6.0*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 6.0*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fyyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=
     >            6.0*(f(1,4,1,i,j,k)
     >            +dx*(f(2,4,1,i,j,k)
     >            +dx*(f(3,4,1,i,j,k)
     >            +dx*f(4,4,1,i,j,k))))
               p2=
     >            6.0*(f(1,4,2,i,j,k)
     >            +dx*(f(2,4,2,i,j,k)
     >            +dx*(f(3,4,2,i,j,k)
     >            +dx*f(4,4,2,i,j,k))))
               p3=
     >            6.0*(f(1,4,3,i,j,k)
     >            +dx*(f(2,4,3,i,j,k)
     >            +dx*(f(3,4,3,i,j,k)
     >            +dx*f(4,4,3,i,j,k))))
               p4=
     >            6.0*(f(1,4,4,i,j,k)
     >            +dx*(f(2,4,4,i,j,k)
     >            +dx*(f(3,4,4,i,j,k)
     >            +dx*f(4,4,4,i,j,k))))
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            f(1,1,4,i,j,k)+dy*(f(1,2,4,i,j,k)+dy*(f(1,3,4,i,j,k)+
     >            dy*f(1,4,4,i,j,k)))
     >        +dx*(f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +dx*(f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +dx*(f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            )))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  4th derivatives (.le.2 in each coordinate)
C
      if(ict(1).eq.4) then
         if(ict(2).eq.1) then
C                               ! fxxyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=
     >            4.0*f(3,3,1,i,j,k)+12.0*dy*f(3,4,1,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,1,i,j,k)+3.0*dy*f(4,4,1,i,j,k))
               p2=
     >            4.0*f(3,3,2,i,j,k)+12.0*dy*f(3,4,2,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,2,i,j,k)+3.0*dy*f(4,4,2,i,j,k))
               p3=
     >            4.0*f(3,3,3,i,j,k)+12.0*dy*f(3,4,3,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,3,i,j,k)+3.0*dy*f(4,4,3,i,j,k))
               p4=
     >            4.0*f(3,3,4,i,j,k)+12.0*dy*f(3,4,4,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2= 2.0*
     >            (f(3,2,2,i,j,k)+2.0*dy*(f(3,3,2,i,j,k)+
     >            1.5*dy*f(3,4,2,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,2,i,j,k)+2.0*dy*(f(4,3,2,i,j,k)+
     >            1.5*dy*f(4,4,2,i,j,k)))
               p3= 2.0*
     >            (f(3,2,3,i,j,k)+2.0*dy*(f(3,3,3,i,j,k)+
     >            1.5*dy*f(3,4,3,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,3,i,j,k)+2.0*dy*(f(4,3,3,i,j,k)+
     >            1.5*dy*f(4,4,3,i,j,k)))
               p4= 2.0*
     >            (f(3,2,4,i,j,k)+2.0*dy*(f(3,3,4,i,j,k)+
     >            1.5*dy*f(3,4,4,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3= 2.0*
     >            (f(3,1,3,i,j,k)+dy*(f(3,2,3,i,j,k)+dy*(f(3,3,3,i,j,k)+
     >            dy*f(3,4,3,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 2.0*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=
     >            2.0*f(2,3,2,i,j,k)+6.0*dy*f(2,4,2,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,2,i,j,k)+6.0*dy*f(3,4,2,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,2,i,j,k)+6.0*dy*f(4,4,2,i,j,k)))
               p3=
     >            2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k)))
               p4=
     >            2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >     f(2,2,3,i,j,k)+
     >         dy*(2.0*f(2,3,3,i,j,k)+dy*3.0*f(2,4,3,i,j,k))
     >            +2.*dx*(
     >     f(3,2,3,i,j,k)+
     >         dy*(2.0*f(3,3,3,i,j,k)+dy*3.0*f(3,4,3,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,3,i,j,k)+
     >         dy*(2.0*f(4,3,3,i,j,k)+dy*3.0*f(4,4,3,i,j,k))
     >            ))
               p4=
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +2.*dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            ))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >            2.0*f(1,3,3,i,j,k)+6.0*dy*f(1,4,3,i,j,k)
     >            +dx*(2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >            +dx*(2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >            +dx*(2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k))))
               p4=
     >            2.0*f(1,3,4,i,j,k)+6.0*dy*f(1,4,4,i,j,k)
     >            +dx*(2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >            +dx*(2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >            +dx*(2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k))))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
      endif
C
C
C----------------------------------
C  4th derivatives (3 in a coordinate)
C
      if(ict(1).eq.-4) then
         if(ict(2).eq.1) then
C                               ! fxxxy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1= 6.0*
     >            (f(4,2,1,i,j,k)+2.0*dy*(f(4,3,1,i,j,k)+
     >            1.5*dy*f(4,4,1,i,j,k)))
               p2= 6.0*
     >            (f(4,2,2,i,j,k)+2.0*dy*(f(4,3,2,i,j,k)+
     >            1.5*dy*f(4,4,2,i,j,k)))
               p3= 6.0*
     >            (f(4,2,3,i,j,k)+2.0*dy*(f(4,3,3,i,j,k)+
     >            1.5*dy*f(4,4,3,i,j,k)))
               p4= 6.0*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2= 6.0*
     >            (f(4,1,2,i,j,k)+dy*(f(4,2,2,i,j,k)+dy*(f(4,3,2,i,j,k)+
     >            dy*f(4,4,2,i,j,k))))
               p3= 6.0*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 6.0*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=
     >            6.0*(f(2,4,1,i,j,k)
     >              +2.0*dx*(f(3,4,1,i,j,k)
     >              +1.5*dx*f(4,4,1,i,j,k)))
               p2=
     >            6.0*(f(2,4,2,i,j,k)
     >              +2.0*dx*(f(3,4,2,i,j,k)
     >              +1.5*dx*f(4,4,2,i,j,k)))
               p3=
     >            6.0*(f(2,4,3,i,j,k)
     >              +2.0*dx*(f(3,4,3,i,j,k)
     >              +1.5*dx*f(4,4,3,i,j,k)))
               p4=
     >            6.0*(f(2,4,4,i,j,k)
     >              +2.0*dx*(f(3,4,4,i,j,k)
     >              +1.5*dx*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            f(2,1,4,i,j,k)+dy*(f(2,2,4,i,j,k)+dy*(f(2,3,4,i,j,k)+
     >            dy*f(2,4,4,i,j,k)))
     >        +2.*dx*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k)))
     >        +1.5*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k)))
     >            ))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fyyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=
     >            6.0*(f(1,4,2,i,j,k)
     >            +dx*(f(2,4,2,i,j,k)
     >            +dx*(f(3,4,2,i,j,k)
     >            +dx*f(4,4,2,i,j,k))))
               p3=
     >            6.0*(f(1,4,3,i,j,k)
     >            +dx*(f(2,4,3,i,j,k)
     >            +dx*(f(3,4,3,i,j,k)
     >            +dx*f(4,4,3,i,j,k))))
               p4=
     >            6.0*(f(1,4,4,i,j,k)
     >            +dx*(f(2,4,4,i,j,k)
     >            +dx*(f(3,4,4,i,j,k)
     >            +dx*f(4,4,4,i,j,k))))
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
C
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >     f(1,2,4,i,j,k)+
     >         dy*(2.0*f(1,3,4,i,j,k)+dy*3.0*f(1,4,4,i,j,k))
     >            +dx*(
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            )))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  5th derivatives (.le.2 in each coordinate)
C
      if(ict(1).eq.5) then
         if(ict(2).eq.1) then
C                               ! fxxyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=
     >            4.0*f(3,3,2,i,j,k)+12.0*dy*f(3,4,2,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,2,i,j,k)+3.0*dy*f(4,4,2,i,j,k))
               p3=
     >            4.0*f(3,3,3,i,j,k)+12.0*dy*f(3,4,3,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,3,i,j,k)+3.0*dy*f(4,4,3,i,j,k))
               p4=
     >            4.0*f(3,3,4,i,j,k)+12.0*dy*f(3,4,4,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3= 2.0*
     >            (f(3,2,3,i,j,k)+2.0*dy*(f(3,3,3,i,j,k)+
     >            1.5*dy*f(3,4,3,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,3,i,j,k)+2.0*dy*(f(4,3,3,i,j,k)+
     >            1.5*dy*f(4,4,3,i,j,k)))
               p4= 2.0*
     >            (f(3,2,4,i,j,k)+2.0*dy*(f(3,3,4,i,j,k)+
     >            1.5*dy*f(3,4,4,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >            2.0*f(2,3,3,i,j,k)+6.0*dy*f(2,4,3,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,3,i,j,k)+6.0*dy*f(3,4,3,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,3,i,j,k)+6.0*dy*f(4,4,3,i,j,k)))
               p4=
     >            2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  5th derivatives (3 in a coordinate)
C
      if(ict(1).eq.-5) then
         if(ict(2).eq.1) then
C                               ! fxxxyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=12.0*(
     >            f(4,3,1,i,j,k)+3.0*dy*f(4,4,1,i,j,k))
               p2=12.0*(
     >            f(4,3,2,i,j,k)+3.0*dy*f(4,4,2,i,j,k))
               p3=12.0*(
     >            f(4,3,3,i,j,k)+3.0*dy*f(4,4,3,i,j,k))
               p4=12.0*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2= 6.0*
     >            (f(4,2,2,i,j,k)+2.0*dy*(f(4,3,2,i,j,k)+
     >            1.5*dy*f(4,4,2,i,j,k)))
               p3= 6.0*
     >            (f(4,2,3,i,j,k)+2.0*dy*(f(4,3,3,i,j,k)+
     >            1.5*dy*f(4,4,3,i,j,k)))
               p4= 6.0*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3= 6.0*
     >            (f(4,1,3,i,j,k)+dy*(f(4,2,3,i,j,k)+dy*(f(4,3,3,i,j,k)+
     >            dy*f(4,4,3,i,j,k))))
               p4= 6.0*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxyyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=12.0*(f(3,4,1,i,j,k)+3.0*dx*f(4,4,1,i,j,k))
               p2=12.0*(f(3,4,2,i,j,k)+3.0*dx*f(4,4,2,i,j,k))
               p3=12.0*(f(3,4,3,i,j,k)+3.0*dx*f(4,4,3,i,j,k))
               p4=12.0*(f(3,4,4,i,j,k)+3.0*dx*f(4,4,4,i,j,k))
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4= 2.0*
     >            (f(3,1,4,i,j,k)+dy*(f(3,2,4,i,j,k)+dy*(f(3,3,4,i,j,k)+
     >            dy*f(3,4,4,i,j,k))))
     >          +6.0*dx*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxyyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=
     >            6.0*(f(2,4,2,i,j,k)
     >              +2.0*dx*(f(3,4,2,i,j,k)
     >              +1.5*dx*f(4,4,2,i,j,k)))
               p3=
     >            6.0*(f(2,4,3,i,j,k)
     >              +2.0*dx*(f(3,4,3,i,j,k)
     >              +1.5*dx*f(4,4,3,i,j,k)))
               p4=
     >            6.0*(f(2,4,4,i,j,k)
     >              +2.0*dx*(f(3,4,4,i,j,k)
     >              +1.5*dx*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fxyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >     f(2,2,4,i,j,k)+
     >         dy*(2.0*f(2,3,4,i,j,k)+dy*3.0*f(2,4,4,i,j,k))
     >            +2.*dx*(
     >     f(3,2,4,i,j,k)+
     >         dy*(2.0*f(3,3,4,i,j,k)+dy*3.0*f(3,4,4,i,j,k))
     >            +1.5*dx*(
     >     f(4,2,4,i,j,k)+
     >         dy*(2.0*f(4,3,4,i,j,k)+dy*3.0*f(4,4,4,i,j,k))
     >            ))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(9).eq.1) then
C                               ! fyyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >            6.0*(f(1,4,3,i,j,k)
     >            +dx*(f(2,4,3,i,j,k)
     >            +dx*(f(3,4,3,i,j,k)
     >            +dx*f(4,4,3,i,j,k))))
               p4=
     >            6.0*(f(1,4,4,i,j,k)
     >            +dx*(f(2,4,4,i,j,k)
     >            +dx*(f(3,4,4,i,j,k)
     >            +dx*f(4,4,4,i,j,k))))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(10).eq.1) then
C                               ! fyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            2.0*f(1,3,4,i,j,k)+6.0*dy*f(1,4,4,i,j,k)
     >            +dx*(2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >            +dx*(2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >            +dx*(2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  6th derivatives (2 in each coordinate)
C
      if(ict(1).eq.6) then
C                               ! fxxyyzz
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            k=kv(v)
            dx=dxv(v)
            dy=dyv(v)
            dz=dzv(v)
C
               p3=
     >            4.0*f(3,3,3,i,j,k)+12.0*dy*f(3,4,3,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,3,i,j,k)+3.0*dy*f(4,4,3,i,j,k))
               p4=
     >            4.0*f(3,3,4,i,j,k)+12.0*dy*f(3,4,4,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
C
            fval(v,iaval)=2.0*p3+6.0*dz*p4
         enddo
      endif
C
C----------------------------------
C  6th derivatives (3 in a coordinate)
C
      if(ict(1).eq.-6) then
         if(ict(2).eq.1) then
C                               ! fxxxyyy
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p1=36.0*f(4,4,1,i,j,k)
               p2=36.0*f(4,4,2,i,j,k)
               p3=36.0*f(4,4,3,i,j,k)
               p4=36.0*f(4,4,4,i,j,k)
C
               fval(v,iaval)=p1+dz*(p2+dz*(p3+dz*p4))
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=12.0*(
     >            f(4,3,2,i,j,k)+3.0*dy*f(4,4,2,i,j,k))
               p3=12.0*(
     >            f(4,3,3,i,j,k)+3.0*dy*f(4,4,3,i,j,k))
               p4=12.0*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3= 6.0*
     >            (f(4,2,3,i,j,k)+2.0*dy*(f(4,3,3,i,j,k)+
     >            1.5*dy*f(4,4,3,i,j,k)))
               p4= 6.0*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxxzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4= 6.0*
     >            (f(4,1,4,i,j,k)+dy*(f(4,2,4,i,j,k)+dy*(f(4,3,4,i,j,k)+
     >            dy*f(4,4,4,i,j,k))))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxyyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p2=12.0*(f(3,4,2,i,j,k)+3.0*dx*f(4,4,2,i,j,k))
               p3=12.0*(f(3,4,3,i,j,k)+3.0*dx*f(4,4,3,i,j,k))
               p4=12.0*(f(3,4,4,i,j,k)+3.0*dx*f(4,4,4,i,j,k))
C
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxxyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4= 2.0*
     >            (f(3,2,4,i,j,k)+2.0*dy*(f(3,3,4,i,j,k)+
     >            1.5*dy*f(3,4,4,i,j,k)))
     >          +6.0*dx*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(8).eq.1) then
C                               ! fxyyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=
     >            6.0*(f(2,4,3,i,j,k)
     >              +2.0*dx*(f(3,4,3,i,j,k)
     >              +1.5*dx*f(4,4,3,i,j,k)))
               p4=
     >            6.0*(f(2,4,4,i,j,k)
     >              +2.0*dx*(f(3,4,4,i,j,k)
     >              +1.5*dx*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(9).eq.1) then
C                               ! fxyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            2.0*f(2,3,4,i,j,k)+6.0*dy*f(2,4,4,i,j,k)
     >              +2.0*dx*(
     >            2.0*f(3,3,4,i,j,k)+6.0*dy*f(3,4,4,i,j,k)
     >              +1.5*dx*(
     >            2.0*f(4,3,4,i,j,k)+6.0*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(10).eq.1) then
C                               ! fyyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            6.0*(f(1,4,4,i,j,k)
     >            +dx*(f(2,4,4,i,j,k)
     >            +dx*(f(3,4,4,i,j,k)
     >            +dx*f(4,4,4,i,j,k))))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  7th derivatives
C
      if(abs(ict(1)).eq.7) then
         if(ict(2).eq.1) then
C                               ! fxxxyyyz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dz=dzv(v)
C
               p2=36.0*f(4,4,2,i,j,k)
               p3=36.0*f(4,4,3,i,j,k)
               p4=36.0*f(4,4,4,i,j,k)
               fval(v,iaval)=p2+dz*(2.0*p3+dz*3.0*p4)
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=12.0*(
     >            f(4,3,3,i,j,k)+3.0*dy*f(4,4,3,i,j,k))
               p4=12.0*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxxyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4= 6.0*
     >            (f(4,2,4,i,j,k)+2.0*dy*(f(4,3,4,i,j,k)+
     >            1.5*dy*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(5).eq.1) then
C                               ! fxxyyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p3=12.0*(f(3,4,3,i,j,k)+3.0*dx*f(4,4,3,i,j,k))
               p4=12.0*(f(3,4,4,i,j,k)+3.0*dx*f(4,4,4,i,j,k))
C
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(6).eq.1) then
C                               ! fxxyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            4.0*f(3,3,4,i,j,k)+12.0*dy*f(3,4,4,i,j,k)
     >              +12.0*dx*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
C
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(7).eq.1) then
C                               ! fxyyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
               dy=dyv(v)
               dz=dzv(v)
C
               p4=
     >            6.0*(f(2,4,4,i,j,k)
     >              +2.0*dx*(f(3,4,4,i,j,k)
     >              +1.5*dx*f(4,4,4,i,j,k)))
C
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  8th derivatives
C
      if(abs(ict(1)).eq.8) then
         if(ict(2).eq.1) then
C                               ! fxxxyyyzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dz=dzv(v)
C
               p3=36.0*f(4,4,3,i,j,k)
               p4=36.0*f(4,4,4,i,j,k)
               fval(v,iaval)=2.0*p3+6.0*dz*p4
            enddo
         endif
C
         if(ict(3).eq.1) then
C                               ! fxxxyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dy=dyv(v)
C
               p4=12.0*(
     >            f(4,3,4,i,j,k)+3.0*dy*f(4,4,4,i,j,k))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
         if(ict(4).eq.1) then
C                               ! fxxyyyzzz
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               j=jv(v)
               k=kv(v)
               dx=dxv(v)
C
               p4=12.0*(f(3,4,4,i,j,k)+3.0*dx*f(4,4,4,i,j,k))
               fval(v,iaval)=6.0*p4
            enddo
         endif
C
      endif
C
C----------------------------------
C  9th derivative
C
      if(abs(ict(1)).eq.9) then
C                               ! fxxxyyyzzz
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            k=kv(v)
C
            p4=36.0*f(4,4,4,i,j,k)
            fval(v,iaval)=6.0*p4
         enddo
      endif
C
      return
      end
 
 
 
