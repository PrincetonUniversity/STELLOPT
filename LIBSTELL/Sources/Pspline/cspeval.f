c  cspeval -- eval cubic spline function and/or derivatives
c
      subroutine cspeval(xget,iselect,fval,x,nx,ilinx,f,ier)
c
      real xget                         ! interpolation target
      real fval(3)                      ! output values
      real x(nx),f(4,nx)                ! spline data
c
      integer iselect(3)                ! output selector
      integer ilinx                     ! =1 if x(...) is evenly spaced
c
c  modification -- dmc 11 Jan 1999 -- remove SAVE stmts; break routine
C    into these parts:
C
C    cspevx -- find grid cell of target pt.
C    cspevfn -- evaluate function using output of bcpsevxy
C
C    in cases where multiple functions are defined on the same grid,
C    time can be saved by using cspevx once and then cspevfn
C    multiple times.
c
c  input:
c     (xget)        location where interpolated value is desired
c                   x(1).le.xget.le.x(nx) expected
c
c     iselect       select desired output
c
c                     iselect(1)=1 -- want function value (f) itself
c                     iselect(2)=1 -- want  df/dx
c                     iselect(3)=1 -- want  d2f/dx2
c
c              example:  iselect(1)=iselect(2)=iselect(3)=1
c                            f, df/dx, and d2f/dx2 are all evaluated
c                            and returned in order in fval(1), fval(2),
c                            and fval(3)
c                        iselect(1)=0, iselect(2)=1, iselect(3)=0
c                            only the 1st derivative is evaluated
c                            and returned in fval(1).
c
c                     set iselect(1)=3 to get d3f/dx3, 1 value only.
c
c                   see fval (output) description.
c
c     x(1...nx)     independent coordinate x, strict ascending
c
c     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
c
c  **CAUTION** actual even spacing of x, is NOT CHECKED HERE!
c
c     f             the function values (at grid points) and spline coefs
c
c  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
c
c      spline value =
c        f(1,i) + dx*f(2,i) + dx**2*f(3,i) + dx**3*f(4,i)
c
c  output:
c      up to 3 elements of fval, ordered as follows:
c        fval(1)=function value or lowest order derivative requested
c        fval(2)=next order derivative
c             etc
c        the ordering is a subset of the sequence given under the "iselect"
c        description.
c
c      ier = 0 -- successful completion; = 1 -- an error occurred.
c
c-------------------------------------------------------------------
c  local
c
      integer :: ia(1) = (/ 0 /)
      real :: dxa(1)
c
c--------------------------
c
      call cspevx(xget,x,nx,ilinx,ia(1),dxa(1),ier)
      if(ier.ne.0) return
c
      call cspevfn(iselect,1,1,fval,ia,dxa,f,nx)
c
      return
      end
c
c-------------------------------------------------------------------------
c  cspevx -- look up x zone
c
c  this is the "first part" of cspeval, see comments, above.
c
      subroutine cspevx(xget,x,nx,ilinx,i,dx,ier)
c
      integer nx                        ! x array dimension
c
      real xget                         ! target point
      real x(nx)                        ! independent coord. array
c
      integer ilinx                     ! =1:  assume x evenly spaced
c
c  output of cspevx
c
      integer i                         ! index to cell containing target pt
      real dx                           ! displacement of target pt w/in cell
                                        ! dx = x-x(i)
c
c  the input argument range is checked...
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
      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0e-7*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
            write(6,1001) xget,x(1),x(nx)
 1001       format(' ?cspeval:  xget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((xget.lt.x(1)-0.5*zxtol).or.
     >         (xget.gt.x(nx)+0.5*zxtol))
     >      write(6,1011) xget,x(1),x(nx)
 1011       format(' %cspeval:  xget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(xget.lt.x(1)) then
               zxget=x(1)
            else
               zxget=x(nx)
            endif
         endif
      endif
      if(ier.ne.0) return
c
c  now find interval in which target point lies..
c
      nxm=nx-1
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
      dx=zxget-x(i)
c
      return
      end
c------------------------------------------------------------------------
c  cspevfn -- OK now evaluate the cubic spline
c
      subroutine cspevfn(ict,ivec,ivd,fval,iv,dxv,f,nx)
c
c  input:
      integer ict(3)                    ! selector:
c        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
c        ict(2)=1 for df/dx   ""
c        ict(3)=1 for d2f/dx2
c
c        set ict(1)=3 to get d3f/dx3 (only)
c
      integer ivec,ivd                  ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c    ivd -- 1st dimension of fval, .ge.ivec
c
c output:
      real fval(ivd,*)                 ! output array
c
c    v = index to element in vector;
c  fval(v,1) = first item requested by ict(...),
c  fval(v,2) = 2nd item requested,  ...etc...
c
c  input:
      integer iv(ivec)                  ! grid cell indices -- vectors
      real dxv(ivec)                    ! displacements w/in cell -- vectors
c
      real f(4,nx)                      ! cubic fcn spline coeffs array
c
c  usage example:
c
c  1.  for each element xx(v) in a vector of x values:
c    find the x zone index and displacement with respect to the
c    lower end point of the zone; store thes in vectors iv and dxv.
c
c  2.  set ict(1)=0, ict(2)=1, ict(3)=0 -- get only the 1st derivative
c
c  3.  ivec is the length of the vector; ivd is the 1st dimension of the
c      array fval to receive the output
c
c      real fval(ivd,3)
c      real xv(ivd)
c      real dxv(ivd)
c      integer iv(ivd)
c      integer ict(3)
c
c      real fspline(4,nx)  ! spline coeffs
c      data ict/0,1,0/     ! this call:  want 1st derivative only
c                          !  this will be output to fval(*,1).
c      ...
c      do iv=1,ivec
c        ...               ! comput indices & displacements
c      enddo
c      call cspevfn(ict,ivec,ivd,fval,iv,dxv,fspline,nx)
c
c--------------------------------------------------------------------
c  local:
c
      integer v                         ! vector element index
c
c  OK can now do evaluations
c
      iaval=0  ! fval addressing
c
      if(ict(1).eq.3) then
c
c  fxxx = d3f/dx3
c
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            fval(v,iaval)=6.0*f(4,i)
         enddo
      else
c
c  normal call...
c
         if(ict(1).gt.0) then
c  evaluate f
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               dx=dxv(v)
               fval(v,iaval)=f(1,i)+dx*(f(2,i)+dx*(f(3,i)+dx*f(4,i)))
            enddo
         endif
c
         if(ict(2).gt.0) then
c  evaluate df/dx
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               dx=dxv(v)
               fval(v,iaval)=f(2,i)+dx*(2.0*f(3,i)+dx*3.0*f(4,i))
            enddo
         endif
c
         if(ict(3).gt.0) then
c  evaluate d2f/dx2
            iaval=iaval+1
            do v=1,ivec
               i=iv(v)
               dx=dxv(v)
               fval(v,iaval)=2.0*f(3,i)+dx*6.0*f(4,i)
            enddo
         endif
      endif
c
      return
      end
c----------------------
