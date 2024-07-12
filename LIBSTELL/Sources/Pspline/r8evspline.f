      subroutine r8evspline(xget,x,nx,ilinx,f,ict,fval,ier)
C
C  use mkspline to set up spline coefficients...
C
C  evaluate a 1d cubic Spline interpolant -- this is C2
C
C  this subroutine calls two subroutines:
C     herm1x  -- find cell containing (xget)
C     fvspline  -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer nx                        ! grid size
      REAL*8 xget                         ! target of this interpolation
      REAL*8 x(nx)                        ! ordered x grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
C
      REAL*8 f(0:1,nx)                    ! function data
C
C  f(0,i) = f @ x(i)
C  f(1,i) = d2f/dx2 @ x(i)
C
C      (these are spline coefficients selected for continuous 2-
C      diffentiability, see mkspline.for)
C
      integer ict(3)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return d2f/dx2  (0, don't)
c
c        set ict(1)=3 to get d3f/dx3 (only)
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
C
C  examples:
C    on input ict = [1,1,1]
C   on output fval= [f,df/dx,d2f/dx2]
C
C    on input ict = [1,0,0]
C   on output fval= [f] ... elements 2 -- 3 never referenced.
C
C    on input ict = [0,0,1]
C   on output fval= [d2f/dx2] ... elements 2 -- 3 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i(1)                         ! cell indices
C
C  normalized displacement from x(i) grid point.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C
      REAL*8 xparam(1)
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i))
C
      REAL*8 hx(1)
      REAL*8 hxi(1)
C
C  0 .le. xparam .le. 1
C
C  ** the interface is very similar to herm2ev.for; can use herm2xy **
C---------------------------------------------------------------------
C
      i(1)=0
      call r8herm1x(xget,x,nx,ilinx,i(1),xparam(1),hx(1),hxi(1),ier)
      if(ier.ne.0) return
c
      call r8fvspline(ict,1,1,fval,i,xparam,hx,hxi,f,nx)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 2d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine r8fvspline(ict,ivec,ivecd,
     >   fval,ii,xparam,hx,hxi,fin,nx)
C
C  use mkspline to set up spline coefficients...
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx,iadr,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xp,xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi
!============
      integer ict(3)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec)                  ! target cells (i)
      REAL*8 xparam(ivec)
                          ! normalized displacements from (i) corners
C
      REAL*8 hx(ivec)                     ! grid spacing, and
      REAL*8 hxi(ivec)                    ! inverse grid spacing 1/(x(i+1)-x(i))
C
      REAL*8 fin(0:1,nx)                  ! interpolant data (cf "evspline")
C
      REAL*8 fval(ivecd,*)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine
C  evspline comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii and parameter inputs
C     xparam,hx,hxi, are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C  Spline evaluation consists of a "mixing" of the interpolant
C  data using the linear functionals xparam, xpi = 1-xparam,
C  xparam**3-xparam, xpi**3-xpi ...
C  and their derivatives as needed.
C
      integer v
C
      REAL*8 sum
      REAL*8 sixth
C
      data sixth/0.166666666666666667_r8/
C
C---------------
C   ...in x direction
C
      iadr=0
C
      if(ict(1).le.2) then
C
C  normal call
C
         if(ict(1).eq.1) then
C
C  get desired values:
C
            iadr = iadr + 1
            do v=1,ivec
               i=ii(v)
C
               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cx=xp*(xp2-1.0_r8)
               cxi=xpi*(xpi2-1.0_r8)
               hx2=hx(v)*hx(v)
C
C  function value:
C
               sum=xpi*fin(0,i) +xp*fin(0,i+1)
               sum=sum+sixth*hx2*(cxi*fin(1,i) +cx*fin(1,i+1))
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
C
               xp=xparam(v)
               xpi=1.0_r8-xp
               xp2=xp*xp
               xpi2=xpi*xpi
C
               cxd=3.0_r8*xp2-1.0_r8
               cxdi=-3.0_r8*xpi2+1.0_r8
 
               sum=hxi(v)*(fin(0,i+1)-fin(0,i))
               sum=sum+sixth*hx(v)*(cxdi*fin(1,i) +cxd*fin(1,i+1))
C
               fval(v,iadr)=sum
            enddo
         endif
C
         if(ict(3).eq.1) then
C
C  d2f/dx2:
C
            iadr=iadr+1
            do v=1,ivec
               i=ii(v)
C
               xp=xparam(v)
               xpi=1.0_r8-xp
C
               sum=xpi*fin(1,i) +xp*fin(1,i+1)
               fval(v,iadr)=sum
            enddo
         endif
C
      else
C
C  return fxxx = d3f/dx3
C
         iadr=1
         do v=1,ivec
            i=ii(v)
            fval(v,iadr)=hxi(v)*(fin(1,i+1)-fin(1,i))
         enddo
C
      endif
C
      return
      end
