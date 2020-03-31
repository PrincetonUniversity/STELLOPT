      subroutine r8mkbicubw(x,nx,y,ny,f,nf2,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcymin,bcymin,ibcymax,bcymax,
     >   wk,nwk,
     >   ilinx,iliny,ier)
C
C  setup a bicubic spline; store coefficients in compact form
C  (as per suggestion of L. Zakharov, PPPL, Feb. 1999)
C     4 coeffs per grid point:  f,fxx,fyy,fxxyy
C
C
C  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER itest,iadfp,isiz1,iadfw,inwk
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 hxlast,hylast
!============
      integer nx                        ! length of x vector
      integer ny                        ! length of y vector
      REAL*8 x(nx)                        ! x vector, strict ascending
      REAL*8 y(ny)                        ! y vector, strict ascending
C
      integer nf2                       ! 2nd dimension of f, nf2.ge.nx
C  input/output:
C
      REAL*8 f(4,nf2,ny)                  ! data & spline coefficients
C
C  on input:  f(1,i,j) = f(x(i),y(j))
C  on output:  f(1,i,j) unchanged
C              f(2,i,j) = d2f/dx2(x(i),y(j))
C              f(3,i,j) = d2f/dy2(x(i),y(j))
C              f(4,i,j) = d4f/dx2dy2(x(i),y(j))
C
C  and the interpolation formula for (x,y) in (x(i),x(i+1))^(y(j),y(j+1))
C  is:
C        hx = x(i+1)-x(i)   hy = y(j+1)-y(j)
C        dxp= (x-x(i))/hx   dxm= 1-dxp     dxp,dxm in (0,1)
C        dyp= (x-x(i))/hx   dym= 1-dyp     dyp,dym in (0,1)
C        dx3p = dxp**3-dxp  dx3m = dxm**3-dxm     dxp3,dxm3 in (0,1)
C
C   finterp = dxm*(dym*f(1,i,j)+dyp*f(1,i,j+1))
C            +dxp*(dym*f(1,i+1,j)+dyp*f(1,i+1,j+1))
C     +1/6*hx**2*
C            dx3m*(dym*f(2,i,j)+dyp*f(2,i,j+1))
C           +dx3p*(dym*f(2,i+1,j)+dyp*f(2,i+1,j+1))
C     +1/6*hy**2*
C            dxm*(dy3m*f(3,i,j)+dy3p*f(3,i,j+1))
C           +dxp*(dy3m*f(3,i+1,j)+dy3p*f(3,i+1,j+1))
C     +1/36*hx**2*hy**2*
C            dx3m*(dym*f(4,i,j)+dyp*f(4,i,j+1))
C           +dx3p*(dym*f(4,i+1,j)+dyp*f(4,i+1,j+1))
C
C  where the f(2:4,*,*) are cleverly evaluated to assure
C  (a) finterp is continuous and twice differentiable across all
C      grid cell boundaries, and
C  (b) all boundary conditions are satisfied.
C
C  the boundary conditions specification options are:
C  inputs:
C
      integer ibcxmin                   ! bc flag for x=xmin
      REAL*8 bcxmin(ny)                   ! bc data vs. y at x=xmin
      integer ibcxmax                   ! bc flag for x=xmax
      REAL*8 bcxmax(ny)                   ! bc data vs. y at x=xmax
C
      integer ibcymin                   ! bc flag for y=ymin
      REAL*8 bcymin(nx)                   ! bc data vs. x at y=ymin
      integer ibcymax                   ! bc flag for y=ymax
      REAL*8 bcymax(nx)                   ! bc data vs. x at y=ymax
C
C  with interpretation:
c   ibcxmin -- indicator for boundary condition at x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- periodic boundary condition
c     =0 -- use "not a knot"
c     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
c     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
c     =5 -- match 1st derivative to 1st divided difference
c     =6 -- match 2nd derivative to 2nd divided difference
c     =7 -- match 3rd derivative to 3rd divided difference
c           (for more detailed definition of BCs 5-7, see the
c           comments of subroutine mkspline)
c   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
c
c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the BC is periodic.
c
C  and analogous interpretation for ibcymin,bcymin,ibcymax,bcymax
C  (df/dy or d2f/dy2 boundary conditions at y=ymin and y=ymax).
C
C  *** workspace ***
C  This FORTRAN-77 routine requires a workspace of
C    nwk = (AT LEAST) 21*nx*ny
C  *** for dynamic allocation of this workspace use subroutine mkbicub
C    which has the same arguments as mkbicubw, except that the workspace
C    does not have to be provided -- mkbicub is FORTRAN-90.
C
C  input:
      integer nwk                       ! size of workspace
C  modified on output:
      REAL*8 wk(nwk)                      ! workspace
C
C  and output arguments
C
      integer ilinx                     ! =1: x grid is "nearly" equally spaced
      integer iliny                     ! =1: y grid is "nearly" equally spaced
C
C  ilinx and iliny are set to zero if corresponding grids are not equally
C  spaced
C
C  and completion code
C
      integer ier                       ! =0 on exit if there is no error.
C
C  if there is an error, ier is set and a message is output on unit 6.
C  these are considered programming errors in the calling routine.
C
C  possible errors:
C    x(...) not strict ascending
C    y(...) not strict ascending
C    nx .lt. 4
C    ny .lt. 4
C    invalid boundary condition flag
C    workspace too small
C
C-----------------------------------------------------
C  local arrays
C
      integer iselect(10)
      integer zvalues(10)
C
      data iselect/-1,0,0,0,0,0,0,0,0,0/
C
C-----------------------------------------------------
C
      itest=21*nx*ny
      if(nwk.lt.itest) then
         write(6,9901) nwk,itest
 9901    format(' ?mkbicubw:  workspace too small:'/
     >          '  user supplied:  nwk=',i7,'; need at least:  ',i7/
     >      '  nwk = at least 21*nx*ny is required.')
         ier=1
         return
      endif
C
      iadfp=1
      isiz1=16*nx*ny
      iadfw=iadfp+isiz1
      inwk = nwk-isiz1
C
      call r8mkbicop(f,nf2,wk(iadfp),nx,ny)
C
C  evaluate 4x4 continuous bicubic spline
C
      call r8bcspline(x,nx,y,ny,wk(iadfp),nx,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcymin,bcymin,ibcymax,bcymax,
     >   wk(iadfw),inwk,
     >   ilinx,iliny,ier)
C
      if(ier.ne.0) return
C
C  convert to 4-coefficient (2nd/4th partials form)
C
      hxlast=x(nx)-x(nx-1)
      hylast=y(ny)-y(ny-1)
      call r8mkbicon(f,nf2,wk(iadfp),nx,ny,hxlast,hylast)
C
      return
      end
C
C----------------------------------------------------------------
C  mkbicop -- copy spline function input data
C
      subroutine r8mkbicop(fin,nf2,fwk,nx,ny)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx,ny,nf2,iy,ix
!============
      REAL*8 fin(4,nf2,ny)
      REAL*8 fwk(4,4,nx,ny)
C
      do iy=1,ny
         do ix=1,nx
            fwk(1,1,ix,iy)=fin(1,ix,iy)
         enddo
      enddo
C
      return
      end
C----------------------------------------------------------------
C  mkbicon -- create compact spline representation from 4x4
C             (bcspline) representation
C
      subroutine r8mkbicon(fin,nf2,fwk,nx,ny,hxlast,hylast)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx,ny,nf2,iy,ix,iflag,ixuse,iyuse,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 hxlast,hylast,dxuse,dyuse
!============
      REAL*8 fin(4,nf2,ny)
      REAL*8 fwk(4,4,nx,ny)
C
C-----------------------------------------------------
C  local arrays
C
      integer iselect(10)
      REAL*8 zvalues(10)
C
      data iselect/-1,0,0,0,0,0,0,0,0,0/
C
C-----------------------------------------------------
C
      do iy=1,ny
         do ix=1,nx
C
C  copy derivatives from result.  Special treatment needed for end zones
C
            iflag=0
            dxuse=0.0_r8
            dyuse=0.0_r8
            ixuse=ix
            iyuse=iy
            if(ix.eq.nx) then
               iflag=1
               dxuse=hxlast
               ixuse=ix-1
            endif
            if(iy.eq.ny) then
               iflag=1
               dyuse=hylast
               iyuse=iy-1
            endif
C
            if(iflag.eq.1) then
               call r8bcspevfn(iselect,1,1,zvalues,
     >            ixuse,iyuse,dxuse,dyuse,
     >            fwk,nx,ny)
               do j=2,4
                  fin(j,ix,iy)=zvalues(j)
               enddo
            else
               fin(2,ix,iy)=2.0_r8*fwk(3,1,ix,iy)
               fin(3,ix,iy)=2.0_r8*fwk(1,3,ix,iy)
               fin(4,ix,iy)=4.0_r8*fwk(3,3,ix,iy)
            endif
C
         enddo
      enddo
C
      return
      end
