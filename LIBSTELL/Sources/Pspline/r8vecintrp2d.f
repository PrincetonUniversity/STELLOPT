      subroutine r8vecintrp2d(ict,ivec,xvec,yvec,ivd,fval,
     >   nx,xpkg,ny,ypkg,jspline,fspl,icoeff,ixdim,iydim,
     >   iwarn,ier)
c
c
c  vectorized hybrid spline evaluation routine -- 2d
c  1.  call vectorized zone lookup routine
c  2.  call vectorized hybrid spline evaluation routine
c
c--------------------------
c  input:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer ict(6)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
c        ict(4)=1 for d2f/dx2 (don't evaluate if ict(4)=0)
c        ict(5)=1 for d2f/dy2 (don't evaluate if ict(5)=0)
c        ict(6)=1 for d2f/dxdy (don't evaluate if ict(6)=0)
c
c        ict(1)=3 followed by ict(2:5) = 1 or 0 allow 3rd derivatives to
c          be selected:  fxxx fxxy fxyy fyyy
c
c        ict(1)=4 followed by ict(2:4) allow 4th derivatives to
c          be selected:  fxxxy fxxyy fxyyy; fxxxx=fyyyy=0
c
c        ict(1)=5 followed by ict(2:4) allow 4th derivatives to
c          be selected:  fxxxyy fxxyyy
c
c        ict(1)=6 specifies 6th derivative:  fxxxyyy (step fcn)
c
c     in hybrid spline evaluation, any derivative along a dimension
c     with zonal interpolation (jspline(i)=-1) gives zero;
c
c     piecewise linear and Hermite interpolation give zero if a 2nd or
c     higher derivative is sought along any dimension.
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c
c  list of (x,y) pairs:
c
      REAL*8 xvec(ivec)                   ! x-locations at which to evaluate
      REAL*8 yvec(ivec)                   ! y-locations at which to evaluate
c
      integer ivd                       ! 1st dimension of output array
c
c    ivd -- 1st dimension of fval, .ge.ivec
c
c output:
      REAL*8 fval(ivd,*)                  ! output array
c
c  fval(1:ivec,1) -- values as per 1st non-zero ict(...) element
c  fval(1:ivec,2) -- values as per 2nd non-zero ict(...) element
c   --etc--
c
c input:
      integer nx,ny                     ! dimension of spline grids
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
c
      integer :: jspline(2)             ! interpolation method on each dim:
c        jspline(1) for x; jspline(2) for y
c        -1: zonal step fcn; 0: pcwise linear; 1: Hermite; 2: compact spline
c
      integer icoeff                    ! #coeffs per data node
      integer ixdim                     ! x dim for fspl
      integer iydim                     ! y dim for fspl
      REAL*8 fspl(icoeff,ixdim,iydim)     ! hybrid spline coefficients
c
c  ixdim=nx-1 for zonal step function along x dimension; o.w. expect ixdim=nx
c  iydim=ny-1 for zonal step function along y dimension; o.w. expect iydim=ny
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local variables and arrays
c
      integer :: iwarn1,iwarn2
c
      integer ix(ivec)                  ! zone indices {j}
      REAL*8 dxn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hx(ivec)                     ! h(j) vector
      REAL*8 hxi(ivec)                    ! 1/h(j) vector
c
      integer iy(ivec)                  ! zone indices {j}
      REAL*8 dyn(ivec)                    ! normalized displacements w/in zones
      REAL*8 hy(ivec)                     ! h(j) vector
      REAL*8 hyi(ivec)                    ! 1/h(j) vector
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
c
      if(nx.lt.2) then
         write(6,*) ' ?vecintrp2d:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ny.lt.2) then
         write(6,*) ' ?vecintrp2d:  ny.lt.2:  ny = ',ny
         ier=1
      endif
c
      call vecin2d_argchk('vecintrp2d',jspline,
     >     icoeff,nx,ny,ixdim,iydim,ier)
 
      if(ivec.le.0) then
         write(6,*) ' ?vecintrp2d:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecintrp2d:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookup
c
      ix=0
      iy=0
      call r8xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
      call r8xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
      iwarn=iwarn1+iwarn2
c
c  vectorized evaluation
c
      call r8fvintrp2d(ict,ivec,ivd,fval,ix,iy,dxn,dyn,
     >   hx,hxi,hy,hyi,jspline,fspl,icoeff,ixdim,iydim)
c
      return
      end
