      subroutine vecintrp3d(ict,ivec,xvec,yvec,zvec,ivd,fval,
     >   nx,xpkg,ny,ypkg,nz,zpkg,jspline,fspl,icoeff,ixdim,iydim,izdim,
     >   iwarn,ier)
c
      implicit NONE
c
c  vectorized hybrid spline evaluation routine -- 3d
c  1.  call vectorized zone lookup routine
c  2.  call vectorized hybrid spline evaluation routine
c
c--------------------------
c  input:
      integer ict(10)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
c        ict(4)=1 for df/dy  (don't evaluate if ict(4)=0)
c        ict(5)=1 for d2f/dx2 (don't evaluate if ict(5)=0)
c        ict(6)=1 for d2f/dy2 (don't evaluate if ict(6)=0)
c        ict(7)=1 for d2f/dz2 (don't evaluate if ict(7)=0)
c        ict(8)=1 for d2f/dxdy (don't evaluate if ict(8)=0)
c        ict(9)=1 for d2f/dxdz (don't evaluate if ict(9)=0)
c        ict(10)=1 for d2f/dydz (don't evaluate if ict(10)=0)
c
c     higher derivatives can be selected by setting ict(1)=3,-3,4,-4,...
c     see fvtricub.for comments.
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
c  list of (x,y,z) triples:
c
      real xvec(ivec)                   ! x-locations at which to evaluate
      real yvec(ivec)                   ! y-locations at which to evaluate
      real zvec(ivec)                   ! z-locations at which to evaluate
c
      integer ivd                       ! 1st dimension of output array
c
c    ivd -- 1st dimension of fval, .ge.ivec
c
c output:
      real fval(ivd,*)                  ! output array
c
c  fval(1:ivec,1) -- values as per 1st non-zero ict(...) element
c  fval(1:ivec,2) -- values as per 2nd non-zero ict(...) element
c   --etc--
c
c input:
      integer nx,ny,nz                  ! dimension of spline grids
      real xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      real ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
      real zpkg(nz,4)                   ! z grid "package" (cf genxpkg)
c
      integer :: jspline(3)             ! interpolation method on each dim:
c        jspline(1) for x; jspline(2) for y; jspline(3) for z
c        -1: zonal step fcn; 0: pcwise linear; 1: Hermite; 2: compact spline
c
      integer icoeff                    ! #coeffs per data node
      integer ixdim                     ! x dim for fspl
      integer iydim                     ! y dim for fspl
      integer izdim                     ! z dim for fspl
      real fspl(icoeff,ixdim,iydim,izdim)     ! hybrid spline coefficients
c
c  ixdim=nx-1 for zonal step function along x dimension; o.w. expect ixdim=nx
c  iydim=ny-1 for zonal step function along y dimension; o.w. expect iydim=ny
c  izdim=nz-1 for zonal step function along z dimension; o.w. expect izdim=nz
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer :: iwarn1,iwarn2,iwarn3
c
      integer, dimension(:), allocatable :: ix,iy,iz
      real, dimension(:), allocatable :: dxn,dyn,dzn
      real, dimension(:), allocatable :: hx,hxi,hy,hyi,hz,hzi
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
c
      if(nx.lt.2) then
         write(6,*) ' ?vecintrp3d:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ny.lt.2) then
         write(6,*) ' ?vecintrp3d:  ny.lt.2:  ny = ',ny
         ier=1
      endif
c
      if(nz.lt.2) then
         write(6,*) ' ?vecintrp3d:  nz.lt.2:  nz = ',nz
         ier=1
      endif
c
      call vecin3d_argchk('vecintrp3d',jspline,
     >     icoeff,nx,ny,nz,ixdim,iydim,izdim,ier)
c
      if(ivec.le.0) then
         write(6,*) ' ?vecintrp3d:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecintrp3d:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
      allocate(ix(ivec), iy(ivec), iz(ivec),
     >   dxn(ivec), dyn(ivec), dzn(ivec),
     >   hx(ivec),  hy(ivec),  hz(ivec),
     >   hxi(ivec), hyi(ivec), hzi(ivec), stat=ier)
c
      if(ier.ne.0) then
         write(6,*) 
     >      ' ?vecintrp3d: memory allocation failure.'
         ier=99
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookup
c
      ix=0; iy=0; iz=0
      call xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
      call xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
      call xlookup(ivec,zvec,nz,zpkg,2,iz,dzn,hz,hzi,iwarn3)
      iwarn=max(iwarn1,iwarn2,iwarn3)
c
c  vectorized evaluation
c
      call fvintrp3d(ict,ivec,ivd,fval,ix,iy,iz,dxn,dyn,dzn,
     >   hx,hxi,hy,hyi,hz,hzi,jspline,fspl,icoeff,ixdim,iydim,izdim)
c
      deallocate(ix,iy,iz,dxn,dyn,dzn,hx,hy,hz,hxi,hyi,hzi)
c
      return
      end
