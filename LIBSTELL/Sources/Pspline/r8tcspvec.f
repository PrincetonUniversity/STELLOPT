      subroutine r8tcspvec(ict,ivec,xvec,yvec,zvec,ivd,fval,
     >   nx,xpkg,ny,ypkg,nz,zpkg,fspl,inf4,inf5,
     >   iwarn,ier)
c
c  vectorized spline evaluation routine -- 3d spline
c  1.  call vectorized zone lookup routine
c  2.  call vectorized spline evaluation routine
c
c--------------------------
c  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iwarn1,iwarn2,iwarn3
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 stat
!============
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
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
c
c  list of (x,y,z) triples:
c
      REAL*8 xvec(ivec)                   ! x-locations at which to evaluate
      REAL*8 yvec(ivec)                   ! y-locations at which to evaluate
      REAL*8 zvec(ivec)                   ! z-locations at which to evaluate
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
      integer nx,ny,nz                  ! dimension of spline grids
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
      REAL*8 zpkg(nz,4)                   ! z grid "package" (cf genxpkg)
      integer inf4                      ! fspl 4th array dimension, .ge.nx
      integer inf5                      ! fspl 5th array dimension, .ge.ny
      REAL*8 fspl(4,4,4,inf4,inf5,nz)     ! (non-compact) spline coefficients
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer, dimension(:), allocatable :: ix,iy,iz
      REAL*8, dimension(:), allocatable :: dxv,dyv,dzv
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
c
      if(nx.lt.2) then
         write(6,*) ' ?tcspvec:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ny.lt.2) then
         write(6,*) ' ?tcspvec:  ny.lt.2:  ny = ',ny
         ier=1
      endif
c
      if(nz.lt.2) then
         write(6,*) ' ?tcspvec:  nz.lt.2:  nz = ',nz
         ier=1
      endif
c
      if(ivec.le.0) then
         write(6,*) ' ?tcspvec:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?tcspvec:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
      allocate(ix(ivec), iy(ivec), iz(ivec),
     >   dxv(ivec), dyv(ivec), dzv(ivec), stat=ier)
c
      if(ier.ne.0) then
         write(6,*)
     >      ' ?tcspvec: memory allocation failure.'
         ier=99
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookups
c
      ix=0; iy=0; iz=0
      call r8xlookup(ivec,xvec,nx,xpkg,1,ix,dxv,dxv,dxv,iwarn1)
      call r8xlookup(ivec,yvec,ny,ypkg,1,iy,dyv,dyv,dyv,iwarn2)
      call r8xlookup(ivec,zvec,nz,zpkg,1,iz,dzv,dzv,dzv,iwarn3)
      iwarn=max(iwarn1,iwarn2,iwarn3)
c
c  vectorized evaluation
c
      call r8tcspevfn(ict,ivec,ivd,fval,ix,iy,iz,dxv,dyv,dzv,
     >   fspl,inf4,inf5,nz)
c
      deallocate(ix,iy,iz,dxv,dyv,dzv)
c
      return
      end
