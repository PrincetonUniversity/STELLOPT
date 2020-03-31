      subroutine r8vecpc2(ict,ivec,xvec,yvec,ivd,fval,
     >   nx,xpkg,ny,ypkg,fspl,inf1,
     >   iwarn,ier)
c
c  vectorized piecewise linear evaluation routine -- 2d
c  1.  call vectorized zone lookup routine
c  2.  call vectorized piecewise linear evaluation routine
c
c--------------------------
c  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iwarn1,iwarn2
!============
      integer ict(4)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for df/dy  (don't evaluate if ict(3)=0)
c        ict(4)=1 for d2f/dxdy (don't evaluate if ict(4)=0)
c
c  note -- derivatives are *not* continuous!
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (piecewise linear values to look up)
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
      integer nx,ny                     ! dimension of axis grids
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
      integer inf1                      ! fspl 3rd array dimension, .ge.nx
      REAL*8 fspl(inf1,ny)                ! Piecewise Linear function
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
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
      if(ivec.le.0) then
         write(6,*) ' ?vecpc2:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecpc2:  output vector dimension less than input ',
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
      call r8pc2fcn(ict,ivec,ivd,fval,ix,iy,dxn,dyn,
     >   hx,hxi,hy,hyi,fspl,inf1,ny)
c
      return
      end
