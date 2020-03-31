      subroutine vecpc2(ict,ivec,xvec,yvec,ivd,fval,
     >   nx,xpkg,ny,ypkg,fspl,inf1,
     >   iwarn,ier)
c
c  vectorized piecewise linear evaluation routine -- 2d
c  1.  call vectorized zone lookup routine
c  2.  call vectorized piecewise linear evaluation routine
c
c--------------------------
c  input:
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
      real xvec(ivec)                   ! x-locations at which to evaluate
      real yvec(ivec)                   ! y-locations at which to evaluate
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
      integer nx,ny                     ! dimension of axis grids
      real xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      real ypkg(ny,4)                   ! y grid "package" (cf genxpkg)
      integer inf1                      ! fspl 3rd array dimension, .ge.nx
      real fspl(inf1,ny)                ! Piecewise Linear function
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
      real dxn(ivec)                    ! normalized displacements w/in zones
      real hx(ivec)                     ! h(j) vector
      real hxi(ivec)                    ! 1/h(j) vector
c
      integer iy(ivec)                  ! zone indices {j}
      real dyn(ivec)                    ! normalized displacements w/in zones
      real hy(ivec)                     ! h(j) vector
      real hyi(ivec)                    ! 1/h(j) vector
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
      call xlookup(ivec,xvec,nx,xpkg,2,ix,dxn,hx,hxi,iwarn1)
      call xlookup(ivec,yvec,ny,ypkg,2,iy,dyn,hy,hyi,iwarn2)
      iwarn=iwarn1+iwarn2
c
c  vectorized evaluation
c
      call pc2fcn(ict,ivec,ivd,fval,ix,iy,dxn,dyn,
     >   hx,hxi,hy,hyi,fspl,inf1,ny)
c
      return
      end
