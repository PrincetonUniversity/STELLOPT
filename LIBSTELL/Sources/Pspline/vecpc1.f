      subroutine vecpc1(ict,ivec,xvec,ivd,fval,nx,xpkg,fspl,
     >   iwarn,ier)
c
c  vectorized piecewise linear evaluation routine -- 1d
c  1.  call vectorized zone lookup routine
c  2.  call vectorized piecewise linear evaluation routine
c
c     note -- derivatives are *not* continuous
c
c--------------------------
c  input:
      integer ict(2)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (piecewise linear values to look up)
c
      real xvec(ivec)                   ! x-locations at which to evaluate
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
      integer nx                        ! dimension of x grid
      real xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      real fspl(nx)                     ! Piecewise Linear function
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer iv(ivec)                  ! zone indices {j}
      real dxn(ivec)                    ! normalized displacements w/in zones
      real h(ivec)                      ! h(j) vector
      real hi(ivec)                     ! 1/h(j) vector
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
c
      if(ivec.le.0) then
         write(6,*) ' ?vecpc1:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecpc1:  output vector dimension less than input ',
     >      'vector dimension.'
         write(6,*) ' ivec=',ivec,' ivd=',ivd
         ier=1
      endif
c
      if(ier.ne.0) return
c
c  vectorized lookup
c
      iv=0
      call xlookup(ivec,xvec,nx,xpkg,2,iv,dxn,h,hi,iwarn)
c
c  vectorized evaluation
c
      call pc1fcn(ict,ivec,ivd,fval,iv,dxn,h,hi,fspl,nx)
c
      return
      end
