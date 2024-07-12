      subroutine spvec(ict,ivec,xvec,ivd,fval,nx,xpkg,fspl,iwarn,ier)
c
c  vectorized spline evaluation routine -- 1d spline
c  1.  call vectorized zone lookup routine
c  2.  call vectorized spline evaluation routine
c
c--------------------------
c  input:
      integer ict(3)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c        ict(3)=1 for d2f/dx2 (don't evaluate if ict(3)=0)
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (spline values to look up)
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
      integer nx                        ! dimension of spline x grid
      real xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      real fspl(4,nx)                   ! (non-compact) spline coefficients
c
c output:
c condition codes, 0 = normal return
      integer iwarn                     ! =1 if an x value was out of range
      integer ier                       ! =1 if argument error detected
c
c---------------------------------------------------------------
c  local arrays
c
      integer iv(ivec)                  ! zone indices
      real dxv(ivec)                    ! displacements w/in zones
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
      if(nx.lt.2) then
         write(6,*) ' ?spvec:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ivec.le.0) then
         write(6,*) ' ?spvec:  vector dimension .le. 0:  ivec = ',ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?spvec:  output vector dimension less than input ',
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
      call xlookup(ivec,xvec,nx,xpkg,1,iv,dxv,dxv,dxv,iwarn)
c
c  vectorized evaluation
c
      call cspevfn(ict,ivec,ivd,fval,iv,dxv,fspl,nx)
c
      return
      end
