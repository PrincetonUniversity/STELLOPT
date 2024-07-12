      subroutine r8vecherm1(ict,ivec,xvec,ivd,fval,nx,xpkg,fspl,
     >   iwarn,ier)
c
c  vectorized hermite evaluation routine -- 1d
c  1.  call vectorized zone lookup routine
c  2.  call vectorized hermite evaluation routine
c
c--------------------------
c  input:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer ict(2)                    ! selector:
c        ict(1)=1 for f      (don't evaluate if ict(1)=0)
c        ict(2)=1 for df/dx  (don't evaluate if ict(2)=0)
c
      integer ivec                      ! vector dimensioning
c
c    ivec-- number of vector pts (hermite values to look up)
c
      REAL*8 xvec(ivec)                   ! x-locations at which to evaluate
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
      integer nx                        ! dimension of hermite x grid
      REAL*8 xpkg(nx,4)                   ! x grid "package" (cf genxpkg)
      REAL*8 fspl(2,nx)                   ! Hermite coefficients (cf herm1ev)
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
      REAL*8 dxn(ivec)                    ! normalized displacements w/in zones
      REAL*8 h(ivec)                      ! h(j) vector
      REAL*8 hi(ivec)                     ! 1/h(j) vector
c
c---------------------------------------------------------------
c
c  error checks
c
      ier=0
      if(nx.lt.2) then
         write(6,*) ' ?vecherm1:  nx.lt.2:  nx = ',nx
         ier=1
      endif
c
      if(ivec.le.0) then
         write(6,*) ' ?vecherm1:  vector dimension .le. 0:  ivec = ',
     >      ivec
         ier=1
      endif
c
      if(ivd.lt.ivec) then
         write(6,*)
     >      ' ?vecherm1:  output vector dimension less than input ',
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
      call r8xlookup(ivec,xvec,nx,xpkg,2,iv,dxn,h,hi,iwarn)
c
c  vectorized evaluation
c
      call r8herm1fcn(ict,ivec,ivd,fval,iv,dxn,h,hi,fspl,nx)
c
      return
      end
