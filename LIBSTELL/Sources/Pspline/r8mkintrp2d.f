      subroutine r8mkintrp2d(x,nx,y,ny,jspline,
     >     f,icoeff,ixdim,iydim,
     >     ibcxmin,bcxmin,ibcxmax,bcxmax,
     >     ibcymin,bcymin,ibcymax,bcymax,
     >     ier)
C
C  setup bicubic spline, or bicubic hermite, or Hybrid linear/zonal with
C  1d cubic Hermite or Spline interpolation in 1 dimension
C
C
C  input:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer nx                        ! length of x vector
      integer ny                        ! length of y vector
      REAL*8 x(nx)                        ! x vector, strict ascending
      REAL*8 y(ny)                        ! y vector, strict ascending
C
      integer :: jspline(2)             ! interpolation method control
C        (1) -- 1st dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
C        (2) -- 2nd dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
C
C    Standard interpolation-- all jspline values match
C      e.g. jspline(1)=jspline(2)=2 for bicubic spline
C    Hybrid interpolation occurs when multiple jspline values are specified.
C
C    RESTRICTION:  jspline = (/ 1, 2 /) or (/ 2, 1 /) not allowed.  I.e.
C    Hermite and Spline interpolation cannot be mixed in a hybrid interpolant.
C    This restriction exists because of technical issues in the
C    implementation (it could be removed in principle but the work to do
C    this has not been scheduled).
C
C  coefficient buffer dimensioning:
      integer :: icoeff                 ! #coefficients per data point
      integer :: ixdim                  ! nx; nx-1 if jspline(1)==-1
      integer :: iydim                  ! ny; ny-1 if jspline(2)==-1
C  input/output:
      REAL*8 f(icoeff,ixdim,iydim)        ! data & spline coefficients
C
C  input bdy condition data:
C   meaning of arguments as described in "mkbicubw.for"
C   *ignored* for dimensions with piecewise linear or zonal interpolation;
C   For Hermite interpolation only the values {-1,0,1} are accepted.
C
      integer ibcxmin                   ! bc flag for x=xmin
      REAL*8 bcxmin(iydim)                ! bc data vs. y at x=xmin
      integer ibcxmax                   ! bc flag for x=xmax
      REAL*8 bcxmax(iydim)                ! bc data vs. y at x=xmax
C
      integer ibcymin                   ! bc flag for y=ymin
      REAL*8 bcymin(ixdim)                ! bc data vs. x at y=ymin
      integer ibcymax                   ! bc flag for y=ymax
      REAL*8 bcymax(ixdim)                ! bc data vs. x at y=ymax
C
      integer ier                       ! =0 on exit if there is no error.
c   ier -- completion code, 0 for normal
C
C-----------------------
      integer :: kspline
      integer :: ii,imul,imin,imax,ickx,icky
      integer :: idum1,idum2
      integer :: ipx,ipy
      REAL*8 :: ztol = 0.0001_r8
      REAL*8, dimension(:,:), allocatable :: wk2
C-----------------------
c
      ier = 0
c
      imin=3
      imax=-2
      imul=1
      do ii=1,2
         imin=min(imin,jspline(ii))
         imax=max(imax,jspline(ii))
         if(jspline(ii).gt.0) imul=imul*2
      enddo
 
      if((imin.lt.-1).or.(imax.gt.2)) then
         ier = 1
         write(6,*)
     >        ' ?mkintrp2d: spline type control out of range -1 to 2: ',
     >        jspline
      endif
      if(ier.ne.0) return
 
      if(imin.eq.imax) then
         kspline=imin   ! same interp type on all dimensions
      else
         kspline=-99    ! hybrid
         if(imin.gt.0) then
            ier = 1
            write(6,*)
     >           ' ?mkintrp2d: spline/Hermite hybrid not supported (',
     >           jspline,')'
         else if(imax.lt.1) then
            ier = 1
            write(6,*)
     >           ' ?mkintrp2d: zonal/linear hybrid not supported (',
     >           jspline,')'
         endif
      endif
      if(ier.ne.0) return
c
      if(imul.ne.icoeff) then
         write(6,*)
     >        ' ?coeff dimension inconsistency for spline type codes ',
     >        jspline
         write(6,*) ' in mkintrp2d: expected: ',imul,' got: ',icoeff
         ier=1
         return
      endif
c
c  check dimensioning consistency
c
      if(jspline(1).eq.-1) then
         ickx=nx-1
      else
         ickx=nx
      endif
c
      if(jspline(2).eq.-1) then
         icky=ny-1
      else
         icky=ny
      endif
c
      if((ickx.ne.ixdim).or.(icky.ne.iydim)) then
         write(6,*)
     >        ' ?mkintrp2d: dimensioning inconsistent with '//
     >        'interpolation controls: ',jspline
         write(6,*) '  expected: ',ickx,icky,'; got: ',ixdim,iydim
         ier=1
         return
      endif
c
      if(jspline(1).le.0) then
         call r8splinck(x,nx,idum1,ztol,ier)
         if(ier.ne.0) then
            write(6,*) ' ?mkintrp2d: x axis not strict ascending.'
            return
         endif
      endif
c
      if(jspline(2).le.0) then
         call r8splinck(y,ny,idum1,ztol,ier)
         if(ier.ne.0) then
            write(6,*) ' ?mkintrp2d: y axis not strict ascending.'
            return
         endif
      endif
c
c  if no work to be done: exit now
      if(imul.eq.1) return
c
c  check Hermite BCs if necessary
c
      if(jspline(1).eq.1) then
         if((min(ibcxmin,ibcxmax).lt.-1).or.
     >        (max(ibcxmin,ibcxmax).gt.1)) then
            write(6,*) ' ?mkintrp2d: Bdy Cond code out of range for'
            write(6,*) '  Hermite interpolation; (-1:1) allowed, '//
     >           'found: ',ibcxmin,ibcxmax
            ier=1
            return
         endif
         ipx=0
         if(ibcxmin.eq.-1) then
            ipx=1
         else if((ibcxmin.eq.1).or.(ibcxmax.eq.1)) then
            ipx=2
         endif
      endif
c
      if(jspline(2).eq.1) then
         if((min(ibcymin,ibcymax).lt.-1).or.
     >        (max(ibcymin,ibcymax).gt.1)) then
            write(6,*) ' ?mkintrp2d: Bdy Cond code out of range for'
            write(6,*) '  Hermite interpolation; (-1:1) allowed, '//
     >           'found: ',ibcymin,ibcymax
            ier=1
            return
         endif
         ipy=0
         if(ibcymin.eq.-1) then
            ipy=1
         else if((ibcymin.eq.1).or.(ibcymax.eq.1)) then
            ipy=2
         endif
      endif
 
      if(kspline.eq.1) then
         ! bicubic Hermite
 
         ! stuff the BCs inside the function data at the right locations:
         call r8util_bcherm2(f, ixdim, iydim,
     >        ibcxmin, ibcxmax, ibcymin, ibcymax,
     >        bcxmin, bcxmax, bcymin, bcymax,
     >        x, y)
 
         call r8akherm2p(x,ixdim,y,iydim,
     >        f,ixdim,idum1,idum2,ipx,ipy,ier)
 
      else if(kspline.eq.2) then
         ! bicubic Spline
 
         call r8mkbicub(x,nx, y,ny, f,ixdim,
     >        ibcxmin,bcxmin,ibcxmax,bcxmax,
     >        ibcymin,bcymin,ibcymax,bcymax,
     >        idum1,idum2,ier)
 
      else
         ! Hybrid
 
         if(jspline(1).gt.0) then
            ! cubic in x direction
            do ii=1,iydim
               if(jspline(1).eq.1) then
                  call r8util_bcherm1(f(1,1,ii), ixdim,
     >                 ibcxmin, ibcxmax, bcxmin(ii), bcxmax(ii), x)
                  call r8akherm1p(x,ixdim,f(1,1,ii),idum1,ipx,ier)
 
               else if(jspline(1).eq.2) then
                  call r8mkspline(x,ixdim,f(1,1,ii),
     >                 ibcxmin,bcxmin(ii), ibcxmax,bcxmax(ii),
     >                 idum1,ier)
               endif
               if(ier.ne.0) exit
            enddo
 
         else
            ! cubic in y direction
            allocate(wk2(2,ny))
 
            do ii=1,ixdim
               wk2(1,1:iydim) = f(1,ii,1:iydim)
               wk2(2,1:iydim) = 0.0_r8
               if(jspline(2).eq.1) then
                  call r8util_bcherm1(wk2, iydim,
     >                 ibcymin, ibcymax, bcymin(ii), bcymax(ii), y)
                  call r8akherm1p(y,iydim,wk2,idum2,ipy,ier)
 
               else if(jspline(2).eq.2) then
                  call r8mkspline(y,iydim,wk2,
     >                 ibcymin,bcymin(ii), ibcymax,bcymax(ii),
     >                 idum2,ier)
 
               endif
               if(ier.ne.0) exit
               f(1:2,ii,1:iydim) = wk2(1:2,1:iydim)
            enddo
 
            deallocate(wk2)
         endif
 
      endif
 
      return
      end
