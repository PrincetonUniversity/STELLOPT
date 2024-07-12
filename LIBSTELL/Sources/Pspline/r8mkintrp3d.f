      subroutine r8mkintrp3d(x,nx,y,ny,z,nz,jspline,
     >     f,icoeff,ixdim,iydim,izdim,
     >     ibcxmin,bcxmin,ibcxmax,bcxmax,
     >     ibcymin,bcymin,ibcymax,bcymax,
     >     ibczmin,bczmin,ibczmax,bczmax,
     >     ier)
c
c  setup a tricubic spline, or tricubic Hermite, or hybrid linear/zonal 2d or
c  1d with 1d or 2d cubic or Hermite spline interpolation
c
C
C  input:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer nx                        ! length of x vector
      integer ny                        ! length of y vector
      integer nz                        ! length of z vector
      REAL*8 x(nx)                        ! x vector, strict ascending
      REAL*8 y(ny)                        ! y vector, strict ascending
      REAL*8 z(nz)                        ! z vector, strict ascending
c
      integer :: jspline(3)             ! interpolation method control
C        (1) -- 1st dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
C        (2) -- 2nd dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
C        (3) -- 3rd dimension: -1: zone, 0: pclin, 1: Hermite, 2: spline
C
C    Standard interpolation-- all jspline values match
C      e.g. jspline(1)=jspline(2)=jspline(3)=2 for tricubic spline
C
C    Hybrid interpolation-- not all jspline values are the same.
C
C    RESTRICTION: if any jspline(...) element has value 1, none can have
C    value 2.  I.e. Spline and Hermite interpolation cannot currently be mixed.
C    This restriction exists because of technical issues in the
C    implementation (it could be removed in principle but the work to do
C    this has not been scheduled).
c
c  coefficient buffer dimensions
      integer :: icoeff                 ! #coefficients per data point
      integer :: ixdim                  ! nx; nx-1 if jspline(1)==-1
      integer :: iydim                  ! ny; ny-1 if jspline(2)==-1
      integer :: izdim                  ! nz; nz-1 if jspline(3)==-1
c  input/output:
      REAL*8 f(icoeff,ixdim,iydim,izdim)  ! data and spline coefficients
c
C  boundary condition data
C    Contiguous storage is assumed-- 1st dimension size must match
C      actual use
C
C
      integer ibcxmin,ibcxmax           ! BC type flag @xmin, xmax
      integer ibcymin,ibcymax           ! BC type flag @ymin, ymax
      integer ibczmin,ibczmax           ! BC type flag @zmin, zmax
C
      REAL*8 bcxmin(iydim,izdim),bcxmax(iydim,izdim) ! xmin & xmax BC data
      REAL*8 bcymin(ixdim,izdim),bcymax(ixdim,izdim) ! ymin & ymax BC data
      REAL*8 bczmin(ixdim,iydim),bczmax(ixdim,iydim) ! zmin & zmax BC data
c
c  where BC data is not required, dummy scalars may be passed.
C  the ibc* flags determine whether BC data isneeded.
c
c  BC data not required for zonal or piecewise linear interpolation
c  for Hermite interpolation ibc* values from set {-1,0,1} are accepted.
c
c  BC data:  bcxmin & bcxmax:  BC vs. y,z @xmin,xmax
C            bcymin & bcymax:  BC vs. x,z @ymin,ymax
C            bczmin & bczmax:  BC vs. x,y @zmin,zmax
C
c   ibcxmin -- indicator for boundary condition at xmin=x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- use periodic boundary condition
c     =0 -- use "not a knot"
c     =1 -- match slope, specified at x(1),y(iy),z(iz) by bcxmin(iy,iz)
c     =2 -- match 2nd derivative, specified at x(1),y(iy),z(iz)
c           by bcxmin(iy,iz
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all y(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all y(j)
c     =5 -- df/dx BC from 1st divided difference
c     =6 -- d2f/dx2 BC from 2nd divided difference (parabolic fit)
c     =7 -- d3f/dx3 BC from 3rd divided difference (cubic fit)
c   ***NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
c
c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
c            and ibcxmax, bcxmax are ignored.
c
c   interpretation of ibcymin,bcymin,ibcymax,bcymax
c     is same as with ibcxmin,...
c
c   interpretation of ibczmin,bczmin,ibczmax,bczmax
c     is same as with ibcxmin,...
c
c   the explicit bdy condition arrays are referenced only if the
c     corresponding "ibc" flag values are set > 0.
c
      integer ier                       ! exit code
c   ier -- completion code, 0 for normal
c
C-----------------------
      integer :: kspline
      integer :: ii,jj,imul,imin,imax,ickx,icky,ickz,inum
      integer :: idum1,idum2,idum3,idimcu
      integer :: ipx,ipy,ipz
      REAL*8 :: ztol = 0.0001_r8
      logical :: ifound(-1:2)
 
      REAL*8, dimension(:,:), allocatable :: wk2
      REAL*8, dimension(:,:,:), allocatable :: wk3
C-----------------------
c
      ier=0
c
      imin=3
      imax=-2
      imul=1
      inum=0
      idimcu=0
      ifound = .FALSE.
 
      do ii=1,3
         imin=min(imin,jspline(ii))
         imax=max(imax,jspline(ii))
         if(jspline(ii).gt.0) then
            idimcu=idimcu+1
            imul=imul*2
         endif
         if(.not.ifound(jspline(ii))) then
            ifound(jspline(ii)) = .TRUE.
            inum = inum+1
         endif
      enddo
 
      if((imin.lt.-1).or.(imax.gt.2)) then
         ier = 1
         write(6,*)
     >        ' ?mkintrp3d: spline type control out of range -1 to 2: ',
     >        jspline
      endif
      if(ier.ne.0) return
 
      if(inum.eq.1) then
         kspline=imin   ! same interp type on all dimensions
      else
         kspline=-99    ! hybrid
         if(ifound(1).and.ifound(2)) then
            ier = 1
            write(6,*)
     >           ' ?mkintrp3d: spline/Hermite hybrid not supported (',
     >           jspline,')'
         endif
      endif
      if(ier.ne.0) return
c
      if(imul.ne.icoeff) then
         write(6,*)
     >        ' ?coeff dimension inconsistency for spline type codes ',
     >        jspline
         write(6,*) ' in mkintrp3d: expected: ',imul,' got: ',icoeff
         ier=1
         return
      endif
c
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
      if(jspline(3).eq.-1) then
         ickz=nz-1
      else
         ickz=nz
      endif
c
      if((ickx.ne.ixdim).or.(icky.ne.iydim).or.(ickz.ne.izdim)) then
         write(6,*)
     >        ' ?mkintrp2d: dimensioning inconsistent with '//
     >        'interpolation controls: ',jspline
         write(6,*) '  expected: ',ickx,icky,ickz,
     >        '; got: ',ixdim,iydim,izdim
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
      if(jspline(3).le.0) then
         call r8splinck(z,nz,idum1,ztol,ier)
         if(ier.ne.0) then
            write(6,*) ' ?mkintrp2d: z axis not strict ascending.'
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
c
      if(jspline(3).eq.1) then
         if((min(ibczmin,ibczmax).lt.-1).or.
     >        (max(ibczmin,ibczmax).gt.1)) then
            write(6,*) ' ?mkintrp2d: Bdy Cond code out of range for'
            write(6,*) '  Hermite interpolation; (-1:1) allowed, '//
     >           'found: ',ibczmin,ibczmax
            ier=1
            return
         endif
         ipz=0
         if(ibczmin.eq.-1) then
            ipz=1
         else if((ibczmin.eq.1).or.(ibczmax.eq.1)) then
            ipz=2
         endif
      endif
c
      if(kspline.eq.1) then
         ! tricubic Hermite
 
         ! put the BCs inside the function data at the right locations...
         call r8util_bcherm3(f, ixdim, iydim, izdim,
     >        ibcxmin,ibcxmax, ibcymin,ibcymax, ibczmin,ibczmax,
     >        bcxmin, bcxmax,  bcymin, bcymax,  bczmin, bczmax,
     >        x, y, z)
 
         call r8akherm3p(x,ixdim, y,iydim, z,izdim, f,ixdim,iydim,
     >        idum1,idum2,idum3, ipx,ipy,ipz, ier)
 
      else if(kspline.eq.2) then
         ! tricubic Spline
 
         call r8mktricub(x,nx,y,ny,z,nz,
     >     f,ixdim,iydim,
     >     ibcxmin,bcxmin,ibcxmax,bcxmax,iydim,
     >     ibcymin,bcymin,ibcymax,bcymax,ixdim,
     >     ibczmin,bczmin,ibczmax,bczmax,ixdim,
     >     idum1,idum2,idum3, ier)
 
      else
         ! Hybrid
         if(idimcu.eq.1) then
            ! cubic along 1 dimension; other dims are step or pclin
 
            if(jspline(1).gt.0) then
                                ! cubic in x direction
               do jj=1,izdim
                  do ii=1,iydim
                     if(jspline(1).eq.1) then
                        call r8util_bcherm1(f(1,1,ii,jj), ixdim,
     >                       ibcxmin, ibcxmax,
     >                       bcxmin(ii,jj), bcxmax(ii,jj), x)
                        call r8akherm1p(x,ixdim,f(1,1,ii,jj),idum1,ipx,
     >                       ier)
 
                     else if(jspline(1).eq.2) then
                        call r8mkspline(x,ixdim,f(1,1,ii,jj),
     >                       ibcxmin,bcxmin(ii,jj),
     >                       ibcxmax,bcxmax(ii,jj),
     >                       idum1,ier)
                     endif
                     if(ier.ne.0) exit
                  enddo
                  if(ier.ne.0) exit
               enddo
 
            else if(jspline(2).gt.0) then
                                ! cubic in y direction
               allocate(wk2(2,ny))
 
               do jj=1,izdim
                  do ii=1,ixdim
                     wk2(1,1:iydim) = f(1,ii,1:iydim,jj)
                     wk2(2,1:iydim) = 0.0_r8
                     if(jspline(2).eq.1) then
                        call r8util_bcherm1(wk2, iydim,
     >                       ibcymin, ibcymax,
     >                       bcymin(ii,jj), bcymax(ii,jj), y)
                        call r8akherm1p(y,iydim,wk2,idum2,ipy,ier)
 
                     else if(jspline(2).eq.2) then
                        call r8mkspline(y,iydim,wk2,
     >                       ibcymin,bcymin(ii,jj),
     >                       ibcymax,bcymax(ii,jj),
     >                       idum2,ier)
 
                     endif
                     if(ier.ne.0) exit
                     f(1:2,ii,1:iydim,jj) = wk2(1:2,1:iydim)
                  enddo
                  if(ier.ne.0) exit
               enddo
 
               deallocate(wk2)
 
            else
                                ! cubic in z direction
               allocate(wk2(2,nz))
 
               do jj=1,iydim
                  do ii=1,ixdim
                     wk2(1,1:izdim) = f(1,ii,jj,1:izdim)
                     wk2(2,1:izdim) = 0.0_r8
                     if(jspline(3).eq.1) then
                        call r8util_bcherm1(wk2, izdim,
     >                       ibczmin, ibczmax,
     >                       bczmin(ii,jj), bczmax(ii,jj), z)
                        call r8akherm1p(z,izdim,wk2,idum2,ipz,ier)
 
                     else if(jspline(3).eq.2) then
                        call r8mkspline(z,izdim,wk2,
     >                       ibczmin,bczmin(ii,jj),
     >                       ibczmax,bczmax(ii,jj),
     >                       idum2,ier)
 
                     endif
                     if(ier.ne.0) exit
                     f(1:2,ii,jj,1:izdim) = wk2(1:2,1:izdim)
                  enddo
                  if(ier.ne.0) exit
               enddo
 
               deallocate(wk2)
 
            endif
         else
            ! cubic along 2 dimensions
            if(jspline(3).le.0) then
               do ii=1,izdim
                  if(jspline(1).eq.1) then
                     call r8util_bcherm2(f(1,1,1,ii),ixdim,iydim,
     >                    ibcxmin, ibcxmax, ibcymin, ibcymax,
     >                    bcxmin(1,ii), bcxmax(1,ii),
     >                    bcymin(1,ii), bcymax(1,ii),
     >                    x,y)
                     call r8akherm2p(x,ixdim,y,iydim,f(1,1,1,ii),ixdim,
     >                    idum1,idum2,ipx,ipy,ier)
                  else
                     call r8mkbicub(x,ixdim,y,iydim,f(1,1,1,ii),ixdim,
     >                    ibcxmin, bcxmin(1,ii), ibcxmax, bcxmax(1,ii),
     >                    ibcymin, bcymin(1,ii), ibcymax, bcymax(1,ii),
     >                    idum1,idum2,ier)
                  endif
                  if(ier.ne.0) exit
               enddo
 
            else if(jspline(2).le.0) then
               allocate(wk3(4,nx,nz))
 
               do ii=1,iydim
                  wk3 = f(1:4,1:nx,ii,1:nz)
                  if(jspline(1).eq.1) then
                     call r8util_bcherm2(wk3,ixdim,izdim,
     >                    ibcxmin, ibcxmax, ibczmin, ibczmax,
     >                    bcxmin(ii,1:nz), bcxmax(ii,1:nz),
     >                    bczmin(1,ii), bczmax(1,ii), x, z)
                     call r8akherm2p(x,ixdim,z,izdim,wk3,ixdim,
     >                    idum1,idum2,ipx,ipz,ier)
                  else
                     call r8mkbicub(x,ixdim,z,izdim,wk3,ixdim,
     >                    ibcxmin,bcxmin(ii,1:nz),
     >                    ibcxmax,bcxmax(ii,1:nz),
     >                    ibczmin,bczmin(1,ii), ibczmax,bczmax(1,ii),
     >                    idum1,idum2,ier)
                  endif
                  if(ier.ne.0) exit
                  f(1:4,1:nx,ii,1:nz) = wk3
               enddo
 
               deallocate(wk3)
 
            else
               allocate(wk3(4,ny,nz))
 
               do ii=1,ixdim
                  wk3 = f(1:4,ii,1:ny,1:nz)
                  if(jspline(2).eq.1) then
                     call r8util_bcherm2(wk3,iydim,izdim,
     >                    ibcymin, ibcymax, ibczmin, ibczmax,
     >                    bcymin(ii,1:nz), bcymax(ii,1:nz),
     >                    bczmin(ii,1:ny), bczmax(ii,1:ny), y, z)
                     call r8akherm2p(y,iydim,z,izdim,wk3,iydim,
     >                    idum1,idum2,ipy,ipz,ier)
                  else
                     call r8mkbicub(y,iydim,z,izdim,wk3,iydim,
     >                    ibcymin,bcymin(ii,1:nz),
     >                    ibcymax,bcymax(ii,1:nz),
     >                    ibczmin,bczmin(ii,1:ny),
     >                    ibczmax,bczmax(ii,1:ny),
     >                    idum1,idum2,ier)
                  endif
                  if(ier.ne.0) exit
                  f(1:4,ii,1:ny,1:nz) = wk3
               enddo
 
               deallocate(wk3)
 
            endif
         endif
 
      endif
 
      return
      end
