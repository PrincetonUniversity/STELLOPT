      subroutine r8herm1ev(xget,x,nx,ilinx,f,ict,fval,ier)
C
C  evaluate a 1d cubic Hermite interpolant -- this is C1.
C
C  this subroutine calls two subroutines:
C     herm1x   -- find cell containing (xget,yget)
C     herm1fcn -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx
!============
      REAL*8 xget                         ! target of this interpolation
      REAL*8 x(nx)                        ! ordered x grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
C
      REAL*8 f(0:1,nx)                    ! function data
C
C       contents of f:
C
C  f(0,i) = f @ x(i)
C  f(1,i) = df/dx @ x(i)
C
      integer ict(2)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C
C output arguments:
C =================
C
      REAL*8 fval(*)                      ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C
C  examples:
C    on input ict = [1,1]
C   on output fval= [f,df/dx]
C
C    on input ict = [1,0]
C   on output fval= [f] ... element 2 never referenced
C
C    on input ict = [0,1]
C   on output fval= [df/dx] ... element 2 never referenced
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i                         ! cell index
C
C  normalized displacement from (x(i)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C
      REAL*8 xparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i))
C
      REAL*8 hx
      REAL*8 hxi
C
C  0 .le. xparam .le. 1
C
C---------------------------------------------------------------------
C
      call r8herm1x(xget,x,nx,ilinx,i,xparam,hx,hxi,ier)
      if(ier.ne.0) return
c
      call r8herm1fcn(ict,1,1,fval,i,xparam,hx,hxi,f,nx)
C
      return
      end
C---------------------------------------------------------------------
c  herm1xy -- look up x zone
c
c  this is the "first part" of herm1ev, see comments, above.
c
      subroutine r8herm1x(xget,x,nx,ilinx,i,xparam,hx,hxi,ier)
c
c  input of herm1x
c  ===============
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nxm,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zxget,zxtol
!============
      integer nx                        ! x array dimension
c
      REAL*8 xget                         ! target point
      REAL*8 x(nx)                        ! indep. coordinate, strict ascending
c
      integer ilinx                     ! =1:  x evenly spaced
c
c  output of herm1x
c  =================
      integer i                         ! index to cell containing target pt
c          on exit:  1.le.i.le.nx-1
c
c  normalized position w/in (i) cell (always btw 0 and 1):
c
      REAL*8 xparam                       ! (xget-x(i))/(x(i+1)-x(i))
c
c  grid spacing
c
      REAL*8 hx                           ! hx = x(i+1)-x(i)
c
c  inverse grid spacing:
c
      REAL*8 hxi                          ! 1/hx = 1/(x(i+1)-x(i))
c
      integer ier                       ! return ier.ne.0 on error
c
c------------------------------------
c
      ier=0
c
c  range check
c
      zxget=xget
      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0E-7_r8*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
            write(6,1001) xget,x(1),x(nx)
 1001       format(' ?herm1ev:  xget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((xget.lt.x(1)-0.5_r8*zxtol).or.
     >         (xget.gt.x(nx)+0.5_r8*zxtol))
     >      write(6,1011) xget,x(1),x(nx)
 1011       format(' %herm1ev:  xget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(xget.lt.x(1)) then
               zxget=x(1)
            else
               zxget=x(nx)
            endif
         endif
      endif
      if(ier.ne.0) return
c
c  now find interval in which target point lies..
c
      nxm=nx-1
c
      if(ilinx.eq.1) then
         ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
         i=min(nxm, ii)
         if(zxget.lt.x(i)) then
            i=i-1
         else if(zxget.gt.x(i+1)) then
            i=i+1
         endif
      else
         if((1.le.i).and.(i.lt.nxm)) then
            if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
               continue  ! already have the zone
            else
               call r8zonfind(x,nx,zxget,i)
            endif
         else
            i=nx/2
            call r8zonfind(x,nx,zxget,i)
         endif
      endif
c
      hx=(x(i+1)-x(i))
c
      hxi=1.0_r8/hx
c
      xparam=(zxget-x(i))*hxi
c
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 1d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine r8herm1fcn(ict,ivec,ivecd,
     >   fval,ii,xparam,hx,hxi,fin,nx)
C
C  input:
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx,i,iadr
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xp,xpi,xp2,xpi2,ax,axbar,bx,bxbar,axp,axbarp,bxp,bxbarp
!============
      integer ict(2)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec)                  ! target cells
      REAL*8 xparam(ivec)
                                  ! normalized displacements in cells
C
      REAL*8 hx(ivec)                     ! grid spacing, and
      REAL*8 hxi(ivec)                    ! inverse grid spacing 1/(x(i+1)-x(i))
C
      REAL*8 fin(0:1,nx)                  ! Hermite dataset
C
C  output:
C
      REAL*8 fval(ivecd,*)                ! interpolation results
C
C  for detailed description of fin, ict and fval see subroutine
C  herm1ev comments.  Note ict is not vectorized -- the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument; ivecd.ge.ivec
C     expected!
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C  Hermite cubic basis functions
C  -->for function value matching
C     a(0)=0 a(1)=1        a'(0)=0 a'(1)=0
C   abar(0)=1 abar(1)=0  abar'(0)=0 abar'(1)=0
C
C   a(x)=-2*x**3 + 3*x**2    = x*x*(-2.0*x+3.0)
C   abar(x)=1-a(x)
C   a'(x)=-abar'(x)          = 6.0*x*(1.0-x)
C
C  -->for derivative matching
C     b(0)=0 b(1)=0          b'(0)=0 b'(1)=1
C   bbar(0)=0 bbar(1)=0  bbar'(0)=1 bbar'(1)=0
C
C     b(x)=x**3-x**2         b'(x)=3*x**2-2*x
C     bbar(x)=x**3-2*x**2+x  bbar'(x)=3*x**2-4*x+1
C
C   ...in x direction
C
      REAL*8 sum
      integer v
C
      do v=1,ivec
         i=ii(v)
         xp=xparam(v)
         xpi=1.0_r8-xp
         xp2=xp*xp
         xpi2=xpi*xpi
C
         iadr=0
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
C
C  Hermite basis functions
C
            ax=xp2*(3.0_r8-2.0_r8*xp)
            axbar=1.0_r8-ax
            bx=-xp2*xpi
            bxbar=xpi2*xp
C
            sum=axbar*fin(0,i) + ax*fin(0,i+1)
C
            sum=sum+hx(v)*(bxbar*fin(1,i) + bx*fin(1,i+1))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
C
            axp=6.0_r8*xp*xpi
            axbarp=-axp
            bxp=xp*(3.0_r8*xp-2.0_r8)
            bxbarp=xpi*(3.0_r8*xpi-2.0_r8)
C
            sum=hxi(v)*(axbarp*fin(0,i) +axp*fin(0,i+1))
C
            sum=sum+ bxbarp*fin(1,i) + bxp*fin(1,i+1)
C
            fval(v,iadr)=sum
         endif
C
      enddo                             ! vector loop
C
      return
      end
