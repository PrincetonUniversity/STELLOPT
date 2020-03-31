      subroutine r8akherm1(x,nx,fherm,ilinx,ier)
C
C  create a data set for Hermite interpolation, based on Akima's method
C  [Hiroshi Akima, Communications of the ACM, Jan 1974, Vol. 17 No. 1]
C
C  1d routine -- default boundary condition (based on divided differences)
C
C  input:
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ier,ilinx
!============
      integer nx                        ! array dimensions
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 fherm(0:1,nx)                ! data/Hermite array
C
C  fherm(0,i) = function value f at x(i)       **on input**
C
C  fherm(1,i) = derivative df/dx at x(i)       **on output**
C
C addl output:
C  ilinx=1 if x is "evenly spaced" ier=0 if no errors
C
C  ** x must be strict ascending **
C
C work by calling akherm1p; no periodic boundary condition
C
      call r8akherm1p(x,nx,fherm,ilinx,0,ier)
C
      return
      end
C----------------------------
      subroutine r8akherm1p(x,nx,fherm,ilinx,ipx,ier)
C
C  Akima subroutine with boundary condition option:
C
C  =>ipx=1 for periodic boundary condition
C
C  =>ipx=2 for user specified boundary condition:
C    in which case fherm(1,1) and fherm(1,nx) must contain
C    the derivatives fx(x(1)) and fx(x(nx)), respectively, and these
C    will not be modified on output.
C
C  =>ipx=0: default boundary conditions
C
C  other arguments as with akherm1.
C
C  input:
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ier,ilinx,ierbc,ix,ixm2,ixm1,ixmm2,ixmm1,ixp2,ixp1
      INTEGER ixpp2,ixpp1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ztol,cxp,cxm,cxpp,cxmm,cxtrap0,cxtrap1
!============
      integer nx                        ! array dimensions
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 fherm(0:1,nx)                ! data/Hermite array
      integer ipx                       ! =1:  f' periodic in x
C
C----------------------------
C
      REAL*8 wx(2)
C
C  error checks...
C
      ztol=1.0E-3_r8
      ier=0
C
      call r8splinck(x,nx,ilinx,ztol,ier)
      if(ier.ne.0) then
         write(6,*) '?akherm1:  x axis not strict ascending.'
         ier=ier+1
      endif
C
      ierbc=0
      call ibc_ck(ipx,'akherm1','Bdy Cond',0,2,ierbc)
      ier=ier+ierbc
      if(ier.gt.0) return
C
C  deal with boundary region.  The boundary derivatives are set
C  and "ghost" points are provided...
C
C  all paths need the 1st div. diffs at the bdy...
C
      cxp=(fherm(0,2)-fherm(0,1))/(x(2)-x(1))
      cxm=(fherm(0,nx)-fherm(0,nx-1))/(x(nx)-x(nx-1))
C
      if(ipx.eq.1) then
C
C  periodic BC
C
C  LHS/RHS
C
         if(nx.gt.2) then
            cxpp=(fherm(0,3)-fherm(0,2))/(x(3)-x(2))
            cxmm=(fherm(0,nx-1)-fherm(0,nx-2))/(x(nx-1)-x(nx-2))
C
            call r8akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,1))
            fherm(1,nx)=fherm(1,1)
         else
            ! nx=2
            fherm(1,1)=cxp  ! =cxm
            fherm(1,nx)=cxm ! =cxp
         endif
c
         cxtrap0=cxm
         cxtrap1=cxp
C
      else if(ipx.eq.0) then
C
C  default BC -- standard numeric extrapolation
C
         if(nx.gt.2) then
            cxpp=(fherm(0,3)-fherm(0,2))/(x(3)-x(2))
            fherm(1,1)=1.5_r8*cxp-0.5_r8*cxpp
C
            cxmm=(fherm(0,nx-1)-fherm(0,nx-2))/(x(nx-1)-x(nx-2))
            fherm(1,nx)=1.5_r8*cxm-0.5_r8*cxmm
C
         else
            ! nx=2
            fherm(1,1)=cxp  ! =cxm
            fherm(1,nx)=cxm ! =cxp
         endif
C
C  extrapolate to slope to ghost points just past bdy...
C
         cxtrap0=2.0_r8*fherm(1,1)-cxp
         cxtrap1=2.0_r8*fherm(1,nx)-cxm
C
      else
C
C  BC supplied by user.  Also use this for extrapolation...
C  extrapolate to slope to ghost points just past bdy...
C
         cxtrap0=2.0_r8*fherm(1,1)-cxp
         cxtrap1=2.0_r8*fherm(1,nx)-cxm
C
      endif
C
C NOTE: this loop is inactive if nx=2
C
      do ix=2,nx-1
c
c  x div. diffs in vicinity
c
         ixm2=ix
         ixm1=ixm2-1
         ixmm2=ixm1
         ixmm1=ixm1-1
c
         ixp2=ix+1
         ixp1=ix
         ixpp2=ix+2
         ixpp1=ixp2
c
         if(ix.eq.2) then
            cxmm=cxtrap0
         else
            cxmm=(fherm(0,ixmm2)-fherm(0,ixmm1))/(x(ixmm2)-x(ixmm1))
         endif
c
         if(ix.eq.nx-1) then
            cxpp=cxtrap1
         else
            cxpp=(fherm(0,ixpp2)-fherm(0,ixpp1))/(x(ixpp2)-x(ixpp1))
         endif
c
         cxm=(fherm(0,ixm2)-fherm(0,ixm1))/(x(ixm2)-x(ixm1))
         cxp=(fherm(0,ixp2)-fherm(0,ixp1))/(x(ixp2)-x(ixp1))
C
         call r8akherm0(cxmm,cxm,cxp,cxpp,wx,fherm(1,ix))
c
      enddo
C
      return
      end
C--------------------------------------
      subroutine r8akherm0(cxmm,cxm,cxp,cxpp,wx,result)
c
c  basic akima formula for 1st derivative at pt p:
c
c     cxmm = numerical slope 2 zones to left
c     cxm  = numerical slope 1 zone to left
c     cxp  = numerical slope 1 zone to right
c     cxpp = numerical slope 2 zones to right
c
c  return slope at point p and weighting factors
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 cxmm,cxm,cxp,cxpp            ! nearby slopes (in)
      REAL*8 wx(2)                        ! weights (out)
      REAL*8 result                       ! Akima nominal slope (out)
c
      wx(2)=abs(cxm-cxmm)
      wx(1)=abs(cxp-cxpp)
      if(wx(1)+wx(2).eq.0.0_r8) then
         wx(1)=1.0_r8
         wx(2)=1.0_r8
      endif
c
c  the derivative -- weighted average of neighbouring slopes
c
      result=(wx(1)*cxm+wx(2)*cxp)/(wx(1)+wx(2))
c
      return
      end
