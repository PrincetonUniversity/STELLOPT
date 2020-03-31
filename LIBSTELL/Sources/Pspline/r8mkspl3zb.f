      subroutine r8mkspl3zb(fun,x,nx,th,nth,ph,nph,fspl,nf2,nf3,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,nb1,
     >   wk,inwk,ilinx,ilinth,ilinph,ier)
C
C  create a tricubic biperiodic spline with knots at grid points and
C  function values from the callable function `fun' passed.
C
C  use tpsplinb to setup the spline coefficients
C
C  periodic boundary condition for theta & phi grids.
C  boundary condition data may be specified at x(1) and x(nx) for each
C  theta & phi:
C     ibcxmin=0 ==> "not a knot" boundary condition, cf cubspl.for
C     ibcxmin=1 ==> match slope, bcxmin(ith,iph) gives df/dx at
C         x(1),th(ith),ph(iph)
C     ibcxmin=2 ==> match 2nd derivative, given at x(1),th(ith),ph(iph)
C          by bcxmin(ith,iph)
C
C     ibcxmax,bcxmax have analogous interpretation -- at x(nx).
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nth,nph,nf2,nf3,nb1,inwk,nx,iph,ith,ix
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fun,zdummy
!============
      external fun                      ! passed real function(x,th)
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 th(nth)                      ! th coordinate array
      REAL*8 ph(nph)                      ! ph coordinate array
C
      REAL*8 fspl(8,nf2,nf3,nph)          ! function data / spline coeff array
      REAL*8 wk(inwk)                     ! workspace for bpsplinb
C
      integer ibcxmin                   ! boundary condition indicator
      REAL*8 bcxmin(nb1,nph)              ! boundary condition data
      integer ibcxmax                   ! boundary condition indicator
      REAL*8 bcxmax(nb1,nph)              ! boundary condition data
C
      integer ilinx                     ! output =1 if x(...) evenly spaced
      integer ilinth                    ! output =1 if th(...) evenly spaced
      integer ilinph                    ! output =1 if ph(...) evenly spaced
C
      integer ier                       ! completion code from bpspline
C
C----------------------------
C
      ier=0
      if(nf2.lt.nx) then
         write(6,'('' ?mkspl3pb -- array dim error, nf2 .lt. nx'')')
         ier=1
      endif
      if(nf3.lt.nth) then
         write(6,'('' ?mkspl3pb -- array dim error, nf3 .lt. nth'')')
         ier=2
      endif
      if(nb1.lt.nth) then
         write(6,'('' ?mkspl3pb -- array dim error, nb1 .lt. nth'')')
         ier=3
      endif
 
      if(ier.ne.0) return
C
      do iph=1,nph
         do ith=1,nth
            do ix=1,nx
               fspl(1,ix,ith,iph)=fun(x(ix),th(ith),ph(iph))
            enddo
         enddo
      enddo
C
      call r8mktricubw(x,nx,th,nth,ph,nph,fspl,nf2,nf3,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,nb1,
     >   -1,zdummy,-1,zdummy,max(nx,nth,nph),
     >   -1,zdummy,-1,zdummy,max(nx,nth,nph),
     >   wk,inwk,ilinx,ilinth,ilinph,ier)
C
      return
      end
