      subroutine r8mkspl2zb(fun,x,nx,th,nth,fspl,nf2,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   wk,inwk,ilinx,ilinth,ier)
C
C  create a bicubic periodic spline with knots at grid points and
C  function values from the callable function `fun' passed.
C
C  use bpsplinb to setup the spline coefficients
C
C  periodic boundary condition for theta grid;
C  boundary condition data may be specified at x(1) and x(nx)
C     ibcxmin=0 ==> "not a knot" boundary condition, cf cubspl.for
C     ibcxmin=1 ==> match slope, bcxmin(ith) gives df/dx at x(1),th(ith).
C     ibcxmin=2 ==> match 2nd derivative, given at x(1),th(ith) by bcxmin(ith)
C
C     ibcxmax,bcxmax have analogous interpretation -- at x(nx)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nth,nf2,inwk,nx,ix,ith
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fun,zdummy
!============
      external fun                      ! passed real function(x,th)
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 th(nth)                      ! th coordinate array
C
      REAL*8 fspl(4,nf2,nth)              ! function data / spline coeff array
      REAL*8 wk(inwk)                     ! workspace for bpsplinb
C
      integer ibcxmin                   ! boundary condition indicator
      REAL*8 bcxmin(nth)                  ! boundary condition data
      integer ibcxmax                   ! boundary condition indicator
      REAL*8 bcxmax(nth)                  ! boundary condition data
C
      integer ilinx                     ! output =1 if x(...) evenly spaced
      integer ilinth                    ! output =1 if th(...) evenly spaced
C
      integer ier                       ! completion code from bpspline
C
C----------------------------
C
      ier=0
      if(nf2.lt.nx) then
         write(6,'('' ?mkspl2pb -- array dim error, nf2.lt.nx'')')
         ier=1
      endif
C
      do ix=1,nx
         do ith=1,nth
            fspl(1,ix,ith)=fun(x(ix),th(ith))
         enddo
      enddo
C
      call r8mkbicubw(x,nx,th,nth,fspl,nf2,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   -1,zdummy,-1,zdummy,
     >   wk,inwk,ilinx,ilinth,ier)
C
      return
      end
