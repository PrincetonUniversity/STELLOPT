      subroutine r8mkspl2p(fun,x,nx,th,nth,fspl,nf3,wk,inwk,
     >   ilinx,ilinth,ier)
C
C  create a bicubic periodic spline with knots at grid points and
C  function values from the callable function `fun' passed.
C
C  use bpspline to setup the spline coefficients
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nth,nf3,inwk,nx,ix,ith
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fun
!============
      external fun                      ! passed real function(x,th)
      REAL*8 x(nx)                        ! x coordinate array
      REAL*8 th(nth)                      ! th coordinate array
C
C  output:
C
      REAL*8 fspl(4,4,nf3,nth)            ! function data / spline coeff array
      REAL*8 wk(inwk)                     ! workspace -- at least 5*max(nx,nth)
C
      integer ilinx                     ! output =1 if x(1...nx) evenly spaced
      integer ilinth                    ! output =1 if th(1..nth) evenly spaced
C
      integer ier                       ! completion code from bpspline 0=OK
C
C----------------------------
C
      ier=0
      if(nf3.lt.nx) then
         write(6,'('' ?mkspl2p -- array dim error, nf3.lt.nx'')')
         ier=1
      endif
      if(inwk.lt.5*max(nx,nth)) then
         write(6,'('' ?mkspl2p -- array dim error, inwk too small'')')
         ier=2
      endif
C
      do ix=1,nx
         do ith=1,nth
            fspl(1,1,ix,ith)=fun(x(ix),th(ith))
         enddo
      enddo
C
      call r8bpspline(x,nx,th,nth,fspl,nf3,wk,inwk,ilinx,ilinth,ier)
C
      return
      end
