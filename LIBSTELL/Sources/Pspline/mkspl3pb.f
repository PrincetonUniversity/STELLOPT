      subroutine mkspl3pb(fun,x,nx,th,nth,ph,nph,fspl,nf4,nf5,
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
      external fun                      ! passed real function(x,th)
      real x(nx)                        ! x coordinate array
      real th(nth)                      ! th coordinate array
      real ph(nph)                      ! ph coordinate array
C
      real fspl(4,4,4,nf4,nf5,nph)      ! function data / spline coeff array
      real wk(inwk)                     ! workspace for bpsplinb
C
      integer ibcxmin                   ! boundary condition indicator
      real bcxmin(nb1,nph)              ! boundary condition data
      integer ibcxmax                   ! boundary condition indicator
      real bcxmax(nb1,nph)              ! boundary condition data
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
      if(nf4.lt.nx) then
         write(6,'('' ?mkspl3pb -- array dim error, nf4 .lt. nx'')')
         ier=1
      endif
      if(nf5.lt.nth) then
         write(6,'('' ?mkspl3pb -- array dim error, nf5 .lt. nth'')')
         ier=2
      endif
      if(nb1.lt.nth) then
         write(6,'('' ?mkspl3pb -- array dim error, nb1 .lt. nth'')')
         ier=3
      endif
      if(inwk.lt.5*max(nx,nth,nph)) then
         write(6,'('' ?mkspl3pb -- array dim error, inwk too small'')')
         ier=4
      endif
 
      if(ier.ne.0) return
C
      do iph=1,nph
         do ith=1,nth
            do ix=1,nx
               fspl(1,1,1,ix,ith,iph)=fun(x(ix),th(ith),ph(iph))
            enddo
         enddo
      enddo
C
      call tpsplinb(x,nx,th,nth,ph,nph,fspl,nf4,nf5,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,nb1,
     >   wk,inwk,ilinx,ilinth,ilinph,ier)
C
      return
      end
