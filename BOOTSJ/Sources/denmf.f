
      subroutine denmf(trigsu, trigsv, ifaxu, ifaxv, irho)
c  Evaluate the coefficients d(m,n) using CRAY fft991, cfft99
c  vectorized 1D FFT routines.  Also evaluate fraction trapped and fraction
c  passing.
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
      integer, dimension(13) :: ifaxu, ifaxv
      real(rprec), dimension(3*nthetah/2 + 1) :: trigsu
      real(rprec), dimension(2*nzetah) :: trigsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, D18 = 1.0e-18_dp,
     1 one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, i
      real(rprec), dimension(nthetah + 2,nzetah) :: a11
C-----------------------------------------------
c
c- First evaluate the complex coefficients d(n,m).
c
c- Load the the array a11 with (Bm/B)**2 * (1-B/Bm)**1.5.
c  Remember that the arrays BFIELD now contains B/Bm.
c  this is a flux surface average

      if (isymm0 .eq. 0) then
         a11(:nthetah,:nzetah) =
     1      (abs(one - bfield(:nthetah,:nzetah)) + D18)**1.5_dp
     2      *b2avg(irho)/(bmax1(irho)*bfield(:nthetah,:nzetah))**2
         a11(nthetah+1,:nzetah) = zero
         a11(nthetah+2,:nzetah) = zero

         call do_fft (a11, dmn, trigsu, trigsv, ifaxu, ifaxv, nthetah,
     1      nzetah, mbuse, nbuse)

      endif

      avgbobm2 = b2avg(irho)/bmax1(irho)**2

c  Now calculate the fraction passing and the fraction trapped.

      call fraction(irho)

      end subroutine denmf
