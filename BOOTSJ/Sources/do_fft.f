
      subroutine do_fft(a11, answer_mn, trigsu, trigsv, ifaxu, ifaxv,
     1   ntheta, nzeta, mbuse, nbuse)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use stel_kinds
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, nzeta, mbuse, nbuse
      integer, dimension(*) :: ifaxu, ifaxv
      real(rprec), dimension(ntheta + 2,nzeta) :: a11
      real(rprec), dimension(3*ntheta/2 + 1) :: trigsu
      real(rprec), dimension(2*nzeta) :: trigsv
      complex(rprec), dimension(-mbuse:mbuse,0:nbuse) :: answer_mn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: inc, jump, isign, jumpv, incv
      real(rprec), dimension(:), allocatable :: work11
      real(rprec) :: factor
C-----------------------------------------------

      allocate (work11(nzeta*(ntheta+2)))
c
c- Forward real to complex transform in the theta direction with index n.
c  i.e., integral of exp(-i*n*theta) * f(theta,zetah).
c
      inc = 1
      jump = ntheta + 2
      isign = -1

      call fft991(a11,work11,trigsu,ifaxu,inc,jump,ntheta,nzeta,isign)

c
c- now forward transform in the zetah direction with index m.
c  i.e., integral of exp(-i*m*zetah) * [theta transform of f(theta,zetah)]
c
      jumpv = 1
      incv = jump/2
      call cfft99(a11,work11,trigsv,ifaxv,incv,jumpv,nzeta,incv,isign)

      deallocate (work11)

c
c- Now reorganize and scale these to get the complex d(n,m)
c  FACTOR = 1 / (number of points used in the forward transform direction).
c  Because of (anti)symmetry, only nonnegative m values are needed.  We also
c  only fill to n = mbuse and m = nbuse since this is all that will be
c  used in the sums over n and m.
c
      factor = one/nzeta
c  store a11 in answer_mn array
      call reorganz (answer_mn, mbuse, nbuse, factor, a11, ntheta,
     1   nzeta)

      end subroutine do_fft
