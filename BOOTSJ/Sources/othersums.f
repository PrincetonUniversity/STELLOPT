
      subroutine othersums(irho)
c--
c     calculates other1
c     of integral of W(lambda) which is not dependant on lambda
c--
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
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, m
      real(rprec) :: denom, qn
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: sumit
C-----------------------------------------------
c
c-  m - poloidal, n - toroidal
c
      if (isymm0 .ne. 0) then
         other1(irho) = zero
         return
      endif
c
c- Form the "other" sums.  The first uses alpha1(m,n), calculated in WOFLAM,
c  and d(m,n) from DENM.
c
c- Load  rfmn with
c
c  (m*R+n*periods*S)/(m-n*periods*q) *
c         exp(m*thetamax-n*zetahmax) * (1.5*alpha1(m,n)+d(m,n))
c
c- for only those harmonics that are going to be used in the sum.
c
      qn = periods*qsafety(irho)*zetasign
      do m = -mbuse, mbuse
         do n = 0, nbuse
            denom = m + n*qn
            if (n.ne.0 .or. m.ne.0) then
               denom = denom/(denom**2 + (damp_bs*m)**2)
               rfmn(m,n) = (m*capr(irho)+n*periods*caps(irho))*denom*(
     1            cos(m*thetamax(irho)-n*zetahmax(irho))*real(1.5_dp*
     2            alpha1mn(m,n)+dmn(m,n))-sin(m*thetamax(irho)-n*
     3            zetahmax(irho))*aimag(1.5_dp*alpha1mn(m,n)+dmn(m,n)))
            else
               rfmn(m,n) = zero
            endif
         end do
      end do

c  First sum,

      other1(irho) = sumit(rfmn,mbuse,nbuse)
c
c- Then multiply by all the stuff out front.
c
      other1(irho) = -other1(irho)*qsafety(irho)/ftrapped(irho)*(one +
     1   aiogar(irho)/qsafety(irho))

      end subroutine othersums
