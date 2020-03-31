
      subroutine fraction(irho)
c--
c  This calculates the fraction passing and the fraction trapped.
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
      integer irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0
      integer, parameter :: n_lambda = 41
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l
      real(rprec), dimension(n_lambda) ::
     1   alambda, avgvpov, antgrand, answer

C-----------------------------------------------

      dlambda = one/(n_lambda - 1)

c  Fill lambda mesh, evaluate V||/V

      do l = 1, n_lambda
         alambda(l) = (l - 1)*dlambda
         avgvpov(l) = sum(sqrt(abs(one - alambda(l)*bfield))*gsqrt_b)
      end do

c  Form integrand

      where (avgvpov .gt. zero) antgrand = alambda/avgvpov

c  Do integral for the passing fraction

      call simpun (alambda, antgrand, n_lambda, answer)

      fpassing(irho) = .75_dp*avgbobm2*answer(n_lambda)*sum_gsqrt_b
      ftrapped(irho) = one - fpassing(irho)

      end subroutine fraction
