

      subroutine tok_fraction(ft_tok,irho)
c--
c  This calculates the fraction passing and the fraction trapped
c  for a tokamak using the Lin-Liu Miller approximation with Houlberg
c  weighting factors.  (Phys. Plas. 2, 1966 (1995).
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
      real(rprec) ft_tok
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer i
      real(rprec) fub, flb, bavg, b2av, sum_gsqrt_tok
      real(rprec) bmax
C-----------------------------------------------

      real(rprec), dimension(:), allocatable ::
     1    bfield_tok,  gsqrt_tok, b2obm_tok, one_b2
      allocate(bfield_tok(nthetah),  gsqrt_tok(nthetah),
     1   b2obm_tok(nthetah), one_b2(nthetah), stat = i)
      if(i .ne. 0) stop 'allocation error of tokamak fields'


c  make symetric copies of field quantites--average over zeta to extract
c  symmetric component.

      do i =1,nthetah
         bfield_tok(i) = sum(bfield(i,:nzetah))/nzetah
      enddo

c  renormalize to max tok field equivalent

      bmax = maxval(bfield_tok)
      bfield_tok = bfield_tok/bmax
      where(bfield_tok .gt. one) bfield_tok = one

c  calculate 1/b**2, b**2

      one_b2 = one/bfield_tok**2
      b2obm_tok = bfield_tok**2

c  jacobian only includes 1/b**2 component and is normalized to bmax tok
c  integrate over jacobian to obtain normalization for averages
      sum_gsqrt_tok= sum(one_b2)


c  find average b, b**2

      bavg = sum(bfield_tok*one_b2)/sum_gsqrt_tok
      b2av = sum(b2obm_tok*one_b2)/sum_gsqrt_tok

c  find upper bound

      fub = one-(one-sqrt(one-bavg)*(one+0.5_dp*bavg))*b2av/bavg**2

c  find lower bound


c  minus <1/b**2>

      flb = - sum(one/b2obm_tok**2)/sum_gsqrt_tok

c  plus <sqrt(1-b)/b**2>

      flb = flb + sum(sqrt(one-bfield_tok)
     1    /b2obm_tok**2)/sum_gsqrt_tok

c  plus <sqrt(1-b)/b)>/2

      flb = flb + 0.5_dp*sum(sqrt(one-bfield_tok)
     1    /bfield_tok/b2obm_tok)/sum_gsqrt_tok

c  1+(previous stuff)*<b**2>


      flb = one + flb*b2av

c finally combine upper and lower limits with Houlberg factors
      ft_tok = 0.25_dp*flb + 0.75_dp*fub

      deallocate(bfield_tok, gsqrt_tok, b2obm_tok, one_b2)

      return
      end subroutine tok_fraction
