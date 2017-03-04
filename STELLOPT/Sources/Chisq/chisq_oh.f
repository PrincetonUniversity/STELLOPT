      SUBROUTINE chisq_oh(sigma, ivar, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, lextcur, oh_coefs, nextcur_vmec
      USE vmec_input, ONLY: lfreeb, extcur
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt
      INTEGER :: num
      REAL(rprec) :: sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: sum0
C-----------------------------------------------

      IF (ABS(sigma) >= bigno .or. .not. lfreeb) RETURN
      num = num + 1

      IF (nopt .gt. 0) THEN

         sum0 = SUM(oh_coefs(:nextcur_vmec)*extcur(:nextcur_vmec))

         index_array(num) = ivar

         wegt(num) = sigma
         chisq_target(num) = zero
         chisq_match(num) = sum0

      ELSE
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
      ENDIF

      END SUBROUTINE chisq_oh
