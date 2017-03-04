      SUBROUTINE chisq_rzbdy(Target_x, x_opt, sigma, ivar, num, nopt,
     1    lmin)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt
      INTEGER :: num
      REAL(rprec), INTENT(in) :: Target_x, x_opt, sigma
      LOGICAL :: lmin
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
      IF (ABS(sigma) .ge. bigno) RETURN

      num = num + 1

      IF (nopt .gt. 0) THEN
         index_array(num) = ivar
         wegt(num) = sigma*Target_x
         IF (wegt(num) .eq. zero) wegt(num) = one
         chisq_target(num) = Target_x
         IF (lmin) THEN
            chisq_match(num) = MAX(Target_x, x_opt)
         ELSE
            chisq_match(num) = MIN(Target_x, x_opt)
         END IF
      ELSE
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
      END IF

      END SUBROUTINE chisq_rzbdy
