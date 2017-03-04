

      SUBROUTINE chisq_jedge (Target, sigma, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno
      USE vmec_input, ONLY: ac
      USE vparams, ONLY: twopi
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, intent(in) :: nopt
      INTEGER :: num
      REAL(rprec) :: target, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: sum1
C-----------------------------------------------
      IF (ABS(sigma) < bigno) THEN
         num = num + 1
         IF (nopt .gt. 0) THEN
            sum1 = sum(ac(0:10))

            index_array(num) = ivar_jedge
            chisq_target(num) = Target
            chisq_match(num) = sum1
            wegt(num) = ABS(sigma)
         ELSE
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_jedge)
         END IF
      END IF

      END SUBROUTINE chisq_jedge
