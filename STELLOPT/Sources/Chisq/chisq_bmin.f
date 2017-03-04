      SUBROUTINE chisq_bmin (sigma, ivar,
     1          bmin_b, num, nopt, mskip, sqrt_nsurf, dtheta)
      USE stel_kinds
      USE chisq_mod
      USE boozer_params, ONLY: nu2_b
      USE optim, ONLY: bigno
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt, mskip
      INTEGER :: num
      REAL(rprec) :: bmin_b(*), sigma, sqrt_nsurf, dtheta
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
      REAL(rprec), PARAMETER :: max_variation = 2.e-2_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n
      REAL(rprec) :: avg_bmin
C-----------------------------------------------
      IF (ABS(sigma) .ge. bigno) RETURN

      IF (nopt .gt. 0) THEN

         avg_bmin = SUM(bmin_b(:nu2_b))
     1            - p5*(bmin_b(1) + bmin_b(nu2_b))
         avg_bmin = avg_bmin * dtheta

         DO n = 1, nu2_b, mskip
            num = num + 1
            index_array(num) = ivar
            wegt(num) = max_variation*sigma*avg_bmin*sqrt_nsurf
            chisq_target(num) = avg_bmin
            chisq_match(num) = bmin_b(n)
         END DO

      ELSE
         DO n = 1, nu2_b, mskip
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
         END DO
      END IF

      END SUBROUTINE chisq_bmin
