      SUBROUTINE chisq_beta (Target_b, sigma, targetW, sigmaW,
     1                       num, nopt)
      USE stel_kinds
      USE chisq_mod
      use optim, ONLY: wp_opt, wb_opt, bigno, lbeta_min,
     1                 vp_opt, pres_opt, nrad
      use vparams, ONLY: twopi
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER, INTENT(inout) :: num
      REAL(rprec), INTENT(in) :: Target_b, sigma, targetW, sigmaW
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      REAL(rprec), PARAMETER :: min_beta = 1.0e-4_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: match, Ekin
C-----------------------------------------------
      IF ( (ABS(sigma) .ge. bigno) .and.
     1     (ABS(sigmaW) .ge. bigno) ) RETURN

      IF (nopt .gt. 0) THEN
         IF (ABS(sigma) < bigno) THEN
            num=num + 1
            index_array(num) = ivar_beta
            match = zero
            IF (wb_opt .gt. zero) match = wp_opt/wb_opt
            chisq_target(num) = Target_b
!
!           matched value of beta can exceed Target_b value with no penalty
!
            IF (lbeta_min) THEN
               chisq_match(num) = MIN(match, Target_b)
            ELSE
               chisq_match(num) = match
            END IF
            wegt(num) = ABS(sigma*MAX(Target_b, min_beta))
         END IF
         
         IF (ABS(sigmaw) < bigno) THEN
            num = num + 1
            index_array(num) = ivar_beta
            ekin = 1.5_rprec * twopi**2 *
     1           sum(vp_opt(2:nrad)*pres_opt(2:nrad))/(nrad-1)
            chisq_match(num) = ekin
            chisq_target(num) = Targetw
            wegt(num) = abs(sigmaw)
         END IF
      ELSE
         IF (ABS(sigma) < bigno) THEN
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_beta)
         ENDIF
         IF (ABS(sigmaw) < bigno) THEN
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_beta)
         ENDIF
      END IF

      END SUBROUTINE chisq_beta
