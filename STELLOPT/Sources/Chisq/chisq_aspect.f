      SUBROUTINE chisq_aspect(match, target_a, sigma, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, laspect_max
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num
      REAL(rprec), INTENT(in) :: match, target_a, sigma
C-----------------------------------------------
      IF (ABS(sigma) .ge. bigno) RETURN

      num = num + 1
      IF (nopt .gt. 0) THEN
         index_array(num) = ivar_aspect
         wegt(num) = sigma
         chisq_target(num) = target_a
!
!        IF (laspect_max), ONLY contribute to chisq IF match > target_a
!       (match < target_a is permitted without penalty)
!
         IF (laspect_max) THEN
            chisq_match(num) = MAX(match, target_a)
         ELSE
            chisq_match(num) = match
         END IF

      ELSE 
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_aspect)
      END IF

      END SUBROUTINE chisq_aspect
