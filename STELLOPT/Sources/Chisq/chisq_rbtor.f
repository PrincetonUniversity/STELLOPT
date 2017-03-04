      SUBROUTINE chisq_rbtor(match, target_rb, sigma, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num
      REAL(rprec) :: match, target_rb, sigma
C-----------------------------------------------
!
!     Use this to match R*Btor in either fixed boundary
!     of free boundary mode (in which CASE, coil currents
!     must SUM correctly to give effective poloidal current)
!
      IF (ABS(sigma) .ge. bigno) RETURN

      num = num + 1
      IF (nopt .gt. 0) THEN
         index_array(num) = ivar_rbtor
         wegt(num) = sigma
         chisq_target(num) = target_rb
         chisq_match(num) = match
      ELSE
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_rbtor)
      END IF

      END SUBROUTINE chisq_rbtor
