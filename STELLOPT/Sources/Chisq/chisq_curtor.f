

      SUBROUTINE chisq_curtor (Target, sigma, num, nrad, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: buco_opt, bigno, curtor_opt
      USE vparams, ONLY: mu0, twopi
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, intent(in) :: nrad, nopt
      INTEGER :: num
      REAL(rprec) :: target, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: curtor
C-----------------------------------------------
      IF (ABS(sigma) .ge. bigno) RETURN
 
      curtor = 0.0
      IF (nopt .gt. 0) THEN
         curtor = twopi * MAXVAL(ABS(buco_opt(2:nrad))) / mu0
         num = num + 1
         index_array(num) = ivar_curtor
         chisq_target(num) = Target
         !chisq_match(num) = curtor
         chisq_match(num) = curtor_opt
         wegt(num) = ABS(sigma)
         chisq_descript(num) = descript(ivar_curtor)
      ELSE
         num = num + 1
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_curtor)
      END IF

      END SUBROUTINE chisq_curtor
