      SUBROUTINE chisq_polflux (Target_in, sigma, hs, nrad, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: iota_opt, phip_opt
      USE vparams, ONLY: twopi
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt
      INTEGER, INTENT(inout) :: num
      REAL(rprec), INTENT(in) :: target_in, sigma, hs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: psi_a
C-----------------------------------------------
      num = num + 1
      IF (nopt .gt. 0)then
         psi_a = ABS(-twopi * hs *
     1               SUM(iota_opt(2:nrad)*phip_opt(2:nrad)) )

         index_array(num) = ivar_fluxp
         chisq_match(num) = psi_a
         chisq_target(num) = MAX(psi_a, target_in)
         wegt(num) = sigma * ABS(target_in)
      ELSE
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_fluxp)
      END IF

      END SUBROUTINE chisq_polflux
