      SUBROUTINE chisq_bpres (pres_opt, sigma, ivar, num, nrad,
     1   nopt, nflag)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, bigno
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nrad, nopt, nflag
      INTEGER :: num
      REAL(rprec) :: pres_opt(*), sigma(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad
      REAL(rprec) :: balpgrad_opt
C-----------------------------------------------
!
!     MATCH EDGE PRESSURE=0
!

      IF ((nflag .eq. 0) .and. (sigma(1) .ge. bigno)) RETURN
      IF ((nflag .gt. 0) .and. ANY(sigma(3:nrad) .ge. bigno)) RETURN
      IF (nflag .eq. 0) THEN
         num = num + 1
         IF (nopt .gt. 0) THEN
            index_array(num) = ivar
            wegt(num) = sigma(1)
            chisq_target(num) = zero
            chisq_match(num) = pres_opt(nrad)
         ELSE
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
         END IF

      ELSE
!
!     MATCH MONOTONICALLY DECREASING PRESSURE PROFILE  (dP/ds<0)
!
         IF (nopt. gt. 0) THEN
            DO jrad = 3, nrad
               num = num + 1
               index_array(num) = ivar
               wegt(num) = sigma(jrad)
               balpgrad_opt = (pres_opt(jrad) - pres_opt(jrad-1))
     1             /REAL(nrad-1, rprec)
               chisq_target(num) = zero
               chisq_match(num) = MAX(balpgrad_opt, zero)
            END DO
         ELSE
            DO jrad = 3, nrad
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
            END DO
         ENDIF
      END IF

      END SUBROUTINE chisq_bpres
