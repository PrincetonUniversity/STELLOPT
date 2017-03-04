      SUBROUTINE chisq_extcur(sigma, target, ivar, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, lextcur, nextcur_vmec
      USE vmec_input, ONLY: lfreeb, extcur
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt
      INTEGER :: num
      REAL(rprec) :: sigma(*), target(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, n
C-----------------------------------------------
!     THIS ATTEMPTS TO MAKE THE COIL CURRENTS SMALL (ZERO) IF POSSIBLE
      n = SIZE(lextcur)

      IF (ALL(ABS(sigma(:n)) >= bigno .or. .not.lextcur(:n))
     1       .or. .not. lfreeb) RETURN

         DO j = 1, n
            IF( ABS(sigma(j)) < bigno .and. lextcur(j)) THEN
               num = num + 1
               IF (nopt .gt. 0) THEN
                  index_array(num) = ivar

                  IF (j <= nextcur_vmec) THEN
                     wegt(num) = sigma(j)
                     chisq_match(num) = extcur(j)
                     chisq_target(num) = target(j)
                  ELSE
                     chisq_match(num) = zero
                  END IF
!                  chisq_target(num) = zero
               ELSE
                  IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
               ENDIF
            END IF
         END DO

      END SUBROUTINE chisq_extcur
