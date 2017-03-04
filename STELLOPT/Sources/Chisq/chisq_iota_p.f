      SUBROUTINE chisq_iota_p(iota, sigma_max, sigma_min, hs,
     1           num, nrad, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt
      INTEGER :: num
      REAL(rprec), INTENT(in) :: iota(nrad), sigma_max(nrad),
     1    sigma_min(nrad), hs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(kind=rprec), PARAMETER :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad, ncount
      REAL(kind=rprec) :: Target_loc, sj, match
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: piota_prime
C-----------------------------------------------
      ncount = count(ABS(sigma_max(2:nrad-1)) < bigno)
     1       + count(ABS(sigma_min(2:nrad-1)) < bigno)

C      PRINT *,'iota_p ncount=',ncount

      IF (ncount == 0) RETURN

      IF (nopt .gt. 0) THEN
         DO jrad = 2, nrad-1
           IF (sigma_max(jrad)<bigno .or. sigma_min(jrad)<bigno) THEN
             sj = hs*(jrad - c1p5)

             Target_loc = piota_prime(sj)
             match = (iota(jrad+1)-iota(jrad-1))/(2*hs)

c             PRINT *,'iota_p, request, sig_mx, sig_mn=',match,Target_loc,
c     1             sigma_max(jrad), sigma_min(jrad)
           ENDIF

           IF (sigma_max(jrad) < bigno) THEN
             num = num + 1
             index_array(num) = ivar_iota_p
             wegt(num) = sigma_max(jrad)
             chisq_target(num) = Target_loc
             chisq_match(num) = MAX(match, Target_loc)
           ENDIF

           IF (sigma_min(jrad) < bigno) THEN
             num = num + 1
             index_array(num) = ivar_iota_p
             wegt(num) = sigma_min(jrad)
             chisq_target(num) = Target_loc
             chisq_match(num) = MIN(match, Target_loc)
           ENDIF

         END DO
      
      ELSE
         DO jrad = 2, nrad-1
            IF (sigma_max(jrad)<bigno) THEN
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) =
     1                           descript(ivar_iota_p)
            END IF
            IF (sigma_min(jrad) < bigno) THEN
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = 
     1                           descript(ivar_iota_p)
            END IF
         END DO
      ENDIF

      END SUBROUTINE chisq_iota_p
