      SUBROUTINE chisq_iota(match, sigma, hs, num, nrad, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno
      USE optim_params, ONLY: target_iota_max, target_iota_min,
     1                        sigma_iota_max, sigma_iota_min,
     2                        target_iota_max_min, sigma_iota_max_min
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt
      INTEGER :: num
      REAL(rprec) :: match(*), sigma(*), hs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad
      REAL(rprec) :: TargetIota, AvgIota, sj
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: piota_t
C-----------------------------------------------
      IF (ANY(ABS(sigma(2:nrad)) .lt. bigno)) THEN
        IF (nopt .gt. 0) THEN
          avgiota = zero
          DO jrad = 2, nrad
            sj = hs*(jrad - c1p5)
            AvgIota = AvgIota + ABS(piota_t(sj))
          END DO
          AvgIota = AvgIota/(nrad - 1)
          DO jrad = 2, nrad
             num = num + 1
             index_array(num) = ivar_iota
             sj = hs*(jrad - c1p5)
             TargetIota = piota_t(sj)
             wegt(num) = sigma(jrad)*AvgIota
             chisq_target(num) = TargetIota
             chisq_match(num) = match(jrad)
          END DO

        ELSE
          DO jrad = 2, nrad
             num = num + 1
             IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_iota)
          END DO
        ENDIF
      ENDIF

      IF (sigma_iota_max < bigno) THEN
         num = num+1
         IF( nopt > 0) THEN
            index_array(num) = ivar_iota_bounds
            wegt(num) = sigma_iota_max
            chisq_target(num) = target_iota_max
            chisq_match(num) = MAX(target_iota_max,
     1                             MAXVAL(match(2:nrad)))
         ELSE
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_iota_bounds)
         ENDIF
      ENDIF

      IF (sigma_iota_max_min < bigno) THEN
         num = num+1
         IF( nopt > 0) THEN
            index_array(num) = ivar_iota_bounds
            wegt(num) = sigma_iota_max_min
            chisq_target(num) = target_iota_max_min
            chisq_match(num) = MAX(target_iota_max_min,
     1                             MAXVAL(match(2:nrad)))
         ELSE
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_iota_bounds)
         ENDIF
      ENDIF

      IF (sigma_iota_min < bigno) THEN
         num = num+1
         IF( nopt > 0) THEN
            index_array(num) = ivar_iota_bounds
            wegt(num) = sigma_iota_min
            chisq_target(num) = target_iota_min
            chisq_match(num) = MIN(target_iota_min,
     1                             MINVAL(match(2:nrad)))
         ELSE
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_iota_bounds)
         ENDIF
      ENDIF

      END SUBROUTINE chisq_iota
