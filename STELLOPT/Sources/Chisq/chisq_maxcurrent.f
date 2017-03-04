      SUBROUTINE chisq_maxcurrent(target_in, sigma, num, nrad, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: buco_opt, bigno
      USE vparams, ONLY: mu0, twopi
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt
      INTEGER, INTENT(inout) :: num
      REAL(rprec), INTENT(in) :: target_in, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad, jcount, i
      REAL(rprec) :: dnorm, dmax_j
C-----------------------------------------------
      IF (ABS(sigma) .ge. bigno) RETURN

      jcount = 0

      IF (nopt .gt. 0) THEN
        DO jrad = 2,nrad
          num = num + 1
          index_array(num) = ivar_maxj
          dmax_j = twopi * buco_opt(jrad) / mu0
          wegt(num) = ABS(sigma*target_in)
          chisq_match(num) = dmax_j
          IF (ABS(dmax_j) .gt. ABS(target_in) + wegt(num)) THEN
            chisq_target(num) = SIGN(target_in, dmax_j)
            jcount = jcount + 1
          ELSE
            chisq_target(num) = chisq_match(num)
          END IF
        END DO
      ELSE
        DO jrad = 2,nrad
          num = num + 1
          IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_maxj)
        END DO
      ENDIF


      IF (jcount .gt. 1) THEN
        i = num - (nrad-2)
        dnorm = SQRT(REAL(jcount,rprec))
        wegt(i:num) = wegt(i:num) * dnorm
      END IF

      END SUBROUTINE chisq_maxcurrent
