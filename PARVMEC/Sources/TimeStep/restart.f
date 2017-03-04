      SUBROUTINE restart_iter(time_step)
      USE vmec_main
      USE xstuff
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: time_step
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p03 = 1.03_dp, cp90 = 0.90_dp
#if defined(SKS)
      REAL(rprec) :: skston, skstoff
#endif
!-----------------------------------------------
#if defined(SKS)
      CALL second0(skston)
#endif
      IF (PARVMEC) THEN
        SELECT CASE (irst)
        CASE DEFAULT
          pxstore(:neqs2) = pxc(:neqs2)
          RETURN
        CASE (2:3)
          pxcdot(:neqs2) = zero
          pxc(:neqs2) = pxstore(:neqs2)
          time_step = time_step*((irst-2)/c1p03 + cp90*(3-irst))
          IF (irst .eq. 2) THEN
            ijacob = ijacob + 1
            iter1 = iter2
          END IF
          irst = 1
          RETURN
        END SELECT
      ELSE
        SELECT CASE (irst)
        CASE DEFAULT
          xstore(:neqs2) = xc(:neqs2)
          RETURN
        CASE (2:3)
          xcdot(:neqs2) = zero
          xc(:neqs2) = xstore(:neqs2)
          time_step = time_step*((irst-2)/c1p03 + cp90*(3-irst))
          IF (irst .eq. 2) THEN
            ijacob = ijacob + 1
            iter1 = iter2
          END IF
          irst = 1
          RETURN
        END SELECT
      END IF
#if defined(SKS)
      CALL second0(skstoff)
      restart_time = restart_time + (skstoff - skston)
#endif
      END SUBROUTINE restart_iter
