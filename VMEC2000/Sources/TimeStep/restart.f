      SUBROUTINE restart_iter(time_step)
      USE vmec_main
      USE xstuff
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: CopyLastNType, ZeroLastNtype
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: time_step
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), PARAMETER :: c1p03 = 1.03_dp, cp90 = 0.90_dp
      REAL(dp)            :: treston, trestoff
!-----------------------------------------------
      CALL second0(treston)

      IF (PARVMEC) THEN
         SELECT CASE (irst)
            CASE DEFAULT
               CALL CopyLastNType(pxc, pxstore)
               RETURN
            CASE (2:3)
               CALL ZeroLastNtype(pxcdot)
               CALL CopyLastNType(pxstore, pxc)
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

      CALL second0(trestoff)
      restart_time = restart_time + (trestoff - treston)

      END SUBROUTINE restart_iter
