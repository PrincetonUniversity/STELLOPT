      FUNCTION pwell_t (x)
      USE optim_params, ONLY: rprec, target_well
      IMPLICIT NONE
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, pwell_t
C-----------------------------------------------
      pwell_t = 0

      DO i = UBOUND(target_well,1), LBOUND(target_well,1), -1
         pwell_t = x*pwell_t + target_well(i)
      END DO

      END FUNCTION pwell_t
