      FUNCTION piota_t (x)
      USE optim_params, ONLY: rprec, target_iota
      IMPLICIT NONE
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, piota_t
C-----------------------------------------------
      piota_t = 0

      DO i = UBOUND(target_iota,1), LBOUND(target_iota,1), -1
         piota_t = x*piota_t + target_iota(i)
      END DO

      END FUNCTION piota_t
