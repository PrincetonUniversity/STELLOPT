      FUNCTION piota_prime (x)
      USE optim_params, ONLY: rprec, target_iota_p
      IMPLICIT NONE
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, piota_prime
C-----------------------------------------------
      piota_prime = 0

      DO i = UBOUND(target_iota_p,1), LBOUND(target_iota_p,1), -1
         piota_prime = x*piota_prime + target_iota_p(i)
      END DO

      END FUNCTION piota_prime
