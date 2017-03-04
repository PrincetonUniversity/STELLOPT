      FUNCTION pseed_t (x)
      USE optim_params, ONLY: rprec, aseedcur
      IMPLICIT NONE
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, pseed_t
C-----------------------------------------------
      pseed_t = 0

      DO i = UBOUND(aseedcur,1), LBOUND(aseedcur,1), -1
         pseed_t = x*pseed_t + aseedcur(i)
      END DO

      END FUNCTION pseed_t
