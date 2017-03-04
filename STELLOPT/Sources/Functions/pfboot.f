      FUNCTION pfboot (x)
      USE optim_params, ONLY: rprec, fboot
      IMPLICIT NONE
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, pfboot
C-----------------------------------------------
      pfboot = 0

      DO i = UBOUND(fboot,1), LBOUND(fboot,1), -1
         pfboot = x*pfboot + fboot(i)
      END DO

      END FUNCTION pfboot
