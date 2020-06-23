      RECURSIVE FUNCTION COMMON_FACTOR(N,m,ifact) RESULT(val)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN) :: ifact
      INTEGER:: val

      val = 0
      IF (mod(N,m)==0) THEN
         val = m
      ELSE
         val = COMMON_FACTOR(N,m+ifact,ifact)
      END IF

      END FUNCTION COMMON_FACTOR