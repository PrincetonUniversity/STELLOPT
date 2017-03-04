      RECURSIVE FUNCTION factorial(num) RESULT(fact)
      USE stel_kinds
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: num
      REAL(rprec) :: fact
      INTEGER:: nfact

      IF(num == 1 .or. num == 0) THEN
         nfact = 1
      ELSE
         nfact = num*factorial(num-1)
      END IF

      fact = nfact

      END FUNCTION factorial
