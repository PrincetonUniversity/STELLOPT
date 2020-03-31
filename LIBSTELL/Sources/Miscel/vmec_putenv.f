      SUBROUTINE vmec_putenv(ename, evalue, ierror)
      IMPLICIT NONE
      INTEGER :: ierror
      CHARACTER(LEN=*), INTENT(in) :: ename, evalue
      CHARACTER(LEN=LEN_TRIM(ename)+LEN_TRIM(evalue)+2) :: temp
      INTEGER, EXTERNAL :: putenv

      temp = TRIM(ename) // "=" // TRIM(evalue) // CHAR(0)

      !ierror = putenv(temp)

      END SUBROUTINE vmec_putenv
