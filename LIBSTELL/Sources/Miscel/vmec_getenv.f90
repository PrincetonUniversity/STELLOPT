      SUBROUTINE vmec_getenv(ename, evalue)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: ename, evalue
#if defined(CRAY)
      INTEGER :: lenname=0, lenval, ierror
      CALL pxfgetenv(ename, lenname, evalue, lenval, ierror)
#else
      CALL getenv(ename, evalue)
#endif
      END SUBROUTINE vmec_getenv
