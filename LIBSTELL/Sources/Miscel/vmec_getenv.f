      SUBROUTINE vmec_getenv(ename, evalue)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: ename, evalue
!DEC$ IF DEFINED (CRAY)
      INTEGER :: lenname=0, lenval, ierror
      CALL pxfgetenv(ename, lenname, evalue, lenval, ierror)
!DEC$ ELSE
      CALL getenv(ename, evalue)
!DEC$ ENDIF
      END SUBROUTINE vmec_getenv
