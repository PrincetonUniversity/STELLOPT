      SUBROUTINE vmec_getpid(ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror
!DEC$ IF .NOT.DEFINED (CRAY) .AND. .NOT.DEFINED(IRIX64) .AND. .NOT.DEFINED(IFORT)
C      INTEGER, EXTERNAL :: getpid
      INTEGER :: getpid

      ierror = 0

      ipid = getpid()

      IF (ipid < 0) ierror = -ipid
!DEC$ ELSE
      CALL pxfgetpid (ipid, ierror)
!DEC$ ENDIF
      END SUBROUTINE vmec_getpid
