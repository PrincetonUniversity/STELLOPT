      SUBROUTINE vmec_getpid(ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror
#if !defined(CRAY) && !defined(IRIX64) && !defined(IFORT)
!      INTEGER, EXTERNAL :: getpid
      INTEGER :: getpid

      ierror = 0

      ipid = getpid()

      IF (ipid < 0) ierror = -ipid
#else
      CALL pxfgetpid (ipid, ierror)
#endif
      END SUBROUTINE vmec_getpid
