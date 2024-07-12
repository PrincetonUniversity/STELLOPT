      SUBROUTINE pxffork_g (ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror
#if defined(WIN32)
      ierror = 0
      ipid = 0
#elif !defined(CRAY) && !defined(IRIX64) && !defined(IFORT)
      INTEGER, EXTERNAL :: fork

      ierror = 0
      !ipid = fork()
      IF (ipid < 0) ierror = -ipid
#else
      CALL pxffork (ipid, ierror)
#endif
      END SUBROUTINE pxffork_g
