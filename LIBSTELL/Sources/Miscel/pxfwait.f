      SUBROUTINE pxfwait_g(istat, iretpid, ierror)
      IMPLICIT NONE
      INTEGER :: istat, iretpid, ierror
#if defined(WIN32)
      istat = 0
      iretpid = 0
      ierror = 0
#elif !defined(CRAY) && !defined(IRIX64) && !defined(IFORT)
      INTEGER, EXTERNAL :: wait

      iretpid = 0
      istat = 0

      !ierror = wait(0)
#else
      CALL pxfwait(istat, iretpid, ierror)
#endif

      END SUBROUTINE pxfwait_g
