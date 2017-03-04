      SUBROUTINE pxfwait_g(istat, iretpid, ierror)
      IMPLICIT NONE
      INTEGER :: istat, iretpid, ierror
!DEC$ IF DEFINED (WIN32)
      istat = 0
      iretpid = 0
      ierror = 0
!DEC$ ELSEIF .NOT.DEFINED (CRAY) .AND. .NOT.DEFINED(IRIX64) .AND. .NOT.DEFINED(IFORT)
!      INTEGER, EXTERNAL :: wait

      iretpid = 0
      istat = 0

!      ierror = wait(0)
!DEC$ ELSE 
!      CALL pxfwait(istat, iretpid, ierror)
!DEC$ ENDIF

      END SUBROUTINE pxfwait_g
