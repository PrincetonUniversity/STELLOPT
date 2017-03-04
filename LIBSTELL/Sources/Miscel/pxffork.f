      SUBROUTINE pxffork_g (ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror
!DEC$ IF DEFINED (WIN32)
      ierror = 0
      ipid = 0
!DEC$ ELSEIF .NOT.DEFINED (CRAY) .AND. .NOT.DEFINED(IRIX64) .AND. .NOT.DEFINED(IFORT)
!      INTEGER, EXTERNAL :: fork

      ierror = 0
!      ipid = fork()
      IF (ipid < 0) ierror = -ipid
!DEC$ ELSE
      CALL pxffork (ipid, ierror)
!DEC$ ENDIF
      END SUBROUTINE pxffork_g
