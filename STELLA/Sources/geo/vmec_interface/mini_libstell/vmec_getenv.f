      SUBROUTINE vmec_getenv(ename, evalue)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: ename, evalue

! MJL 2016-09-22 The rest of this subroutine is commented out because there was an error on my laptop
! associated with pxfgetenv not being found, and this subroutine is not needed for regcoil anyway.

! !DEC$ IF DEFINED (CRAY)
!       INTEGER :: lenname=0, lenval, ierror
!       CALL pxfgetenv(ename, lenname, evalue, lenval, ierror)
! !DEC$ ELSE
!       CALL getenv(ename, evalue)
! !DEC$ ENDIF
      END SUBROUTINE vmec_getenv
