      SUBROUTINE myfork(i, maxprocess, wrapper, fcn)
      USE system_mod
      IMPLICIT NONE      
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER i, maxprocess
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL wrapper, fcn
!DEC$ IF .NOT.DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: pid, status, iretpid, ierror, werror
      INTEGER, SAVE :: nprocess = 0
C-----------------------------------------------
 
      IF (i .eq. 1) nprocess = 0
      ierror = -1
 
!     Child process: limit number to max_process to avoid potential system hang-up
 
      DO WHILE(ierror .ne. 0)
 
         IF (nprocess .lt. maxprocess) CALL pxffork (pid, ierror)
 
         IF (ierror .ne. 0) THEN
!           wait for next available processor
            CALL pxfwait (status, iretpid, werror)
!           IF (status.gt.0 .and. nprocess.ge.1) THEN
            IF (nprocess .ge. 1) nprocess = nprocess - 1
!           ELSE
!              nprocess = 0
!           END IF
         END IF
      END DO
 
      IF (pid .eq. 0) THEN
         CALL wrapper (i, fcn)
!DEC$ IF DEFINED (CRAY)
         CALL EXIT(1)
!DEC$ ELSE
         STOP
!DEC$ ENDIF
      END IF
 
      nprocess = nprocess + 1
!DEC$ ENDIF
      END SUBROUTINE myfork
