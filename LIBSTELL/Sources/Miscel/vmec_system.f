      SUBROUTINE vmec_system(cmd, ierror)
      INTEGER, OPTIONAL :: ierror
      INTEGER :: ireturn
      CHARACTER(LEN=*), INTENT(in) :: cmd

!DEC$ IF DEFINED (CRAY)
      INTEGER, EXTERNAL :: ishell
      ireturn = ishell(TRIM(cmd))
!DEC$ ELSEIF DEFINED (RISC)
      CALL system(TRIM(cmd), ireturn)
!DEC$ ELSEIF DEFINED (IRIX64)
      CALL system(TRIM(cmd))
      ireturn = 0
!DEC$ ELSEIF DEFINED (LINUX) .OR. DEFINED(OSF1) .OR. DEFINED(MACOSX)
C      INTEGER, EXTERNAL :: system
      INTEGER :: system
      ireturn = system(TRIM(cmd))
!DEC$ ELSEIF DEFINED(WIN32) .OR. DEFINED(SUNOS)
      INTEGER, EXTERNAL :: system
      ireturn = system(TRIM(cmd))
!DEC$ ELSE
      INTEGER, EXTERNAL :: system
      ireturn = system(TRIM(cmd) // CHAR(0))
!DEC$ ENDIF
      IF (PRESENT(ierror)) ierror = ireturn

      END SUBROUTINE vmec_system
