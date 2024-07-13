      SUBROUTINE vmec_system(cmd, ierror)
      INTEGER, OPTIONAL :: ierror
      INTEGER :: ireturn
      CHARACTER(LEN=*), INTENT(in) :: cmd

#if defined(CRAY)
      INTEGER, EXTERNAL :: ishell
      !ireturn = ishell(TRIM(cmd))
#elif defined(RISC)
      !CALL system(TRIM(cmd), ireturn)
#elif defined(IRIX64)
      !CALL system(TRIM(cmd))
      ireturn = 0
#elif defined(LINUX) || defined(OSF1) || defined(DARWIN)
!      INTEGER, EXTERNAL :: system
!      INTEGER :: system
!      ireturn = system(TRIM(cmd))
#elif defined(WIN64)
      CALL SYSTEM(TRIM(cmd), ireturn)
#elif defined(WIN32) || defined(SUNOS)
      INTEGER, EXTERNAL :: system
      ireturn = system(TRIM(cmd))
#else
!      INTEGER, EXTERNAL :: system
      CHARACTER(LEN=LEN_TRIM(cmd)+1) :: cmd1
      cmd1 = TRIM(cmd) // CHAR(0)
!      ireturn = system(TRIM(cmd1))
#endif
      IF (PRESENT(ierror)) ierror = ireturn

      END SUBROUTINE vmec_system
