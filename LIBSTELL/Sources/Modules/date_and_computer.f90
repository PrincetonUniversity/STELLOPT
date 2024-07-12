      MODULE date_and_computer
      USE system_mod, ONLY: getenv
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=3), DIMENSION(12), PARAMETER :: months =            &
        (/ 'Jan','Feb','Mar','Apr','May','Jun',                         &
           'Jul','Aug','Sep','Oct','Nov','Dec' /)
      CHARACTER(LEN=9), PARAMETER :: finfo='VInfo.txt'
      CHARACTER(LEN=100) :: computer, os, os_release

#if defined(DARWIN)
      CHARACTER(LEN=2), PARAMETER :: os_flag = '-s'
#else
      CHARACTER(LEN=2), PARAMETER :: os_flag = '-o'
#endif
     
      CONTAINS

      SUBROUTINE GetComputerInfo
      CHARACTER(LEN=100) :: temp
#if defined(WIN32)
      computer = ' Window_NT'
      os       = ' MS Windows 2000'
      os_release  = ' 5.00'
#else
      OPEN(unit=10101,file=FINFO,status='replace')
      temp = "hostname >> " // FINFO // char(0)
      CALL system(temp)
      READ(10101,'(a)') computer
      CLOSE(10101,status='delete')
      OPEN(unit=10101,file=FINFO,status='replace')
      temp = "uname " // os_flag // " >> " // FINFO // char(0)
      CALL system(temp)
      READ(10101,'(a)') os
      CLOSE(10101,status='delete')
      OPEN(unit=10101,file=FINFO,status='replace')
      temp = "uname -r >> " // FINFO // char(0)
      CALL system(temp)
      READ(10101,'(a)') os_release
      CLOSE(10101,status='delete')
#endif
      END SUBROUTINE GetComputerInfo

      END MODULE date_and_computer
