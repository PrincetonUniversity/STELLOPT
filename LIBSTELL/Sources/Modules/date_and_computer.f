      MODULE date_and_computer
      USE system_mod, ONLY: system, getpid
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=3), DIMENSION(12), PARAMETER :: months =
     1  ( / 'Jan','Feb','Mar','Apr','May','Jun',
     2      'Jul','Aug','Sep','Oct','Nov','Dec' / )
      CHARACTER(LEN=100) :: computer, os, os_release
      CONTAINS

      SUBROUTINE GetComputerInfo
!DEC$ IF DEFINED (WIN32)
      computer = ' Window_NT'
      os       = ' MS Windows 2000'
      os_release  = ' 5.00'
!DEC$ ELSE
      INTEGER :: ierror, ipid, iunit=10
      CHARACTER(LEN=200) :: fileId

!     Get unique unit number from pid
!      CALL getpid(ipid, ierror)
!      iunit = iunit + ipid
!      IF (iunit .gt. 100000) iunit = iunit-100000
!      WRITE (fileId,'(a,i6.6)') "fort.", iunit

!      CALL system('uname -a > ' // TRIM(fileId), ierror)
!      IF (ierror .eq. 0) THEN
!         CALL safe_open(iunit, ierror, fileId, 'old','formatted')
!         IF (ierror .eq. 0) THEN
!            READ (iunit, *) os, computer, os_release
!            CLOSE(iunit, status='delete', iostat=ierror)
!         END IF
!      END IF
!      IF (ierror .ne. 0) THEN
!         os = 'Unknown'
!         computer = 'Unknown'
!         os_release = '0'
!      END IF
!     SAL - USE GNU F77 Intrisic Procedures
      CALL GETENV('HOST',computer)
      CALL GETENV('OSTYPE',os)
      CALL GETENV('HOSTTYPE',os_release)
!      os_release = '0'
!DEC$ ENDIF

      END SUBROUTINE GetComputerInfo

      END MODULE date_and_computer
