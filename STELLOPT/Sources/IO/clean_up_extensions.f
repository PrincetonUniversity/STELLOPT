      SUBROUTINE clean_up_extensions (extension, minext)
      USE optim, ONLY: list, move
      USE system_mod
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*) :: extension, minext
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit=20, istat, indx
      CHARACTER(LEN=300) :: line, file_name, temp, temp2
C-----------------------------------------------
!
!     AT END OF RUN, REPLACES ALL THE _OPT FILE EXTENSIONS WITH THE FINAL .MIN EXTENSION
!     FIRST FIND (POSSIBLY) INCORRECT FILE EXTENSION IN THREED1 FILE
!
      file_name = 'threed1.' // TRIM(extension)
      CALL safe_open(iunit, istat, file_name, 'old', 'formatted')
      IF (istat .ne. 0) RETURN

      indx = 0

      DO WHILE (istat.eq.0)
         READ (iunit, '(a)', iostat=istat, END=100) line
         indx = INDEX(line, 'SHOT ID')
         IF (indx .ne. 0) EXIT
      END DO

 100  CONTINUE
      CLOSE (unit=iunit)
      IF (indx .eq. 0) RETURN

      temp = line(indx+10:)
      indx = INDEX(temp, '_opt', back=.true.)
      IF (indx .eq. 0) RETURN

      indx = indx - 1
      line = temp(indx:)
      istat = INDEX(temp,' ')
      IF (istat .ne. 0) line = temp(indx:istat-1)

!     Make list of (text) file names to check and replace _opt with .min extension
!     Skip boozmn file, wout*.nc, which are binary files

      temp = list // "*." // TRIM(extension) // " > ls_file"
      CALL system(temp, istat)
      CALL safe_open(iunit, istat, "ls_file", "old", "formatted")
      DO WHILE (istat.eq.0 .or. istat.ge.127)
         READ (iunit, '(a)', END=150) file_name
!        skip binary files
         temp2 = ADJUSTL(file_name)
         IF (temp2(1:7) == "boozmn.".or.
     1       temp2(1:4) == "wout") CYCLE
!DEC$ IF .NOT.DEFINED(WIN32)
         temp = '/bin/sed -e "s/' // TRIM(line) // '/'
     1          // line(1:1) // TRIM(minext) // '/g" '
     2          // TRIM(file_name) // ' > tempxyz; ' 
     3          // move // 'tempxyz ' // TRIM(file_name)
         CALL system(temp, istat)
!DEC$ ENDIF
      END DO
 150  CONTINUE

      CLOSE (iunit, status='delete')

      END SUBROUTINE clean_up_extensions
