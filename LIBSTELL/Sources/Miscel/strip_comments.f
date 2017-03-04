      SUBROUTINE strip_comments(input_file,output_file)
      USE safe_open_mod
      IMPLICIT NONE
!
!     strips comment lines (starting with '!') from input_file
!     renames clean file input_file // '.stripped'
!
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*) :: input_file
      CHARACTER(LEN=*) :: output_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_strip = 20
      INTEGER :: istat, iustrip, iunew, itry
      CHARACTER(LEN=1024) :: line
      LOGICAL :: lex
      INTEGER, EXTERNAL :: getcwd
      

C-----------------------------------------------
      itry = 1
 99   INQUIRE (file=input_file, exist=lex, iostat=istat)
      IF (.not.lex .and. itry .lt. 200) THEN
         itry = itry + 1
         GOTO 99
      END IF
      IF (istat.ne.0 .or. .not. lex) THEN
!         istat = getcwd(line)
         PRINT *,TRIM(input_file),' does not exist in directory: ',
     1           TRIM(line)
         CALL FLUSH(6)
         STOP
      END IF

      iustrip = unit_strip
      CALL safe_open(iustrip, istat, input_file, 'old','formatted')
      IF (istat .ne. 0) THEN
        line = 'error opening ' // TRIM(input_file) //
     1         ' in STRIP_COMMENTS'
        PRINT *, line
        PRINT *,'istat = ', istat
        STOP
      END IF
      iunew = iustrip + 1
      CALL safe_open(iunew, istat, TRIM(output_file),
     1     'replace', 'formatted')
!      CALL safe_open(iunew, istat, TRIM(input_file) // '.stripped',
!     1     'replace', 'formatted')
      IF (istat .ne. 0) THEN
        line = 'error opening ' // TRIM(output_file) //
     1      ' in STRIP_COMMENTS'
!        line = 'error opening ' // TRIM(input_file) //
!     1      '.stripped  in STRIP_COMMENTS'
        PRINT *, line
        PRINT *,'istat = ', istat
        STOP
      END IF
      DO
         READ(iustrip, '(a)', END=100) line
         line = ADJUSTL(line)
         IF (line(1:1) == '!') CYCLE
         WRITE(iunew, '(a)') TRIM(line)
      END DO
      
      CALL FLUSH(iunew)
          

 100  CONTINUE

      CLOSE(iustrip)
      CLOSE(iunew)

      END SUBROUTINE strip_comments
