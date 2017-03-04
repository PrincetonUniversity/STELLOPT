      SUBROUTINE chisq_diagnostics (num, nopt, iflag, extension, 
     1                              lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, exe_suffix, remove, v3post_in,
     1                 ndiagno, sigma_diagno, data_diagno, 
     2                 remove
      USE read_v3post_mod
      USE system_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num, iflag
      CHARACTER(LEN=*), INTENT(in) :: extension
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, ierr, kdiag
      CHARACTER(LEN=200) :: cmd, filename
C-----------------------------------------------
      iflag = 0

      IF (nopt > 0) THEN
!
!     THE v3post_in FILE IS THE SAME FOR ALL PROCESSORS, SO NO COPY NEEDS TO BE MADE
!
!DEC$ IF DEFINED(OLDSTYLE_DIAGNO)

         cmd = TRIM(home_dir) // 'xv3post' // exe_suffix // 
     1         TRIM(v3post_in) // ' ' // TRIM(extension) // 
     2         ' > v3dump.' // extension
         CALL system(cmd, ierr)
         IF (ierr.lt.127 .and. ierr.ne.0) THEN
            IF (lscreen) PRINT *,
     1                   'XV3POST failed in chisq_diagnostics call: ',
     2                   'ierr = ',ierr
            iflag = -34
            RETURN
         END IF
         
         cmd = remove // 'v3dump.' // extension
         CALL system(cmd, ierr)

!DEC$ ELSE
         CALL v3post_sub(TRIM(v3post_in), TRIM(extension), ierr)
         IF (ierr .ne. 0) THEN
            PRINT *, 'ierr = ',ierr,' v3post_in file: ',TRIM(extension)
            iflag = -35
            RETURN
         END IF
!DEC$ ENDIF

!
!  Open and read v3post.extension file written out by xv3post (may be .nc file)
!
         CALL read_v3post_file (extension, filename, istat)

         IF (istat .ne. 0) THEN
            IF (lscreen) PRINT *,
     1         'Could not read v3post output file, istat = ', istat
            iflag = -36
            RETURN
         ENDIF

!
!  Match chisq for each diagnostic
!
         IF (SIZE(signal_diag_cext) .ne. ndiagno) THEN
            PRINT *,' ndiagno=',ndiagno,' != SIZE(signal array)=',
     1      SIZE(signal_diag_cext)
            STOP
         END IF
         DO kdiag = 1, ndiagno
            num = num + 1
            wegt(num) = ABS(sigma_diagno(kdiag))
            chisq_match(num) = signal_diag_cext(kdiag) 
     1                       + signal_diag_plasma(kdiag)
            chisq_target(num) = data_diagno(kdiag)
            index_array(num) = ivar_diagnostics
         END DO

         CALL read_v3post_deallocate

      ELSE

         DO kdiag = 1, ndiagno
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_diagnostics)
         END DO
       
      END IF

      END SUBROUTINE chisq_diagnostics
