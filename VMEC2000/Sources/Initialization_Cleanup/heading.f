      SUBROUTINE heading(extension, time_slice, iseq_count, lmac,
     1                   lscreen, lwrite)
      USE vmec_main, ONLY: rprec
      USE vparams, ONLY: nthreed, nmac
      USE vmec_params, ONLY: version_
      USE date_and_computer
      USE parallel_include_module, ONLY: grank
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iseq_count
      REAL(rprec) :: time_slice
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lmac, lscreen, lwrite
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=100), PARAMETER ::
     &   banner = ' THIS IS PARVMEC (PARALLEL VMEC), VERSION '
      CHARACTER(LEN=*), PARAMETER :: VersionID1 =
     &   ' Lambda: Full Radial Mesh. L-Force: hybrid full/half.'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: imon, nout
      CHARACTER(LEN=10) :: date0, time0, zone0
      CHARACTER(LEN=50) :: dateloc, Version
      LOGICAL :: lfirst
C-----------------------------------------------
!
!     Open output files
!
      IF (grank .NE. 0) THEN
         lscreen = .FALSE.
      END IF

      CALL open_output_files(extension, iseq_count, lmac, lscreen,
     &                       lfirst, lwrite)

      IF (.NOT.lfirst .OR. .NOT.lwrite) RETURN

!     FORTRAN-90 ROUTINE
      CALL DATE_AND_TIME(date0, time0, zone0)
      READ(date0(5:6),'(i2)') imon
      WRITE(dateloc,1000) months(imon), date0(7:8), date0(1:4),
     &   time0(1:2), time0(3:4), time0(5:6)

      IF (lscreen) THEN
         CALL GetComputerInfo

         Version = TRIM(ADJUSTL(version_))
         WRITE (nthreed,1002) TRIM(banner), TRIM(Version),
     &                        TRIM(VersionID1), TRIM(computer),
     &                        TRIM(os), TRIM(os_release), TRIM(dateloc)

         IF (lfirst) THEN
            WRITE (*,1001) iseq_count + 1, time_slice, TRIM(extension)
            WRITE (*,1003) TRIM(banner), TRIM(Version),
     &                     TRIM(VersionID1), TRIM(computer), TRIM(os),
     &                     TRIM(os_release), TRIM(dateloc)
         END IF
      ENDIF

      DO nout = nthreed, nthreed + 1
         imon = nout
         IF (imon .eq. nthreed + 1) THEN
            imon = nmac
         END IF
         IF (imon .eq. nmac .and. .not.lmac) CYCLE
         WRITE (imon,1004) TRIM(extension), iseq_count, time_slice
      END DO

1000  FORMAT('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)
1001  FORMAT('  SEQ = ',i4,' TIME SLICE',1p,e12.4/
     &       '  PROCESSING INPUT.',a)
1002  FORMAT(a,1x,a,/a,//,' COMPUTER: ',a,2x,' OS: ',a,2x,
     &       ' RELEASE: ',a,2x,a)
1003  FORMAT(1x,a,1x,a,/1x,a,//,'  COMPUTER: ',a,2x,' OS: ',a,2x,
     &       ' RELEASE: ',a,2x,a)
1004    FORMAT(' SHOT ID.: ',a,2x,'SEQ. NO.:',i4,/,
     &         ' TIME SLICE = ',f5.0,' ms')

      END SUBROUTINE heading
