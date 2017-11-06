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
#if defined(SKS)
     1   banner = 
     2   ' THIS IS PARVMEC (PARALLEL VMEC), VERSION ' 
#else
     1   banner = ' THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION '
#endif
      CHARACTER(LEN=*), PARAMETER :: VersionID1 =
     1   ' Lambda: Full Radial Mesh. L-Force: hybrid full/half.'
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
      IF (grank .NE. 0) lscreen=.FALSE.

      CALL open_output_files (extension, iseq_count, lmac, lscreen,
     1                        lfirst, lwrite)

      IF (.NOT.lfirst .OR. .NOT.lwrite) RETURN

!     FORTRAN-90 ROUTINE
      CALL DATE_AND_TIME(date0,time0,zone0)
      READ(date0(5:6),'(i2)')imon
      WRITE(dateloc,100)months(imon),date0(7:8),date0(1:4),
     1  time0(1:2),time0(3:4),time0(5:6)
 100  FORMAT('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

      IF (lscreen) THEN
         CALL GetComputerInfo

         IF (lfirst) WRITE (*,'(a,i4,a,1p,e12.4/2a)')
     1   '  SEQ = ', iseq_count+1,
     2   ' TIME SLICE',time_slice,'  PROCESSING INPUT.', TRIM(extension)

         Version = TRIM(ADJUSTL(version_))
         WRITE(nthreed,'(a,1x,a,/,a,//,3(2a,2x),a)') TRIM(banner),
     1     TRIM(Version), TRIM(VersionID1), 
     2     ' COMPUTER: ', TRIM(computer), ' OS: ', TRIM(os),
     3     ' RELEASE: ', TRIM(os_release), TRIM(dateloc)
         IF (lfirst)
     1   WRITE (*,'(1x,a,1x,a,/,1x,a,//,1x,3(2a,2x),a)') TRIM(banner),
     2   TRIM(Version), TRIM(VersionID1), 
     3     ' COMPUTER: ', TRIM(computer), ' OS: ', TRIM(os),
     4     ' RELEASE: ', TRIM(os_release), TRIM(dateloc)
      ENDIF

      DO nout = nthreed, nthreed+1
        imon = nout
        IF (imon .eq. nthreed+1) imon = nmac
        IF (imon.eq.nmac .and. .not.lmac) CYCLE
        WRITE (imon,3) TRIM(extension),iseq_count,time_slice
      ENDDO

 3    FORMAT(' SHOT ID.: ',a,2x,'SEQ. NO.:',i4,/,
     1       ' TIME SLICE = ',f5.0,' ms')

      END SUBROUTINE heading
