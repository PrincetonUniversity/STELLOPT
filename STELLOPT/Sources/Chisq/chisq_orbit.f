      SUBROUTINE chisq_orbit (sigma, ivar, num,
     1                        nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, link
      USE optim_params, ONLY: lorbit_opt
      USE system_mod, ONLY: system
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ivar, nopt, iflag, num
      REAL(rprec)  :: sigma
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_orbit = 79
      INTEGER :: iunit, ierr, istat
      REAL(rprec) :: transit_completed, transit_asked,
     1               average_loss_time
      INTEGER :: nzterm, np_start, np_escape
      REAL(rprec), PARAMETER :: zero = 0
      CHARACTER(LEN=200) :: version
      LOGICAL :: ex
C-----------------------------------------------

      IF (nopt > 0) THEN

         iunit = unit_orbit

         IF (lscreen) THEN
           WRITE(6,*)
     1           'Running ORBIT for Fast Ion Loss Calculations'
         ENDIF

         ierr = 0
         version = TRIM(home_dir) // 'xmkjmc'
         version = TRIM(version) // ' ' // extension
         CALL system(version,ierr)
         IF( ierr .lt. 127 .and. ierr .ne. 0 ) THEN
            iflag = -28
            RETURN
         ENDIF

         version = TRIM(home_dir) // 'xeq3d'
         version = TRIM(version) // ' -b ' // extension
         CALL system(version,ierr)
         IF( ierr .lt. 127 .and. ierr .ne. 0 ) THEN
            iflag = -29
            RETURN
         ENDIF

         version = TRIM(home_dir) // 'xorbit3d'
         version = TRIM(version) // ' -b ' // extension
         CALL system(version,ierr)
         IF( ierr .lt. 127 .and. ierr .ne. 0 ) THEN
            iflag = -30
            RETURN
         ENDIF

         CALL safe_open(iunit, istat,'orbsum.'//TRIM(extension),
     1                  'unknown', 'formatted')
         IF( istat .ne. 0 ) THEN
            iflag = -31
            RETURN
         ELSE
            READ(unit_orbit,*) nzterm, np_start, np_escape,
     >                   transit_completed, transit_asked,
     >                   average_loss_time

            CLOSE(unit_orbit)
         ENDIF
!
         num = num + 1
         index_array(num) = ivar
         chisq_match(num) =  average_loss_time
c            chisq_match(num) = float(np_escape)
c     >                         /float(np_start)
c     >                         /transit_completed
c     >                         *transit_asked
         chisq_target(num) =  zero
         wegt(num) = sigma

      ELSE

         version = TRIM(home_dir) // 'xmkjmc'
         INQUIRE(file=version, exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            PRINT *, 'xmkjmc not found in ' // TRIM(home_dir)
            lorbit_opt = .false.
         ENDIF

         ierr=0
         version = link // "../mkjmc_param.in ./mkjmc_param.in"
         CALL system(TRIM(version), ierr)
         IF (ierr .lt. 127 .and. ierr .ne. 0) THEN
            PRINT *,' MKJMC control file mkjmc_param.in not linked'
            lorbit_opt = .false.
         ENDIF

         version = TRIM(home_dir) // 'xeq3d'
         INQUIRE(file=version, exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            PRINT *, 'xeq3d not found in ' // TRIM(home_dir)
            lorbit_opt = .false.
         ENDIF

         version = TRIM(home_dir) // 'xorbit3d'
         INQUIRE(file=version, exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            PRINT *, 'xorbit3d not found in ' // TRIM(home_dir)
            lorbit_opt = .false.
         ENDIF

         version = link // "../orbit_param.in ./orbit_param.in"
         CALL system(TRIM(version), ierr)
         IF (ierr .lt. 127 .and. ierr .ne. 0) THEN
         PRINT *,' ORBIT control file orbit_param.in not linked'
         lorbit_opt = .false.
         ENDIF

         IF (lorbit_opt) THEN
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
         ENDIF
      ENDIF

      END SUBROUTINE chisq_orbit
