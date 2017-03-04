      SUBROUTINE chisq_nescoil (number, nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, bigno, exe_suffix
      USE optim_params, ONLY: lnescoil_opt, coil_separation,
     1     target_coil_jmax, sigma_coil_complex, sigma_coil_jmax,
     2     sigma_berr_ave, target_coil_complex
      USE system_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: number, nopt, iflag
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nescoil0 = 28
      REAL(rprec), PARAMETER :: zero = 0, one = 1
      INTEGER, PARAMETER :: ivals = 3
      CHARACTER(LEN=10), PARAMETER :: search_string(ivals) =
     1   (/ 'Complexity', 'J Surf MAX', 'Berr ave, ' /)
      INTEGER, PARAMETER :: search_index(ivals) =  (/ 14, 24, 22 /)
      REAL(rprec) :: search_value(ivals)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit, istat, istat1, index1, itype
      REAL(rprec) :: coil_complex_opt, jsheet_max, berr_ave, berr_max
      CHARACTER(len=LEN_TRIM(home_dir) + LEN_TRIM(extension)+30)
     1    :: version
      CHARACTER(LEN=200) :: line
      LOGICAL :: ex, ex1
C-----------------------------------------------
!
!     COMPUTE COIL COMPLEXITY FROM NESCOIL POTCOEFFS FILE
!
      version = TRIM(home_dir) // 'xbnorm' // TRIM(exe_suffix)

      IF (nopt .gt. 0) THEN
         iunit = nescoil0


!     Write coil_separation for USE by BNORM (NESCOIL)
         version = TRIM(version) // ' wout.' // TRIM(extension)
         WRITE (version,'(a,1p,e12.3)') TRIM(version), coil_separation
         IF (lscreen) THEN
           PRINT '(/,a)',' Running BNORM code...'
           PRINT 100,' Coil separation = ', coil_separation
         END IF
         CALL system(version)                              !!Produce bnorm, nescin files

         version = TRIM(home_dir) // 'xnescoil' // TRIM(exe_suffix)
         IF (lscreen)
     1      PRINT '(/,a)',' Running NESCOIL code ...'
         CALL load_physics_codes (version, 'nescin', ' ',
     1        'nescout', extension, iunit, iflag)
         IF (iflag .ne. 0) RETURN
!
!     Read in current complexity,amaximum current density, and berr from nescout file (iunit)
!
         DO itype = 1, ivals
            index1 = -1
            DO WHILE (index1 .le. 0)
               READ(iunit, '(a)', iostat=istat) line
               IF (istat .ne. 0) STOP ' ISTAT != 0 in chisq_nescoil'
               index1 = index(line, search_string(itype))
               IF (index1 .gt. 0) THEN
                  IF (itype .eq. 3) THEN
                     READ (line(search_index(itype):), *)
     1               search_value(itype), berr_max
                  ELSE
                     READ (line(search_index(itype):), *)
     1               search_value(itype)
                  END IF
               END IF
            END DO
         END DO

         coil_complex_opt = search_value(1)
         jsheet_max       = search_value(2)
         berr_ave         = search_value(3)

         IF (lscreen) THEN
            PRINT 100,' Coil Complexity = ', coil_complex_opt
            PRINT 100,' Maximum Sheet Current Density = ', jsheet_max
            PRINT 100,' Berr-ave = ', berr_ave
         END IF

         CLOSE (iunit)

!        COIL COMPLEXITY (NO PENALTY IF ACTUAL COMPLEXITY IS BELOW TARGET)
         IF (sigma_coil_complex .lt. bigno) THEN
            number = number + 1
            index_array(number) = ivar_coil_complexity
            wegt(number) = sigma_coil_complex
            chisq_target(number) = target_coil_complex
            chisq_match(number)  =
     1         MAX(target_coil_complex, coil_complex_opt)
         END IF

!        MAXIMUM CURRENT DENSITY
         IF (sigma_coil_jmax .lt. bigno) THEN
            number = number + 1
            index_array(number) = ivar_coil_jmax
            wegt(number) = sigma_coil_jmax
            chisq_target(number) = MIN(ABS(jsheet_max),target_coil_jmax)
            chisq_match(number)  = ABS(jsheet_max)
         END IF

!        AVERAGE BERR
         IF (sigma_berr_ave .lt. bigno) THEN
            number = number + 1
            index_array(number) = ivar_berr_ave
            wegt(number) = sigma_berr_ave
            chisq_target(number) = zero
            chisq_match(number)  = ABS(berr_ave)
         END IF

      ELSE
         INQUIRE(file=TRIM(version), exist=ex, iostat=istat)
         INQUIRE(file=TRIM(home_dir) // 'xnescoil' // TRIM(exe_suffix), 
     1           exist=ex1, iostat=istat1)
         IF (istat.ne.0 .or. .not.ex) THEN
            IF (lscreen)
     1          PRINT *, 'xbnorm file not found in ' // TRIM(home_dir)
            lnescoil_opt = .false.
         ELSE IF (istat1.ne.0 .or. .not.ex1) THEN
            IF (lscreen)
     1          PRINT *, 'xnescoil file not found in ' // TRIM(home_dir)
            lnescoil_opt = .false.
         ELSE
            IF (sigma_coil_complex .lt. bigno) THEN
               number = number+1
               IF (nopt .eq. -2) chisq_descript(number) = 
     1                           descript(ivar_coil_complexity)
            ENDIF
            IF (sigma_coil_jmax .lt. bigno) THEN
               number = number+1
               IF (nopt .eq. -2) chisq_descript(number) = 
     1                           descript(ivar_coil_jmax)
            ENDIF
            IF (sigma_berr_ave .lt. bigno) THEN
               number = number+1
               IF (nopt .eq. -2) chisq_descript(number) = 
     1            descript(ivar_berr_ave)
            ENDIF
         ENDIF
      ENDIF

 100  FORMAT(a,1p,e12.3)

      END SUBROUTINE chisq_nescoil
