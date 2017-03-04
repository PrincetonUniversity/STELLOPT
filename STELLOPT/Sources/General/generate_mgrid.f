      SUBROUTINE generate_mgrid (extension, lscreen, iflag)
c-----------------------------------------------
c   M o d u l e s
c-----------------------------------------------
      USE stel_kinds
      USE vmec_input, ONLY: nzeta, extcur
      USE optim_params, ONLY: rgrid_min, rgrid_max, zgrid_min, 
     1    zgrid_max, lcoil_geom, lvac_opt
      USE optim, ONLY: rmax_opt, rmin_opt, zmax_opt, aminor_opt,
     1    home_dir, exe_suffix, remove, move
      USE safe_open_mod
      USE system_mod
      USE coilsnamin, ONLY: lmodular, lsaddle, lbcoil, lvf, lsurfv,
     1    lsadsfv
      IMPLICIT NONE
c-----------------------------------------------
c   D u m m y   A r g u m e n t s
c-----------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: extension
      INTEGER, INTENT(out) :: iflag
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      INTEGER :: ierr, iunit=10, i1, j1, istat
      REAL(rprec) :: rmin, rmax, zmin, zmax, delta
      CHARACTER(LEN=200) :: parse_string, output_file, coils_file,
     1                      line, cmd
      CHARACTER(LEN=4) :: cmd_line
      LOGICAL :: lexist
c-----------------------------------------------
!
!     Generates an mgrid file (for free-boundary calculation) corresponding to PRESENT state of coils
!     Parse the extcur file to read in extcur array
!

!     Step 1: generate coils.extension file 
!
      IF (lcoil_geom) THEN
!     Call xcoilgeom (the reduced xcoilopt executable)
!     It reads the COILSIN NAMELIST in the INPUT.EXTENSION file. Redirect screen output from
!     xcoilgeom to a file, COIL_TARGETS.EXTENSION, for future parsing to get engineering target
!     chi-sq values
         output_file = 'coil_targets.' // TRIM(extension)
         parse_string = TRIM(home_dir) // 'xcoilgeom' // exe_suffix
     1               // TRIM(extension) // ' > ' // TRIM(output_file)
         IF (lscreen) PRINT *, 'Creating MGRID file'
         CALL system (parse_string, ierr)
         IF (ierr.lt.127 .and. ierr.ne.0) THEN
            IF (lscreen) 
     1      PRINT *, 'XCOILGEOM failed in generate_mgrid: ','ierr = ',
     2      ierr
            iflag = -22
            RETURN
         END IF
         coils_file = "coils." // TRIM(extension)

      ELSE IF (lvac_opt) THEN
!     Call biotsavart function write_coils_file(extension) to generate rotated/shifted coils file
!     write out coils file "coils.extension_rot"
         cmd = TRIM(home_dir) // 'xvacopt' // TRIM(exe_suffix)
         cmd_line = ' F T'
         CALL load_physics_codes (cmd, 'input', cmd_line, ' ',
     1        extension, iunit, ierr)
         IF (ierr .ne. 0) THEN
            iflag = -22
            RETURN
         END IF
         coils_file = "coils." // TRIM(extension) // "_rot"
      END IF

!     b. Make sure coils.extension file was successfully created

      INQUIRE (file=coils_file, exist=lexist, iostat=istat)
      IF (istat.ne.0 .or. .not.lexist) THEN
         iflag = -23
         RETURN
      END IF

!     Step 2: generate actual mgrid.extension file from coils.extension
!     a. Estimate grid box dimensions from wout file dimension (previous run of vmec)
!        or rgrid_min,max, zgrid_min,max in namelist
!        Use plasma max, min +- minor radius, and adjust to make delr = delz
      IF (rmax_opt .gt. rmin_opt) THEN
         rmin = MIN(rmin_opt, rgrid_min)
         zmax = MAX(zmax_opt, zgrid_max)
         rmax = MAX(rmax_opt, rgrid_max)
         zmin =-zmax
         delta = ((rmax - rmin) - (zmax - zmin))/2
         IF (delta .gt. zero) THEN
            zmax = zmax + delta
            zmin = zmin - delta
         ELSE IF (delta .lt. zero) THEN
            rmax = rmax - delta
            rmin = rmin + delta
            IF (rmin .lt. EPSILON(aminor_opt)) THEN
               rmax = rmax - rmin
               rmin = EPSILON(aminor_opt)
            END IF
         END IF
         rgrid_min = rmin;  rgrid_max = rmax
         zgrid_min = zmin;  zgrid_max = zmax
      END IF
      IF (rgrid_min.ge.rgrid_max .or. zgrid_min.ge.zgrid_max) THEN
            PRINT *, 'rgrid_min MUST BE < rgrid_max'
            PRINT *, 'zgrid_min MUST BE < zgrid_max'
            STOP
      END IF

!     b. Parse command line for xgrid file [produces mgrid_extension.nc (or .bin)]

      WRITE (parse_string, '(3a,2(1x,a1),1x,4(1p,e20.12,1x),i4)')
     1   TRIM(home_dir),'xgrid ', TRIM(coils_file(7:)), 'S','Y',
     2   rgrid_min, rgrid_max, zgrid_min, zgrid_max, nzeta

!DEC$ IF .NOT.DEFINED (WIN32)
      IF (.not.lscreen) parse_string = TRIM(parse_string) //
     1   ' > /dev/null'
!DEC$ ENDIF


!     c. Make system call to run xgrid, which generates mgrid file
      CALL system (parse_string, ierr)

!     d. Clean up by removing coils file
      CALL system (remove // TRIM(coils_file))

!     e. Parse extcur file if NOT in lvac_opt=T mode
      IF (.not.lvac_opt) THEN
         parse_string = "extcur." // TRIM(coils_file(7:))
         CALL safe_open(iunit, ierr, parse_string, 'old', 'formatted')
         IF (ierr .ne. 0) THEN
            PRINT *,
     1      ' Error opening extcur file in STELLOPT generate_mgrid: ',
     2      ' ierr = ', ierr
            iflag = -24
            RETURN
         END IF

         extcur = 0

         DO WHILE (ierr .eq. 0)
            READ (iunit,'(a)', iostat=ierr) line
            IF (ierr .ne. 0) EXIT
            i1 = INDEX(line, '(')
            IF (i1 .eq. 0) EXIT
            j1 = INDEX(line,')') - 1
            READ (line(i1+1:j1), *) j1
            IF (j1.le.0 .or. j1.gt.SIZE(extcur)) CYCLE
            i1 = INDEX(line,'=')
            IF (i1 .eq. 0) EXIT
            READ (line(i1+1:), *) extcur(j1)
         END DO

         CLOSE (iunit)

      ELSE
         coils_file = "mgrid_" // TRIM(extension)
         INQUIRE (file=TRIM(coils_file) // "_rot", exist=lexist, 
     1            iostat=istat)
         IF (lexist) THEN
            cmd =move // TRIM(coils_file) // "_rot " // TRIM(coils_file)
         ELSE
            INQUIRE (file=TRIM(coils_file) // "_rot.nc", exist=lexist, 
     1               iostat=istat)
            IF (lexist) THEN
               cmd =move // TRIM(coils_file) // "_rot.nc " // 
     1                      TRIM(coils_file) // ".nc"
            ELSE
               INQUIRE (file=TRIM(coils_file) // "_rot.bin", 
     1                  exist=lexist, iostat=istat)
               IF (lexist) THEN
                  cmd =move // TRIM(coils_file) // "_rot.bin " // 
     1                         TRIM(coils_file) // ".bin"
               END IF
            END IF
         END IF 
         CALL system(cmd, ierr)
         IF (ierr .ne. 0) THEN
            PRINT *,' CMD: ', TRIM(cmd)
            STOP 'CMD FAILURE IN GENERATE_MGRID'
         END IF
         cmd = remove // "extcur." // TRIM(coils_file(7:)) // "_rot"
         CALL system(cmd, ierr)
      END IF

      iflag = 0

      END SUBROUTINE generate_mgrid
