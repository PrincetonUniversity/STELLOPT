      SUBROUTINE chisq_dkes(sigma, num, nopt, iflag, nsval,
     1   extension, command)
      USE stel_kinds
      USE chisq_mod
      USE optim_params, ONLY: ldkes_opt, ndkes_mask, dkes_nu,
     1                        dkes_efield, nsurf_mask
      USE optim, ONLY: home_dir, iunit_dkes, exe_suffix, remove
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
c---------------------------------------------------------
c    D u m m y  A r g u m e n t s
c---------------------------------------------------------
      INTEGER, INTENT(in) :: nopt, nsval
      INTEGER, INTENT(inout) :: num, iflag
      REAL(rprec), INTENT(in) :: sigma
      CHARACTER(LEN=*), INTENT(in) :: extension, command
c
c---------------------------------------------------------
c    L o c a l s  P a r a m e t e r s
c---------------------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      CHARACTER(LEN=*), PARAMETER :: input_file = 'boozer',
     1   output_file = 'opt_dkes'
c---------------------------------------------------------
c    L o c a l  V a r i a b l e s
c---------------------------------------------------------
      CHARACTER(LEN=30) :: dkes_args
      CHARACTER(len=LEN_TRIM(home_dir)+20) :: version
      CHARACTER(len=LEN_TRIM(command)+ 40) :: comdline
      INTEGER :: iunit, ierror
      REAL (rprec) ::  L11p, L33p, L31p
      REAL (rprec) ::  L11m, L33m, L31m
      REAL (rprec) ::  scal11, scal33, scal13
      LOGICAL, PARAMETER :: ldebug = .false.
c---------------------------------------------------------

      IF (ldebug .and. myid.eq.master) THEN

        PRINT *, "********** Beginning chisq_dkes **************"
        PRINT *, "nopt = ", nopt
        PRINT *, "num = ", num
        PRINT *, "ldkes_opt = ", ldkes_opt
        PRINT *, "surf_dkes = ", ndkes_mask(nsval)
        PRINT *, "dkes_nu = ", dkes_nu(nsval)
        PRINT *, "dkes_efield = ", dkes_efield(nsval)
        PRINT *, "sigma_dkes = ", sigma
        PRINT *, "extension = ", TRIM(extension)
        PRINT *, "command = ", TRIM(command)
        PRINT *, "unit_dkes = ", unit_dkes

      END IF

      iflag = 0

      IF (nopt .gt. 0) THEN

         iunit = unit_dkes

!        VERSION: EXECUTABLE FOR XDKES
         version = TRIM(home_dir) // 'xdkes' // TRIM(exe_suffix)

!        USE INTERNAL FILE TO WRITE DKES NAMELIST INPUT INTO A CHARACTER VARIABLE
!        NAMELIST INPUT = SURFACE NUMBER, COLLISIONALITY, EFIELD
         WRITE(dkes_args,'(1x,i3,1x,f10.5,1x,f10.5)')
     1        nsval, dkes_nu(nsval), dkes_efield(nsval)

!        ADD DKES NAMELIST INPUT TO COMMAND LINE
!        ADD LOGICAL SCREEN OUTPUT VARIABLE TO COMMAND LINE
!        CALL TO DKES WILL AUTOMATICALLY ACTIVATE THE INTERNAL DKES_INPUT_PREPARE CALL
         comdline = TRIM(dkes_args) // ' ' // TRIM(command)

         CALL load_physics_codes (version, input_file, comdline,
     1         output_file, extension, iunit, iflag)
         IF (iflag .ne. 0) RETURN

!        READ OPT_DKES FILE AND STORE DATA IN DKES_OPT FILE FOR ALL SURFACES
         READ(iunit,'(3(2x,e24.13))', iostat=ierror) L11p, L33p, L31p
         READ(iunit,'(3(2x,e24.13))', iostat=ierror) L11m, L33m, L31m
         READ(iunit,'(3(2x,e24.13))', iostat=ierror)
     1        scal11, scal33, scal13

         WRITE(iunit_dkes,'(/,a,i5)') 'RADIAL SURFACE JS =', nsval
         WRITE(iunit_dkes,'(3(2x,e24.13))') L11p, L33p, L31p
         WRITE(iunit_dkes,'(3(2x,e24.13))') L11m, L33m, L31m
         WRITE(iunit_dkes,'(3(2x,e24.13))') scal11, scal33, scal13

         CLOSE (iunit, status='delete')
         CALL system(remove // "dkesout." // TRIM(extension))

         IF (ldebug .and. myid.eq.master) THEN
            PRINT *, "OPT_DKES FILE"
            PRINT *,  L11p, L33p, L31p
            PRINT *,  L11m, L33m, L31m
            PRINT *,  scal11, scal33, scal13
         ENDIF
         IF (ierror .ne. 0) THEN
            iflag = -20
            RETURN
         END IF

!        LOAD CHISQ ARRAYS: NOTE THERE IS NO NORMALIZATION FOR L11
!
         num = num + 1
         index_array(num) = ivar_dkes
         wegt(num) = sigma
         chisq_target(num) = zero
         chisq_match(num) = 0.5_dp*(L11p + L11m)

      ELSE
         num = num + 1
         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_dkes)
      END IF

      IF (ldebug .and. myid.eq.master) THEN
         PRINT *, "at END of chisq_dkes"
         PRINT *, "nopt = ", nopt
         PRINT *, "num = ", num
         PRINT *, "********** Ending chisq_dkes **************"
      END IF

      END SUBROUTINE chisq_dkes
