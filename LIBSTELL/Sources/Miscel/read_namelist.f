      SUBROUTINE read_namelist(iunit, io_stat, lc_name)
      USE vmec_input, ONLY: read_indata_namelist,
     1   read_mse_namelist
      USE vmec_seq, ONLY: vseq
      USE bootsj_input, ONLY: read_boot_namelist
      USE optim_params, ONLY: read_optimum_namelist, lprof_opt,
     2    lcurprof_opt, numjstar, nsd
      USE coilsnamin, ONLY: read_coils_namelist
      USE gade_mod, ONLY: read_gade_namelist
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, io_stat
      CHARACTER(LEN=*) :: lc_name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ifind
      CHARACTER(LEN=1), PARAMETER :: lead = '&'
      CHARACTER(LEN=LEN_TRIM(lc_name)+1) :: namelc
!-----------------------------------------------

      io_stat = -1
      REWIND (iunit)
      namelc = lead // TRIM(ADJUSTL(lc_name))

      ifind = MIN(len_trim(namelc), 132)
      IF (namelc(1:ifind) .eq. '&indata') THEN
         CALL read_indata_namelist (iunit, io_stat)

      ELSE IF (namelc(1:ifind) .eq. '&optimum') THEN
         CALL read_optimum_namelist (iunit, io_stat)
         IF (io_stat .gt. 0) RETURN
!
!        obsolete assignments
!
         IF (lcurprof_opt) lprof_opt = .true.

      ELSE IF (namelc(1:ifind) .eq. '&bootin') THEN
         CALL read_boot_namelist (iunit, io_stat)
      ELSE IF (namelc(1:ifind) .eq. '&mseprofile') THEN
         CALL read_mse_namelist (iunit, io_stat)
      ELSE IF (namelc(1:ifind) .eq. '&vseq') THEN
         READ (iunit, nml=vseq, iostat=io_stat)
      ELSE IF (namelc(1:ifind) .eq. '&coilsin') THEN
         CALL read_coils_namelist (iunit, io_stat)
      ELSE IF (namelc(1:ifind) .eq. '&ga_de') THEN
         CALL read_gade_namelist (iunit, io_stat)
      END IF

!DEC$ IF DEFINED (IRIX64)
      IF (io_stat .eq. -1) io_stat = 0                  !Catches EOF on ORIG2000 Machine
!DEC$ ENDIF
      END SUBROUTINE read_namelist
