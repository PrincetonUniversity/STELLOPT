!-----------------------------------------------------------------------
!     Module:        thrift_input_mod
!     Authors:       C. van Ham
!     Date:          11/XX/2022
!     Description:   This module contains the THRIFT input namelist and
!                    subroutine which initializes and reads the
!                    THRIFT input namelist.
!-----------------------------------------------------------------------
      MODULE thrift_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE thrift_vars
      USE thrift_runtime
      USE safe_open_mod

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Namelists
!         &thrift_input
!            nr             Number of radial gridpoints
!-----------------------------------------------------------------------
      NAMELIST /thrift_input/ nparallel_runs,bootstrap_type,mboz,nboz, &
                              nrho, nssize, tstart, tend, ntimesteps, jtol, &
                              picard_factor, npicard, lverbj, &
                              eccd_type, &
                              vessel_ecrh, mirror_ecrh, wmode_ecrh, &
                              targettype_ecrh, antennatype_ecrh, &
                              antennaposition_ecrh, &
                              targetposition_ecrh, rbeam_ecrh, &
                              rfocus_ecrh, nra_ecrh, nphi_ecrh, &
                              freq_ecrh, power_ecrh
      
!-----------------------------------------------------------------------
!     Subroutines
!         init_thrift_input:   Initializes the namelist
!         read_thrift_input:   Reads the namelist
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE init_thrift_input
      IMPLICIT NONE
      bootstrap_type     = 'bootsj'
      eccd_type          = ''
      nparallel_runs     = 1
      mboz               = 32
      nboz               = 16
      nrho               = 16
      nssize             = 16
      ntimesteps         = 32
      npicard            = 5
      tstart             = 0.1
      tend               = 1.0
      jtol               = 0.01
      picard_factor      = 0.5
      leccd              = .FALSE.
      lnbcd              = .FALSE.
      lohmic             = .FALSE.
      lverbj             = .FALSE.
      ! TRAVIS vars
      vessel_ecrh     = ''
      mirror_ecrh     = ''
      freq_ecrh       = -1
      wmode_ecrh      = -1
      targettype_ecrh  = 'cyl'
      antennatype_ecrh = 'cyl'
      antennaposition_ecrh = 0
      targetposition_ecrh = 0
      rbeam_ecrh = 0
      rfocus_ecrh = 0
      nra_ecrh = 0
      nphi_ecrh = 8
      RETURN
      END SUBROUTINE init_thrift_input
      
      SUBROUTINE read_thrift_input(filename, istat)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL :: lexist
      INTEGER :: iunit
      CHARACTER(LEN=1000) :: line

      ! Read namelist
      IF (filename /= 'IMAS') THEN
         istat=0
         iunit=12
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) stop 'Could not find input file'
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'thrift_input in: '//TRIM(filename),istat)
         READ(iunit,NML=thrift_input,IOSTAT=istat)
         IF (istat /= 0) THEN
            backspace(iunit)
            read(iunit,fmt='(A)') line
            write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
            CALL handle_err(NAMELIST_READ_ERR,'thrift_input in: '//TRIM(filename),istat)
         END IF
         CLOSE(iunit)
      END IF
      CALL tolower(bootstrap_type)
      CALL tolower(eccd_type)
      leccd = eccd_type .ne. ''
      RETURN
      END SUBROUTINE read_thrift_input

      SUBROUTINE write_thrift_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      INTEGER, INTENT(out) :: istat
      INTEGER :: ik, n
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: vecvar2  = "(2X,A,'(',I3.3,',',I3.3,')',1X,'=',1X,ES22.12E3)"
      istat = 0
      WRITE(iunit_out,'(A)') '&THRIFT_INPUT'
      WRITE(iunit_out,'(A)') '!---------- GENERAL PARAMETERS ------------'
      WRITE(iunit_out,outint) 'NPARALLEL_RUNS',nparallel_runs
      WRITE(iunit_out,outstr) 'BOOTSTRAP_TYPE',bootstrap_type
      WRITE(iunit_out,outflt) 'JTOL',jtol
      WRITE(iunit_out,outint) 'NPICARD',npicard
      WRITE(iunit_out,outflt) 'PICARD_FACTOR',picard_factor
      WRITE(iunit_out,'(A)') '!---------- GRID PARAMETERS ------------'
      WRITE(iunit_out,outint) 'NRHO',nrho
      WRITE(iunit_out,outint) 'NTIMESTEPS',ntimesteps
      WRITE(iunit_out,outflt) 'TEND',tend
      WRITE(iunit_out,'(A)') '!---------- BOOZER TRANSFORMATION ------------'
      WRITE(iunit_out,outint) 'MBOZ',mboz
      WRITE(iunit_out,outint) 'MBOZ',nboz
      WRITE(iunit_out,'(A)') '/'
      RETURN
      END SUBROUTINE write_thrift_namelist

      END MODULE THRIFT_INPUT_MOD