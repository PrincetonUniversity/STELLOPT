!-----------------------------------------------------------------------
!     Module:        stelltran_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/21/2015
!     Description:   This module contains the STELLTRAN input namelist
!                    and subroutine which initializes and reads the
!                    STELLTRAN input namelist.
!-----------------------------------------------------------------------
      MODULE stelltran_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_vars
      USE stelltran_data
!       USE stelltran_rho_vars
      USE stelltran_equilutils
      USE safe_open_mod, ONLY: safe_open
      USE vmec_params, ONLY: version_vmec=> version_
      USE mpi_params
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF        
!-----------------------------------------------------------------------
!     Input Namelists
!         &optimum
!            equil_type         Equilibrium Code:
!                                  'VMEC2000' (default)
!                                  'SPEC'
!            db_file            Path to timeslice database file
!
!   NOTE:  All varaibles which start with target have an similar
!          varaible starting with sigma which defines the error bars.
!-----------------------------------------------------------------------
      REAL(rprec) :: rad_param
      NAMELIST /stelltran_input/ equil_type, db_file, fit_tol, mboz,&
                                 nboz, run_type, dt, ntimesteps, nprof,&
                                 vessel_ecrh, mirror_ecrh, antennatype_ecrh,&
                                 targettype_ecrh, antennaposition_ecrh,&
                                 targetposition_ecrh, rbeam_ecrh, rfocus_ecrh,&
                                 nra_ecrh, nphi_ecrh, wmode_ecrh, freq_ecrh,&
                                 rad_param
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_stellopt_input:   Reads optimum namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_stelltran_input(filename, istat, ithread)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      INTEGER, INTENT(in) :: ithread
      LOGICAL :: lexist
      INTEGER :: i, iunit, local_master
      ! Initializations
      equil_type      = 'VMEC2000'
      db_file         = 'empty'
      run_type        = 'analysis'
      fit_tol         = 0.01
      ntimesteps      = 100 ! Only needed for predictive case.
      nprof           = 50
      vessel_ecrh     = ''
      mirror_ecrh     = ''
      freq_ecrh       = -1
      wmode_ecrh      = ''
      targettype_ecrh  = 'cyl'
      antennatype_ecrh = 'cyl'
      antennaposition_ecrh = 0
      targetposition_ecrh = 0
      rbeam_ecrh = 0
      rfocus_ecrh = 0
      nra_ecrh = 0
      nphi_ecrh = 8
      ! Read name list
      lexist            = .false.
      istat=0
      iunit=12
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) CALL handle_err(FILE_EXIST_ERR,TRIM(filename),istat)
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat /= 0) CALL handle_err(FILE_OPEN_ERR,TRIM(filename),istat)
      READ(iunit,NML=stelltran_input,IOSTAT=istat)
      IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'STELLTRAN_INPUT in: '//TRIM(filename),istat)
      CALL FLUSH(iunit)
      CLOSE(iunit)
      ! Fix String vars
      equil_type=TRIM(equil_type)
      equil_type=ADJUSTL(equil_type)
      db_file=TRIM(db_file)
      db_file=ADJUSTL(db_file)
      run_type=TRIM(run_type)
      run_type=ADJUSTL(run_type)
      ! Handle Some things
      CALL tolower(equil_type)
      IF ((myid == master) .and. (TRIM(equil_type(1:4)) == 'vmec') ) THEN
         WRITE(6,*)        " Equilibrium calculation provided by: "
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,"(2X,A)") "=========       Variational Moments Equilibrium Code (v "//TRIM(version_vmec)//")           ========="
         WRITE(6,"(2X,A)") "=========                (S. Hirshman, J. Whitson)                      ========="
         WRITE(6,"(2X,A)") "=========         http://vmecwiki.pppl.wikispaces.net/VMEC              ========="
         WRITE(6,"(2X,A)") "================================================================================="
         WRITE(6,*)        "    "
      END IF
      END SUBROUTINE read_stelltran_input
      
      SUBROUTINE write_stelltran_namelist(iunit,istat)
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(in) :: istat
      INTEGER :: ik, n, m, u, v
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,L1,2(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3))"
      WRITE(iunit,'(A)') '&STELLTRAN_INPUT'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(A)') '!       Transport Run Control Parameters'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,outstr) 'EQUIL_TYPE',TRIM(equil_type)
      WRITE(iunit,outstr) 'EQUIL_TYPE',TRIM(db_file)
      WRITE(iunit,outstr) 'RUN_TYPE',TRIM(run_type)
      WRITE(iunit,outflt) 'FIT_TOL',fit_tol
      WRITE(iunit,outint) 'NTIMESTEPS',ntimesteps
      WRITE(iunit,outint) 'NPROF',nprof
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(A)') '!          ECE Reflectometry OPTIMIZATION'
      WRITE(iunit,'(A)') '!----------------------------------------------------------------------'
      WRITE(iunit,'(2X,A,I3.3)') 'NRA_ECRH = ',nra_ecrh
      WRITE(iunit,'(2X,A,I3.3)') 'NPHI_ECRH = ',nphi_ecrh
      IF (LEN_TRIM(vessel_ecrh) > 1) WRITE(iunit,outstr) 'VESSEL_ECRH',TRIM(vessel_ecrh)
      IF (LEN_TRIM(mirror_ecrh) > 1) WRITE(iunit,outstr) 'MIRROR_ECRH',TRIM(mirror_ecrh)
      IF (LEN_TRIM(targettype_ecrh) > 1) WRITE(iunit,outstr) 'TARGETTYPE_ECRH',TRIM(targettype_ecrh)
      IF (LEN_TRIM(antennatype_ecrh) > 1) WRITE(iunit,outstr) 'ANTENNATYPE_ECRH',TRIM(antennatype_ecrh)
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_stelltran_namelist
      
      END MODULE stelltran_input_mod
