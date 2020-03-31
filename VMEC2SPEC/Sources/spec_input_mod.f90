!-----------------------------------------------------------------------
!     Module:        spec_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/09/2012
!     Description:   This module contains the vmec2spec input namelist and
!                    subroutine which initializes and reads the
!                    vmec2spec input namelist.
!-----------------------------------------------------------------------
      MODULE spec_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE spec_runtime
      USE spec_background
      USE safe_open_mod, ONLY: safe_open
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Namelists
!         fieldlines_input
!-----------------------------------------------------------------------
      NAMELIST /vmec2spec_input/ nvol, ni, tflux
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_fieldlines_input:   Reads fieldlines_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_vmec2spec_input(filename, istat)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL :: lexist
      INTEGER :: i, iunit
      ! Initializations
      nvol   = 8
      ni     = 4
      tflux  = 0
      pflux  = 0
      pl     = 0
      ql     = 0
      pr     = 0
      qr     = 0
      ! Read namelist
      istat=0
      iunit=12
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) stop 'Could not find input file'
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'vmec2spec_input in: input.'//TRIM(id_string),istat)
      READ(iunit,NML=vmec2spec_input,IOSTAT=istat)
      IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'vmec2spec_input in: input.'//TRIM(id_string),istat)
      CLOSE(iunit)
      END SUBROUTINE read_vmec2spec_input

      END MODULE spec_input_mod
