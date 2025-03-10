!-----------------------------------------------------------------------
!     Subroutine:    spec_init
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/09/2012
!     Description:   This subroutine initializes a run of the SPEC code.
!                    The SPEC code can be run in three modes:
!                    1) From a VMEC input file.  (default)
!                    2) From a VMEC wout file.   (-wout)
!                    The first argument to the executable is always the
!                    name of a file.  The optional flags control how the
!                    code is executed.
!-----------------------------------------------------------------------
      SUBROUTINE spec_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE spec_runtime
      USE spec_background, ONLY: nvol
      USE spec_input_mod
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error Flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ierr
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      CALL read_vmec2spec_input(TRIM(id_string),ierr)
      IF (lverb) THEN
         WRITE(6,*) '-----SPEC File Parameters-----'
         write(6,*)                '      nvol: ',nvol
         !IF (m_new > 0) write(6,*)            '         m: ',m_new
         !IF (n_new > 0) write(6,*)            '         n: ',n_new
      END IF      
      IF (lwout) THEN
         lwout = .true.
         CALL spec_init_wout
      ELSE IF (INDEX(TRIM(id_string(1:5)),'input') .gt. 0) THEN
         lwout = .false.
         CALL spec_init_input
      ELSE
         CALL handle_error(BAD_INPUT_ERR,TRIM(id_string),-1)
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE spec_init
