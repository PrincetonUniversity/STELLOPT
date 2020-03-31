!-----------------------------------------------------------------------
!     Subroutine:    pies_init
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine initializes a run of the PIES code.
!                    The pies code can be run in three modes:
!                    1) From a VMEC input file.  (default)
!                    2) From a VMEC wout file.   (-wout)
!                    3) From a PIES netCDF file. (-lrestart)
!                    The first argument to the executable is always the
!                    name of a file.  The optional flags control how the
!                    code is executed.
!-----------------------------------------------------------------------
      SUBROUTINE pies_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_runtime
      USE pies_background, ONLY: k
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error Flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ierr
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,*) '-----PIES File Parameters-----'
         WRITE(6,*)                           '  extsurfs: ',extsurfs
         IF (k > 0) write(6,*)                '         k: ',k
         IF (m_new > 0) write(6,*)            '         m: ',m_new
         IF (n_new > 0) write(6,*)            '         n: ',n_new
         IF (free_override .eq. 0) write(6,*) '  Fixed Boundary Enforced!'
         IF (free_override .eq. 1) write(6,*) '  Free  Boundary Enforced!'
         IF (lmake_coils) WRITE(6,*)          '  Create coil_data file!' 
      END IF      
      IF (INDEX(TRIM(id_string(1:4)),'wout') .gt. 0) THEN
         lwout = .true.
         CALL pies_init_wout
      ELSE IF (INDEX(TRIM(id_string(1:5)),'input') .gt. 0) THEN
         lwout = .false.
         CALL pies_init_input
      ELSE
         CALL handle_error(BAD_INPUT_ERR,TRIM(id_string),-1)
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE pies_init
