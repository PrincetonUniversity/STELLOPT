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
      USE pies_input_mod
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, ierr
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      OPEN(iunit,FILE='input.' // TRIM(id_string))
      CALL read_pies_input(iunit,ierr)
      
      write(6,*) '-----PIES File Parameters-----'
      write(6,*) '  extsurfs: ',extsurfs
      write(6,*) '    lideal: ',lideal
      IF (lrestart) THEN
      !   CALL pies_init_restart
      ELSE IF (lwout) THEN
         CALL pies_init_wout
      ELSE
         CALL pies_init_input
      END IF
      CLOSE(iunit)
      !CALL write_pies_netcdf
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE pies_init
