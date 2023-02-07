!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_diagnostics
!     Authors:       S. Lazerson
!     Date:          11/24/2022
!     Description:   This subroutine calls the relevant diagnostic codes.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_run_diagnostics
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_equil
      USE thrift_vars
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) RETURN

      ! Check and run the diagno code
      IF (ldiagno) CALL thrift_paraexe('diagno',proc_string,lscreen_subcodes)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_diagnostics

