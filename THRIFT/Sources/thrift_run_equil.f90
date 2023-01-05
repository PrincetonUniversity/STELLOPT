!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_equil
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine calls the equilibrium code.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_run_equil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_equil
      USE thrift_vars
      USE read_wout_mod, ONLY: read_wout_deallocate, read_wout_file, &
                               betatot, Aminor, Rmajor,phi, Volume
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      IF (lvmec) THEN
         CALL thrift_paraexe('paravmec_run',proc_string,lscreen_subcodes)
         ! Read the VMEC output
         CALL read_wout_deallocate; ier = 0
         CALL read_wout_file(TRIM(proc_string),ier)
         eq_beta      = betatot
         eq_Aminor    = Aminor
         eq_Rmajor    = Rmajor
         eq_phiedge   = MAXVAL(phi)
         eq_volume    = Volume
         ! Load equil_utils for helpers later on
         CALL thrift_load_vmec
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_equil

