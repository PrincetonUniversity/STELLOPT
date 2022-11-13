!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_bootstrap
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine calls the relevant bootstrap
!                    current code.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_run_bootstrap
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_profiles_mod
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         THRIFT_JBOOT(:,mytimestep) = 0
         RETURN
      END IF

      SELECT CASE(TRIM(bootstrap_type))
         CASE ('bootsj')
            CALL thrift_paraexe('booz_xform',proc_string,lscreen_subcodes)
            CALL thrift_paraexe('bootsj',proc_string,lscreen_subcodes)
         CASE ('sfincs')
         CASE ('test')
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_bootstrap

