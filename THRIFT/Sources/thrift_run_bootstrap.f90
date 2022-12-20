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
      USE thrift_equil
      USE thrift_profiles_mod
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i
      REAL(rprec) :: pprime, rho, sflux
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         THRIFT_JBOOT(:,mytimestep) = 0
         RETURN
      END IF

      SELECT CASE(TRIM(bootstrap_type))
         CASE ('model','simple','test')
            ! j_BS = sqrt(epsilon) Rmajor *dp/dPsi
            ! epsilon = r/R (inverse aspect ratio)
            ! dp/dPsi : Pa/Wb (toroidal flux derivative)
            DO i = 1, nrho
               rho = THRIFT_RHO(i)
               sflux = rho*rho
               CALL get_prof_pprime(rho,THRIFT_T(mytimestep),pprime) ! dpdrho
               pprime = pprime * 0.5 / rho ! dpds = dpdrho * drho/ds = dpdrho / (ds/drho) = dpdrho / (2*rho)
               THRIFT_JBOOT(i,mytimestep) = SQRT(eq_Aminor*eq_Rmajor*rho) * pprime / eq_phiedge
            END DO
         CASE ('bootsj')
            CALL thrift_paraexe('booz_xform',proc_string,lscreen_subcodes)
            CALL thrift_paraexe('bootsj',proc_string,lscreen_subcodes)
         CASE ('sfincs')
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_bootstrap

