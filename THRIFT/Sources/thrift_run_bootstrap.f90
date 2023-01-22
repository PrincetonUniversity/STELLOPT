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
      INTEGER :: i, ier
      REAL(rprec) :: pprime, rho, epsilon, phip
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ier = 0

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
            ! dp/dPsi = dp/drho * drho/dPhi
            !         = dp/drho / (dPhi/drho)
            DO i = 1, nrho
               rho = THRIFT_RHO(i)
               epsilon = rho*eq_Aminor/eq_Rmajor ! Inverse Aspect Ratio
               CALL get_prof_pprime(rho,THRIFT_T(mytimestep),pprime) ! dpdrho
               CALL EZspline_interp(phip_spl,rho,phip,ier) ! dPsi/drho
               pprime = pprime / phip
               THRIFT_JBOOT(i,mytimestep) = SQRT(epsilon) * pprime * eq_Rmajor
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

