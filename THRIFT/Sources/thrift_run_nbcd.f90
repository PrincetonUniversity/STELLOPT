!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_eccd
!     Authors:       S. Lazerson
!     Date:          02/04/2023
!     Description:   This subroutine calls the relevant NBCD code.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_run_nbcd
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_equil
      USE thrift_vars
      USE thrift_profiles_mod
!-----------------------------------------------------------------------
!     Local Variables
!         i           Loop index
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), PARAMETER :: amu = 1.66053906660D-27
      INTEGER :: i
      REAL(rprec) :: mass, NBI_A, NBI_E, NBI_Z, NBI_alpha, NBI_theta, NBI_Pdense ! probably should be inputs
      REAL(rprec) :: ne, Te, Ti, Zeff, InvAspect, G, NBI_V, rho, timenow, jopb
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         THRIFT_JNBCD(:,mytimestep) = 0
         RETURN
      END IF

      SELECT CASE(TRIM(eccd_type))
         CASE ('model','okano')
            timenow = THRIFT_T(mytimestep)
            ! Use the Okano model for NBCD
            mass = amu
            NBI_A = 1
            NBI_E = 55D3 * 1.60217663D-19
            NBI_Z = 1
            NBI_alpha = 0.3
            NBI_theta = 1.22173
            DO i = 1, nrho
               rho = THRIFT_RHO(i)
               !! This part should be replaced by namelist var
               NBI_Pdense = 1D6*(1-rho**0.25)
               !!!!!!!!!!!!!!!!!
               NBI_V = SQRT(2*NBI_E/(NBI_A*amu))
               InvAspect = rho*eq_Aminor/eq_Rmajor
               CALL get_prof_ne(rho,   timenow, ne)
               CALL get_prof_te(rho,   timenow, Te)
               CALL get_prof_te(rho,   timenow, Ti)
               CALL get_prof_zeffne(rho, timenow, Zeff)
               CALL vmec_ohkawa(rho*rho,zeff,G) ! (1-l31)/Zeff
               G = 1.0D0-NBI_Z*G
               CALL nbcd_okano(ne, Te, Ti, mass, Zeff, InvAspect, G, &
                               NBI_A, NBI_V, NBI_Z, NBI_alpha, NBI_theta, &
                               jopb)
               THRIFT_JNBCD(i,mytimestep) = jopb * NBI_Pdense
            END DO
         CASE ('beams3d')
            !CALL thrift_paraexe('beams3d',proc_string,lscreen_subcodes)
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_nbcd

