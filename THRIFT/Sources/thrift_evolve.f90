!-----------------------------------------------------------------------
!     Subroutine:    thrift_evolve
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine evolves the current profile
!-----------------------------------------------------------------------
      SUBROUTINE thrift_evolve
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: deltaj, jold
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Initialize the current density
      THRIFT_J       = 0
      THRIFT_JBOOT   = 0
      THRIFT_JPLASMA = 0
      THRIFT_JECCD   = 0
      THRIFT_JNBCD   = 0
      THRIFT_JOHMIC  = 0

      ! Allocate the convergence helper
      ALLOCATE(deltaj(nrho), jold(nrho))
      jold   = 1E3 ! so on loop 1 we don't divide by zero
      
      ! Loop over timesteps
      DO mytimestep = 1, ntimesteps
      
         ! Setup the profiles
         CALL thrift_equil_p

         ! Converge equilibrium currents
         deltaj = 10*jtol
         DO WHILE (ANY(deltaj > jtol))
            ! Update equilbrium current
            CALL thrift_equil_j
            ! Run equilibrium
            CALL thrift_run_equil
            ! Calculate Bootstrap
            CALL thrift_run_bootstrap
            ! Calculate Current Drive
            IF (leccd)  CALL thrift_run_ECCD
            IF (lnbcd)  CALL thrift_run_NBI
            IF (lohmic) CALL thrift_run_ohmic
            ! Calc the inductive chagne in current
            CALL thrift_jinductive
            ! Update total current
            THRIFT_J(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep) &
                                   + THRIFT_JBOOT(:,mytimestep) &
                                   + THRIFT_JECCD(:,mytimestep) &
                                   + THRIFT_JNBCD(:,mytimestep) &
                                   + THRIFT_JOHMIC(:,mytimestep)
            ! Check the convergence
            deltaj = ABS( THRIFT_J(:,mytimestep) - jold) / ABS(jold)
         END DO
         jold = THRIFT_J(:,mytimestep)
      END DO

      ! Deallocate helpers
      DEALLOCATE(deltaj,jold)
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_evolve

