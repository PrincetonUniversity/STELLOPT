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
      LOGICAL :: lfirst_pass
      INTEGER :: nsubsteps
      REAL(rprec) :: alpha
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: deltaj, jold
      CHARACTER(len = 16)     :: temp1_str, temp2_str
      
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
      alpha = 0.05
      jold   = 1E3 ! so on loop 1 we don't divide by zero
      lscreen_subcodes = .TRUE.
      lfirst_pass = .TRUE.
      
      ! Loop over timesteps
      DO mytimestep = 2, ntimesteps
      
         ! Setup the profiles
         CALL thrift_equil_p

         ! Converge equilibrium currents
         deltaj = 10*jtol; nsubsteps = 0
         DO WHILE (ANY(deltaj > jtol))

            ! Update Substeps
            nsubsteps = nsubsteps + 1

            ! Get out if too many substeps
            IF (nsubsteps > 5) EXIT

            ! We need minor_radius to calc Itor
            IF (lfirst_pass) minor_radius = 1.0

            ! Create filename
            WRITE(temp1_str,'(i5)') mytimestep
            WRITE(temp2_str,'(i3)') nsubsteps
            proc_string = TRIM(TRIM(id_string) // '.' //  &
                  TRIM(ADJUSTL(temp1_str)) // '_' // &
                  TRIM(ADJUSTL(temp2_str)))

            ! Update equilbrium current
            CALL thrift_equil_j

            ! Run equilibrium
            CALL thrift_run_equil

            ! Calculate Bootstrap
            CALL thrift_run_bootstrap

            ! Calculate Current Drive
            !IF (leccd)  CALL thrift_run_ECCD
            !IF (lnbcd)  CALL thrift_run_NBI
            !IF (lohmic) CALL thrift_run_ohmic

            ! Calc the inductive chagne in current
            !CALL thrift_jinductive

            ! Update total current
            THRIFT_J(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep) &
                                   + THRIFT_JBOOT(:,mytimestep) &
                                   + THRIFT_JECCD(:,mytimestep) &
                                   + THRIFT_JNBCD(:,mytimestep) &
                                   + THRIFT_JOHMIC(:,mytimestep)

            ! Check the convergence
            deltaj = 0
            IF (lfirst_pass) THEN
               WHERE(ABS(THRIFT_J(:,mytimestep))>0) deltaj = ABS(THRIFT_J(:,mytimestep))
               jold = alpha*THRIFT_J(:,mytimestep)
            ELSE
               WHERE(ABS(jold)>0) deltaj = ABS( THRIFT_J(:,mytimestep) - jold) / ABS(jold)
               THRIFT_J(:,mytimestep) = jold + alpha*(THRIFT_J(:,mytimestep)-jold)
               jold = THRIFT_J(:,mytimestep)
            END IF

            ! Print Header
            IF (lverb .and. lfirst_pass) THEN
               WRITE(6,*)'    T     NSUB     ITOR     MAX(deltaj)'
               WRITE(6,*)'---------------------------------------'
            END IF

            ! Print progress
            WRITE(6,'(2X,F5.3,1X,I5,1X,ES10.2,1X,ES10.2)') &
                THRIFT_T(mytimestep),nsubsteps,SUM(THRIFT_J(:,mytimestep)),MAXVAL(deltaj)

            ! Turn off screen output after one run
            lscreen_subcodes = .FALSE.

            ! End of first pass
            lfirst_pass = .FALSE.

         END DO
      END DO

      ! Deallocate helpers
      DEALLOCATE(deltaj,jold)
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_evolve

