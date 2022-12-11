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
      USE thrift_equil
      USE thrift_vars
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lfirst_pass, lfirst_sub_pass
      INTEGER :: nsubsteps
      REAL(rprec) :: alpha
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: deltaj, jold
      CHARACTER(len = 16)     :: temp1_str, temp2_str
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! For when we collect the workers
      IF (myworkid .ne. master) RETURN

      ! Initialize the current density
      THRIFT_J        = 0
      THRIFT_JBOOT    = 0
      THRIFT_JPLASMA  = 0
      THRIFT_JECCD    = 0
      THRIFT_JNBCD    = 0
      THRIFT_JOHMIC   = 0
      THRIFT_JSOURCE  = 0

      ! Allocate the convergence helper
      ALLOCATE(deltaj(nrho), jold(nrho))
      alpha = 0.05
      jold   = 1E3 ! so on loop 1 we don't divide by zero
      lscreen_subcodes = .TRUE.
      lfirst_pass = .TRUE.
      
      ! Loop over timesteps
      DO mytimestep = 1, ntimesteps
      
         ! Setup the profiles
         CALL thrift_equil_p

         ! Converge Source Currents
         deltaj = 10*jtol; nsubsteps = 0; eq_beta = 1E-9
         lfirst_sub_pass = .TRUE.
         DO WHILE (ANY(deltaj > jtol))

            ! Update Substeps
            nsubsteps = nsubsteps + 1

            ! Get out if too many substeps
            IF (nsubsteps > npicard .or. eq_beta==0) EXIT

            ! Create filename
            WRITE(temp1_str,'(i5)') mytimestep
            WRITE(temp2_str,'(i3)') nsubsteps
            proc_string = TRIM(TRIM(id_string) // '.' //  &
                  TRIM(ADJUSTL(temp1_str)) // '_' // &
                  TRIM(ADJUSTL(temp2_str)))

            ! Update equilbrium current
            CALL thrift_equil_j(lfirst_sub_pass)

            ! Run equilibrium
            CALL thrift_run_equil

            ! Calculate Bootstrap
            CALL thrift_run_bootstrap

            ! Calculate Current Drive
            IF (leccd)  CALL thrift_run_ECCD
            !IF (lnbcd)  CALL thrift_run_NBI
            !IF (lohmic) CALL thrift_run_ohmic

            ! Update total source current (picard iteration)
            THRIFT_JSOURCE(:,mytimestep) = (1-picard_factor)*THRIFT_JSOURCE(:,mytimestep) &
                                           + picard_factor*(  THRIFT_JBOOT(:,mytimestep) &
                                                            + THRIFT_JECCD(:,mytimestep) &
                                                            + THRIFT_JNBCD(:,mytimestep) &
                                                            + THRIFT_JOHMIC(:,mytimestep))
            ! Update total current
            THRIFT_J(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep) &
                                   + THRIFT_JSOURCE(:,mytimestep)


            ! Check the convergence
            deltaj = 0
            IF (nsubsteps==1) THEN
               WHERE(ABS(THRIFT_J(:,mytimestep))>0) deltaj = ABS(THRIFT_J(:,mytimestep))
               jold = alpha*THRIFT_J(:,mytimestep)
            ELSE
               WHERE(ABS(jold)>0) deltaj = ABS( THRIFT_J(:,mytimestep) - jold) / ABS(jold)
               THRIFT_J(:,mytimestep) = jold + alpha*(THRIFT_J(:,mytimestep)-jold)
               jold = THRIFT_J(:,mytimestep)
            END IF

            ! Print Header
            IF (lverb .and. lfirst_pass) THEN
               WRITE(6,*)'    T     NSUB     BETA     ITOR     MAX(deltaj)'
               WRITE(6,*)'------------------------------------------------'
            END IF

            ! Print progress
            IF (lverb) WRITE(6,'(2X,F5.3,3X,I5,4X,F5.2,2(1X,ES10.2))') &
                THRIFT_T(mytimestep),nsubsteps,eq_beta*100,SUM(THRIFT_J(:,mytimestep))/nrho,MAXVAL(deltaj)

            ! Turn off screen output after one run
            lscreen_subcodes = .FALSE.

            ! End of first pass
            lfirst_pass = .FALSE.
            lfirst_sub_pass = .FALSE.

         END DO

         ! Calc the inductive chagne in current
         ! This should change only J_PLASMA given J_SOURCE
         CALL thrift_jinductive

      END DO

      ! Deallocate helpers
      DEALLOCATE(deltaj,jold)

      ! Get back workers
      CALL thrift_paraexe('exit','test',lscreen_subcodes)
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_evolve

