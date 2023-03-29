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
      USE thrift_funcs
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lfirst_pass, lfirst_sub_pass
      INTEGER :: i, ier
      REAL(rprec) :: alpha, rho, s, aminor1, aminor2, temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: deltaj, jold
      CHARACTER(len = 16)     :: temp1_str, temp2_str
      CHARACTER(len = 79)     :: header_str,progress_str
      CHARACTER(len = 79)     :: temp_prog_str
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! For when we collect the workers
      IF (myworkid .ne. master) RETURN

      ! Initialize the current density
      THRIFT_J        = 0; THRIFT_JPLASMA  = 0; THRIFT_JSOURCE  = 0
      THRIFT_JBOOT    = 0; THRIFT_JECCD    = 0; THRIFT_JNBCD    = 0
      THRIFT_JOHMIC   = 0; THRIFT_JBOOT_S  = 0;
      ! Initialize enclosed currents
      THRIFT_I        = 0; THRIFT_IPLASMA  = 0; THRIFT_ISOURCE  = 0
      THRIFT_IBOOT    = 0; THRIFT_IECCD    = 0; THRIFT_INBCD    = 0
      THRIFT_IOHMIC   = 0; THRIFT_UGRID    = 0;   
      ! Initialize profile variables
      THRIFT_P        = 0; THRIFT_PPRIME   = 0; THRIFT_ETAPARA  = 0
      ! Initialize magnetic variables
      THRIFT_VP       = 0; THRIFT_PHIEDGE  = 0; THRIFT_S11      = 0
      THRIFT_BAV      = 0; THRIFT_BSQAV    = 0
      THRIFT_AMINOR   = 0; THRIFT_RMAJOR   = 0
      ! Initialize coefficients
      THRIFT_COEFF_A = 0; THRIFT_COEFF_B = 0; THRIFT_COEFF_C = 0; THRIFT_COEFF_D = 0;
      THRIFT_COEFF_BP= 0; THRIFT_COEFF_CP= 0; THRIFT_COEFF_DP= 0;
      THRIFT_ALPHA1  = 0; THRIFT_ALPHA2  = 0; THRIFT_ALPHA3  = 0; THRIFT_ALPHA4  = 0;
      THRIFT_MATLD   = 0; THRIFT_MATMD   = 0; THRIFT_MATUD   = 0; THRIFT_MATRHS  = 0;


      ! Allocate the convergence helper
      ALLOCATE(deltaj(nrho), jold(nrho))
      alpha = 0.05
      jold   = 1E3 ! so on loop 1 we don't divide by zero
      lscreen_subcodes = .TRUE.
      lfirst_pass = .TRUE.
      
      ! Loop over timesteps
      write(6,*) 'start loop'
      DO mytimestep = 1, ntimesteps
      
            write(6,*) 'setup p'

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
            WRITE(temp1_str,'(i3.3)') mytimestep
            WRITE(temp2_str,'(i3.3)') nsubsteps
            proc_string = TRIM(TRIM(id_string) // '.' //  &
                  TRIM(ADJUSTL(temp1_str)) // '_' // &
                  TRIM(ADJUSTL(temp2_str)))
                  write(6,*) 'setup j'

            ! Update equilbrium current
            CALL thrift_equil_j(lfirst_sub_pass)

            write(6,*) 'run equil'

            ! Run equilibrium
            CALL thrift_run_equil

            ! Update equilibrium/profile variables
            write(6,*) 'update vars'
            CALL update_vars

            ! Calculate Bootstrap
            write(6,*) 'bootstrap'
            CALL thrift_run_bootstrap
            write(6,*) 'eccd'

            ! Calculate Current Drive
            IF (leccd)  CALL thrift_run_ECCD
            !IF (lnbcd)  CALL thrift_run_NBI
            !IF (lohmic) CALL thrift_run_ohmic

            ! Update total source current (picard iteration)
            !THRIFT_JSOURCE(:,mytimestep) = (1-picard_factor)*THRIFT_JSOURCE(:,mytimestep) &
            !                               + picard_factor*(  THRIFT_JBOOT(:,mytimestep) &
            !                                                + THRIFT_JECCD(:,mytimestep) &
            !                                                + THRIFT_JNBCD(:,mytimestep) &
            !                                                + THRIFT_JOHMIC(:,mytimestep))
            THRIFT_JSOURCE(:,mytimestep) =   THRIFT_JBOOT(:,mytimestep) &
                                           + THRIFT_JECCD(:,mytimestep) &
                                           + THRIFT_JNBCD(:,mytimestep) &
                                           + THRIFT_JOHMIC(:,mytimestep)
                                           write(6,*) 'jinductive'

            ! Update the plasma current  
            CALL thrift_jinductive

            
            ! Update total current
            IF (nsubsteps==1) THEN
              THRIFT_J(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep) &
                                     + THRIFT_JSOURCE(:,mytimestep)
            ELSE
              THRIFT_J(:,mytimestep) = THRIFT_J(:,mytimestep)*(1-picard_factor) &
                     +  picard_factor*(THRIFT_JPLASMA(:,mytimestep) &
                                     + THRIFT_JSOURCE(:,mytimestep))
            END IF

            ! Check the convergence
            deltaj = 0
            IF (nsubsteps==1) THEN
            deltaj = ABS(THRIFT_J(:,mytimestep))
               jold = THRIFT_J(:,mytimestep)
            ELSE
               WHERE(ABS(jold)>0) deltaj = ABS( THRIFT_J(:,mytimestep) - jold) / ABS(jold)
               jold = THRIFT_J(:,mytimestep)
            END IF
            ! Calculate total currents
            CALL curden_to_curtot(THRIFT_JBOOT(:,  mytimestep),THRIFT_IBOOT(:,  mytimestep))
            CALL curden_to_curtot(THRIFT_JECCD(:,  mytimestep),THRIFT_IECCD(:,  mytimestep))
            CALL curden_to_curtot(THRIFT_JNBCD(:,  mytimestep),THRIFT_INBCD(:,  mytimestep))
            CALL curden_to_curtot(THRIFT_JOHMIC(:, mytimestep),THRIFT_IOHMIC(:, mytimestep))
            CALL curden_to_curtot(THRIFT_JPLASMA(:,mytimestep),THRIFT_IPLASMA(:,mytimestep))
            THRIFT_ISOURCE(:,mytimestep)  = THRIFT_IBOOT(:,  mytimestep) &
                                          + THRIFT_IECCD(:,  mytimestep) &
                                          + THRIFT_INBCD(:,  mytimestep) &
                                          + THRIFT_IOHMIC(:, mytimestep)
            THRIFT_I(:,mytimestep)        = THRIFT_IPLASMA(:,mytimestep) &
                                          + THRIFT_ISOURCE(:,mytimestep)

            ! Print Header
            IF (lverb .and. lfirst_pass) THEN
            
               header_str = '   T  NSUB  BETA        ITOR     IPLASMA       IBOOT'
               IF (leccd)  header_str = TRIM(header_str)//'       IECCD'
               IF (lnbcd)  header_str = TRIM(header_str)//'       INBCD'
               IF (lohmic) header_str = TRIM(header_str)//'      IOHMIC'
               header_str = TRIM(header_str)//' MAX dJ/JOLD'
               
               WRITE(6,*)''
               WRITE(6,*) header_str
               WRITE(6,*)'==============================================================================='
            END IF

            ! Print progress
            IF (lverb) THEN
                  WRITE(progress_str,'(1X,F6.3,1X,I2,1X,F5.2,3(1X,ES11.3))') &
                  THRIFT_T(mytimestep),nsubsteps,eq_beta*100,THRIFT_I(nsj,mytimestep),&
                  THRIFT_IPLASMA(nsj,mytimestep), THRIFT_IBOOT(nsj,mytimestep)
                  IF (leccd) THEN
                        WRITE(temp_prog_str,'(ES11.3)') THRIFT_IECCD(nsj,mytimestep)
                        progress_str = TRIM(progress_str)//" "//temp_prog_str
                  END IF
                  IF (lnbcd) THEN
                        WRITE(temp_prog_str,'(ES11.3)') THRIFT_INBCD(nsj,mytimestep)
                        progress_str = TRIM(progress_str)//" "//temp_prog_str
                  END IF
                  IF (lohmic) THEN
                        WRITE(temp_prog_str,'(ES11.3)') THRIFT_IOHMIC(nsj,mytimestep)
                        progress_str = TRIM(progress_str)//" "//temp_prog_str
                  END IF
                  WRITE(temp_prog_str,'(ES11.3)') MAXVAL(deltaj)
                  progress_str = TRIM(progress_str)//" "//temp_prog_str
                WRITE(6,*) progress_str
            END IF
            ! Turn off screen output after one run
            lscreen_subcodes = .FALSE.

            ! End of first pass
            lfirst_pass = .FALSE.
            lfirst_sub_pass = .FALSE.

         END DO

      END DO

      ! Deallocate helpers
      DEALLOCATE(deltaj,jold)

      ! Get back workers
      CALL thrift_paraexe('exit','test',lscreen_subcodes)
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_evolve

