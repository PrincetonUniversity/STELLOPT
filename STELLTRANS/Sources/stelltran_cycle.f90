!-----------------------------------------------------------------------
!     Subroutine:    stelltran_cycle
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/21/2015
!     Description:   This is the main cycle routine for the transport
!                    simulation.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_cycle
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_data
!-----------------------------------------------------------------------
!     Local Variables
!          lconsistent    Logical to determine equilibrium consistency
!          itime          Timeslice index
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lconsistent
      INTEGER :: itime
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      lveryfirst_pass = .TRUE.
      DO itime=1, ntimesteps
         ! Update the equilibrium
         CALL stelltran_updateeq(itime)
         
         ! Set some pass parameters
         lconsistent = .FALSE.
         lfirst_pass = .TRUE.
         
         ! Iterate till the profiles are consistent
         DO WHILE (.not. lconsistent)
            ! Calculate an equilibrium
            CALL stelltran_runeq(itime)
            
            ! Calculate basic stuff (volume, etc.)
            CALL stelltran_loadeq
            
            ! Calculate and check current
            CALL stelltran_current(itime)
            	  
            ! Calculate and check species profiles
            ! Don't need this handled by stelltran_updateeq
            !CALL stelltran_thermal(itime)
         
            ! Print a default message at the top
            IF (lverb .and. lveryfirst_pass) &
               WRITE(6,'(A)') '===iter=====Te[keV]==Ti[keV]==Ne[10^20]==Zeff====Johm=====Jboot====Jrf======Jbeam==='
            	    
            ! Check the consistency of the equilibirum
            CALL stelltran_eqconsistency(itime,lconsistent)
            lveryfirst_pass = .FALSE.
            lfirst_pass = .FALSE.
         END DO

         SELECT CASE(TRIM(run_type))
            CASE('analysis')
                  ! If analytic, use profiles to find transport coefficients
                  CALL stelltran_balance
            CASE('predictive')
                  ! If predictive, use equilibrium to calculate neoclassical transport coefficients
                  ! For initial testing, just load in the things from sfincs we already have.
                  ! Load the next coefficients
                  CALL stelltran_find_coeffs(itime) ! Form of this is pretty... specific currently
                  CALL stelltran_bal_update(itime)
                  print*, 'here after bal update'
         END SELECT

         ! Calculate the particle/power balance if requested
!          IF (.true.) THEN
!             CALL stelltran_part_balance
!             CALL stelltran_pow_balance
!          END IF 
         ! Save the run to use in the next timestep
         proc_string_old = proc_string
         
         ! Save the timestep
         ! CALL stelltran_save
         	 
      END DO
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_cycle
