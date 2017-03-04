!-----------------------------------------------------------------------
!     Subroutine:    stelltran_runeq
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/21/2015
!     Description:   This routine handels execution of the equilibrium.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_runeq(itime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE vmec_input
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
      USE vmec_main, ONLY:  multi_ns_grid
      USE mpi_params
      
      USE vmec_input
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          ier        Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: itime
      
      LOGICAL :: lscreen
      INTEGER :: ier
      INTEGER :: vctrl_array(5)
      REAL(rprec) :: curtor_first
      CHARACTER(256) :: reset_string, temp_str
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------
      WRITE(temp_str,'(i5.5)') itime
      proc_string = TRIM(TRIM(id_string) // '_time' // TRIM(ADJUSTL(temp_str)))
	
      ! Run the equilibrium
      SELECT CASE (TRIM(equil_type))
         CASE('vmec2000')
            ! Update profiles
            IF (lveryfirst_pass) curtor_first = curtor
            CALL stelltran_profiles_vmec
            IF (lveryfirst_pass) curtor = curtor_first
            ! Update the internal vars
            CALL stelltran_reinit_vmec
            ! Setup run
            vctrl_array(1) = restart_flag+timestep_flag+output_flag+reset_jacdt_flag ! Need restart to get profile variations
            vctrl_array(2) = 0     ! vmec error flag  
            vctrl_array(3) = -1    ! Use multigrid
            vctrl_array(4) = MAXLOC(ns_array,DIM=1)
            vctrl_array(5) = myid ! Output file sequence number
            reset_string = proc_string
            lscreen = .false.
            IF (lfirst_pass) THEN
               vctrl_array(4) = -1
               reset_string = proc_string_old
            END IF
            IF (lveryfirst_pass) THEN
               reset_string = ''
               vctrl_array(1) = timestep_flag+output_flag+reset_jacdt_flag
               vctrl_array(4) = 0
               lscreen = .true.
            END IF
            ! Run VMEC
            CALL runvmec(vctrl_array,proc_string,lscreen,0,reset_string)
            ier=vctrl_array(2)
            IF (ier /= successful_term_flag  .or. &
                   ier /= norm_term_flag) THEN
                WRITE(6,'(A)') !!!!! VMEC failed to converge  !!!!!'
            END IF
      END SELECT
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_runeq
