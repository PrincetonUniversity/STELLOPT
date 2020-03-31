!-----------------------------------------------------------------------
!     Subroutine:    stellopt_write_eqfile
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          11/12/2019
!     Description:   This subroutine handles writing out of equilibrium
!                    files.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_write_eqfile
      USE stellopt_runtime
      USE stellopt_input_mod
      USE mpi_inc
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
      USE read_wout_mod, ONLY: read_wout_file, write_wout_file, read_wout_deallocate
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        NONE
!----------------------------------------------------------------------
      IMPLICIT NONE

!----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        vcrtl_array VMEC Control Array
!----------------------------------------------------------------------
      INTEGER                :: ier
      INTEGER                ::  vctrl_array(5)

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ier = 0
      SELECT CASE(TRIM(equil_type))
         CASE ('vmec2000_old','animec','flow','satire')
            vctrl_array(1) = output_flag ! Output to file
            vctrl_array(2) = 0     ! vmec error flag  
            vctrl_array(3) = 0    ! Use multigrid
            vctrl_array(4) = 0     ! Iterative 
            vctrl_array(5) = myid ! Output file sequence number
            CALL runvmec(vctrl_array,proc_string,.false.,MPI_COMM_SELF,'')
         CASE('parvmec','paravmec','vmec2000','vboot')
            CALL stellopt_paraexe('paravmec_write',proc_string,.false.)
         CASE('vmec2000_oneeq')
            CALL read_wout_deallocate
            CALL read_wout_file(TRIM(proc_string_old),ier)
            CALL write_wout_file('wout_'//TRIM(proc_string)//'.nc',ier)
         CASE('test')
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_write_eqfile