!-----------------------------------------------------------------------
!     Subroutine:    stellopt_coil_to_vac
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/15/2016
!     Description:   This subroutine reads a coils file and generates
!                    the appropriate vacuum grid file.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_coil_to_vac(lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY: proc_string
      USE stellopt_vars, ONLY: equil_type
      USE write_mgrid, only: mgrid_ext, lstell_sym
      USE makegrid_global, only: task, lscreen_mgrid => lscreen
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        lscreen   Terminal output
!----------------------------------------------------------------------
      LOGICAL :: lscreen
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      CALL tolower(equil_type)
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec')
            ! Adjust the namelist variables
            task='MGRID'
            mgrid_ext=TRIM(proc_string)
            lscreen_mgrid = lscreen
            ! Call the Routine
            
      END SELECT
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_coil_to_vac
