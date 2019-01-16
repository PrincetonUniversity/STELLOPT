!-----------------------------------------------------------------------
!     Subroutine:    stellopt_mercier
!     Authors:       J.Schmitt (jcschmitt at auburn dot edu)
!     Date:          2019
!     Description:   This subroutine calculates the Mercier growth
!                    rate.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_mercier(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE stellopt_targets, ONLY: sigma_mercier_criterion, target_mercier_criterion
      USE equil_vals, ONLY: mercier_criterion
      ! VMEC
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
!      INTEGER ::  ier, ik
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  MERCIER STABLITY CALCULATION  ------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('animec','flow','satire','spec')
               print *, 'K============= UNTESTED SECTION OF THE CODE'
         CASE('vmec2000','parvmec','paravmec','vboot','vmec2000_oneeq')
            IF (iflag .ne. 0) RETURN

            IF (.not. ALLOCATED(mercier_criterion)) PRINT *, "<----mercier_criterion not allocated in stellopt_mecier"
            IF (.not. ALLOCATED(mercier_criterion)) STOP
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' -------------------------  MERCIER CALCULATION DONE  ---------------------'
      CALL FLUSH(6)
      RETURN
  48  format('====================================================')
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_mercier

