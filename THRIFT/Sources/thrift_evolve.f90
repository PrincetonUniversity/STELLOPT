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
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Loop over timesteps
      DO mytimestep = 1, ntimesteps
         
         ! Initialize J
         IF (mytimestep == 1) THEN
            THRIFT_J(:,mytimestep) = 0
         END IF

         ! Iterate over J_Source till convereged
         ! CALL thrift_sourceterm

         ! Calc the DT

      END DO
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_evolve

