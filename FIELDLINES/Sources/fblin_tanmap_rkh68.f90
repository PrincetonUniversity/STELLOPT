!-----------------------------------------------------------------------
!     Function:      fblin_tanmap_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following (NAG) with tangent map
!                    integration.
!-----------------------------------------------------------------------
      SUBROUTINE fblin_tanmap_rkh68(phi,q,qdot,istat)
!-----------------------------------------------------------------------
!     Libraries NONE
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2),q(3),q(4),q(5),q(6)) = (R,Z,T11,T12,T21,T22)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: istat
      DOUBLE PRECISION :: phi, q(6), qdot(6)
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      CALL fblin_tanmap_nag(phi,q,qdot)
      istat = 0
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin_tanmap_rkh68
