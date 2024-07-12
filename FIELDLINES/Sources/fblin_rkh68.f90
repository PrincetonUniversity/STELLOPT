!-----------------------------------------------------------------------
!     Function:      fblin_rkh68
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/16/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following (Runge-Kutta-Huta 6th order 8-stage).
!-----------------------------------------------------------------------
      SUBROUTINE fblin_rkh68(phi,q,qdot,istat)
!-----------------------------------------------------------------------
!     Libraries None
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: istat
      DOUBLE PRECISION :: phi, q(2), qdot(2)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      istat    = 0
      CALL fblin_nag(phi,q,qdot)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin_rkh68
