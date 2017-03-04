!-----------------------------------------------------------------------
!     Function:      fxmn_monit
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine is the monitoring funciton for the
!                    iterative xmn calculation
!-----------------------------------------------------------------------
      SUBROUTINE fxmn_monit(fmin,fmax,sim,n,n1,ncall)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Variables
!          fmin      Smallest function value in simplex
!          fmax      Largest function value in simplex
!          sim       Position of the current simplex
!          n         Number of variables
!          nvert     Number of vertices in simplex
!          ncall     Number of function evaluations
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: n, n1, ncall
      DOUBLE PRECISION :: fmin, fmax
      DOUBLE PRECISION :: sim(n1,n)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      WRITE(6,'(1X,A,I5,A,F10.8)') 'After ',ncall,' function calls, the value is',fmin
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fxmn_monit
