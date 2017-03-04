!-----------------------------------------------------------------------
!     Function:      monit_find_axis
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/30/2012
!     Description:   This subroutine monitors the find axis routine.
!-----------------------------------------------------------------------
      SUBROUTINE monit_find_axis(fmin, fmax, sim, n, is, ncall, iuser, ruser)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Variables
!          r0         Initial R position
!          phi0       Initial phi position
!          z0         Initial z position
!          phi1       Final phi position
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: n, is, ncall, iuser(*)
      REAL    :: fmin, fmax, sim(is,n), ruser(*)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      PRINT *,'           ',ncall,fmin,fmax,sim(1,1),sim(2,1),sim(3,1),sim(1,2),sim(2,2),sim(2,3)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE monit_find_axis
