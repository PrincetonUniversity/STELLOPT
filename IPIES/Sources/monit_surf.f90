!-----------------------------------------------------------------------
!     Subroutine:    monit_surf
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/12/2011
!     Description:   This is the monitoring function for linetomn.
!-----------------------------------------------------------------------
      SUBROUTINE monit_surf(fmin,fmax,sim,n,is,ncall,iuser,ruser)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Parameters
!          n           Number of Fourier coefficients
!          xc          Fourier coefficients (1:n)=(1:mnmax)
!          fc          Total distance to points of field-line
!          iuser(1)    NAG work array - Number of points in field-line
!          iuser(2)    NAG work array - Signs
!          iuser(3)    NAG work array - Flag to perform pre-calc (0=No,1=YES)
!          ruser       NAG work array - xu,xv,Fln
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER             :: n, is, ncall, iuser(*)
      DOUBLE PRECISION    :: fmin, fmax, sim(is,n), ruser(*)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      PRINT *,'MONIT_SURF',fmin,fmax,n,is,ncall
      CALL FLUSH(6)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE monit_surf
