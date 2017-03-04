!-----------------------------------------------------------------------
!     Function:      fphif
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/2/2011
!     Description:   This function determines the error due to following
!                    field lines a finite length.
!                    Note that phgaus is hardcoded as 1e-8.
!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION fphif(phif)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
!-----------------------------------------------------------------------
!     Input Parameters
!          phif          Length of field line in phi
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: phif
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      DOUBLE PRECISION    :: pi2
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------
      pi2 = 4 * ATAN(1._rprec)
      fphif = exp(-(phif/phgaus)*(phif/phgaus)) * phgaus / (phif * sqrt(pi2)) - 0.5 * follow_tol
      RETURN
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------
      END FUNCTION fphif
