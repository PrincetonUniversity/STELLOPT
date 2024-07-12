      SUBROUTINE lmdif1_mp(fcn, m, n, x, fvec, tol, epsfcn,
     1           nfev_end, diag, mode, info, lwa)
!
!!    THIS IS THE MULTIPLE PROCESSOR VERSION OF LEV_OPT USING MPI CALLS
!!          (D. A. Spong 10/27/00)

!!!   epsfcn was added as a dummy argument, eliminating the circular
!!!   MODULE reference.
!!    diag and mode were added as arguments, so that different choices
!!    for these can be tried in the future (as suggested by L. A. Berry)
!!    Currently mode is set = 1 in the calling routine (stellopt) so that
!!    diag is internally generated.  See definitions of mode and diag below.
!!
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, lwa, nfev_end, mode
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(in) :: tol, epsfcn
      REAL(rprec), DIMENSION(n), INTENT(inout) :: x, diag
      REAL(rprec), DIMENSION(m), INTENT(out) :: fvec
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn
C-----------------------------------------------
      CALL lmdif1(fcn, m, n, x, fvec, tol, epsfcn, nfev_end,
     1           diag, mode, info, lwa, 0, 0)

      END SUBROUTINE lmdif1_mp
