      SUBROUTINE getfsq(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: rprec, ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: medge
      REAL(rprec), INTENT(out) :: gnormr, gnormz
      REAL(rprec), INTENT(in)  :: gnorm
      REAL(rprec), DIMENSION(ns,mnsize*ntmax), INTENT(in) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jsmax
!-----------------------------------------------
      jsmax = ns1 + medge
      gnormr = gnorm * SUM(gcr(:jsmax,:)**2)
      gnormz = gnorm * SUM(gcz(:jsmax,:)**2)

      END SUBROUTINE getfsq
