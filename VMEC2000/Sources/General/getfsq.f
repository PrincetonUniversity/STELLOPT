      SUBROUTINE getfsq_par(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: rprec, ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: medge
      REAL(dp), INTENT(out) :: gnormr, gnormz
      REAL(dp), INTENT(in)  :: gnorm
      REAL(dp), DIMENSION(mnsize,ns,ntmax), INTENT(IN) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jsmax, nsmin, nsmax, l
      REAL(dp) :: tmpgcx(ns,2), totalgcx(2)
!-----------------------------------------------
      IF (.NOT. lactive) RETURN

      jsmax = ns1 + medge
      nsmin=tlglob; nsmax=MIN(trglob,jsmax)
      IF (trglob .GT. jsmax) tmpgcx(jsmax+1:trglob,1:2) = 0

      DO l = nsmin, nsmax
         tmpgcx(l,1) = SUM(gcr(:,l,:)**2)
         tmpgcx(l,2) = SUM(gcz(:,l,:)**2)
      END DO
      DO l = 1, 2
        CALL Gather1XArray(tmpgcx(:,l))
        totalgcx(l) = SUM(tmpgcx(:,l))
      END DO

      gnormr = gnorm * totalgcx(1)
      gnormz = gnorm * totalgcx(2)

      END SUBROUTINE getfsq_par

      SUBROUTINE getfsq(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: rprec, ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: medge
      REAL(dp), INTENT(out) :: gnormr, gnormz
      REAL(dp), INTENT(in)  :: gnorm
      REAL(dp), DIMENSION(ns,mnsize*ntmax), INTENT(in) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jsmax
!-----------------------------------------------
      jsmax = ns1 + medge
      gnormr = gnorm * SUM(gcr(:jsmax,:)**2)
      gnormz = gnorm * SUM(gcz(:jsmax,:)**2)

      END SUBROUTINE getfsq
