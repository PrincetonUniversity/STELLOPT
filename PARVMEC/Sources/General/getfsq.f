#if defined(SKS)
      SUBROUTINE getfsq_par(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: rprec, ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: medge
      REAL(rprec), INTENT(out) :: gnormr, gnormz
      REAL(rprec), INTENT(in)  :: gnorm
      REAL(rprec), DIMENSION(mnsize,ns,ntmax), INTENT(in) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jsmax, nsmin, nsmax
      REAL(rprec) :: totalgcr, totalgcz, tmpgcr, tmpgcz
      REAL(rprec) :: skston, skstoff
      REAL(rprec), DIMENSION(2) :: tmpgcx, totalgcx
!-----------------------------------------------

      jsmax = ns1 + medge
      nsmin=tlglob; nsmax=MIN(trglob,jsmax)


      IF (lactive) THEN
        CALL second0(skston)
        tmpgcr=SUM(gcr(:,nsmin:nsmax,:)**2)
        tmpgcz=SUM(gcz(:,nsmin:nsmax,:)**2)
        tmpgcx(1)=tmpgcr; tmpgcx(2)=tmpgcz
        CALL MPI_Allreduce(tmpgcx,totalgcx,2,MPI_REAL8,MPI_SUM,
     1                     NS_COMM,MPI_ERR)
        totalgcr=totalgcx(1); totalgcz=totalgcx(2)
        CALL second0(skstoff)
        allreduce_time = allreduce_time + (skstoff - skston)
      END IF

      gnormr = gnorm * totalgcr
      gnormz = gnorm * totalgcz

      END SUBROUTINE getfsq_par
#endif

      SUBROUTINE getfsq(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: rprec, ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module
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
