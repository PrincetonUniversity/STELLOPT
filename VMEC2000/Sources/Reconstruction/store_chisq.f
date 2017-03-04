      SUBROUTINE store_chisq
      USE vmec_main
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
!
!       COMPUTES CHI**2 FROM DIFFERENT SUBROUTINES AT VARIOUS TIME-STEPS
!       WRITTEN BY D.K. LEE (3/93)
!

      IF (MOD(iter2,nstep).ne.10 .and. iequi.eq.0) RETURN
         chisqerr(jchix1) = SUM(chisqerr(:jchix))
         IF (.not.lpprof) chisqerr(jchix1) = chisqerr(jchix1)
     1      - chisqerr(ithom0)
         nchistp = nchistp + 1
         IF (nchistp .gt. mstp) RETURN
         chi2(:,nchistp) = chisqerr
         nchi2(nchistp) = iter2 - 10
         IF (iequi .eq. 1) nchi2(nchistp) = iter2
         IF (iter2 .eq. 10) nchi2(nchistp) = 1

      END SUBROUTINE store_chisq
