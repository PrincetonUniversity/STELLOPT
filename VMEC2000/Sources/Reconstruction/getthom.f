      SUBROUTINE getthom(amat_i, amat_p, data_array, re, ro, kcthom)
      USE vmec_main
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER kcthom
      REAL(rprec), DIMENSION(isnodes,*) :: amat_i
      REAL(rprec), DIMENSION(ipnodes,*) :: amat_p
      REAL(rprec), DIMENSION(*) :: data_array
      REAL(rprec), DIMENSION(ns,nzeta,*) :: re, ro
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, ks
      INTEGER, DIMENSION(itse + 2) :: isortp
      REAL(rprec), DIMENSION(itse) :: datalsq_p
      REAL(rprec), DIMENSION(itse + 2) :: rgrid
      REAL(rprec):: sig_fac, t1
      REAL(rprec) :: rmat1u(itse)
      LOGICAL :: l1v(itse)
C-----------------------------------------------

!
!       THIS SUBROUTINE COMPUTES THE LEAST-SQUARES AMATRIX AND DATA
!       ARRAYS FOR MATCHING TO THE PRESSURE DATA AT EQUALLY-SPACED KNOTS
!       IN SQRT(PHI-FLUX) SPACE (IPNODES, INCLUDING S=0 AND S=1)
!
!       COMING INTO THIS ROUTINE, THE RTHOM, DATATHOM HAVE BEEN
!       PREVIOUSLY ORDERED SO RTHOM(1) < RTHOM(2) < ....
!

!       IF (LPOFR), user has input P(R,Z)
!       If (.NOT.LPOFR),  Then user has input P(s), not P(R)
!

      IF (.not.lpofr) THEN                       !p(R) or p(s) ?

!       CONSTRUCT RTHOM BASED ON P(s)
         rthom(:itse) = sthom(:itse)

         CALL pofs (re, ro, ns, rthom, itse)
         rthommax = rthom(itse)
         rthommin = rthom(1)

      ENDIF

!
!       IF NO PRESSURE DATA, MAKE SURE CHISQ-THOM <= CHI_target
!
      sig_fac = one

!
!       CONSTRUCT EVERYTHING BASED ON P(R)
!       FOR POINTS OUTSIDE GRID, SET R = either rmin,rmax
!
      kcthom = 0
      IF (itse .gt. 0) THEN
         l1v(:itse) = .false.
         datalsq_p(:itse) = datathom(:itse)*pfac   !sorted data array
         sigma_thom(:itse) = sigma_thom(:itse)/sig_fac
         pcalc(:itse) = 1.0/sigma_thom(:itse)      !pcalc = sorted sigma array
         rmat1u(:itse) = rthom(:itse)
         WHERE (rmat1u(:itse) .lt. rinner)
            rmat1u(:itse) = rinner
         ELSEWHERE
            l1v(:itse) = rmat1u(:itse) .gt. router
         END WHERE
         WHERE (l1v(:itse)) rmat1u(:itse) = router
         rgrid(kcthom+1:itse+kcthom) = rmat1u(:itse)
         kcthom = itse + kcthom
      ENDIF


!
!       FIND S,THETA INDICES FOR GRIDDED R-THOM ARRAY (RGRID)
!
      CALL findphi (re, ro, rgrid, delse2, delso2, rmid, indexs2,
     1   indexu2, indexr, kcthom)

!
!       COMPUTE MATRIX ELEMENTS FOR PRESSURE SPLINE NODES CORRESPONDING
!       TO ORDERED RGRID ARRAY
!
      CALL getspline (amat_p, pknots, hthom, delse2, hs, indexs2,
     1   isortp, kcthom, ipnodes)

!
!       MATCH PRESSURE SPLINE KNOTS TO THOMSON SCATTERING DATA
!
!       ON ENTRY INTO THIS LOOP, PCALC = 1/SIGMA_THOM
!

      DO i = 1, kcthom
         ks = isortp(i)         !Index BEFORE sorting on pknots (INDEXed

         data_array(i) = datalsq_p(ks)*pcalc(ks)
         t1 = pthommax*pcalc(ks)
         amat_p(:ipnodes,i) = t1*amat_p(:ipnodes,i)
      END DO
      IF (.not.lpofr) rthompeak = rgrid(isortp(1))

      itse2 = kcthom
      amat_i(:,:itse2) = zero

      END SUBROUTINE getthom
