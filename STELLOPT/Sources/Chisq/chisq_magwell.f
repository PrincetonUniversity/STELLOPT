      SUBROUTINE chisq_magwell(hs, sigma, num, nrad, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, vp_opt, phip_opt
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt
      INTEGER :: num
      REAL(rprec) :: hs, sigma(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad
      REAL(rprec) :: TargetWell, AvgWell, vpf(nrad), sj, vpp
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: pwell_t
C-----------------------------------------------
      IF (nopt .gt. 0) THEN

         vpf(1) = 1.5_dp*vp_opt(2)/phip_opt(2)
     1          - 0.5_dp*vp_opt(3)/phip_opt(3)
         vpf(nrad) = 1.5_dp*vp_opt(nrad)/phip_opt(nrad) -
     1               0.5_dp*vp_opt(nrad-1)/phip_opt(nrad-1)
         vpf(2:nrad-1) = 0.5_dp*(vp_opt(3:nrad)/phip_opt(3:nrad) +
     1                           vp_opt(2:nrad-1)/phip_opt(2:nrad-1))

         AvgWell = 0
         DO jrad = 2, nrad
           sj = hs * (jrad - one)
           AvgWell = AvgWell + ABS(pwell_t(sj)) * hs
         END DO
         DO jrad = 2, nrad
            num = num + 1
            index_array(num) = ivar_well
            vpp = (vpf(jrad)-vpf(1))/vpf(1)
            IF( AvgWell .lt. 1.0e-20_dp )
     1          vpp = (vpf(jrad)-vpf(jrad-1))/vpf(1)
            sj = hs*(jrad - one)
            TargetWell = pwell_t(sj)
            IF( AvgWell .lt. 1.0e-20_dp .and. vpp .lt. 0.0_dp)
     1         TargetWell = vpp
            wegt(num) = sigma(jrad)
            IF(AvgWell .ge. 1.0e-20_dp)
     1            wegt(num) = sigma(jrad)*avgwell
            chisq_target(num) = targetwell
            chisq_match(num) = vpp
         END DO
      ELSE
         DO jrad = 2, nrad
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_well)
         END DO
      ENDIF

      END SUBROUTINE chisq_magwell
