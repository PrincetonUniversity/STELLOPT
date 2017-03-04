      SUBROUTINE chisq_curvature (sigma, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE boozer_params
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num
      REAL(rprec) :: sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: test_pts = 4
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lzth
      REAL(rprec) :: test_fcn
      REAL(rprec), DIMENSION(4) :: curv_kur
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL flux_surf_curv
C-----------------------------------------------

      IF (nopt .gt. 0) THEN
         CALL flux_surf_curv (curv_kur, rmnc_bdy, zmns_bdy,
     1        xm_bdy, xn_bdy)
         DO lzth = 1, test_pts
            num = num + 1
            index_array(num) = ivar_curve
c
c       As a target to prevent flux surface cusping, use a tanh functional of
c       the curvature kurtosis around 0,90,180, and 270 degree phi planes.
c       This is designed to allow moderate elliptical and triangular elongations,
c       but not strong cusps. note: the numbers in this functional may need
c       readjustment for particular situations.
c
            test_fcn = 0.44_dp + p5*
     1        TANH((curv_kur(lzth)-20.0_dp)/15.0_dp)
            wegt(num) = sigma
            chisq_target(num) = 1.e-3_dp
            chisq_match(num) = test_fcn
         END DO
      ELSE
         DO lzth = 1, test_pts
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_curve)
         END DO
      END IF


      END SUBROUTINE chisq_curvature
