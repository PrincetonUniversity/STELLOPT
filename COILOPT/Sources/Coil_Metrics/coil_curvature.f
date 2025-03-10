      SUBROUTINE coil_curvature
      USE boundary
      USE modular_coils
      USE Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, ncpts
      REAL(rprec) :: crv_max, crv, csum
      REAL(rprec) :: u0, v0, du
      REAL(rprec) :: tx, ty, tz
!-----------------------------------------------

      ncpts = 4*nwire
      du = 1.0_dp/ncpts
      DO i=1,nmid
        crv_max = 0
        csum = 0
        DO j = 1,ncpts
          u0 = (j-1)*du
          CALL modular_curve (i, u0, v0, tx, ty, tz, crv)
          crv_max = MAX(crv_max, crv)
          csum = csum + crv
        END DO
        rc_min(i) = 10000.0_dp
        IF (crv_max .gt. 0.0_dp) rc_min(i) = 1.0_dp/crv_max
        cu_sum(i) = csum*du
      END DO

      END SUBROUTINE coil_curvature
