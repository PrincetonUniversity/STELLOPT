      SUBROUTINE init_tf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE tf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nvariables, nx_max, ny_max, nx, ny, ntf
      REAL(rprec) :: xvariables(*)
      REAL(rprec) :: tfc_xmin, tfc_xmax
      REAL(rprec) :: tfc_ymin, tfc_ymax
      REAL(rprec) :: tfc_zmin, tfc_zmax
      REAL(rprec) :: tfc_dphi, tfc_phi, tfc_rad
      REAL(rprec) :: tfc_dx, tfc_dy, tfc_scl
!-----------------------------------------------
      mtfcoil = 1
      mtfwire = 2

!     Model QPS tf coil with RETURN legs ...
      nx_max = 6
      ny_max = 2
      tfc_dphi = twopi/(nx_max*ny_max)
      tfc_phi = 0.5_dp*tfc_dphi
!     tfc_scl = 0.9_dp ! since the following are engineering values
      tfc_scl = 1.0_dp ! for vmec input based on actual plasma SIZE
      tfc_rad = 2.15_dp/tfc_scl
      tfc_xmin = -0.375_dp/tfc_scl
      tfc_xmax =  0.375_dp/tfc_scl
      tfc_dx = 0.150_dp/tfc_scl
      tfc_ymin = 0._dp
      tfc_ymax = 0.013_dp/tfc_scl
      tfc_dy = 0.026_dp/tfc_scl
      tfc_zmin = -2.48_dp/tfc_scl
      tfc_zmax =  2.48_dp/tfc_scl
      IF (lqos) THEN
         mtfwire = 5
!        distribute straight tf filaments in a rectangle in x-y plane
         mtfcoil = nx_max*ny_max
         ntf = 0
         DO ny = 1, ny_max
            DO nx = 1, nx_max
               ntf = ntf + 1
               IF (ny .eq. 1) THEN
                  tfc_x(ntf,1) = (tfc_xmax - (nx-1)*tfc_dx)
               ELSE
                  tfc_x(ntf,1) = (tfc_xmin + (nx-1)*tfc_dx)
               END IF
               tfc_x(ntf,2) = tfc_x(ntf,1)
               tfc_x(ntf,3) = tfc_rad*COS(tfc_phi)
               tfc_x(ntf,4) = tfc_x(ntf,3)
               tfc_x(ntf,5) = tfc_x(ntf,1)
               tfc_y(ntf,1) = (tfc_ymax - (ny-1)*tfc_dy)
               tfc_y(ntf,2) = tfc_y(ntf,1)
               tfc_y(ntf,3) = tfc_rad*SIN(tfc_phi)
               tfc_y(ntf,4) = tfc_y(ntf,3)
               tfc_y(ntf,5) = tfc_y(ntf,1)
               tfc_z(ntf,1) = tfc_zmax
               tfc_z(ntf,2) = tfc_zmin
               tfc_z(ntf,3) = tfc_zmin
               tfc_z(ntf,4) = tfc_zmax
               tfc_z(ntf,5) = tfc_zmax
               tfc_cur(ntf) = i_tfc/mtfcoil
               tfc_phi = tfc_phi + tfc_dphi
            END DO
         END DO
      ELSE
!        single filament at x=0, y=0
         tfc_zmin = -1000
         tfc_zmax =  1000
         tfc_x(1,1) = 0
         tfc_x(1,2) = 0
         tfc_y(1,1) = 0
         tfc_y(1,2) = 0
         tfc_z(1,1) = tfc_zmax
         tfc_z(1,2) = tfc_zmin
         tfc_cur(1) = i_tfc
      END IF

      xvariables (1) = i_tfc
      nvariables = 1

      END SUBROUTINE init_tf_coils
