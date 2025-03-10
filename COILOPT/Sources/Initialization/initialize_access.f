      SUBROUTINE initialize_access
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE bcoils_mod
      USE Vcoilpts
      USE Vwire
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, status
      REAL (rprec) :: dx_access, dy_access, dz_access
!-----------------------------------------------

      n_access_pts = nwdim1

      ALLOCATE(x_access(ncdim, nwdim1), y_access(ncdim, nwdim1),
     1         z_access(ncdim, nwdim1), stat=status)
      IF(status /= 0) STOP "Cannot ALLOCATE access points"

      DO n = 1, n_access
         dx_access = (x1_access(n) - x0_access(n))/(n_access_pts - 1)
         dy_access = (y1_access(n) - y0_access(n))/(n_access_pts - 1)
         dz_access = (z1_access(n) - z0_access(n))/(n_access_pts - 1)
         DO i = 1, n_access_pts
            x_access (n,i) = x0_access (n) + (i-1)*dx_access
            y_access (n,i) = y0_access (n) + (i-1)*dy_access
            z_access (n,i) = z0_access (n) + (i-1)*dz_access
         END DO
      END DO

      END SUBROUTINE initialize_access
