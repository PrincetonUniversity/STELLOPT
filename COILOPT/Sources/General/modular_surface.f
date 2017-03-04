      SUBROUTINE modular_surface (u, v, x, y, z, nx, ny, nz)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE modular_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: k
      REAL(rprec) :: u, v
      REAL(rprec) :: x, y
      REAL(rprec) :: nx, ny, nz, nmag
      REAL(rprec) :: th0, th1, ph0, ph1
      REAL(rprec) :: r, ru, rv
      REAL(rprec) :: z, zu, zv
      REAL(rprec) :: ck, sk
!-----------------------------------------------

!     theta and derivative wrt u

      th0 = twopi*u
      th1 = twopi

!     phi and derivative wrt v

      ph0 = twopi*v/nfper
      ph1 = twopi/nfper

!     compute R, Z and derivatives wrt u and v

      r   = 0
      ru  = 0
      rv  = 0
      z   = 0
      zu  = 0
      zv  = 0
      DO k = 1, numsurf
         ck = COS(m_num(k)*th0 + nfper*n_num(k)*ph0)
         sk = SIN(m_num(k)*th0 + nfper*n_num(k)*ph0)
!     R ...
         r  = r  + rmn_sf(k)*ck
         ru = ru - rmn_sf(k)*m_num(k)*sk
         rv = rv - rmn_sf(k)*n_num(k)*sk
!     Z ...
         z  = z  + zmn_sf(k)*sk
         zu = zu + zmn_sf(k)*m_num(k)*ck
         zv = zv + zmn_sf(k)*n_num(k)*ck
      END DO
      ru = ru*th1
      zu = zu*th1
      rv = rv*ph1*nfper
      zv = zv*ph1*nfper

      sk = SIN(ph0)
      ck = COS(ph0)

!     RETURN X and Y

      x = r*ck
      y = r*sk

!     compute Nx, Ny, Nz

      nx = ru*zv*sk - zu*(rv*sk + r*ck)
      ny = zu*(rv*ck - r*sk) - ru*zv*ck
      nz = r*ru
      nmag = SQRT (nx**2 + ny**2 + nz**2)
      nx = -nx/nmag
      ny = -ny/nmag
      nz = -nz/nmag

      END SUBROUTINE modular_surface
