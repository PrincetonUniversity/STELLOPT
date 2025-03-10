      SUBROUTINE modular_curve (n, u, v, tx, ty, tz, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE modular_coils
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: k, n, ncl
      REAL(rprec) :: tx, ty, tz
      REAL(rprec) :: ck, sk, u, v,
     1   ph0, ph1, ph2, r0, r1, r2, z0, z1, z2,
     2   x0, x1, x2, y0, y1, y2, s1, c
      REAL(rprec) :: th0, th1 
      REAL(rprec) :: fth0, fth1, fth2
!-----------------------------------------------

!     compute theta and derivative wrt u. Note: u=s in Eq. 2 of Strickler, et al

      ncl = n
      th0 = twopi*u
      th1 = twopi

!     compute phi and derivatives wrt u

      IF (nf_rho .GT. 0) THEN
!     Case with possible windbacks - call theta_f to modify argument th0
!     put theta and derivatives in temporary variables fth0, fth1, fth2 so
!     the winding law calculation will not interfere with the winding
!     surface evaluation below
         fth0 = th0
         CALL theta_f (ncl, fth0, fth1, fth2)
         ph0 = phi_full(ncl)
         ph1 = 0
         ph2 = 0
         DO k = 0,nf_phi
            ck = COS(k*fth0)
            sk = SIN(k*fth0)
            ph0 = ph0 +  modular(ncl)%phic(k)*ck
     1                +  modular(ncl)%phis(k)*sk
            ph1 = ph1 + (modular(ncl)%phis(k)*ck
     1                -  modular(ncl)%phic(k)*sk)*k
            ph2 = ph2 + (modular(ncl)%phis(k)*sk
     1                +  modular(ncl)%phic(k)*ck)*k**2
         END DO
!     note the order of the next two lines is important
         ph2 = twopi**2*(ph2*fth1**2 + ph1*fth2)
         ph1 = twopi*fth1*ph1
      ELSE
!     Case without windbacks
         ph0 = phi_full(ncl)
         ph1 = 0
         ph2 = 0
         DO k = 0,nf_phi
            ck = COS(k*th0)
            sk = SIN(k*th0)
            ph0 = ph0 +  modular(ncl)%phic(k)*ck
     1                +  modular(ncl)%phis(k)*sk
            ph1 = ph1 + (modular(ncl)%phis(k)*ck
     1                -  modular(ncl)%phic(k)*sk)*k
            ph2 = ph2 + (modular(ncl)%phis(k)*sk
     1                +  modular(ncl)%phic(k)*ck)*k**2
         END DO
         ph1 = th1*ph1
         ph2 = -th1**2*ph2
      END IF

!     compute R, Z and derivatives wrt u 

      r0 = 0
      r1 = 0
      r2 = 0
      z0 = 0
      z1 = 0
      z2 = 0
      DO k = 1, numsurf
         ck = COS(m_num(k)*th0 + nfper*n_num(k)*ph0)
         sk = SIN(m_num(k)*th0 + nfper*n_num(k)*ph0)
!     R ...
         r0 = r0 + rmn_sf(k)*ck
         r1 = r1 - rmn_sf(k)*(m_num(k)*th1 + nfper*n_num(k)*ph1)*sk
         r2 = r2 - rmn_sf(k)*(nfper*n_num(k)*ph2*sk
     1      + (m_num(k)*th1 + nfper*n_num(k)*ph1)**2*ck)
!     Z ...
         z0 = z0 + zmn_sf(k)*sk
         z1 = z1 + zmn_sf(k)*(m_num(k)*th1 + nfper*n_num(k)*ph1)*ck
         z2 = z2 + zmn_sf(k)*(nfper*n_num(k)*ph2*ck
     1      - (m_num(k)*th1 + nfper*n_num(k)*ph1)**2*sk)
      END DO

!     compute X and derivatives wrt u

      x0 = r0*COS(ph0)
      x1 = r1*COS(ph0) - r0*ph1*SIN(ph0)
      x2 = r2*COS(ph0) - r1*ph1*SIN(ph0)
     1   - (r0*ph1**2*COS(ph0) + (r0*ph2 + r1*ph1)*SIN(ph0))

!     compute Y and derivatives wrt u

      y0 = r0*SIN(ph0)
      y1 = r1*SIN(ph0) + r0*ph1*COS(ph0)
      y2 = r2*SIN(ph0) + r1*ph1*COS(ph0)
     1   + (-r0*ph1**2*SIN(ph0) + (r0*ph2 + r1*ph1)*COS(ph0))

!     compute ds/du

      s1 = SQRT (x1*x1 + y1*y1 + z1*z1)

!     TANgent vector dx/ds, dy/ds, dz/ds

      tx = x1/s1
      ty = y1/s1
      tz = z1/s1

!     compute curvature

      c = SQRT ((y1*z2 - y2*z1)**2
     1        + (z1*x2 - z2*x1)**2
     2        + (x1*y2 - x2*y1)**2)/s1**3

!     RETURN toroidal v
      v = nfper*ph0/twopi

      END SUBROUTINE modular_curve
