      SUBROUTINE eval_modular_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      USE boundary, ONLY: nfp, rmnc_b, zmns_b, xm_b, xn_b, mnmax
      USE Vwire, ONLY: nwire, nwire1
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nc, icoil, i, j, k, n, itheta
      REAL(rprec) ::
     1   x(nwdim), y(nwdim), z(nwdim),
     2   di0, r, phi0, phi2, phi1,
     3   phi_tot, zw, rw, dtheta, theta, theta1, theta2,
     4   dssq, dr_ext
!-----------------------------------------------

!        LSYMM      NODD       modular Configuration
!        _____      ____       _____________________
!          F          0        coils on both v=0 and v=1/2 symmetry planes
!          T          0        no coils on either v=0 or v=1/2 planes
!          T          1        coils on v=0 planes
!          F          1        coils on v=1/2 planes

      nc = nmod_coils / nfp                  !coils per field period
      IF (nmod_coils .le. 0) RETURN

!     offset in secular part of phi
      IF (lsymm) THEN
!        symmetry coil at phi = 0 for nc odd
         di0 = (nc + 1)*.5_dp
      ELSE
         IF (nodd .eq. 1) THEN
!           symmetry coil at phi = pi/3 for nc odd
            di0 = 0.5_dp
         ELSE
!           symmetry coils at phi = 0 and phi = pi/3 for nc even
            di0 = 1
         END IF
      END IF

      DO n = 1, nmod_coils
        phi_full(n) = (twopi*REAL(n-di0,rprec))/nmod_coils
      END DO

!     Set coil segment intervals in poloidal angle

      dtheta = twopi/nwire

!     Zero the phi segment points in the modular coil structure

      DO n = 1,nmid
        DO i = 1,nwire1
          modular(n)%phi(i) = zero
        END DO
      END DO

!     First consider the CASE with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp. This implies that the number of coils
!     per field period must be even (nodd = 0).

      IF (nodd.EQ.0 .AND. .NOT.lsymm) THEN

!     For each modular coil, compute phi at nwire1 points equally spaced in theta
!     NOTE: theta=s in Eq.2 of Strickler, et al

!     Symmetry coil at phi = 0.
        n = 1
        DO itheta = 1, nwire1
          theta = dtheta*(itheta-1)
          IF (nf_rho .GT. 0) CALL theta_f (n, theta, theta1, theta2)
          DO k = 1,nf_phi
             modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                              + (modular(n)%phis(k)*SIN(k*theta))
          END DO
        END DO

!     Coils 2 through nmid - 1
        DO n = 2, nmid-1
          DO itheta = 1, nwire1
            theta = dtheta*(itheta-1)
            IF (nf_rho .GT. 0) CALL theta_f (n, theta, theta1, theta2)
            modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                             + modular(n)%phic(0)
            DO k = 1,nf_phi
               modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                                +(modular(n)%phic(k)*COS(k*theta)
     2                                + modular(n)%phis(k)*SIN(k*theta))
            END DO
          END DO
        END DO

!     Symmetry coil at phi = pi/nfp.
        n = nmid
        DO itheta = 1, nwire1
          theta = dtheta*(itheta-1)
          IF (nf_rho .GT. 0) CALL theta_f (n, theta, theta1, theta2)
          DO k = 1,nf_phi
             modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                              + (modular(n)%phis(k)*SIN(k*theta))
          END DO
        END DO

!     Compute remaining phi for the field period using stellarator symmetry

        DO i = 1, nmid
           curcon(i) = modular(i)%current
        END DO

        DO i = 1, nmid-2
        curcon(nc+1-i) = modular(i+1)%current
        DO itheta = 1,nwire+1
        modular(nc+1-i)%phi(nwire1-itheta+1) =-modular(i+1)%phi(itheta)
        END DO
        END DO

!     Compute R and Z for first field period

        DO n = 1, nc
          phi1 = modular(n)%phi(1)
          phi2 = phi1
          DO itheta = 1,nwire
             theta = dtheta*(itheta-1)
             phi0 = modular(n)%phi(itheta)
             phi_tot = phi0 + phi_full(n)
             CALL rz_surf(theta,phi_tot,rw,zw,numsurf,
     1          rmn_sf,zmn_sf,m_num,n_num,nfper)
!            Extend symmetry and adjacent coils for NCSX NBI access
             IF (lncsx) THEN
               CALL radial_ext (n, theta, dr_ext)
               rw = rw + dr_ext
             END IF
             rcoil(n,itheta) = rw
             zcoil(n,itheta) = zw
             rcoil(n,itheta) = ABS(rcoil(n,itheta))
             phi0 = modular(n)%phi(itheta)
             phi1 = MIN(phi1,phi0)
             phi2 = MAX(phi2,phi0)
          END DO
          phimin(n) = phi1
          phimax(n) = phi2
        END DO

!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), or on neither. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).

      ELSE

!     For each modular, compute phi at nwire1 points equally spaced in theta

      DO n = 1,nmid-nodd
        DO itheta = 1,nwire1
          theta = dtheta*(itheta-1)
          IF (nf_rho .GT. 0) CALL theta_f (n, theta, theta1, theta2)
          modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                           + modular(n)%phic(0)
          DO k = 1,nf_phi
             modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                              +(modular(n)%phic(k)*COS(k*theta)
     2                              + modular(n)%phis(k)*SIN(k*theta))
          END DO
        END DO
      END DO

      IF (nodd .EQ. 1) THEN               !central coil for odd no coils
         n = nmid
         DO itheta = 1,nwire1
           theta = dtheta*(itheta-1)
           IF (nf_rho .gt. 0) CALL theta_f (n, theta, theta1, theta2)
           DO k = 1,nf_phi
              modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                               + (modular(n)%phis(k)*SIN(k*theta))
           END DO
         END DO
      END IF

!     Compute remaining phi for the field period using stellarator symmetry

      DO i = 1, nmid
         curcon(i) = modular(i)%current
      END DO
      DO i = 1,nmid-nodd
         curcon(nc+1-i) = modular(i)%current
         DO itheta = 1,nwire+1
           modular(nc+1-i)%phi(nwire1-itheta+1) =-modular(i)%phi(itheta)
         END DO
      END DO

!     Compute R and Z for first field period

      DO n = 1, nc
        phi1 = modular(n)%phi(1)
        phi2 = phi1
        DO itheta = 1,nwire
           theta = dtheta*(itheta-1)
           phi0 = modular(n)%phi(itheta)
           phi_tot = phi0 + phi_full(n)
           CALL rz_surf(theta,phi_tot,rw,zw,numsurf,
     1        rmn_sf,zmn_sf,m_num,n_num,nfper)
!          ExtEND symmetry and adjacent coils for NCSX NBI access
           IF (lncsx) THEN
              CALL radial_ext (n, theta, dr_ext)
              rw = rw + dr_ext
           END IF

           rcoil(n,itheta) = rw
           zcoil(n,itheta) = zw
           rcoil(n,itheta) = ABS(rcoil(n,itheta))
           phi0 = modular(n)%phi(itheta)
           phi1 = MIN(phi1,phi0)
           phi2 = MAX(phi2,phi0)
        END DO
        phimin(n) = phi1
        phimax(n) = phi2
      END DO

      ENDIF   ! END IF ((nodd .eq. 0) .and. (lsymm .eqv. .false.))

!     Store currents for ALL field periods

      DO j = 1,nfp-1
         DO i = 1,nc
           curcon(i + j*nc) = curcon(i)
         END DO
      END DO

!     Store x,y,z coordinates for coils as x_mod(i,j,n), ..., where
!     i = segment number, j = filament number (=1 for central filament),
!     and n = coil number

      DO 100 n = 1, nmod_coils
         icoil = 1 + MOD(n-1,nc)              !coil INDEX, MOD nc
         ymin_cls(n) = HUGE(ymin_cls(n))
         DO i = 1,nwire
           r = rcoil(icoil,i)
           z(i) = zcoil(icoil,i)
           x(i) = r*COS(modular(icoil)%phi(i) + phi_full(n))
           y(i) = r*SIN(modular(icoil)%phi(i) + phi_full(n))
           x_mod(i,1,n) = x(i)
           y_mod(i,1,n) = y(i)
           z_mod(i,1,n) = z(i)
           ymin_cls(n) = MIN(ymin_cls(n),ABS(y(i)))
         END DO

         x_mod(nwire1,1,n) = x_mod(1,1,n)
         y_mod(nwire1,1,n) = y_mod(1,1,n)
         z_mod(nwire1,1,n) = z_mod(1,1,n)

!     Compute length of coil n

         mod_length(n) = zero
         DO i = 1,nwire-1
           dssq = (x(i+1) - x(i))**2 + (y(i+1) - y(i))**2
     1          + (z(i+1) - z(i))**2
           mod_length(n) = mod_length(n) + SQRT(dssq)
         END DO
         dssq = (x(1) - x(nwire))**2 + (y(1) - y(nwire))**2
     1        + (z(1) - z(nwire))**2
         mod_length(n) = mod_length(n) + SQRT(dssq)

 100  CONTINUE

      RETURN
!CHECK THAT WINDING SURFACE IS CORRECT
      phi_tot=twopi/8
      xn_b = -xn_b/nfper
      DO i = 1, nwire
         theta=dtheta*(i-1)      
         CALL rz_surf(theta,phi_tot,rw,zw,numsurf,
     1                rmn_sf,zmn_sf,m_num,n_num,nfper)
         WRITE(33, 200) rw, zw
         CALL rz_surf(theta,phi_tot,rw,zw,mnmax,
     1                rmnc_b,zmns_b,int(xm_b),int(xn_b),nfper)
         WRITE(36, 200) rw, zw
      END DO

      xn_b = -xn_b*nfper
 200  FORMAT(1p,2e12.4)

      END SUBROUTINE eval_modular_coils
