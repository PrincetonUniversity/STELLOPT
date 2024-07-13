      SUBROUTINE write_surfaces
      USE stel_constants
      USE boundary
      USE modular_coils
      USE saddle_coils
      USE saddle_surface
      USE Vwire, ONLY: nwire
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: nphi = 9
      INTEGER :: i, k, m, n, mn, itheta, iphi, nscmax
      INTEGER :: npt
      INTEGER, DIMENSION(ncdim) :: npts
      REAL(rprec) :: v0, phi0, theta, dtheta, dphi, rws, zws
      REAL(rprec) :: ck, sk, cosmu, sinmu, cosnv, sinnv, rps, zps
      REAL(rprec), DIMENSION (nwdim1,nphi) :: rw, zw, rs, zs, rp, zp
      REAL(rprec), DIMENSION (nwdim1,2) :: xw, yw, xp, yp
      REAL(rprec), DIMENSION (ncdim, 5*ncdim) :: rpts, zpts
      REAL(rprec), DIMENSION (5*ncdim, nphi) :: rplt, zplt
!-----------------------------------------------
      OPEN(unit=22,file='surfaces.dat',status='unknown')
      OPEN(unit=25,file='topview.dat',status='unknown')
      OPEN(unit=26,file='rzpts.dat',status='unknown')

      nscMAX = ncdim/2
      dtheta = twopi/(nwire-1)
      dphi = twopi/(16*nfp)
      theta = zero

      DO n = 1, nphi
        phi0 = -pi/nfp + (n-1)*dphi

        IF (lsaddle) THEN
           v0 = nfp*phi0/twopi
           CALL eval_saddle_rz (v0, npts, rpts, zpts)
           m = 0
           DO k = 1, nsmid
              npt = npts(k)
              DO i = 1, nscmax
                 m = m + 1
                 IF (i .le. npt) THEN
                    rplt(m,n) = rpts(k,i)
                    zplt(m,n) = zpts(k,i)
                 ELSE
                    rplt(m,n) = 1000
                    zplt(m,n) = 1000
                 END IF
              END DO
           END DO
        END IF

        DO itheta = 1,nwire
           theta = dtheta*REAL(itheta-1)

!     modular winding surface at phi
           CALL rz_surf(theta,phi0,rws,zws,numsurf,
     1                   rmn_sf,zmn_sf,m_num,n_num,nfp)
           rw(itheta,n) = rws
           zw(itheta,n) = zws

!     Saddle winding surface edge at phi
           rws = 0.0
           zws = 0.0
           DO k = 1, numsurf_sad
              ck = COS(m_sad(k)*theta + nfp*n_sad(k)*phi0)
              sk = SIN(m_sad(k)*theta + nfp*n_sad(k)*phi0)
              rws = rws + rmn_sad(k)*ck
              zws = zws + zmn_sad(k)*sk
           END DO
           rs(itheta,n) = rws
           zs(itheta,n) = zws

!     Plasma surface at phi
           rps = 0.0
           zps = 0.0
           DO mn = 1, mnmax
              cosmu = COS(xm_b(mn)*theta)
              sinmu = SIN(xm_b(mn)*theta)
              cosnv = COS(xn_b(mn)*phi0)
              sinnv = SIN(xn_b(mn)*phi0)
              rps = rps + rmnc_b(mn)*(cosmu*cosnv + sinmu*sinnv)
              zps = zps + zmns_b(mn)*(sinmu*cosnv - cosmu*sinnv)
           END DO
           rp(itheta,n) = rps
           zp(itheta,n) = zps

        END DO

      END DO
      DO i = 1,nwire
        WRITE(22,1000)(rw(i,mn),zw(i,mn), mn=1,9)
      END DO
      WRITE(22,'(/)')
      DO i = 1,nwire
        WRITE(22,1000)(rs(i,mn),zs(i,mn), mn=1,9)
      END DO
      WRITE(22,'(/)')
      DO i = 1,nwire
        WRITE(22,1000)(rp(i,mn),zp(i,mn), mn=1,9)
      END DO
      IF (lsaddle) THEN
         DO m = 1,nsmid*nscmax
            WRITE(26,1000) (rplt(m,n), zplt(m,n), n=1,9)
         END DO
      END IF
 1000 FORMAT(1p,18e12.4)

!     Now WRITE data to plot topview

      dphi = twopi/(nwire-1)
      DO n = 1,2
        theta = 0.5*(n-1)*twopi
        DO iphi = 1,nwire
           phi0 = dphi*(iphi-1)

!     modular winding surface edge at phi
           CALL rz_surf(theta,phi0,rws,zws,numsurf,
     1                   rmn_sf,zmn_sf,m_num,n_num,nfp)
           xw(iphi,n) = rws*COS(phi0)
           yw(iphi,n) = rws*SIN(phi0)

!     Plasma surface edge at phi
           rps = 0.0
           zps = 0.0
           DO mn = 1, mnmax
              cosmu = COS(xm_b(mn)*theta)
              sinmu = SIN(xm_b(mn)*theta)
              cosnv = COS(xn_b(mn)*phi0)
              sinnv = SIN(xn_b(mn)*phi0)
              rps = rps + rmnc_b(mn)*(cosmu*cosnv + sinmu*sinnv)
              zps = zps + zmns_b(mn)*(sinmu*cosnv - cosmu*sinnv)
           END DO
           xp(iphi,n) = rps*COS(phi0)
           yp(iphi,n) = rps*SIN(phi0)

        END DO
      END DO

!     DO i = 1,nwire
!       WRITE(25,1000) xw(i,1), yw(i,1)
!     END DO
!     WRITE(25,'(/)')
      DO i = 1,nwire
        WRITE(25,1000) xp(i,1), yp(i,1)
      END DO
      WRITE(25,'(/)')
!     DO i = 1,nwire
!       WRITE(25,1000) xw(i,2), yw(i,2)
!     END DO
!     WRITE(25,'(/)')
      DO i = 1,nwire
        WRITE(25,1000) xp(i,2), yp(i,2)
      END DO

      CLOSE(22)
      CLOSE(25)
      CLOSE(26)

!     WRITE winding surface coefficients
      OPEN (unit=27,file='rz_coeff.dat',status='unknown')
      OPEN (unit=28,file='rz_coeff2.dat',status='unknown')
      WRITE (27,115) numsurf
      WRITE (27,120) (m_num(i), n_num(i), rmn_sf(i), zmn_sf(i),
     1                i=1, numsurf)
      WRITE (28,115) numsurf_sad
      WRITE (28,120) (m_sad(i), n_sad(i), rmn_sad(i), zmn_sad(i),
     1                i=1, numsurf_sad)
  115 FORMAT(i4)
  120 FORMAT(2i4,1p,2e16.8)
      CLOSE (27)
      CLOSE (28)

!     WRITE plasma surface coefficients
      OPEN (unit=31,file='rz_plasma.dat',status='unknown')
      OPEN (unit=32,file='rz_plasma2.dat',status='unknown')
      WRITE (31,135) mnmax
      WRITE (31,140) (xm_b(i), -xn_b(i)/nfp, rmnc_b(i), zmns_b(i),
     1                i=1, mnmax)
      WRITE (32,135) mnmax
      WRITE (32,145) (xm_b(i), -xn_b(i)/nfp, rmnc_b(i), zmns_b(i),
     1                i=1, mnmax)
  135 FORMAT(i4)
  140 FORMAT(2f8.0,2e22.12)
  145 FORMAT(2f10.0,2f20.12)
      CLOSE (31)
      CLOSE (32)

      END SUBROUTINE write_surfaces
