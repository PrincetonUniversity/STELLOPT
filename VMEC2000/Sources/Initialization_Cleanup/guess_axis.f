      SUBROUTINE guess_axis(r1, z1, ru0, zu0)
      USE vmec_main
      USE vmec_params, ONLY: nscale, signgs
      USE realspace, ONLY: sqrts
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nzeta,ntheta3,0:1),
     1     INTENT(in) :: r1, z1
      REAL(rprec), DIMENSION(ns,nzeta,ntheta3), INTENT(in) :: ru0, zu0
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: limpts = 61
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
      REAL(rprec), DIMENSION(nzeta) :: rcom, zcom
      REAL(rprec), DIMENSION(ntheta1) :: r1b, z1b, rub, zub
      REAL(rprec), DIMENSION(ntheta1) :: r12, z12
      REAL(rprec), DIMENSION(ntheta1) :: rs, zs, tau, ru12, zu12, tau0
      REAL(rprec) :: rlim, zlim
      REAL(rprec) :: rmax, rmin, zmax, zmin, dzeta
      REAL(rprec) :: ds, mintau, mintemp
C-----------------------------------------------
!
!     COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
!     LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
!     SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
!     YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
!     CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED
!
      ns12 = (ns+1)/2

      planes: DO iv = 1, nzeta
         IF (.not.lasym .and. iv.gt.nzeta/2+1) THEN
            rcom(iv) = rcom(nzeta+2-iv)
            zcom(iv) =-zcom(nzeta+2-iv)
            CYCLE
         END IF
         r1b(:ntheta3) = r1(ns,iv,:,0) + r1(ns,iv,:,1)
         z1b(:ntheta3) = z1(ns,iv,:,0) + z1(ns,iv,:,1)
         r12(:ntheta3) = r1(ns12,iv,:,0) + r1(ns12,iv,:,1)*sqrts(ns12)
         z12(:ntheta3) = z1(ns12,iv,:,0) + z1(ns12,iv,:,1)*sqrts(ns12)
         rub(:ntheta3) = ru0(ns,iv,:)
         zub(:ntheta3) = zu0(ns,iv,:)
         ru12(:ntheta3) =  p5*(ru0(ns,iv,:) + ru0(ns12,iv,:))
         zu12(:ntheta3) =  p5*(zu0(ns,iv,:) + zu0(ns12,iv,:))

         IF (.not.lasym) THEN
!
!     USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
!     TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)
!
         ivminus = MOD(nzeta + 1 - iv,nzeta) + 1           !!(twopi-v)
         DO iu = 1+ntheta2, ntheta1
            iu_r = ntheta1 + 2 - iu
            r1b(iu) = r1(ns,ivminus,iu_r,0) + r1(ns,ivminus,iu_r,1)
            z1b(iu) =-(z1(ns,ivminus,iu_r,0) + z1(ns,ivminus,iu_r,1))
            r12(iu) = r1(ns12,ivminus,iu_r,0) +
     1                r1(ns12,ivminus,iu_r,1)*sqrts(ns12)
            z12(iu) =-(z1(ns12,ivminus,iu_r,0) +
     1                z1(ns12,ivminus,iu_r,1)*sqrts(ns12))
            rub(iu) =-ru0(ns,ivminus,iu_r)
            zub(iu) = zu0(ns,ivminus,iu_r)
            ru12(iu)=-p5*(ru0(ns,ivminus,iu_r) + ru0(ns12,ivminus,iu_r))
            zu12(iu)= p5*(zu0(ns,ivminus,iu_r) + zu0(ns12,ivminus,iu_r))
         END DO

         END IF
!
!        Scan over r-z grid for interior point
!
         rmin = MINVAL(r1b);  rmax = MAXVAL(r1b)
         zmin = MINVAL(z1b);  zmax = MAXVAL(z1b)
         rcom(iv) = (rmax + rmin)/2; zcom(iv) = (zmax + zmin)/2

!
!        Estimate jacobian based on boundary and 1/2 surface
!
         ds = (ns - ns12)*hs
         DO iu = 1, ntheta1
            rs(iu) = (r1b(iu) - r12(iu))/ds + r1(1,iv,1,0)
            zs(iu) = (z1b(iu) - z12(iu))/ds + z1(1,iv,1,0)
            tau0(iu) = ru12(iu)*zs(iu) - zu12(iu)*rs(iu)
         END DO

         mintau = 0

         DO nlim = 1, limpts
            zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)
            IF (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) THEN
               zlim = 0
               IF (nlim .gt. 1) EXIT
            END IF
!
!           Find value of magnetic axis that maximizes the minimum jacobian value
!
            DO klim = 1, limpts
               rlim = rmin + ((rmax - rmin)*(klim-1))/(limpts-1)
               tau = signgs*(tau0 - ru12(:)*zlim + zu12(:)*rlim)
               mintemp = MINVAL(tau)
               IF (mintemp .gt. mintau) THEN
                  mintau = mintemp
                  rcom(iv) = rlim
                  zcom(iv) = zlim
!           If up-down symmetric and lasym=T, need this to pick z = 0
               ELSE IF (mintemp .eq. mintau) THEN
                  IF (ABS(zcom(iv)).gt.ABS(zlim)) zcom(iv) = zlim
               END IF
            END DO
         END DO

      END DO planes

!
!     FOURIER TRANSFORM RCOM, ZCOM
!
      dzeta = two/nzeta
      DO n = 0, ntor
         raxis_cc(n) = dzeta*SUM(cosnv(:,n)*rcom(:))/nscale(n)
         zaxis_cs(n) =-dzeta*SUM(sinnv(:,n)*zcom(:))/nscale(n)
         raxis_cs(n) =-dzeta*SUM(sinnv(:,n)*rcom(:))/nscale(n)
         zaxis_cc(n) = dzeta*SUM(cosnv(:,n)*zcom(:))/nscale(n)
         IF (n.eq.0 .or. n.eq.nzeta/2) THEN
            raxis_cc(n) = p5*raxis_cc(n)
            zaxis_cc(n) = p5*zaxis_cc(n)
         END IF
      END DO

!  100 FORMAT(' n = ',i4,' raxis = ',1pe10.3,' zaxis = ',1pe10.3)

      END SUBROUTINE guess_axis
