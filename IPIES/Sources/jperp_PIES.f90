!-----------------------------------------------------------------------
!     Subroutine:    jperp_PIES
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/13/2012
!     Description:   This subroutine calculates the perpendicular
!                    current density using the old PIES methodology.
!-----------------------------------------------------------------------
      SUBROUTINE jperp_PIES
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error flag
!          ik           Radial dummy index
!          dpds         Pressure gradient dp/ds
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, ik, mn
      REAL(rprec) :: pi2, mu0, I_norm
      REAL(rprec) :: dpds(0:k), denom(0:k)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      mu0 = 2*pi2*1e-7
      DO ik = 0, k
         CALL EZspline_interp(p_spl,rho(ik),press(ik),ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_interp/jperp',ier)
         CALL EZspline_interp(ip_spl,rho(ik),iprime(ik),ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_interp/jperp',ier)
         CALL EZspline_derivative(p_spl,1,rho(ik),dpds(ik),ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_derivative/jperp',ier)
      END DO
      WHERE(press == 0.0) dpds = 0.0
      ! Now we calculate jperp
      DO mn = 1, mnmax
         jsmns(mn,:) = 0.0
         IF ((xm(mn) /= 0) .and. (xn(mn) /=0)) THEN
            denom(:) = 1._rprec/(nfp*xn(mn)-iota(:)*xm(mn))
            jumnc(mn,:) = nfp*xn(mn)*dpds(:)*denom(:)
            jvmnc(mn,:) = xm(mn)*dpds(:)*denom(:)
         ELSE
            jumnc(mn,:) = iota(:) * iprime(:) + dpdpsi(:)
            jvmnc(mn,:) = iprime(:)
         END IF
      END DO
      PRINT *,jumnc(:,k)
      PRINT *,jvmnc(:,k)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jperp_PIES
