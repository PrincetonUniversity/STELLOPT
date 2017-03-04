!-----------------------------------------------------------------------
!     Function:      fblin
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following
!-----------------------------------------------------------------------
      SUBROUTINE fblin(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_fieldlines
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (xsi,eta)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(2), qdot(2)
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier, ik
      REAL(rprec) :: beta_temp, bxsi_temp
      REAL(rprec) :: rho_local, theta_local, phi_local
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      rho_local = q(1)
      theta_local = q(2)
      phi_local = phi
      IF (rho_local > 1) rho_local = 1
      theta_local = MOD(theta_local,thmx)
      IF (theta_local < 0) theta_local = theta_local + thmx
      phi_local = MOD(phi,phmx)
      IF (phi_local < 0) phi_local = phi_local + phmx
      ! Spline Bxsi
      CALL EZspline_interp(bxsi_spl,rho_local,theta_local,phi_local,bxsi_temp,ier)
      ! Spline Beta
      CALL EZspline_interp(beta_spl,rho_local,theta_local,phi_local,beta_temp,ier)
      IF (rho_local == 1) THEN
         qdot(1) = 0
         qdot(2) = 1
      ELSE
         qdot(1) = bxsi_temp
         qdot(2) = beta_temp
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin
