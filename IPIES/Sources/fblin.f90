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
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (xsi,eta)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(2*last_line+2), qdot(2*last_line+2)
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier, ik
      REAL(rprec) :: eta_temp, xsi_temp, pi2
      REAL(rprec) :: rho_local(0:last_line), theta_local(0:last_line), phi_array(0:last_line)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      theta_local = 0.0
      DO ik = 0, last_line
         xsi_temp = q(2*ik+1)
         eta_temp = q(2*ik+2)
         rho_local(ik)   = sqrt(xsi_temp*xsi_temp+eta_temp*eta_temp)
         IF (rho_local(ik) > 0.0) theta_local(ik) = atan2(eta_temp,xsi_temp)
      END DO
      WHERE (rho_local > 1.0) rho_local = 0.99
      ! Initialize values of phi
      phi_array(0:last_line) = MOD(phi,pi2)
      ! Spline Bxsi
      CALL EZspline_interp(bxsi_spl,last_line+1,rho_local,theta_local,phi_array(0:last_line),bxsi_array(0:last_line),ier)
      ! Spline Beta
      CALL EZspline_interp(beta_spl,last_line+1,rho_local,theta_local,phi_array(0:last_line),beta_array(0:last_line),ier)
      ! Handle Wall hits
      WHERE(rho_local > 1.0) bxsi_array = 0.0
      WHERE(rho_local > 1.0) beta_array = 0.0
      ! Output Qdot
      DO ik = 0, last_line
         qdot(2*ik+1) = bxsi_array(ik)
         qdot(2*ik+2) = beta_array(ik)
      END DO
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin
