!-----------------------------------------------------------------------
!     Function:      fblin_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following (NAG)
!-----------------------------------------------------------------------
      SUBROUTINE fblin_nag(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(2), qdot(2)
!-----------------------------------------------------------------------
!     Local Variables
!          ier        Error flag
!          r_temp     R
!          z_temp     Z
!          phi_temp   phi (radians)
!          br_temp    br/bphi evaluated from splines
!          bz_temp    bz/bphi evaluated from splines
!          br_dot     Random force for diffusion
!-----------------------------------------------------------------------
      INTEGER :: ier
      DOUBLE PRECISION :: r_temp, phi_temp, z_temp, br_temp, bz_temp
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      r_temp = q(1)
      z_temp = q(2)
      phi_temp = MOD(phi,delta_phi)
      IF (phi_temp < 0) phi_temp = delta_phi + phi_temp
      br_temp = 0.0; bz_temp = 0.0
      CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
      IF (ier == 0) THEN
         CALL EZspline_interp(BR_spl,r_temp,phi_temp,z_temp,br_temp,ier)
         CALL EZspline_interp(BZ_spl,r_temp,phi_temp,z_temp,bz_temp,ier)
      END IF
      qdot(1) = br_temp
      qdot(2) = bz_temp
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin_nag
