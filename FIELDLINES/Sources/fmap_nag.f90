!-----------------------------------------------------------------------
!     Function:      fmap_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          09/01/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    return map calculation.
!-----------------------------------------------------------------------
      SUBROUTINE fmap_nag(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE fieldlines_runtime, ONLY: lmu, mu, ladvanced
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(6), qdot(6)
!-----------------------------------------------------------------------
!     Local Variables
!          ier        Error flag
!          r_temp     R
!          z_temp     Z
!          phi_temp   phi (radians)
!          br_temp    br/bphi evaluated from splines
!          bz_temp    bz/bphi evaluated from splines
!          br_r       dR_dot/dR for return map
!          br_z       dR_dot/dZ for return map
!          bz_r       dZ_dot/dR for return map
!          bz_z       dZ_dot/dZ for return map
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, br_temp, bz_temp
      REAL(rprec) :: br_r,br_z,bz_r,bz_z
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      r_temp = q(1)
      z_temp = q(2)
      phi_temp = MOD(phi,delta_phi)
      IF (phi_temp < 0) phi_temp = delta_phi + phi_temp
      IF (ladvanced) THEN
         ! Not implemented
      ELSE
         CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
         IF (ier == 0) THEN
            CALL EZspline_interp(BR_spl,r_temp,phi_temp,z_temp,br_temp,ier)
            CALL EZspline_interp(BZ_spl,r_temp,phi_temp,z_temp,bz_temp,ier)
            CALL EZspline_derivative(BR_spl,1,0,0,r_temp,phi_temp,z_temp,br_r,ier)
            CALL EZspline_derivative(BR_spl,0,0,1,r_temp,phi_temp,z_temp,br_z,ier)
            CALL EZspline_derivative(BZ_spl,1,0,0,r_temp,phi_temp,z_temp,bz_r,ier)
            CALL EZspline_derivative(BZ_spl,0,0,1,r_temp,phi_temp,z_temp,bz_z,ier)
         ELSE
            br_temp = 0.0
            bz_temp = 0.0
         END IF
      END IF
      qdot(1) = br_temp
      qdot(2) = bz_temp
      qdot(3) = br_r
      qdot(4) = br_z
      qdot(5) = bz_r
      qdot(6) = bz_z
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fmap_nag
