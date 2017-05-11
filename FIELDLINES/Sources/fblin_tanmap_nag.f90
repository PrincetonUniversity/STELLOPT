!-----------------------------------------------------------------------
!     Function:      fblin_tanmap_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following (NAG) with tangent map
!                    integration.
!-----------------------------------------------------------------------
      SUBROUTINE fblin_tanmap_nag(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE fieldlines_runtime, ONLY: lmu, mu, ladvanced
      USE sheppack
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2),q(3),q(4),q(5),q(6)) = (R,Z,T11,T12,T21,T22)
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
!          br_dot     Random force for diffusion
!-----------------------------------------------------------------------
      INTEGER :: ier
      DOUBLE PRECISION :: r_temp, phi_temp, z_temp, br_temp, bz_temp
      DOUBLE PRECISION :: drdr, drdz, dzdr,dzdz
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      r_temp = q(1)
      z_temp = q(2)
      phi_temp = MOD(phi,delta_phi)
      IF (phi_temp < 0) phi_temp = delta_phi + phi_temp
      CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
      IF (ier == 0) THEN
         ! Trajectory
         CALL EZspline_interp(BR_spl,r_temp,phi_temp,z_temp,br_temp,ier)
         CALL EZspline_interp(BZ_spl,r_temp,phi_temp,z_temp,bz_temp,ier)
         ! Map
         CALL EZspline_derivative(BR_spl,1,0,0,r_temp,phi_temp,z_temp,drdr,ier)
         CALL EZspline_derivative(BR_spl,0,0,1,r_temp,phi_temp,z_temp,drdz,ier)
         CALL EZspline_derivative(BZ_spl,1,0,0,r_temp,phi_temp,z_temp,dzdr,ier)
         CALL EZspline_derivative(BZ_spl,0,0,1,r_temp,phi_temp,z_temp,dzdz,ier)
      ELSE
         br_temp = 0.0
         bz_temp = 0.0
         drdr    = 0.0
         drdz    = 0.0
         dzdr    = 0.0
         dzdz    = 0.0
      END IF
      qdot(1) = br_temp
      qdot(2) = bz_temp
      qdot(3) = drdr*q(3) + drdz*q(5)
      qdot(4) = drdr*q(4) + drdz*q(6)
      qdot(5) = dzdr*q(3) + dzdz*q(5)
      qdot(6) = dzdr*q(4) + dzdz*q(6)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin_tanmap_nag
