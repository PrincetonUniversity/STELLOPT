!-----------------------------------------------------------------------
!     Function:      jacobian_lsode
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/16/2012
!     Description:   This subroutine calculates the Jacobian of the ODE for 
!                    field line following (LSODE).  Currently just a
!                    dummy routine.
!-----------------------------------------------------------------------
      SUBROUTINE jacobian_lsode(neq,phi,q,ml,mp,pd,nrpd)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid, ONLY: BR_spl, BZ_spl, delta_phi
      USE EZspline_obj
      USE EZspline                                            ! MPI
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: neq, nrpd, ml, mp
      DOUBLE PRECISION :: phi, q(2), pd(2,2)
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, br_temp, bz_temp, br_dot, bz_dot
      REAL(rprec) :: x,y,z,vx, vy, vz
      DOUBLE PRECISION :: R_grad(3), Z_grad(3)
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
         CALL EZspline_gradient(BR_spl,r_temp,phi_temp,z_temp,R_grad,ier)
         CALL EZspline_gradient(BZ_spl,r_temp,phi_temp,z_temp,Z_grad,ier)
         pd(1,1) = R_grad(1) !dBR/dR
         pd(1,2) = R_grad(3) ! dBR/dZ
         pd(2,1) = Z_grad(1) ! dBZ/dR
         pd(2,2) = Z_grad(3) ! dBZ/dZ
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jacobian_lsode
