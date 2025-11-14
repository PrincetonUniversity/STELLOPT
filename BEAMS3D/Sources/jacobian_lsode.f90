!-----------------------------------------------------------------------
!     Function:      jacobian_lsode
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/16/2012
!     Description:   This subroutine calculates the Jacobian of the ODE for 
!                    field line following (LSODE).  Currently just a
!                    dummy routine.
!-----------------------------------------------------------------------
      SUBROUTINE jacobian_lsode(neq,t,q,ml,mp,pd,nrpd)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_grid
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
      DOUBLE PRECISION :: t, q(2), pd(2,2)
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, br_temp, bz_temp, br_dot, bz_dot
      REAL(rprec) :: x,y,z,vx, vy, vz
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      pd = 0.0
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jacobian_lsode
