!-----------------------------------------------------------------
!     Function:      fieldlines_BCART
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          01/07/2026
!     Description:   Returns Bx, By, Bz
!-----------------------------------------------------------------
SUBROUTINE fieldlines_BCART(x,y,z,Bx,By,Bz)
   !--------------------------------------------------------------
   !     Input Parameters
   !         r, phi, z     Cylindrical coordiantes
   !     Output Parameters
   !         Br, Bphi, Bz  Magnetic field components
   !--------------------------------------------------------------
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: x, y, z
   DOUBLE PRECISION, INTENT(OUT) :: Bx, By, Bz

   !--------------------------------------------------------------
   !     Local Variables
   !        phi_temp     Helpers (r,phi,z)
   !--------------------------------------------------------------
   DOUBLE PRECISION :: r_temp, phi_temp, Br_temp, Bphi_temp

   !--------------------------------------------------------------
   !     Begin Subroutine
   !--------------------------------------------------------------

   ! Cartesian coordiantes to cylindrical
   r_temp = sqrt(x*x+y*y)
   phi_temp = atan2(y,x)

   ! Call cylindrical routine
   CALL fieldlines_BCYL(r_temp,phi_temp,z,Br_temp,Bphi_temp,Bz)

   ! Cyl vectors to cartesian
   Bx = Br_temp*cos(phi_temp)-Bphi_temp*sin(phi_temp)
   By = Br_temp*sin(phi_temp)+Bphi_temp*cos(phi_temp)

   RETURN

END SUBROUTINE fieldlines_BCART