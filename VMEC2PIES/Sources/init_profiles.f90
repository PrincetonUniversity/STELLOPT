!-----------------------------------------------------------------------
!     Subroutine:    init_profiles
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/15/2011
!     Description:   This subroutine calculates the total current
!                    density from the pressure and toroidal current.
!                    We define the profiles as a function of the inital
!                    background coordinates so we must normalize our
!                    knots to the value of phiedge we wish to define
!                    the coords.  However, when we evaluate the splines
!                    later we define our PIES normalized flux in terms
!                    of phiedge.  We then throw out values of the
!                    normalized flux greater than 1.0.
!-----------------------------------------------------------------------
      SUBROUTINE init_profiles
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_background, ONLY: k, rho
      USE pies_runtime
      USE pies_profile
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          bcs1        Spline Boundary Conditions (not-a-knot)
!          ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, bcs1(2)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/0,0/) 
      CALL EZspline_init(p_spl,k+1,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/init_profiles',ier)
      CALL EZspline_init(ip_spl,k+1,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/init_profiles',ier)
      p_spl%isHermite = 1
      ip_spl%isHermite = 1
!      p_spl%x1 = torflux/MAXVAL(torflux)
!      ip_spl%x1 = torflux/MAXVAL(torflux)
      CALL EZspline_setup(p_spl,press,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
      CALL EZspline_setup(ip_spl,iprime,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE init_profiles
