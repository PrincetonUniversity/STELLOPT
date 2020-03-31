
SUBROUTINE neo_bderiv(theta, phi, f, g, dfdx, dfdy, dgdx, dgdy)
!
! Calculates first and second derivatives of b using the 2d-splines
!
! Input:  theta, phi
! Output: f      db/dt              (t = theta)
!         g      db/dp              (p = phi)
!         dfdx   d^2b/dt^2
!         dfdy   d^2b/(dt dp)
!         dgdx   d^2b/(dt dp)
!         dgdy   d^2b/dp^2
!
! Input/output consistent for neo_zeros2d
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp), INTENT(in)    ::   theta, phi
  REAL(kind=dp), INTENT(out)   ::   f, g, dfdx, dfdy, dgdx, dgdy

  REAL(kind=dp), DIMENSION(2)  ::   fderiv
  REAL(kind=dp), DIMENSION(3)  ::   sderiv
! **********************************************************************
! Evaluation of pointer
! **********************************************************************
  CALL poi2d(theta_int,phi_int,mt,mp,                                  &
             theta_start,theta_end,phi_start,phi_end,                  &
             theta,phi,theta_ind,phi_ind,theta_d,phi_d,ierr)
! **********************************************************************
! Evaluation of 2d-splines (first and second derivatives)
! **********************************************************************
!  print *, 'before eva2d_fd'
  CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,fderiv)
!  print *, 'before eva2d_sd'
  CALL eva2d_sd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,sderiv)
!  print *, 'after eva2d_sd'
! **********************************************************************
! Outputvalues (for neo_zeros2d)
! **********************************************************************
  f    = fderiv(1)
  g    = fderiv(2)
  dfdx = sderiv(1)
  dfdy = sderiv(2)
  dgdx = sderiv(2)
  dgdy = sderiv(3)
END SUBROUTINE neo_bderiv
