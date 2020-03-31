
SUBROUTINE neo_eval(theta,phi,bval,gval,kval,pval,qval)
! Evaluation of Splines
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp), INTENT(in)  ::   theta, phi
  REAL(kind=dp), INTENT(out) ::   bval, gval, kval, pval, qval
! **********************************************************************
! Evaluation of pointer
! **********************************************************************
  CALL poi2d(theta_int,phi_int,mt,mp,                                  &
             theta_start,theta_end,phi_start,phi_end,                  &
             theta,phi,theta_ind,phi_ind,theta_d,phi_d,ierr)
! **********************************************************************
! Evaluation of 2d-splines
! **********************************************************************
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             b_spl,bval)
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             g_spl,gval)
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             k_spl,kval)
  CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,            &
             p_spl,pval)
  IF (calc_cur .EQ. 1) THEN
     CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
                q_spl,qval)
  END IF
  RETURN
END SUBROUTINE neo_eval
