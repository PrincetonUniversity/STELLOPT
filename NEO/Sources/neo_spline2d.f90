
SUBROUTINE neo_spline2d
! Creation of Spline Arrays
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
! **********************************************************************
! Allocation of spline arrays is done in neo_prep
! **********************************************************************
! **********************************************************************
! Double periodic splines (parameter mt=1 and mp=1)
! **********************************************************************
! Spline for mod b
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,b,b_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for sqrg11
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,sqrg11,g_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for geodesic curviture
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,kg,k_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for parallel derivative
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,pard,p_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for quasi-toroidal phi component of b
  IF (calc_cur .EQ. 1) THEN
!    IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
     CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,bqtphi,q_spl)
!    IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
  END IF
! **********************************************************************
! Spline test
! **********************************************************************
  IF (spline_test .GT. 0) THEN
!     IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_spline_test'
     CALL neo_spline_test
!     IF (write_progress .NE. 0) WRITE (w_us,*) 'after neo_spline_test'
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_spline2d
