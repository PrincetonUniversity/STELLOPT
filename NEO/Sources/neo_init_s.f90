
SUBROUTINE neo_init_s(psi,dpsi)
! Initialization for Specific Magnetic Surface
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

  REAL(kind=dp),                INTENT(out)     :: psi, dpsi
  INTEGER,       DIMENSION(2)                   :: b_minpos, b_maxpos

  REAL(kind=dp), PARAMETER                      :: eps_newt = 1.0e-10_dp
  INTEGER                                       :: iter, error
  INTEGER,       PARAMETER                      :: iterma_newt = 100
  REAL(kind=dp)                                 :: gval_bmin
  REAL(kind=dp)                                 :: kval_bmin,pval_bmin
  REAL(kind=dp)                                 :: gval_bmax
  REAL(kind=dp)                                 :: kval_bmax,pval_bmax
  REAL(kind=dp)                                 :: qval_bmin,qval_bmax

  REAL(kind=dp)  :: tht, pht, f, g, dfdx, dfdy, dgdx, dgdy
  REAL(kind=dp)  :: thi, phi
  REAL(kind=dp)  :: iot,bval,gval,kval,pval
  INTEGER        :: i
! **********************************************************************
! Calculate Fourier sums and derived quantities
! **********************************************************************
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_fourier'
  CALL neo_fourier
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_fourier'
! **********************************************************************
! Calculation of dpsi
! **********************************************************************
  psi = es(psi_ind)
  IF (psi_ind .EQ. 1) THEN
     dpsi = 0
  ELSE
     dpsi = psi - es(psi_ind-1)
  ENDIF
! **********************************************************************
! Initilaze spline arrays
! **********************************************************************
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_spline2d'
  CALL neo_spline2d
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_spline2d'
! **********************************************************************
! Calculate absolute minimum and maximum of b and its location (theta, phi)
! **********************************************************************
  b_minpos   = MINLOC(b)
  b_min      = b(b_minpos(1),b_minpos(2))
  theta_bmin = theta_arr(b_minpos(1))
  phi_bmin   = phi_arr(b_minpos(2))

  b_maxpos   = MAXLOC(b)
  b_max      = b(b_maxpos(1),b_maxpos(2))
  theta_bmax = theta_arr(b_maxpos(1))
  phi_bmax   = phi_arr(b_maxpos(2))

! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_zeros2d (min)'
  CALL neo_zeros2d(theta_bmin, phi_bmin, eps_newt, iterma_newt, iter, error)
  CALL neo_eval(theta_bmin,phi_bmin,b_min,gval_bmin,kval_bmin,&
                pval_bmin,qval_bmin)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_zeros2d'

! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_zeros2d (max)'
  CALL neo_zeros2d(theta_bmax, phi_bmax, eps_newt, iterma_newt, iter, error)
  CALL neo_eval(theta_bmax,phi_bmax,b_max,gval_bmax,kval_bmax,&
                pval_bmax,qval_bmax)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_zeros2d'
! **********************************************************************
! ATTENTION
! Set bmref to the absolute maximum of b on flux surface
! This is absolutely necessary for internal routines
! Rescaling is done at the end of the main program neo.f90
! **********************************************************************
  bmref = b_max
  RETURN
END SUBROUTINE neo_init_s
