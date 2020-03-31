
MODULE neo_exchange
! Working parameters
  USE neo_precision

  REAL(kind=dp)              :: b_min, b_max
  INTEGER                    :: nper
  REAL(kind=dp)              :: rt0, rt0_g
  REAL(kind=dp)              :: bmref, bmref_g
  INTEGER                    :: nstep_per, nstep_min, nstep_max
  INTEGER                    :: write_integrate
  INTEGER                    :: write_diagnostic
  INTEGER                    :: write_cur_inte
  REAL(kind=dp)              :: acc_req
  INTEGER                    :: no_bins
  INTEGER                    :: psi_ind
  INTEGER                    :: calc_nstep_max
  REAL(kind=dp)              :: theta_bmin, phi_bmin
  REAL(kind=dp)              :: theta_bmax, phi_bmax
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: iota
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: curr_pol
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: curr_tor
  REAL(kind=dp)              :: fac
  INTEGER                    :: calc_cur
  INTEGER                    :: hit_rat, nfp_rat, nfl_rat
  REAL(kind=dp)              :: delta_theta_rat
END MODULE neo_exchange
