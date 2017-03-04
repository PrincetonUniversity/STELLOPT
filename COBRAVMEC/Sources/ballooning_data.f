      MODULE ballooning_data
      USE stel_kinds
      IMPLICIT NONE
      INTEGER:: np0
      INTEGER :: tsymm = 1                                                                   ! By default, assume no stellarator symmetry
      INTEGER :: np0_in = 201, kth, k_w 
      INTEGER,PARAMETER :: lmax=10, krich=3, km = krich-1
      INTEGER, PARAMETER :: jnewton_max=5000
      REAL(rprec), PARAMETER :: tole = 1.e-3_dp
      REAL(rprec), PARAMETER :: newton_tol = tole
      REAL(rprec) :: init_zeta, init_theta
      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::
     1    init_zeta_v, init_theta_v
      REAL(rprec) :: init_zetak_st, init_alpha_st, init_thetak_tok,
     1    init_alpha_tok
      LOGICAL:: l_geom_input, l_tokamak_input
      LOGICAL:: lfail_balloon

!      lmax = #MAX iterations in Richardson's; krich= MIN.#samples to interpolate
!      x0 = radial wave number (default=0.)
!      kth = Mode label (=1, most unstable; 2 = second most ....)
!      k_w: number of potential wells included in integration domain
!      tole= tolerance for eigenvalue fit error
!      np0_in=coarsest grid used in richardson's; kth-largest eigenvalue obtained
!      jnewton_max= maximum number of iterations in Newton-Raphson to get theta
!      newton_tol= tolerance in Newton_Raphson
!
!      Initial location of field line can be given in several formats.
!      If L_GEOM_INPUT = T, then it is given by (init_theta, init_zeta)
!      If L_GEOM_INPUT = F, then it is given in (ALPHA, ANGLE) format,
!        being ALPHA the field line label, and ANGLE either a toroidal/poloidal initial location.
!        Then, if L_TOKAMAK_INPUT = T, (ALPHA, ANGLE) = (init_alpha_tok, init_thetak_tok), where
!          the meaning of ALPHA_TOK and THETAK are a la Connor & Hastie
!        IF L_TOKAMAK_INPUT = F, then (ALPHA, ANGLE) = (init_alpha_st, init_zetak_tok), where
!          the meaning of ALPHA_ST and ZETAK are given a la Sanchez & Hirshman
!
      END MODULE ballooning_data
