      MODULE bnorm_mod
      USE stel_constants
      USE vcoilpts
      INTEGER, PARAMETER :: bnorm_dim=1000
      INTEGER :: mnbn_max
      REAL(rprec) :: vacfld_wgt
      REAL(rprec), DIMENSION(bnorm_dim) :: xbn_m, xbn_n, bn_coef
      INTEGER :: mmax_bmn, nmax_bmn, mnmax_bmn
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xm_bmn, xn_bmn
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bmn, bmn_error
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: luv, bsl_error
      CHARACTER(LEN=200) :: bnorm_file
      LOGICAL :: lbnorm
      END MODULE bnorm_mod
