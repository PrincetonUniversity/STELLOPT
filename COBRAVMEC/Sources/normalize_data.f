      MODULE normalize_data
      USE stel_kinds
      IMPLICIT NONE
      LOGICAL:: lasym_v                                                  ! 110909 RS: logical (=T, for ASYMMETRIC input; =F, stellarator symmetry)
      INTEGER :: nfp_v
      REAL(rprec), PARAMETER :: mu_0=1.2566368e-6_dp
      REAL(rprec) :: r0, b0_v, amin, beta0, twopi
      END MODULE normalize_data
