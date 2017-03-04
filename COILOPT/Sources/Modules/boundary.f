      MODULE boundary
      USE stel_constants
      INTEGER :: mpol, ntor, nfp, mnmax, ns, nu, nv, nuv, nedge,
     1           mnmax_nyq
      REAL(rprec) :: TotalArea, sum_d_area, rbphi_avg
      REAL(rprec) :: iota_b, phip_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: phib, thetab
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rb, zb
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rb_th, rb_ph
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: zb_th, zb_ph
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: n_r, n_phi, n_z
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: d_area
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: x_p, y_p, z_p
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: phi_d, theta_d
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xm_b, xn_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rmnc_b, zmns_b, lmns_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rmnc_a, zmns_a
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bnormal_match
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: b_error, b_mod
      CHARACTER boundary_mn_file*30
      END MODULE boundary   
