      MODULE boozer_params
      USE stel_kinds
      INTEGER :: mboz, nboz, nu_boz, nv_boz, nunv, mnboz, nu2_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1  rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  rmnc_opt, zmns_opt,lmns_opt
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  rmns_opt, zmnc_opt,lmnc_opt
      END MODULE boozer_params
