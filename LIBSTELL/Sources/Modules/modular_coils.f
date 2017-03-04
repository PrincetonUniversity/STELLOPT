      MODULE modular_coils
      USE Vcoilpts
      INTEGER, PARAMETER :: nfourier=20, nsurf_p=200
      INTEGER :: nf_rho, nf_phi, nstep, niter_opt, nmod_coils,
     1    nmid, nodd, nfper, nmod_coils_per_period,
     2    m_in, n_in, nmod_unique_coils, nmod_coeffs, nmod_currents
      INTEGER :: numsurf
      REAL(rprec), TARGET :: dcp_wgt
      REAL(rprec) :: epsfcn, dcp_exp, dcp_tgt
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  dcc_wgt, dcc_exp, dcc_tgt, cc_min
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  rc_wgt, rc_exp, rc_tgt, rc_min
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  lmod_wgt, lmod_tgt, cmod_scl
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  ymin_wgt, ymin_tgt, ymin_cls
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  cu_wgt, cu_tgt, cu_sum
      REAL(rprec), DIMENSION(ncdim), TARGET :: r_ext
      REAL(rprec), DIMENSION(ncdim), TARGET :: curmod
      REAL(rprec), DIMENSION(ncdim) :: phimin, phimax
      REAL(rprec), DIMENSION(ncdim) :: mod_length
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: x_mod, y_mod, z_mod
      REAL(rprec), DIMENSION(ncdim) :: curcon
      REAL(rprec), DIMENSION(ncdim,0:nfourier) :: phic, phis, rhoc, rhos
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rho , phi, 
     1                                            rcoil, zcoil
      INTEGER, DIMENSION(nsurf_p) :: m_num,n_num
      REAL(rprec), DIMENSION(nsurf_p) :: rmn_sf, zmn_sf
      REAL(rprec), DIMENSION(ncdim) :: phi_full
      REAL(rprec) :: p_d_min, p_d_max
      REAL (rprec), DIMENSION (ncdim) :: b_max
      LOGICAL :: lmodular
      LOGICAL :: lmodcur
      LOGICAL :: lsurfv
      LOGICAL :: lncsx
      LOGICAL :: lsymm
      END MODULE modular_coils
