      MODULE saddle_coils
      USE Vcoilpts
      INTEGER, PARAMETER :: mfourier = 100
      INTEGER :: nsad_v, nsad_u, nsad_coeffs, nfils
      INTEGER :: nsad_coils_per_period, nsad_coils, nsmid, nsodd
      INTEGER :: nsad_unique_coils
      INTEGER :: nsad_currents, num_cursad
      INTEGER, DIMENSION(ncdim) :: nsad_group
      LOGICAL :: ls_cur(ncdim)
      REAL(rprec), DIMENSION(ncdim) :: csad_scl
      REAL(rprec), DIMENSION(ncdim,0:mfourier) :: sad_v_c,
     1                  sad_v_s, sad_u_c, sad_u_s
      REAL(rprec), DIMENSION(ncdim) :: sad_u0, sad_v0,
     1                  sad_phi0, sad_theta0
      REAL(rprec), DIMENSION(ncdim) :: cursad, c_sad, sad_length
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  dsc_wgt, dsc_exp, dsc_tgt, sc_min
      REAL(rprec), DIMENSION(ncdim,ncdim), TARGET :: 
     1                  sc_dmin_tgt, sc_dmin_wgt, sc_dmin
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  rs_wgt, rs_exp, rs_tgt, rs_min
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  cs_wgt, cs_tgt, cs_sum
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  dscxp_wgt, dscxp_exp, dscxp_tgt, scxp_min
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  lsad_wgt, lsad_tgt
      REAL(rprec), DIMENSION(ncdim), TARGET ::
     1                  rmax_sad, rmax_wgt, rmax_tgt
      REAL(rprec), DIMENSION(ncdim), TARGET :: ymin_sad
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: u_sad, v_sad
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: x_sad, y_sad,
     1                  z_sad
      REAL(rprec) :: deln, delt, p_s_min, bkp_min, bkp_wgt, bkp_tgt
      REAL(rprec), DIMENSION(ncdim), TARGET :: csc_wgt, scd_wgt
      REAL(rprec), DIMENSION(ncdim) :: csc_tgt, scd_tgt
      LOGICAL :: lsaddle, lsadsfv, lsadcur, lsmod, lsadshape,
     1           lspline, lctrlpt, lsplbkp
      INTEGER, DIMENSION(ncdim,0:mfourier) :: nvar_vc, nvar_uc
      END MODULE saddle_coils
