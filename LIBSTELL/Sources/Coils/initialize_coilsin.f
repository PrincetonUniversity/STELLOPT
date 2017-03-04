      SUBROUTINE initialize_coilsin
      USE stel_kinds
      USE coilsnamin

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
!-----------------------------------------------
      bcoil_file = "bcoilxyz.dat"
      lvf = .false.
      lvfc = .true.
      lcc_vf = .true.
      lvfvar = .false.
      lvfr = .false.
      lvfz = .false.
      ltfc = .false.
      ltfcv = .false.
      lbcoil = .false.
      lbcoil_cur = .true.
      lbnorm = .false.
      lmodcur = .false.
      lsadcur = .false.
      lsadshape = .true.
      lmodular = .false.
      lsaddle = .true.
      lsmod = .false.
      lspline = .false.
      lctrlpt = .false.
      lsplbkp = .false.
      lncsx = .false.
      lqos = .false.
      lsymm = .true.
      laccess = .false.
      laxis = .false.
      lpolcur = .false.
      lp_bg = .true.
      lrestart = .false.
      ls_cur = .true.

      nwdim = 400
      nvf_fix = 2               !(=2 for QPS)
      n_access = 0
      nf_phi = 0
      nf_rho = 0
      nmod_coils_per_period = 0
      nsad_coils_per_period = 0
      nsad_u = 0
      nsad_v = 0
      numsurf = 0
      numsurf_sad = 0
      nstep  = 10
      niter_opt = 1
      nopt_alg = 0
      nopt_wsurf = -1
      nfils = 1
      nvar_vc = 1
      nvar_uc = 1

      i_pol = 0
      i_tfc = 0
      rhoc = 0
      rhos = 0
      phic = 0
      phis = 0
      curmod = 0
      dcc_wgt = 0
      dcc_exp = 0
      dcc_tgt = 0
      rc_wgt = 0
      rc_exp = 0
      rc_tgt = 0
      lmod_wgt = 0
      lmod_tgt = 0
      lsad_wgt = 0
      lsad_tgt = 0
      sad_v_c = 0
      sad_v_s = 0
      sad_u_c = 0
      sad_u_s = 0
      sad_v0 = 0
      sad_u0 = 0
      cursad = 0
      cc_vf = 0
      r_ext = 0
      ymin_wgt = 0
      ymin_tgt = 0
      rmax_wgt = 0
      rmax_tgt = 0
      rs_wgt = 0
      rs_exp = 0
      rs_tgt = 0
      dac_wgt = 0
      dac_exp = 0
      dac_tgt = 0
      cvf_wgt = 0
      cvf_tgt = 0
      rvf_wgt = 0
      rvf_tgt = 0
      rc_vf = 0
      zc_vf = 0
      rcfc_vf = 0
      rcfs_vf = 0
      cs_wgt = 0
      cs_tgt = 0
      csc_wgt = 0
      csc_tgt = 0
      scd_wgt = 0
      scd_tgt = 0
      dpc_wgt = 0
      vacfld_wgt = 0
      mc_bg = 0
      bcoil_cur = 0
      dscxp_wgt = 0
      dscxp_exp = 0
      dscxp_tgt = 0
      delt = 0.032_dp
      deln = 0.055_dp           ! NCSX value
      bkp_wgt = 10000
      bkp_tgt = 0.02_dp
      mxb_wgt = 0
      csad_scl = 1
      nmid = 0
      nsmid = 0
      n_access = 0
      num_vf = 0
      nrvf_c = 0
      dsc_wgt = 0
      sc_dmin_wgt = 0

      DO i = 1, ncdim
         nsad_group(i) = i
      END DO

      END SUBROUTINE initialize_coilsin
