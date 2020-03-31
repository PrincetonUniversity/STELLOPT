      MODULE coilsnamin
      USE modular_coils
      USE saddle_coils
      USE saddle_surface
      USE vf_coils
      USE tf_coils
      USE bcoils_mod
      USE bnorm_mod
      USE control_mod
      USE vcoilpts
      IMPLICIT NONE

      NAMELIST /coilsin/ nmod_coils_per_period, nf_phi, nf_rho, epsfcn,
     1   lvf, lmodcur, lsurfv, lbnorm, rhoc, rhos, phic, phis, curmod,
     2   dcc_wgt, dcc_exp, dcc_tgt, dcp_wgt, dcp_exp, dcp_tgt, rc_wgt,
     3   rc_exp, rc_tgt, lmod_wgt, lmod_tgt, niter_opt, nstep, i_pol,
     4   numsurf, m_num, n_num, rmn_sf, zmn_sf, num_vf, rc_vf, zc_vf,
     5   cc_vf, lbcoil, lncsx, lsymm, lsaddle, nsad_coils_per_period,
     6   ltfc, ltfcv, i_tfc, lsadsfv, nsad_u, nsad_v, nfils, sad_v_c,
     7   sad_v_s, sad_u_c, sad_u_s, sad_v0, sad_u0, lsadcur, lpolcur,
     8   cursad, numsurf_sad, m_sad, n_sad, rmn_sad, zmn_sad, dsc_wgt,
     9   dsc_exp, dsc_tgt, r_ext, ymin_wgt, ymin_tgt, lmodular, lqos,
     A   rs_wgt, rs_exp, rs_tgt, lsmod, laccess, n_access, x0_access,
     B   y0_access, z0_access, x1_access, y1_access, z1_access,
     C   dac_wgt, dac_exp, dac_tgt, cvf_wgt, cvf_tgt, cs_wgt, cs_tgt,
     D   cu_wgt, cu_tgt, dpc_wgt, mc_bg, lp_bg, bcoil_cur, dscxp_wgt,
     E   dscxp_exp, dscxp_tgt, deln, delt, lspline, lsplbkp, nvar_vc,
     F   nvar_uc, bkp_wgt, bkp_tgt, mxb_wgt, lvfvar, nrvf_c, rcfc_vf,
     G   rcfs_vf, lvfr, lvfz, nopt_alg, lrestart, nopt_wsurf, rvf_wgt,
     H   rvf_tgt, nsad_group, ls_cur, csad_scl, lsad_wgt, lsad_tgt,
     I   rmax_wgt, rmax_tgt, bcoil_file, lbcoil_cur, lsadshape,
     J   lvfc, lcc_vf, csc_wgt, csc_tgt, scd_wgt, scd_tgt, lctrlpt,
     K   nwdim, nvf_fix, laxis, vacfld_wgt, sc_dmin_tgt, sc_dmin_wgt
!
!        VARIABLE DESCRIPTIONS
!
!        SCALARS
!
!        nmod_coils_per_period - no. of coils per fp, modular representation
!        nsad_coils_per_period - no. of coils per fp, saddle representation
!        nf_phi - no. of phi Fourier modes, modular rep.
!        nf_rho - no. of rho Fourier modes, modular rep.
!        epsfcn - step size for derivative approximation
!        niter_opt - max. number of function evaluations
!        nstep 
!        i_pol - total poloidal current per fp
!        num_vf  - no. of vf coil pairs
!        nsad_u  - no. of u-coefficients (Fourier or spline) in saddle rep. 
!        nsad_v  - no. of v-coefficients (Fourier or spline) in saddle rep.
!        nfils - no. of filaments per modular/saddle coil (1, 3, or 5)
!        n_access  - number of access zone constraints
!        deln  -  dist. from coil centerline to filament, normal to winding surface
!        delt  - dist. from coil centerline to filament, tangent to winding surface
!        nvar_vc - nvar_vc(i,j) = 0/1 to fix/vary spline coefficient i, coil type j
!        nvar_uc 
!        nrvf_c 
!        nopt_alg 
!        nopt_wsurf 
!        nsad_group = array of integers defining current group for 'saddle' coils
!        nwdim 
!        nvf_fix 
!       
!        LOGICAL CONTROL VARIABLES
!       
!        lvf
!        lmodcur = T to vary currents in modular representation
!        lsurfv = T to vary winding surface coefficients in modular representation
!        lbnorm = T to read bnorm coeffs. and match B-normal at plasma boundary
!        lbcoil  = T to read file containing background coils
!        lncsx  = T to extend coils in v=0 plane for NBI access
!        lsymm  = T for coil in v=0 plane with modular representation
!        lsaddle  = T to use saddle OR modular coil representation
!        ltfc 
!        ltfcv 
!        lsadsfv = T to vary winding surface in saddle representation 
!        lsadcur  = T to vary currents in saddle representation
!        lpolcur = T to implement constraint on total poloidal currents
!        lmodular  = T to use modular coil representation (defunct, superceded by lsaddle)
!        lqos = T for code to generate qos TF winding 
!        lsmod  = T to implement modular mode in saddle representation
!        laccess  = T to implement access zone constraints
!        lp_bg
!        lspline = T for spline option in saddle representation 
!        lsplbkp = T to vary spline breakpoints (not used yet)
!        lvfvar  = T to vary VF coil geometry
!        lvfr  = T to vary r-coefficients in VF coil representation
!        lvfz  = T to vary z-coordinates in VF coil representation
!        lrestart 
!        ls_cur 
!        lbcoil_cur 
!        lsadshape = T to vary coil geometry
!        lvfc  = T to vary VF coil currents
!        lcc_vf 
!        lctrlpt = T for control point spline representation (also need lspline=T)
!        laxis = T to include magnetic axis in coils.ext file
!       
!        COILS FOR MODULAR REPRESENTATION
!       
!        rhoc = coeffs. of cosine terms in modular representation of poloidal angle (u)
!        rhos = coeffs. of sine terms in modular representation of poloidal angle (u)
!        phic = coeffs. of cosine terms in modular representation of toroidal angle (v)
!        phis = coeffs. of sine terms in modular representation of of toroidal angle (v)
!        curmod = coil currents for modulars
!        numsurf  - number of winding surface modes, modular rep.
!        m_num - vector of poloidal mode numbers, modular rep.
!        n_num  - vector of toroidal mode numbers, modular rep.
!       
!        COILS FOR SADDLE REPRESENTATION
!
!        sad_v_c = coeffs. of cosine terms for v-coordinate, saddle representation
!        sad_v_s  = coeffs. of sine terms for v-coordinate, saddle representation
!        sad_u_c  = coeffs. of cosine terms for u-coordinate, saddle representation
!        sad_u_s  = coeffs. of sine terms for u-coordinate, saddle representation
!        sad_v0 = const. term used in original b-spline (non-control-pt.) repr.
!        sad_u0 = const. term used in original b-spline (non-control-pt.) repr.
!        cursad = coil currents for saddles
!        numsurf_sad  - no. of winding surface modes, saddle rep.
!        m_sad  - vector of poloidal mode numbers, saddle rep.
!        n_sad  - vector of toroidal mode numbers, saddle rep.
!       
!        SURFACE ARRAYS
!
!        rmn_sf  = R-coefficients, winding surface for modular representation
!        zmn_sf  = Z-coefficients, winding surface for modular representation
!        rmn_sad = R-coefficients, winding surface for saddle representation
!        zmn_sad  = R-coefficients, winding surface for saddle representation
!        
!        TARGET/WEIGHTS
!
!        dcc_wgt = weights for coil-coil spacing penalties, modular rep.
!        dcc_exp = exponents for coil-coil spacing, = -1 to use linear constraints
!        dcc_tgt = targets for coil-coil spacing penalties
!        dcp_wgt = weight for coil-plasma spacing penalties
 
!        dcp_exp = exponent for coil-plasma spacing, = -1 to use linear constraint
!        dcp_tgt = target tor coil-plasma spacing penalty
!        rc_wgt = weights for coil radius of curvature penalties
!        rc_exp = exp. for coil radius of curvature penalty, = -1 for linear constraints
!        rc_tgt  = targets for coil radius of curvature penalties
!        dsc_wgt = weights for coil-coil spacing penalties, 'saddle' rep.
!        dsc_exp = exponent for 'saddle' coil-coil spacing, = -1 for linear con.
!        dsc_tgt = targets for 'saddle' coil-coil spacing penalties
!        sc_dmin_tgt = 2D array of saddle coil-coil spacing penalties for each
!                      coil pair (diagonal is excluded)
!        sc_dmin_wgt = 2D array of saddle coil-coil weights
!        r_ext 
!        ymin_wgt = weights for y-min constraints
!        ymin_tgt = targets for y-min constraints
!        rs_wgt = weights for 'saddle' coil radius of curvature penalties
!        rs_exp = exp. For 'saddle' coil radius of curvature penalty, = -1 for linear
!        rs_tgt  = targets for 'saddle' coil radius of curvature penalties
!        lsad_wgt = weights for 'saddle' coil length constraints
!        lsad_tgt = target values for 'saddle' coil length constraints
!        lmod_wgt = weights for modular coil length constraints
!        lmod_tgt = target values for modular coil length constraints
!        dac_wgt = weights for access penalties
!        dac_exp = exponent for access penalties
!        dac_tgt = target distance for access penalties
!        cvf_wgt 
!        cvf_tgt 
!        cs_wgt 
!        cs_tgt
!        cu_wgt 
!        cu_tgt 
!        dpc_wgt 
!        dscxp_wgt = weights for 'saddle' coil linear current density penalties
!        dscxp_exp = exponents for 'saddle' coil linear current density penalties
!        dscxp_tgt = targets for 'saddle' coil linear current density penalties
!        bkp_wgt 
!        bkp_tgt 
!        mxb_wgt 
!        rmax_wgt 
!        rmax_tgt 
!        rvf_wgt
!        rvf_tgt 
!        csc_wgt 
!        csc_tgt 
!        scd_wgt 
!        scd_tgt 

!        rc_vf 
!        zc_vf
!        cc_vf 
!        i_tfc 
!        x0_access = coordinates of endpoints for lines defining access zones
!        y0_access  = 
!        z0_access = 
!        x1_access = 
!        y1_access = 
!        z1_access= 
!        mc_bg 
!        bcoil_cur = 'background' coil currents
!        rcfc_vf
!        rcfs_vf 
!        csad_scl 
!        bcoil_file = name of 'background' coil file

      CONTAINS

      SUBROUTINE read_coils_namelist (iunit, istat)
      INTEGER :: iunit, istat

      READ (iunit, nml=coilsin, iostat=istat)

      END SUBROUTINE read_coils_namelist

      END MODULE coilsnamin
