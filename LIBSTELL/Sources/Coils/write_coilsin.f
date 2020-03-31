      SUBROUTINE write_coilsin (iunit, istat)
      USE stel_constants
      USE coilsnamin
      USE write_array_generic
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i,j
      REAL(rprec) :: tempc(nmid), temps(nmid)
      REAL(rprec) :: temp2c(MAX(nsmid,num_vf)), 
     1               temp2s(MAX(nsmid,num_vf))
!-----------------------------------------------
      istat = 0

!     Write out two-dimensional arrays with indices explicitly
!     so that IF NFPDIM changes, the NAMELIST can still be READ
!     correctly

      WRITE(iunit,'(1x,a8)')'&COILSIN'
      WRITE(iunit,'(1x,3a)') "BCOIL_FILE = '", TRIM(bcoil_file),"'"
      WRITE(iunit,100) 'LRESTART = ', lrestart
      WRITE(iunit,100) 'LSURFV = ', lsurfv
      WRITE(iunit,100) 'LSADSFV = ', lsadsfv
      WRITE(iunit,100) 'LVF = ', lvf
      WRITE(iunit,100) 'LVFC = ', lvfc
      WRITE(iunit,100) 'LVFVAR = ', lvfvar
      WRITE(iunit,100) 'LVFR = ', lvfr
      WRITE(iunit,100) 'LVFZ = ', lvfz
      WRITE(iunit,100) 'LTFC = ', ltfc
      WRITE(iunit,100) 'LTFCV = ', ltfcv
      WRITE(iunit,100) 'LSADDLE = ', lsaddle
      WRITE(iunit,100) 'LSMOD = ', lsmod
      WRITE(iunit,100) 'LSPLINE = ', lspline
      WRITE(iunit,100) 'LCTRLPT = ', lctrlpt
      WRITE(iunit,100) 'LSPLBKP = ', lsplbkp
      WRITE(iunit,100) 'LMODULAR = ', lmodular
      WRITE(iunit,100) 'LMODCUR = ', lmodcur
      WRITE(iunit,100) 'LSADCUR = ', lsadcur
      WRITE(iunit,100) 'LSADSHAPE = ', lsadshape
      WRITE(iunit,100) 'LPOLCUR = ', lpolcur
      WRITE(iunit,100) 'LBNORM = ', lbnorm
      WRITE(iunit,100) 'LBCOIL = ', lbcoil
      WRITE(iunit,100) 'LBCOIL_CUR = ', lbcoil_cur
      WRITE(iunit,100) 'LNCSX = ', lncsx
      WRITE(iunit,100) 'LQOS = ', lqos
      WRITE(iunit,100) 'LSYMM = ', lsymm
      WRITE(iunit,100) 'LACCESS = ', laccess
      WRITE(iunit,100) 'LAXIS = ', laxis
      WRITE(iunit,200) 'NOPT_ALG = ', nopt_alg
      WRITE(iunit,200) 'NOPT_WSURF = ', nopt_wsurf
      WRITE(iunit,200) 'NITER_OPT = ',niter_opt, 'NSTEP = ', nstep
      WRITE(iunit,200) 'NWDIM = ', nwdim
      WRITE(iunit,200) 'NVF_FIX = ', nvf_fix
      WRITE(iunit,200) 'NF_PHI = ',nf_phi
      IF (nf_rho .gt. 0) THEN
         WRITE(iunit,200) 'NF_RHO = ',nf_rho
      END IF
      WRITE(iunit,300) 'EPSFCN = ', epsfcn
      WRITE(iunit,300) 'I_POL = ', i_pol
      WRITE(iunit,300) 'I_TFC = ', i_tfc
      WRITE(iunit,200) 'NMOD_COILS_PER_PERIOD = ',
     1                  nmod_coils_per_period
      WRITE(iunit,200) 'NUM_VF = ', num_vf
!     Use following weights, etc., for both modular and saddle minimum
!     coil-plasma distance penalties
      WRITE(iunit,300) 'DCP_WGT = ', dcp_wgt, 'DCP_exp = ', dcp_exp,
     1                 'DCP_TGT = ', dcp_tgt
      WRITE(iunit,450) 'DPC_WGT = ', dpc_wgt
      WRITE(iunit,450) 'VACFLD_WGT = ', vacfld_wgt
      WRITE(iunit,450) 'MXB_WGT = ', mxb_wgt
      IF (nmod_coils_per_period .GT. 0) THEN
         CALL write_array(iunit,'DCC_WGT', dcc_wgt, nmid)
         CALL write_array(iunit,'DCC_EXP', dcc_exp, nmid)
         CALL write_array(iunit,'DCC_TGT', dcc_tgt, nmid)
         CALL write_array(iunit,'LMOD_WGT',lmod_wgt, nmid)
         CALL write_array(iunit,'LMOD_TGT',lmod_tgt, nmid)
         CALL write_array(iunit,'RC_WGT', rc_wgt, nmid)
         CALL write_array(iunit,'RC_EXP', rc_exp, nmid)
         CALL write_array(iunit,'RC_TGT', rc_tgt, nmid)
         CALL write_array(iunit,'CU_WGT', cu_wgt, nmid)
         CALL write_array(iunit,'CU_TGT', cu_tgt, nmid)
      END IF
      IF (lmodular .AND. (.NOT.lsaddle)) THEN
         CALL write_array(iunit,'YMIN_WGT', ymin_wgt, nmid)
         CALL write_array(iunit,'YMIN_TGT', ymin_tgt, nmid)
      END IF
      IF (lsaddle) THEN
         CALL write_array(iunit,'YMIN_WGT', ymin_wgt, nsmid)
         CALL write_array(iunit,'YMIN_TGT', ymin_tgt, nsmid)
      END IF
      IF (lncsx) THEN
         CALL write_array(iunit,'R_EXT', r_ext, nmid)
      END IF
      IF (laccess) THEN
         WRITE(iunit,200) 'N_ACCESS = ', n_access
         CALL write_array(iunit,'X0_ACCESS', x0_access, n_access)
         CALL write_array(iunit,'Y0_ACCESS', y0_access, n_access)
         CALL write_array(iunit,'Z0_ACCESS', z0_access, n_access)
         CALL write_array(iunit,'X1_ACCESS', x1_access, n_access)
         CALL write_array(iunit,'Y1_ACCESS', y1_access, n_access)
         CALL write_array(iunit,'Z1_ACCESS', z1_access, n_access)
         CALL write_array(iunit,'DAC_WGT', dac_wgt, n_access)
         CALL write_array(iunit,'DAC_EXP', dac_exp, n_access)
         CALL write_array(iunit,'DAC_TGT', dac_tgt, n_access)
      END IF
      IF (lbcoil) THEN
         CALL write_array(iunit,'MC_BG', mc_bg, mbcoils)
         CALL write_array(iunit,'LP_BG',lp_bg, mbcoils)
         CALL write_array(iunit,'BCOIL_CUR', bcoil_cur, mbcoils)
      END IF
      IF (num_vf .GT. 0) THEN
         CALL write_array(iunit,'LCC_VF', lcc_vf, num_vf)
         CALL write_array(iunit,'CC_VF', cc_vf, num_vf)
         CALL write_array(iunit,'RC_VF', rc_vf, num_vf)
         CALL write_array(iunit,'ZC_VF', zc_vf, num_vf)
         CALL write_array(iunit,'CVF_WGT', cvf_wgt, num_vf)
         CALL write_array(iunit,'CVF_TGT', cvf_tgt, num_vf)
         CALL write_array(iunit,'RVF_WGT', rvf_wgt, num_vf)
         CALL write_array(iunit,'RVF_TGT', rvf_tgt, num_vf)
         WRITE(iunit,200) 'NRVF_C = ', nrvf_c
         IF (nrvf_c .GT. 0) THEN
            DO j = 1, nrvf_c
              DO i = 1, num_vf
                temp2c(i) = rcfc_vf(i,j)
                IF (ABS(temp2c(i)) .lt. 1.e-10_dp) temp2c(i) = 0
                temp2s(i) = rcfs_vf(i,j)
                IF (ABS(temp2s(i)) .lt. 1.e-10_dp) temp2s(i) = 0
              END DO
              CALL write_array(iunit,'RCFC_VF', temp2c, num_vf,j)
              CALL write_array(iunit,'RCFS_VF', temp2s, num_vf,j)
            END DO
         END IF
      END IF
      IF (nmod_coils_per_period .gt. 0) THEN
         CALL write_array(iunit,'CURMOD', curmod, nmid)
      END IF
      IF (nf_phi .GT. 0) THEN
         DO j = 0,nf_phi
           DO i = 1,nmid
             tempc(i) = phic(i,j)                        !modular(i)%phic(j)
             temps(i) = phis(i,j)                        !modular(i)%phis(j)
             IF (ABS(tempc(i)) .lt. 1.e-10_dp) tempc(i) = 0
             IF (ABS(temps(i)) .lt. 1.e-10_dp) temps(i) = 0
           END DO
           CALL write_array(iunit,'PHIC', tempc, nmid,j)
           CALL write_array(iunit,'PHIS', temps, nmid,j)
         END DO
      END IF
      IF (nf_rho .GT. 0) THEN
         DO j = 0,nf_rho
           DO i = 1,nmid
             tempc(i) = rhoc(i,j)                        !modular(i)%rhoc(j)
             temps(i) = rhos(i,j)                        !modular(i)%rhos(j)
             IF (ABS(tempc(i)) .lt. 1.e-10_dp) tempc(i) = 0
             IF (ABS(temps(i)) .lt. 1.e-10_dp) temps(i) = 0
           END DO
           CALL write_array(iunit,'RHOC', tempc, nmid,j)
           CALL write_array(iunit,'RHOS', temps, nmid,j)
         END DO
      END IF

! 2/20/98 WHM Write out surface coefficients

      IF (numsurf .GT. 0) THEN
         WRITE(iunit,200)'NUMSURF = ',numsurf
         CALL write_array(iunit,'M_NUM', m_num, numsurf)
         CALL write_array(iunit,'N_NUM', n_num, numsurf)
         CALL write_array(iunit,'RMN_SF', rmn_sf, numsurf)
         CALL write_array(iunit,'ZMN_SF', zmn_sf, numsurf)
      END IF

!     Saddle coil input parameters

      WRITE(iunit,200) 'NSAD_COILS_PER_PERIOD = ',nsad_coils_per_period
      WRITE(iunit,200) 'NSAD_U = ',nsad_u
      WRITE(iunit,200) 'NSAD_V = ',nsad_v
      IF (nsad_coils_per_period .gt. 0) THEN
         WRITE(iunit,200) 'NFILS = ',nfils
         WRITE(iunit,300) 'DELN = ', deln, 'DELT = ', delt
         CALL write_array(iunit,'NSAD_GROUP', nsad_group, nsmid)
         CALL write_array(iunit,'CSAD_SCL', csad_scl, nsmid)
         CALL write_array(iunit,'CURSAD', cursad, nsmid)
         CALL write_array(iunit,'LS_CUR', ls_cur, nsmid)
         CALL write_array(iunit,'CSC_WGT', csc_wgt, nsmid)
         CALL write_array(iunit,'CSC_TGT', csc_tgt, nsmid)
         CALL write_array(iunit,'DSC_WGT', dsc_wgt, nsmid)
         CALL write_array(iunit,'DSC_EXP', dsc_exp, nsmid)
         CALL write_array(iunit,'DSC_TGT', dsc_tgt, nsmid)

         DO j = 1, nsad_coils
           CALL write_array(iunit, 'SC_DMIN_TGT', 
     1          sc_dmin_tgt(1:nsad_unique_coils,j),nsad_unique_coils,j)
           CALL write_array(iunit, 'SC_DMIN_WGT', 
     1          sc_dmin_wgt(1:nsad_unique_coils,j),nsad_unique_coils,j)
         END DO

         CALL write_array(iunit,'RS_WGT', rs_wgt, nsmid)
         CALL write_array(iunit,'RS_EXP', rs_exp, nsmid)
         CALL write_array(iunit,'RS_TGT', rs_tgt, nsmid)
         CALL write_array(iunit,'CS_WGT', cs_wgt, nsmid)
         CALL write_array(iunit,'CS_TGT', cs_tgt, nsmid)
         CALL write_array(iunit,'LSAD_WGT', lsad_wgt, nsmid)
         CALL write_array(iunit,'LSAD_TGT', lsad_tgt, nsmid)
         CALL write_array(iunit,'RMAX_WGT', rmax_wgt, nsmid)
         CALL write_array(iunit,'RMAX_TGT', rmax_tgt, nsmid)
         CALL write_array(iunit,'SAD_V0', sad_v0, nsmid)
         CALL write_array(iunit,'SAD_U0', sad_u0, nsmid)
         CALL write_array(iunit,'DSCXP_WGT', dscxp_wgt, nsmid)
         CALL write_array(iunit,'DSCXP_exp', dscxp_exp, nsmid)
         CALL write_array(iunit,'DSCXP_TGT', dscxp_tgt, nsmid)
         CALL write_array(iunit,'SCD_WGT', scd_wgt, nsmid)
         CALL write_array(iunit,'SCD_TGT', scd_tgt, nsmid)
         IF (lspline .and. lsplbkp) THEN
            WRITE(iunit,350) 'BKP_WGT = ',bkp_wgt,'BKP_TGT = ',bkp_tgt
         END IF
      END IF
      IF (nsad_v .GT. 0) THEN
         DO j = 0,nsad_v
           DO i = 1,nsmid
             temp2c(i) = sad_v_c(i,j)                    !saddle(i)%v_c(j)
             temp2s(i) = sad_v_s(i,j)                    !saddle(i)%v_s(j)
             IF (ABS(temp2c(i)) .lt. 1.e-10_dp) temp2c(i) = 0
             IF (ABS(temp2s(i)) .lt. 1.e-10_dp) temp2s(i) = 0
           END DO
           IF (lspline) THEN
              CALL write_array(iunit,'NVAR_VC', nvar_vc(1:nsmid,j), 
     1                         nsmid, j)
           END IF
           CALL write_array(iunit,'SAD_V_C', temp2c, nsmid, j)
           CALL write_array(iunit,'SAD_V_S', temp2s, nsmid, j)
         END DO
      END IF
      IF (nsad_u .GT. 0) THEN
         DO j = 0,nsad_u
           DO i = 1,nsmid
             temp2c(i) = sad_u_c(i,j)                   !saddle(i)%u_c(j)
             temp2s(i) = sad_u_s(i,j)                   !saddle(i)%u_s(j)
             IF (ABS(temp2c(i)) .lt. 1.e-10_dp) temp2c(i) = 0
             IF (ABS(temp2s(i)) .lt. 1.e-10_dp) temp2s(i) = 0
           END DO
           IF (lspline) THEN
              CALL write_array(iunit,'NVAR_UC', nvar_uc(1:nsmid,j), 
     1                         nsmid, j)
           END IF
           CALL write_array(iunit,'SAD_U_C', temp2c, nsmid, j)
           CALL write_array(iunit,'SAD_U_S', temp2s, nsmid, j)
         END DO
      END IF

!     Saddle surface coefficients

      IF (numsurf_sad .gt. 0) THEN
         WRITE(iunit,200)'NUMSURF_SAD = ', numsurf_sad
         CALL write_array(iunit,'M_SAD', m_sad, numsurf_sad)
         CALL write_array(iunit,'N_SAD', n_sad, numsurf_sad)
         CALL write_array(iunit,'RMN_SAD', rmn_sad, numsurf_sad)
         CALL write_array(iunit,'ZMN_SAD', zmn_sad, numsurf_sad)
      END IF

      WRITE(iunit,'(a)')'/'

 100  FORMAT(4(1x,a,l2,','))
 200  FORMAT(4(1x,a,i6,','))
 300  FORMAT(3(1x,a,1pe12.4,','))
 350  FORMAT(2(1x,a,1pe12.4,','))
 450  FORMAT(1x,a,1pe12.4,',')


      END SUBROUTINE write_coilsin
