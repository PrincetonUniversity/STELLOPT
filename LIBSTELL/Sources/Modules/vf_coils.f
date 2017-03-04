      MODULE vf_coils
      USE Vcoilpts
      INTEGER :: nvf_fix
      INTEGER :: num_vf, nvf, nrvf_c, nvf_coeffs
      REAL(rprec), DIMENSION(ncdim) :: rc_vf, zc_vf, cc_vf
      REAL(rprec), DIMENSION(ncdim,ncdim) :: rcfc_vf, rcfs_vf
      REAL(rprec), DIMENSION(ncdim) :: rvf, zvf, cvf, cvf_tgt
      REAL(rprec), DIMENSION(ncdim), TARGET :: cvf_wgt, rvf_wgt
      REAL(rprec), DIMENSION(ncdim) :: rvf_tgt, rvf_max
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: x_vf, y_vf, z_vf
      REAL(rprec) :: cvf_ssq
      LOGICAL :: lvf, lvfc, lvfvar, lvfr, lvfz
      LOGICAL, DIMENSION(ncdim) :: lcc_vf
      END MODULE vf_coils
