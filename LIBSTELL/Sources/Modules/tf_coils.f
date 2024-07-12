      MODULE tf_coils
      USE stel_constants
      USE vcoilpts
      REAL(rprec) :: i_pol, i_tfc, pol_cur
      REAL(rprec), TARGET :: dpc_wgt
      INTEGER :: mtfcoil, mtfwire
      REAL(rprec), DIMENSION(ncdim) :: tfc_cur
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: tfc_x, tfc_y, tfc_z
      LOGICAL :: ltfc, ltfcv
      LOGICAL :: lqos
      LOGICAL :: lpolcur
      END MODULE tf_coils
