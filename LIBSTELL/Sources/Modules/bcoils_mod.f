      MODULE bcoils_mod
      USE vcoilpts
      INTEGER :: mbcoils, mc_max
      INTEGER, DIMENSION(ncdim) :: mbwires, mc_bg
      REAL(rprec) :: mxb_wgt
      REAL(rprec), DIMENSION(ncdim) :: bcoil_cur, cc_bg
      REAL(rprec), DIMENSION(:, :), ALLOCATABLE ::  
     1             bcoil_x, bcoil_y, bcoil_z
      LOGICAL :: lbcoil, lbcoil_cur
      LOGICAL :: lp_bg(ncdim)
      INTEGER :: n_access, n_access_pts
      REAL(rprec), DIMENSION(ncdim) :: x0_access, y0_access, z0_access
      REAL(rprec), DIMENSION(ncdim) :: x1_access, y1_access, z1_access
      REAL(rprec), DIMENSION(ncdim) :: dac_exp, dac_tgt, acc_min
      REAL(rprec), DIMENSION(ncdim), TARGET :: dac_wgt
      REAL(rprec), DIMENSION(:, :), ALLOCATABLE ::
     1  x_access, y_access, z_access
      LOGICAL :: laccess, laxis
      CHARACTER(LEN=200) :: bcoil_file
      END MODULE bcoils_mod
