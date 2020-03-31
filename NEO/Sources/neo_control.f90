
MODULE neo_control
! Control parameters from inpt file
  USE neo_precision
  CHARACTER(120)                     :: in_file                             !SPH
  CHARACTER(120)                     :: out_file                            !SPH
  CHARACTER(120)                      :: cur_file
  INTEGER                            :: theta_n
  INTEGER                            :: phi_n
  INTEGER                            :: s_ind_in
  INTEGER                            :: write_progress
  INTEGER                            :: write_output_files
  INTEGER                            :: spline_test
  INTEGER                            :: max_m_mode, max_n_mode
  INTEGER                            :: lab_swi, inp_swi, ref_swi, eout_swi
  INTEGER                            :: no_fluxs
  INTEGER, DIMENSION(:), ALLOCATABLE :: fluxs_arr
END MODULE neo_control
