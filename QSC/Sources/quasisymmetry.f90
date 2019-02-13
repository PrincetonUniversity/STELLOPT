! Main program.  10/12/18.wrapper for QSC().

program quasisymmetry

  use quasisymmetry_variables, only: total_time, general_option, general_option_single, general_option_scan, &
       N_procs, mpi_rank, proc0

  implicit none

  include 'mpif.h'

  call QSC()  !10/12/18.(6b)bulk of orig quasisymmetry.f90.

end program quasisymmetry
