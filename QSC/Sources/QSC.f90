subroutine QSC(xt)   !10/12/18,3/1/19.(6b,9g)bulk of orig quasisymmetry.f90.

  use quasisymmetry_variables, only: start_time, total_time, general_option, general_option_single, general_option_scan, &
       N_procs, mpi_rank, proc0, &
       xtqsc,max_axis_nmax, max_n, R0c,R0s,Z0c,Z0s              !1/29/19.(8o7)=(8l15)
  use vmec_input, only: raxis_cc, raxis_cs, zaxis_cc, zaxis_cs  !1/29/19.(8o6)=(8l14)
  use mpi_params  !10/17/18.

  implicit none

  include 'mpif.h'

  integer :: tic, toc, countrate, ierr
!  real :: start_time, end_time  !out-3/2/19
  real :: end_time
  character(len=120) :: xt  !hm-10/14/18.

  xtqsc=xt
!  write(0,*)'hm-10/16.bgn QSC. xt,xtqsc=',trim(xt),' ',trim(xtqsc)
!  call mpi_init(ierr)
!  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
!  call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
  mpi_rank=0 ; N_procs=1         !hm-10/15/18.(6e6,6e20)
  proc0 = (mpi_rank==0)

!  if (proc0) then
!     print "(a)"," -------------------------------------------------------------"
!     print *,"Quasisymmetry solver"
!  end if
   !call system_clock(tic,countrate)
  call cpu_time(start_time)

!  call quasisymmetry_read_input()   !hm-10/28/18.(6e21d)mv to bfr QSC().

!1/29/19.(8o6)initialize [R0c,..] fra [raxis_cc,..], inverse of qs_rd_input.
!  max_n = min(ntord,max_axis_nmax)  !as qs_wr_vm_input,qs_rd_input. (8o7)c-out
  R0c(1:max_n+1) = raxis_cc(0:max_n)
  R0s(1:max_n+1) = -raxis_cs(0:max_n)
  Z0c(1:max_n+1) = zaxis_cc(0:max_n)
  Z0s(1:max_n+1) = -zaxis_cs(0:max_n)

  call quasisymmetry_validate_input()
  general_option = general_option_single  !10/10/18. = "single", fra qs_vars.f90
  select case (trim(general_option))
  case (general_option_single)
     call quasisymmetry_single_solve()
     call quasisymmetry_write_vmec_input()
!  case (general_option_scan)
!     call quasisymmetry_scan()
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  !call system_clock(toc)
  !total_time = real(toc-tic)/countrate
  call cpu_time(end_time)
  total_time = end_time - start_time

!  call quasisymmetry_write_output()
!  if (proc0) then          !1/2/19.(8b)c-out
!     print "(a)"," -------------------------------------------------------------"
!     print "(a,es10.3,a)","Quasisymmetry solver is complete. Total time=",total_time," sec."
!  end if

!  call mpi_finalize(ierr)

end subroutine QSC
