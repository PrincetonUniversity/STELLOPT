subroutine quasisymmetry_read_input(xt)

  use quasisymmetry_variables
!  use mpi_params, only: myid,master
  use vmec_input, only: raxis_cc, raxis_cs, zaxis_cc, zaxis_cs  !12/21/18.(7m12c)
  use vparams, only: ntord  !12/22/18.(7m12d)

  implicit none

  integer :: numargs
  character(len=200) :: input_filename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999
  real(dp) :: threshold
  character(len=120) :: xt  !hm-10/24/18.
!  integer :: max_n  !hm-12/21/18.(7m12c)as qs_wr_vm_input. 1/29/19.(8o7)xfer to qs_vars.

  xtqsc=xt   !10/28/18.(6e21d)
  max_n = min(ntord,max_axis_nmax) !12/21/18.(7m12c).xferred fra qs_wr_vm_input.

  R0s = 0
  R0c = 0
  Z0s = 0
  Z0c = 0

  ! getcarg is in LIBSTELL
  ! call getcarg(1, input_filename, numargs)
  numargs=1   !hm-10/15/18.
  input_filename="input." // xtqsc  !10/28/18.(6e21d)

!  if (numargs<1) then  !c-out-10/28/18.
! ...
  new_vmec_filename = input_filename  !10/28/18.(6e21d)  &quasisymmetry now in input.<xt>.
  vmec_template_filename = input_filename !11/15/18.(7l7).

  fileUnit=11
  open(unit=fileUnit, file=input_filename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(input_filename)
     stop
  else
     read(fileUnit, nml=quasisymmetry, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error!  I was able to open the file ", trim(input_filename), &
               " but not read data from the quasisymmetry namelist in it."
        if (didFileAccessWork==-1) then
           print *,"Make sure there is a carriage return after the / at the end of the namelist!"
        end if
        stop
     end if
     if (proc0) print *,"Successfully read parameters from quasisymmetry namelist in ", trim(input_filename), "."
  end if
  close(unit = fileUnit)

!  verbose = (trim(verbose_option)==verbose_option_all .or. (proc0 .and. trim(verbose_option)==verbose_option_proc0))       -out-12/21/18.(7m12c)
  verbose=.false.       !12/21/18.(7m12c)
  proc0=.true.       !12/22/18.(7m12d)

  N_phi_original = N_phi   !12/4/18.(7l23a)xferred in fra qs_validate_input.

!!$  ! Count the number of nonzero entries at the beginning of N_phis:
!!$  N_N_phis = 0
!!$  do i=1,max_N_N_phis
!!$     if (N_phis(i) == 0) then
!!$        N_N_phis = i-1
!!$        exit
!!$     end if
!!$  end do

  ! Find how many Fourier modes we are keeping in the axis shape:
  threshold = 1.0d-14
  if (trim(general_option)==general_option_single) then
     ! Set axis_nmax using the single-run options.
     do i = max_axis_nmax+1,1,-1
        !print *,'********* i=',i
        !if (R0c(i) .ne. 0 .or. R0s(i) .ne. 0 .or. Z0s(i) .ne. 0 .or. Z0c(i) .ne. 0) then
        if (abs(R0c(i)) > threshold .or. abs(R0s(i)) > threshold .or. abs(Z0s(i)) > threshold .or. abs(Z0c(i)) > threshold) then
           !print *,"R0c(i):",R0c(i)
           !print *,"R0s(i):",R0s(i)
           !print *,"Z0c(i):",Z0c(i)
           !print *,"Z0s(i):",Z0s(i)
           axis_nmax = i-1
           exit
        end if
     end do

     if (proc0) then
        print *,"R0c:", R0c(1:axis_nmax+1)
        print *,"R0s:", R0s(1:axis_nmax+1)
        print *,"Z0c:", Z0c(1:axis_nmax+1)
        print *,"Z0s:", Z0s(1:axis_nmax+1)
     end if
  else
     ! Set axis_nmax using the scan options
     do i = max_axis_nmax+1,1,-1
        if (abs(R0c_min(i)) > threshold .or. abs(R0s_min(i)) > threshold .or. abs(Z0s_min(i)) > threshold .or. abs(Z0c_min(i)) > threshold &
             .or. abs(R0c_max(i)) > threshold .or. abs(R0s_max(i)) > threshold .or. abs(Z0s_max(i)) > threshold .or. abs(Z0c_max(i)) > threshold) then
           axis_nmax = i-1
           exit
        end if
     end do
  end if

  if (proc0) print *,"axis_nmax:",axis_nmax

!12/21/18.(7m12c)xferred in [raxis_cc,..]=[R0c,..] fra qs_wr_vm_input.
  ! To convert sin(...) modes to vmec, we introduce a minus sign. This is because in vmec,
  ! R and Z ~ sin(m theta - n phi), which for m=0 is sin(-n phi) = -sin(n phi).
  raxis_cc(0:max_n) = R0c(1:max_n+1)
  raxis_cs(0:max_n) = -R0s(1:max_n+1)
  zaxis_cc(0:max_n) = Z0c(1:max_n+1)
  zaxis_cs(0:max_n) = -Z0s(1:max_n+1)

end subroutine quasisymmetry_read_input
