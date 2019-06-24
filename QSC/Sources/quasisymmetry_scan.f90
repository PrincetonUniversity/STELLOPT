subroutine quasisymmetry_scan

  use quasisymmetry_variables

  implicit none

  include 'mpif.h'

  integer :: j, k, j_eta_bar, j_sigma_initial
  integer*8 :: j_scan_local, N_Fourier_scan, N_scan_local, j_scan, j_Fourier_scan, index
  integer, dimension(max_axis_nmax+1, 4) :: scan_state
  integer :: ierr, tag
  integer*8 :: dummy_long(1), scan_index_min, scan_index_max, print_summary_stride
  integer :: mpi_status(MPI_STATUS_SIZE)
  integer, parameter :: buffer_length = 100
  character(len=buffer_length) :: proc_assignments_string
  real(dp), dimension(:), allocatable :: iotas_local, max_elongations_local, rms_curvatures_local, max_curvatures_local, axis_lengths_local
  real(dp), dimension(:), allocatable :: standard_deviations_of_R_local, standard_deviations_of_Z_local, max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local
  integer, dimension(:), allocatable :: axis_helicities_local, B_helicities_local, effective_nfps_local
  logical, dimension(:), allocatable :: iota_tolerance_achieveds_local, elongation_tolerance_achieveds_local, Newton_tolerance_achieveds_local
  integer*8, dimension(:), allocatable :: N_solves_kept
  real(dp), dimension(:), allocatable :: scan_eta_bar_local, scan_sigma_initial_local
  real(dp), dimension(:,:), allocatable :: scan_R0c_local, scan_R0s_local, scan_Z0c_local, scan_Z0s_local
  logical, dimension(:), allocatable :: nonzero_modes
  logical :: keep_going
  real(dp) :: thresh, time1, time2
  integer*8, dimension(max_axis_nmax+1, 4) :: N_scan_array_long

  allocate(nonzero_modes(axis_nmax))

  ! Clean up scan arrays
  do j = 1,max_axis_nmax+1
     if (R0s_N_scan(j)<1) R0s_N_scan(j) = 1
     if (R0c_N_scan(j)<1) R0c_N_scan(j) = 1
     if (Z0s_N_scan(j)<1) Z0s_N_scan(j) = 1
     if (Z0c_N_scan(j)<1) Z0c_N_scan(j) = 1

     N_scan_array(j,1) = max(R0s_N_scan(j),1)
     N_scan_array(j,2) = max(R0c_N_scan(j),1)
     N_scan_array(j,3) = max(Z0s_N_scan(j),1)
     N_scan_array(j,4) = max(Z0c_N_scan(j),1)
  end do
  if (eta_bar_N_scan<1) eta_bar_N_scan = 1
  if (sigma_initial_N_scan<1) sigma_initial_N_scan = 1

  ! Set up eta_bar values for scan:
  allocate(eta_bar_values(eta_bar_N_scan))
  if (trim(eta_bar_scan_option)==eta_bar_scan_option_2_sided_log) then
     if (mod(eta_bar_N_scan,2)==1) stop "Error! When eta_bar_scan_type='2_sided_log', eta_bar_N_scan must be even."
     do j_eta_bar = 1, eta_bar_N_scan/2
        eta_bar_values(j_eta_bar + eta_bar_N_scan/2)   =  exp(log(eta_bar_min) + (log(eta_bar_max) - log(eta_bar_min)) * (j_eta_bar - 1) / max(eta_bar_N_scan/2 - 1, 1))
        eta_bar_values(eta_bar_N_scan/2 - j_eta_bar+1) = -exp(log(eta_bar_min) + (log(eta_bar_max) - log(eta_bar_min)) * (j_eta_bar - 1) / max(eta_bar_N_scan/2 - 1, 1))
     end do
  else
     do j_eta_bar = 1, eta_bar_N_scan
        select case (trim(eta_bar_scan_option))
        case (eta_bar_scan_option_linear)
           eta_bar_values(j_eta_bar) = eta_bar_min + (eta_bar_max - eta_bar_min) * (j_eta_bar - 1) / max(eta_bar_N_scan - 1, 1)
        case (eta_bar_scan_option_log)
           eta_bar_values(j_eta_bar) = exp(log(eta_bar_min) + (log(eta_bar_max) - log(eta_bar_min)) * (j_eta_bar - 1) / max(eta_bar_N_scan - 1, 1))
        case default
           print *,"Error! Invalid eta_bar_scan_option:", eta_bar_scan_option
           stop
        end select
     end do
  end if
  if (proc0) print *,"eta_bar_values:",eta_bar_values

  ! Set up sigma_initial values for scan:
  allocate(sigma_initial_values(sigma_initial_N_scan))
  if (trim(sigma_initial_scan_option)==sigma_initial_scan_option_2_sided_log) then
     if (mod(sigma_initial_N_scan,2)==0) stop "Error! When sigma_initial_scan_type='2_sided_log', sigma_initial_N_scan must be odd."
     sigma_initial_values = 0
     do j_sigma_initial = 1, sigma_initial_N_scan/2
        sigma_initial_values(j_sigma_initial + sigma_initial_N_scan/2 + 1) =  exp(log(sigma_initial_min) + (log(sigma_initial_max) - log(sigma_initial_min)) * (j_sigma_initial - 1) / max(sigma_initial_N_scan/2 - 1, 1))
        sigma_initial_values(sigma_initial_N_scan/2 - j_sigma_initial + 1) = -exp(log(sigma_initial_min) + (log(sigma_initial_max) - log(sigma_initial_min)) * (j_sigma_initial - 1) / max(sigma_initial_N_scan/2 - 1, 1))
     end do
  else
     do j_sigma_initial = 1, sigma_initial_N_scan
        select case (trim(sigma_initial_scan_option))
        case (sigma_initial_scan_option_linear)
           sigma_initial_values(j_sigma_initial) = sigma_initial_min + (sigma_initial_max - sigma_initial_min) * (j_sigma_initial - 1) / max(sigma_initial_N_scan - 1, 1)
        case (sigma_initial_scan_option_log)
           sigma_initial_values(j_sigma_initial) = exp(log(sigma_initial_min) + (log(sigma_initial_max) - log(sigma_initial_min)) * (j_sigma_initial - 1) / max(sigma_initial_N_scan - 1, 1))
        case default
           print *,"Error! Invalid sigma_initial_scan_option:", sigma_initial_scan_option
           stop
        end select
     end do
  end if
  if (proc0) print *,"sigma_initial_values:",sigma_initial_values

  !N_scan = product(R0s_N_scan)*product(R0c_N_scan)*product(Z0s_N_scan)*product(Z0c_N_scan)*product(eta_bar_N_scan)
  N_scan_array_long = N_scan_array ! Cast to long int
  N_Fourier_scan = product(N_scan_array_long)
  N_scan = N_Fourier_scan * eta_bar_N_scan * sigma_initial_N_scan
  if (proc0) print *,"N_Fourier_scan:",N_Fourier_scan,"N_scan:",N_scan

  !print *,"N_procs=",N_procs, "  my_rank=",my_rank

  if (N_scan >= N_procs) then
     scan_index_min = 1 + (mpi_rank * N_scan) / N_procs
     scan_index_max = ((mpi_rank+1) * N_scan) / N_procs
  else
     ! An uncommon case: There are more procs than solves to do
     scan_index_min = min(mpi_rank+1, N_scan)
     if (mpi_rank < N_scan) then
        scan_index_max = scan_index_min
     else
        scan_index_max = scan_index_min - 1
     end if
  end if
  N_scan_local = scan_index_max - scan_index_min + 1
  print_summary_stride = max(N_scan/N_procs/20, 1000)

  write (proc_assignments_string,fmt="(a,i5,a,i5,a,i20,a,i20)") " Proc ",mpi_rank," of",N_procs," will handle solves",scan_index_min," to",scan_index_max

  ! Print the processor/radius assignments in a coordinated manner.
  !dummy = 0
  !tag = 0
  if (proc0) then
     print "(a)",trim(proc_assignments_string)
     do j = 1,N_procs - 1
        ! To avoid a disordered flood of messages to the masterProc,
        ! ping each proc 1 at a time by sending a dummy value:
        !call MPI_SEND(dummy,1,MPI_INT,j,tag,MPI_COMM_WORLD,ierr)
        ! Now receive the message from proc i:
        !call MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,j,MPI_ANY_TAG,MPI_COMM_WORLD,mpi_status,ierr)
        tag = j
        call MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,j,tag,MPI_COMM_WORLD,mpi_status,ierr)
        print "(a)",trim(proc_assignments_string)
     end do
  else
     !! First, wait for the dummy message from proc 0:
     !call MPI_RECV(dummy,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,mpi_status,ierr)
     ! Now send the message to proc 0:
     tag = mpi_rank 
     call MPI_SEND(proc_assignments_string,buffer_length,MPI_CHAR,0,tag,MPI_COMM_WORLD,ierr)
  end if

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  allocate(iotas_local(N_scan_local))
  allocate(max_elongations_local(N_scan_local))
  allocate(rms_curvatures_local(N_scan_local))
  allocate(max_curvatures_local(N_scan_local))
  allocate(max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(N_scan_local))
  allocate(axis_lengths_local(N_scan_local))
  allocate(standard_deviations_of_R_local(N_scan_local))
  allocate(standard_deviations_of_Z_local(N_scan_local))
  allocate(axis_helicities_local(N_scan_local))
  allocate(B_helicities_local(N_scan_local))
  allocate(effective_nfps_local(N_scan_local))
  allocate(iota_tolerance_achieveds_local(N_scan_local))
  allocate(elongation_tolerance_achieveds_local(N_scan_local))
  allocate(Newton_tolerance_achieveds_local(N_scan_local))

  allocate(scan_eta_bar_local(N_scan_local))
  allocate(scan_sigma_initial_local(N_scan_local))
  allocate(scan_R0c_local(N_scan_local,axis_nmax+1))
  allocate(scan_R0s_local(N_scan_local,axis_nmax+1))
  allocate(scan_Z0c_local(N_scan_local,axis_nmax+1))
  allocate(scan_Z0s_local(N_scan_local,axis_nmax+1))

  scan_state = 1
  j_scan = 0
  j_scan_local = 0
  do j_Fourier_scan = 1, N_Fourier_scan
!!$     print *,"scan_state:"
!!$     do k = 1,4
!!$        print *,scan_state(:,k)
!!$     end do

     ! Set Fourier amplitudes for axis:
     select case (trim(Fourier_scan_option))
     case (Fourier_scan_option_linear)
        do j = 1, axis_nmax+1
           R0s(j) = R0s_min(j) + (R0s_max(j) - R0s_min(j)) * (scan_state(j,1) - 1) / max(N_scan_array(j,1) - 1, 1)
           R0c(j) = R0c_min(j) + (R0c_max(j) - R0c_min(j)) * (scan_state(j,2) - 1) / max(N_scan_array(j,2) - 1, 1)
           Z0s(j) = Z0s_min(j) + (Z0s_max(j) - Z0s_min(j)) * (scan_state(j,3) - 1) / max(N_scan_array(j,3) - 1, 1)
           Z0c(j) = Z0c_min(j) + (Z0c_max(j) - Z0c_min(j)) * (scan_state(j,4) - 1) / max(N_scan_array(j,4) - 1, 1)
        end do
     case (Fourier_scan_option_2_sided_log)
        do j = 1, axis_nmax+1
           R0s(j) = evaluate_2_sided_log(R0s_min(j), R0s_max(j), scan_state(j,1), N_scan_array(j,1))
           R0c(j) = evaluate_2_sided_log(R0c_min(j), R0c_max(j), scan_state(j,2), N_scan_array(j,2))
           Z0s(j) = evaluate_2_sided_log(Z0s_min(j), Z0s_max(j), scan_state(j,3), N_scan_array(j,3))
           Z0c(j) = evaluate_2_sided_log(Z0c_min(j), Z0c_max(j), scan_state(j,4), N_scan_array(j,4))
        end do
     case (Fourier_scan_option_2_sided_log_except_Z0s1)
        do j = 1, axis_nmax+1
           R0s(j) = evaluate_2_sided_log(R0s_min(j), R0s_max(j), scan_state(j,1), N_scan_array(j,1))
           R0c(j) = evaluate_2_sided_log(R0c_min(j), R0c_max(j), scan_state(j,2), N_scan_array(j,2))
           Z0s(j) = evaluate_2_sided_log(Z0s_min(j), Z0s_max(j), scan_state(j,3), N_scan_array(j,3))
           Z0c(j) = evaluate_2_sided_log(Z0c_min(j), Z0c_max(j), scan_state(j,4), N_scan_array(j,4))
        end do
        ! Correct Z0s(2) (i.e. n=1)
        if (scan_state(2,3)==1) then
           Z0s(2) = 0
        else
           !Z0s(2) = Z0s_min(2) + (Z0s_max(2) - Z0s_min(2)) * (scan_state(2,3) - 1) / max(N_scan_array(2,3) - 1, 1)
           Z0s(2) = exp(log(Z0s_min(2)) + (log(Z0s_max(2)) - log(Z0s_min(2))) * (scan_state(2,3) - 2) / max(N_scan_array(2,3) - 2, 1))
        end if
     case default
        stop "Error! Unrecognized Fourier_scan_option"
     end select

     ! Compute the effective number of field periods
     thresh = 1.0d-12
     do j = 1, axis_nmax
        nonzero_modes(j) = abs(R0s(j+1))>thresh .or. abs(R0c(j+1))>thresh .or. abs(Z0s(j+1))>thresh .or. abs(Z0c(j+1))>thresh
     end do
     !print *,"nonzero_modes:",nonzero_modes
     if (any(nonzero_modes)) then
        effective_nfp = nfp
        do j=axis_nmax, 2, -1 ! Let's see if the effective nfp is nfp*j.
           !print *,"Trying effective_nfp=nfp*",j
           keep_going = .true.
           do k=1, axis_nmax
              if (mod(k,j).ne.0 .and. nonzero_modes(k)) then
                 keep_going = .false.
                 exit
              end if
           end do
           if (keep_going) then
              effective_nfp = nfp*j
              exit
           end if
        end do
     else
        effective_nfp = 1
     end if
     !print *,"Found effective_nfp=",effective_nfp

     do j_sigma_initial = 1, sigma_initial_N_scan
        sigma_initial = sigma_initial_values(j_sigma_initial)
        do j_eta_bar = 1, eta_bar_N_scan
           eta_bar = eta_bar_values(j_eta_bar)
           
           j_scan = j_scan + 1
           
           ! Only proceed on the appropriate proc
           if (j_scan < scan_index_min) cycle
           if (j_scan > scan_index_max) cycle  ! Could replace with goto?
           
           if (verbose) then
              print "(a)"," ###################################################################################"
              print "(a,i20,a,i20)"," Scan case",j_scan," of",N_scan
              print *,"eta_bar =",eta_bar
              print *,"sigma_initial = ",sigma_initial
              print *,"R0s:",R0s(1:axis_nmax+1)
              print *,"R0c:",R0c(1:axis_nmax+1)
              print *,"Z0s:",Z0s(1:axis_nmax+1)
              print *,"Z0c:",Z0c(1:axis_nmax+1)
              print *,"nonzero_modes:",nonzero_modes
              print *,"effective_nfp:",effective_nfp
           end if
           
           if (trim(verbose_option)==verbose_option_summary .and. mod(j_scan - scan_index_min + 1,print_summary_stride)==1) then
              print "(a,i4,a,i20,a,i20)", " Proc",mpi_rank,": solve",j_scan - scan_index_min + 1," of",N_scan_local
           end if
           
           if (consider_only_nfp .and. (nfp .ne. effective_nfp)) then
              if (verbose) print *,"Skipping this case since effective_nfp=",effective_nfp
              cycle ! Skip this set of Fourier amplitudes
           end if

!           call quasisymmetry_single_solve()
           call QSC()      !6/4/19.(9u6)
           if (skipped_solve) cycle ! In case R0 <= 0 or some other reason caused quasisymmetry_single_solve to exit prematurely.
           if (max_elongation > max_elongation_to_keep) cycle
           if (abs(iota) < min_iota_to_keep) cycle
           if (max_modBinv_sqrt_half_grad_B_colon_grad_B > max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep) cycle
           
           ! If we made it this far, then record the results
           j_scan_local = j_scan_local + 1
           
           scan_eta_bar_local(j_scan_local) = eta_bar
           scan_sigma_initial_local(j_scan_local) = sigma_initial
           scan_R0c_local(j_scan_local,:) = R0c(1:axis_nmax+1)
           scan_R0s_local(j_scan_local,:) = R0s(1:axis_nmax+1)
           scan_Z0c_local(j_scan_local,:) = Z0c(1:axis_nmax+1)
           scan_Z0s_local(j_scan_local,:) = Z0s(1:axis_nmax+1)
           
           iotas_local(j_scan_local) = iota
           max_elongations_local(j_scan_local) = max_elongation
           rms_curvatures_local(j_scan_local) = rms_curvature
           max_curvatures_local(j_scan_local) = max_curvature
           max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(j_scan_local) = max_modBinv_sqrt_half_grad_B_colon_grad_B
           axis_lengths_local(j_scan_local) = axis_length
           standard_deviations_of_R_local(j_scan_local) = standard_deviation_of_R
           standard_deviations_of_Z_local(j_scan_local) = standard_deviation_of_Z
           axis_helicities_local(j_scan_local) = axis_helicity
           B_helicities_local(j_scan_local) = B_helicity
           effective_nfps_local(j_scan_local) = effective_nfp
           iota_tolerance_achieveds_local(j_scan_local) = iota_tolerance_achieved
           elongation_tolerance_achieveds_local(j_scan_local) = elongation_tolerance_achieved
           Newton_tolerance_achieveds_local(j_scan_local) = Newton_tolerance_achieved
        end do
     end do
  
     ! Update scan state for the next solve
     k = 1
     !do k = 1,4
     do while (k <= 4)
        j = 1
        !do j = 1,max_axis_nmax+1
        do while (j <= max_axis_nmax+1)
           if (scan_state(j,k) < N_scan_array(j,k)) then
              scan_state(j,k) = scan_state(j,k) + 1
              ! exit both j and k loops
              k = 4
              j = max_axis_nmax+1
           else
              scan_state(j,k) = 1
           end if
           j = j + 1
        end do
        k = k + 1
     end do
  end do

  call cpu_time(time2)
  print "(a,i5,a,es10.3,a)"," Proc",mpi_rank," finished after",time2 - start_time," sec."

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! So proc0 does not start printing results until all procs have finished.

  call cpu_time(time1)

  ! Finally, send all results to proc 0.
  if (proc0) then
     if (N_procs > 1) then
        print "(a)"," ###################################################################################"
        print *,"Transferring results to proc 0"
     end if

     allocate(N_solves_kept(N_procs))
     N_solves_kept(1) = j_scan_local
     do j = 1, N_procs-1
        ! Use mpi_rank as the tag
        call mpi_recv(dummy_long,1,MPI_INTEGER8,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        N_solves_kept(j+1) = dummy_long(1)
     end do
     print *,"# solves kept on each proc:",N_solves_kept
     N_scan = sum(N_solves_kept)
     print *,"Total # of solves kept:",N_scan

     ! Now that we know the total number of runs that were kept, we can allocate the arrays for the final results:
     allocate(iotas(N_scan))
     allocate(max_elongations(N_scan))
     allocate(rms_curvatures(N_scan))
     allocate(max_curvatures(N_scan))
     allocate(max_modBinv_sqrt_half_grad_B_colon_grad_Bs(N_scan))
     allocate(axis_lengths(N_scan))
     allocate(standard_deviations_of_R(N_scan))
     allocate(standard_deviations_of_Z(N_scan))
     allocate(axis_helicities(N_scan))
     allocate(B_helicities(N_scan))
     allocate(effective_nfps(N_scan))
     allocate(iota_tolerance_achieveds(N_scan))
     allocate(elongation_tolerance_achieveds(N_scan))
     allocate(Newton_tolerance_achieveds(N_scan))

     allocate(scan_eta_bar(N_scan))
     allocate(scan_sigma_initial(N_scan))
     allocate(scan_R0c(N_scan,axis_nmax+1))
     allocate(scan_R0s(N_scan,axis_nmax+1))
     allocate(scan_Z0c(N_scan,axis_nmax+1))
     allocate(scan_Z0s(N_scan,axis_nmax+1))

     ! Store results from proc0 in the final arrays:
     iotas(1:N_solves_kept(1)) = iotas_local(1:N_solves_kept(1))
     max_elongations(1:N_solves_kept(1)) = max_elongations_local(1:N_solves_kept(1))
     rms_curvatures(1:N_solves_kept(1)) = rms_curvatures_local(1:N_solves_kept(1))
     max_curvatures(1:N_solves_kept(1)) = max_curvatures_local(1:N_solves_kept(1))
     max_modBinv_sqrt_half_grad_B_colon_grad_Bs(1:N_solves_kept(1)) = max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(1:N_solves_kept(1))
     axis_lengths(1:N_solves_kept(1)) = axis_lengths_local(1:N_solves_kept(1))
     standard_deviations_of_R(1:N_solves_kept(1)) = standard_deviations_of_R_local(1:N_solves_kept(1))
     standard_deviations_of_Z(1:N_solves_kept(1)) = standard_deviations_of_Z_local(1:N_solves_kept(1))
     axis_helicities(1:N_solves_kept(1)) = axis_helicities_local(1:N_solves_kept(1))
     B_helicities(1:N_solves_kept(1)) = B_helicities_local(1:N_solves_kept(1))
     effective_nfps(1:N_solves_kept(1)) = effective_nfps_local(1:N_solves_kept(1))
     iota_tolerance_achieveds(1:N_solves_kept(1)) = iota_tolerance_achieveds_local(1:N_solves_kept(1))
     elongation_tolerance_achieveds(1:N_solves_kept(1)) = elongation_tolerance_achieveds_local(1:N_solves_kept(1))
     Newton_tolerance_achieveds(1:N_solves_kept(1)) = Newton_tolerance_achieveds_local(1:N_solves_kept(1))

     scan_eta_bar(1:N_solves_kept(1)) = scan_eta_bar_local(1:N_solves_kept(1))
     scan_sigma_initial(1:N_solves_kept(1)) = scan_sigma_initial_local(1:N_solves_kept(1))
     scan_R0c(1:N_solves_kept(1),:) = scan_R0c_local(1:N_solves_kept(1),:)
     scan_R0s(1:N_solves_kept(1),:) = scan_R0s_local(1:N_solves_kept(1),:)
     scan_Z0c(1:N_solves_kept(1),:) = scan_Z0c_local(1:N_solves_kept(1),:)
     scan_Z0s(1:N_solves_kept(1),:) = scan_Z0s_local(1:N_solves_kept(1),:)

     index = N_solves_kept(1) + 1
     do j = 1, N_procs-1
        print "(a,i20,a,i4)"," Proc 0 is receiving results from ",N_solves_kept(j+1)," solves on proc",j
        call mpi_recv(iotas(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_elongations(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(rms_curvatures(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_curvatures(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_modBinv_sqrt_half_grad_B_colon_grad_Bs(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(axis_lengths(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(standard_deviations_of_R(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(standard_deviations_of_Z(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(axis_helicities(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(B_helicities(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(effective_nfps(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(Newton_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(iota_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(elongation_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)

        call mpi_recv(scan_eta_bar(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_sigma_initial(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_R0c(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_R0s(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_Z0c(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_Z0s(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        !do k = 1,max_axis_nmax+1
        !   call mpi_recv(scan_R0c(index:index+N_solves_kept(j+1)-1,k),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        !end do

        index = index + N_solves_kept(j+1)
     end do


  else
     ! Send the number of runs this proc kept:
     dummy_long(1) = j_scan_local
     ! Use mpi_rank as the tag
     call mpi_send(dummy_long,1,MPI_INTEGER8,0,mpi_rank,MPI_COMM_WORLD,ierr)

     ! Send the other results:
     call mpi_send(iotas_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_elongations_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(rms_curvatures_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_curvatures_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(axis_lengths_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(standard_deviations_of_R_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(standard_deviations_of_Z_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(axis_helicities_local(1:j_scan_local),j_scan_local,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(B_helicities_local(1:j_scan_local),j_scan_local,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(effective_nfps_local(1:j_scan_local),j_scan_local,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(Newton_tolerance_achieveds_local(1:j_scan_local),j_scan_local,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(iota_tolerance_achieveds_local(1:j_scan_local),j_scan_local,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(elongation_tolerance_achieveds_local(1:j_scan_local),j_scan_local,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     
     call mpi_send(scan_eta_bar_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_sigma_initial_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_R0c_local(1:j_scan_local,:),j_scan_local*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_R0s_local(1:j_scan_local,:),j_scan_local*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_Z0c_local(1:j_scan_local,:),j_scan_local*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_Z0s_local(1:j_scan_local,:),j_scan_local*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     !do k = 1,max_axis_nmax+1
     !   call mpi_send(scan_R0c_local(1:j_scan_local,k),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     !end do
  end if

  if (proc0) then
     call cpu_time(time2)
     print "(a,es10.3,a)"," Time for communication:",time2 - time1," sec."
     print "(a)"," ###################################################################################"
     print *,"Scan complete."
     if (N_scan < 5000) then
        print "(a,99999(f8.2))"," iotas:",iotas
        print *," "
        print "(a,99999(f8.1))"," elongations:",max_elongations
        print *," "
        print "(a,99999(f8.2))"," rms_curvatures:",rms_curvatures
        print *," "
        print "(a,99999(f8.2))"," max_curvatures:",max_curvatures
        print *," "
        print "(a,99999(f8.2))"," max_modBinv_sqrt_half_grad_B_colon_grad_Bs:",max_modBinv_sqrt_half_grad_B_colon_grad_Bs
        print *," "
        print "(a,99999(f8.2))"," axis_lengths:",axis_lengths
        print *," "
        print "(a,99999(i2))","                 effective_nfps:",effective_nfps
        print *," "
        print "(a,99999(i2))","                axis_helicities:",axis_helicities
        print "(a,99999(i2))","                   B_helicities:",B_helicities
        print *," axis_helicities==B_helicities:",B_helicities==axis_helicities
        print *,"    Newton_tolerance_achieveds:",Newton_tolerance_achieveds
        print *,"      iota_tolerance_achieveds:",iota_tolerance_achieveds
        print *,"elongation_tolerance_achieveds:",elongation_tolerance_achieveds
     end if
  end if

contains

  real(dp) function evaluate_2_sided_log(xmin,xmax,j,N)
    real(dp), intent(in) :: xmin, xmax
    integer, intent(in) :: j, N
    
    if (N<2) then
       evaluate_2_sided_log = xmin ! Note this is not 0. This is so R0c(1) is not 0.
       return
    end if
    if (mod(N,2)==0) stop "Error! When Fourier_scan_option='Fourier_scan_option_2_sided_log', all axis N_scan parameters must be odd or 0."

    if (j == (N+1)/2) then
       evaluate_2_sided_log = 0
    elseif (j > (N+1)/2) then
       evaluate_2_sided_log =  exp(log(xmin) + (log(xmax) - log(xmin)) * (j - (N+1)/2 - 1) / max((N-1)/2 - 1, 1))
    else
       evaluate_2_sided_log = -exp(log(xmin) + (log(xmax) - log(xmin)) * ((N+1)/2 - j - 1) / max((N-1)/2 - 1, 1))
    end if

  end function evaluate_2_sided_log

end subroutine quasisymmetry_scan

