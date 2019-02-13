subroutine quasisymmetry_scan

  use quasisymmetry_variables

  implicit none

  include 'mpif.h'

  integer :: j_scan, j_scan_local, j_Fourier_scan, j, k, N_Fourier_scan, j_B1s, j_B1c, j_sigma_initial, index
  integer, dimension(max_axis_nmax+1, 4) :: scan_state
  integer :: ierr, tag, dummy(1), scan_index_min, scan_index_max, N_scan_local
  integer :: mpi_status(MPI_STATUS_SIZE)
  integer, parameter :: buffer_length = 100
  character(len=buffer_length) :: proc_assignments_string
  real(dp), dimension(:), allocatable :: iotas_local, max_elongations_local
  integer, dimension(:), allocatable :: helicities_local
  logical, dimension(:), allocatable :: iota_tolerance_achieveds_local, elongation_tolerance_achieveds_local, Newton_tolerance_achieveds_local
  integer, dimension(:), allocatable :: N_solves_kept
  real(dp), dimension(:), allocatable :: scan_B1c_local, scan_B1s_local, scan_sigma_initial_local
  real(dp), dimension(:,:), allocatable :: scan_R0c_local, scan_R0s_local, scan_Z0c_local, scan_Z0s_local

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
  if (B1s_N_scan<1) B1s_N_scan = 1
  if (B1c_N_scan<1) B1c_N_scan = 1
  if (sigma_initial_N_scan<1) sigma_initial_N_scan = 1

  allocate(B1c_values(B1c_N_scan))
  do j_B1c = 1, B1c_N_scan
     select case (trim(B1c_scan_option))
     case (B1c_scan_option_linear)
        B1c_values(j_B1c) = B1c_min + (B1c_max - B1c_min) * (j_B1c - 1) / max(B1c_N_scan - 1, 1)
     case (B1c_scan_option_log)
        B1c_values(j_B1c) = exp(log(B1c_min) + (log(B1c_max) - log(B1c_min)) * (j_B1c - 1) / max(B1c_N_scan - 1, 1))
     case default
        print *,"Error! Invalid B1c_scan_option:", B1c_scan_option
        stop
     end select
  end do
  if (proc0) print *,"B1c_values:",B1c_values

  !N_scan = product(R0s_N_scan)*product(R0c_N_scan)*product(Z0s_N_scan)*product(Z0c_N_scan)*product(B1s_N_scan)*product(B1c_N_scan)
  N_Fourier_scan = product(N_scan_array)
  N_scan = N_Fourier_scan * B1s_N_scan * B1c_N_scan * sigma_initial_N_scan
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

  write (proc_assignments_string,fmt="(a,i5,a,i5,a,i9,a,i9)") "Proc ",mpi_rank," of",N_procs," will handle solves",scan_index_min," to",scan_index_max

  ! Print the processor/radius assignments in a coordinated manner.
  dummy = 0
  tag = 0
  if (proc0) then
     print *,trim(proc_assignments_string)
     do j = 1,N_procs - 1
        ! To avoid a disordered flood of messages to the masterProc,
        ! ping each proc 1 at a time by sending a dummy value:
        call MPI_SEND(dummy,1,MPI_INT,j,tag,MPI_COMM_WORLD,ierr)
        ! Now receive the message from proc i:
        call MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,j,MPI_ANY_TAG,MPI_COMM_WORLD,mpi_status,ierr)
        print *,trim(proc_assignments_string)
     end do
  else
     ! First, wait for the dummy message from proc 0:
     call MPI_RECV(dummy,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,mpi_status,ierr)
     ! Now send the message to proc 0:
     call MPI_SEND(proc_assignments_string,buffer_length,MPI_CHAR,0,tag,MPI_COMM_WORLD,ierr)
  end if

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  allocate(iotas_local(N_scan_local))
  allocate(max_elongations_local(N_scan_local))
  allocate(helicities_local(N_scan_local))
  allocate(iota_tolerance_achieveds_local(N_scan_local))
  allocate(elongation_tolerance_achieveds_local(N_scan_local))
  allocate(Newton_tolerance_achieveds_local(N_scan_local))

  allocate(scan_B1c_local(N_scan_local))
  allocate(scan_B1s_local(N_scan_local))
  allocate(scan_sigma_initial_local(N_scan_local))
  allocate(scan_R0c_local(N_scan_local,max_axis_nmax+1))
  allocate(scan_R0s_local(N_scan_local,max_axis_nmax+1))
  allocate(scan_Z0c_local(N_scan_local,max_axis_nmax+1))
  allocate(scan_Z0s_local(N_scan_local,max_axis_nmax+1))

  scan_state = 1
  j_scan = 0
  j_scan_local = 0
  do j_Fourier_scan = 1, N_Fourier_scan
!!$     print *,"scan_state:"
!!$     do k = 1,4
!!$        print *,scan_state(:,k)
!!$     end do

     ! Set Fourier amplitudes for axis:
     do j = 1, max_axis_nmax+1
        R0s(j) = R0s_min(j) + (R0s_max(j) - R0s_min(j)) * (scan_state(j,1) - 1) / max(N_scan_array(j,1) - 1, 1)
        R0c(j) = R0c_min(j) + (R0c_max(j) - R0c_min(j)) * (scan_state(j,2) - 1) / max(N_scan_array(j,2) - 1, 1)
        Z0s(j) = Z0s_min(j) + (Z0s_max(j) - Z0s_min(j)) * (scan_state(j,3) - 1) / max(N_scan_array(j,3) - 1, 1)
        Z0c(j) = Z0c_min(j) + (Z0c_max(j) - Z0c_min(j)) * (scan_state(j,4) - 1) / max(N_scan_array(j,4) - 1, 1)
     end do

     do j_sigma_initial = 1, sigma_initial_N_scan
        sigma_initial = sigma_initial_min + (sigma_initial_max - sigma_initial_min) * (j_sigma_initial - 1) / max(sigma_initial_N_scan - 1, 1)
        do j_B1s = 1, B1s_N_scan
           B1s_over_B0 = B1s_min + (B1s_max - B1s_min) * (j_B1s - 1) / max(B1s_N_scan - 1, 1)
           do j_B1c = 1, B1c_N_scan
              B1c_over_B0 = B1c_values(j_B1c)

              j_scan = j_scan + 1

              ! Only proceed on the appropriate proc
              if (j_scan < scan_index_min) cycle
              if (j_scan > scan_index_max) cycle  ! Could replace with goto?

              if (verbose) then
                 print "(a)"," ###################################################################################"
                 print "(a,i7,a,i7)"," Scan case",j_scan," of",N_scan
                 print *,"B1s =", B1s_over_B0, "   B1c =",B1c_over_B0
                 print *,"sigma_initial = ",sigma_initial
                 print *,"R0s:",R0s(1:axis_nmax+1)
                 print *,"R0c:",R0c(1:axis_nmax+1)
                 print *,"Z0s:",Z0s(1:axis_nmax+1)
                 print *,"Z0c:",Z0c(1:axis_nmax+1)
              end if
              
              if (trim(verbose_option)==verbose_option_summary .and. mod(j_scan - scan_index_min + 1,1000)==1) then
                 print "(a,i4,a,i9,a,i9)", " Proc",mpi_rank,": solve",j_scan - scan_index_min + 1," of",N_scan_local
              end if
              
!              call quasisymmetry_single_solve()
              call QSC()
              if (max_elongation > max_elongation_to_keep) cycle
              !if (abs(iota) < 1e-3) cycle

              ! If we made it this far, then record the results
              j_scan_local = j_scan_local + 1

              scan_B1c_local(j_scan_local) = B1c_over_B0
              scan_B1s_local(j_scan_local) = B1s_over_B0
              scan_sigma_initial_local(j_scan_local) = sigma_initial
              scan_R0c_local(j_scan_local,:) = R0c
              scan_R0s_local(j_scan_local,:) = R0s
              scan_Z0c_local(j_scan_local,:) = Z0c
              scan_Z0s_local(j_scan_local,:) = Z0s

              iotas_local(j_scan_local) = iota
              max_elongations_local(j_scan_local) = max_elongation
              helicities_local(j_scan_local) = helicity
              iota_tolerance_achieveds_local(j_scan_local) = iota_tolerance_achieved
              elongation_tolerance_achieveds_local(j_scan_local) = elongation_tolerance_achieved
              Newton_tolerance_achieveds_local(j_scan_local) = Newton_tolerance_achieved
           end do
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

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! So proc0 does not start printing results until all procs have finished.

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
        call mpi_recv(dummy,1,MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        N_solves_kept(j+1) = dummy(1)
     end do
     print *,"# solves kept on each proc:",N_solves_kept
     N_scan = sum(N_solves_kept)
     print *,"Total # of solves kept:",N_scan

     ! Now that we know the total number of runs that were kept, we can allocate the arrays for the final results:
     allocate(iotas(N_scan))
     allocate(max_elongations(N_scan))
     allocate(helicities(N_scan))
     allocate(iota_tolerance_achieveds(N_scan))
     allocate(elongation_tolerance_achieveds(N_scan))
     allocate(Newton_tolerance_achieveds(N_scan))

     allocate(scan_B1c(N_scan))
     allocate(scan_B1s(N_scan))
     allocate(scan_sigma_initial(N_scan))
     allocate(scan_R0c(N_scan,max_axis_nmax+1))
     allocate(scan_R0s(N_scan,max_axis_nmax+1))
     allocate(scan_Z0c(N_scan,max_axis_nmax+1))
     allocate(scan_Z0s(N_scan,max_axis_nmax+1))

     ! Store results from proc0 in the final arrays:
     iotas(1:N_solves_kept(1)) = iotas_local(1:N_solves_kept(1))
     max_elongations(1:N_solves_kept(1)) = max_elongations_local(1:N_solves_kept(1))
     helicities(1:N_solves_kept(1)) = helicities_local(1:N_solves_kept(1))
     iota_tolerance_achieveds(1:N_solves_kept(1)) = iota_tolerance_achieveds_local(1:N_solves_kept(1))
     elongation_tolerance_achieveds(1:N_solves_kept(1)) = elongation_tolerance_achieveds_local(1:N_solves_kept(1))
     Newton_tolerance_achieveds(1:N_solves_kept(1)) = Newton_tolerance_achieveds_local(1:N_solves_kept(1))

     scan_B1c(1:N_solves_kept(1)) = scan_B1c_local(1:N_solves_kept(1))
     scan_B1s(1:N_solves_kept(1)) = scan_B1s_local(1:N_solves_kept(1))
     scan_sigma_initial(1:N_solves_kept(1)) = scan_sigma_initial_local(1:N_solves_kept(1))
     scan_R0c(1:N_solves_kept(1),:) = scan_R0c_local(1:N_solves_kept(1),:)
     scan_R0s(1:N_solves_kept(1),:) = scan_R0s_local(1:N_solves_kept(1),:)
     scan_Z0c(1:N_solves_kept(1),:) = scan_Z0c_local(1:N_solves_kept(1),:)
     scan_Z0s(1:N_solves_kept(1),:) = scan_Z0s_local(1:N_solves_kept(1),:)

     index = N_solves_kept(1) + 1
     do j = 1, N_procs-1
        print "(a,i8,a,i4)"," Proc 0 is receiving results from ",N_solves_kept(j+1)," solves from proc",j
        call mpi_recv(iotas(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_elongations(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(helicities(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(Newton_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(iota_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(elongation_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)

        call mpi_recv(scan_B1c(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_B1s(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_sigma_initial(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_R0c(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(max_axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_R0s(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(max_axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_Z0c(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(max_axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_Z0s(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(max_axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        !do k = 1,max_axis_nmax+1
        !   call mpi_recv(scan_R0c(index:index+N_solves_kept(j+1)-1,k),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        !end do

        index = index + N_solves_kept(j+1)
     end do


  else
     ! Send the number of runs this proc kept:
     dummy(1) = j_scan_local
     ! Use mpi_rank as the tag
     call mpi_send(dummy,1,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)

     ! Send the other results:
     call mpi_send(iotas_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_elongations_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(helicities_local(1:j_scan_local),j_scan_local,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(Newton_tolerance_achieveds_local(1:j_scan_local),j_scan_local,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(iota_tolerance_achieveds_local(1:j_scan_local),j_scan_local,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(elongation_tolerance_achieveds_local(1:j_scan_local),j_scan_local,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     
     call mpi_send(scan_B1c_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_B1s_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_sigma_initial_local(1:j_scan_local),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_R0c_local(1:j_scan_local,:),j_scan_local*(max_axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_R0s_local(1:j_scan_local,:),j_scan_local*(max_axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_Z0c_local(1:j_scan_local,:),j_scan_local*(max_axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_Z0s_local(1:j_scan_local,:),j_scan_local*(max_axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     !do k = 1,max_axis_nmax+1
     !   call mpi_send(scan_R0c_local(1:j_scan_local,k),j_scan_local,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     !end do
  end if

  if (proc0) then
     print "(a)"," ###################################################################################"
     print *,"Scan complete."
     if (N_scan < 5000) then
        print "(a,99999(f8.2))"," iotas:",iotas
        print *," "
        print "(a,99999(f8.1))"," elongations:",max_elongations
        print *," "
        print "(a,99999(i2))","                     helicities:",helicities
        print *,"    Newton_tolerance_achieveds:",Newton_tolerance_achieveds
        print *,"      iota_tolerance_achieveds:",iota_tolerance_achieveds
        print *,"elongation_tolerance_achieveds:",elongation_tolerance_achieveds
     end if
  end if

end subroutine quasisymmetry_scan
