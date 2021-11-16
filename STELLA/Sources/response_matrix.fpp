module response_matrix

  use netcdf
  use mpi

  implicit none

  public :: init_response_matrix, finish_response_matrix
  public :: read_response_matrix
  public :: response_matrix_initialized

  private

  logical :: response_matrix_initialized = .false.
  integer, parameter :: mat_unit = 70

#if defined MPI && defined ISO_C_BINDING
  integer :: window = MPI_WIN_NULL
#endif

contains

  subroutine init_response_matrix
  
    use linear_solve, only: lu_decomposition
    use fields_arrays, only: response_matrix
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, is_idx
    use kt_grids, only: naky
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: neigen, ikxmod
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use extended_zgrid, only: periodic
    use job_manage, only: time_message
    use mp, only: proc0, job, mp_abort
    use run_parameters, only: mat_gen, lu_option_switch
    use run_parameters, only: lu_option_none, lu_option_local, lu_option_global
    use system_fortran, only: systemf
#if defined MPI && defined ISO_C_BINDING
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use mp, only: curr_focus, sgproc0, mp_comm, sharedsubprocs, scope, barrier
    use mp, only: real_size, nbytes_real
    use mpi
#endif

    implicit none

    integer :: iky, ie, iseg, iz
    integer :: ikx
    integer :: nz_ext, nresponse
    integer :: idx
    integer :: izl_offset, izup
#if defined MPI && defined ISO_C_BINDING
    integer :: prior_focus, ierr
    integer :: disp_unit = 1
    integer (kind=MPI_ADDRESS_KIND) :: win_size
    integer*8 :: cur_pos
    type(c_ptr) :: bptr, cptr
#endif
    real :: dum
    complex, dimension (:), allocatable :: phiext
    complex, dimension (:,:), allocatable :: gext
    logical :: debug = .false.
    character(100) :: message_dgdphi, message_QN, message_lu
    real, dimension (2) :: time_response_matrix_dgdphi
    real, dimension (2) :: time_response_matrix_QN
    real, dimension (2) :: time_response_matrix_lu

    ! Related to the saving of the the matrices in netcdf format
    character(len=15) :: fmt, job_str
    character(len=100) :: file_name
    integer :: istatus
    istatus = 0

    if (proc0.and.debug) then
      write (*,*) " "
      write (*,'(A)') "    ############################################################"
      write (*,'(A)') "                         RESPONSE MATRIX"
      write (*,'(A)') "    ############################################################"
    end if
    message_dgdphi = '     calculate dgdphi: '
    message_QN     = '     calculate QN:     '
    message_lu     = '     calculate LU:     ' 
    time_response_matrix_dgdphi = 0
    time_response_matrix_QN     = 0
    time_response_matrix_lu     = 0


!   All matrices handled by processor i_proc and job are stored
!   on a single file named: response_mat_job.iproc
    fmt = '(I5.5)'
    if (proc0.and.mat_gen) then

       call systemf('mkdir -p mat')

       write (job_str, '(I1.1)') job
       file_name = './mat/response_mat_'//trim(job_str)

       open(unit=mat_unit, status='replace', file=file_name, &
            position='rewind', action='write', form='unformatted')
       write(unit=mat_unit) naky
    end if

    if (response_matrix_initialized) return
    response_matrix_initialized = .true.
  
    if (.not.allocated(response_matrix)) allocate (response_matrix(naky))

#if defined ISO_C_BINDING && defined MPI

!   Create a single shared memory window for all the response matrices and 
!   permutation arrays.
!   Creating a window for each matrix/array would lead to performance
!   degradation on some clusters
    if(window.eq.MPI_WIN_NULL) then
      prior_focus = curr_focus
      call scope(sharedsubprocs)
      win_size = 0
      if(sgproc0) then
        do iky = 1, naky
          do ie = 1, neigen(iky)
            if (periodic(iky)) then
               nresponse = nsegments(ie,iky)*nzed_segment
            else
               nresponse = nsegments(ie,iky)*nzed_segment+1
            end if
            win_size =   win_size &
                       + int(nresponse,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND &
                       + int(nresponse**2,MPI_ADDRESS_KIND)*2*real_size 
          enddo
        enddo
      endif
      call mpi_win_allocate_shared(win_size,disp_unit,MPI_INFO_NULL,mp_comm, &
                                   bptr,window,ierr)

      if(.not.sgproc0) then
        !make sure all the procs have the right memory address
        call mpi_win_shared_query(window, 0, win_size, disp_unit, bptr, ierr)
      endif
      call mpi_win_fence(0,window,ierr)
      
      !the following is a hack that allows us to perform pointer arithmetic in Fortran
      cur_pos = transfer(bptr,cur_pos) 

      call scope(prior_focus)
    endif
#endif
  
    ! for a given ky and set of connected kx values
    ! give a unit impulse to phi at each zed location
    ! in the extended domain and solve for h(zed_extended,(vpa,mu,s))

  
    do iky = 1, naky
     
       if (proc0.and.mat_gen) THEN
          write(unit=mat_unit) iky, neigen(iky)
       end if
     
       ! the response matrix for each ky has neigen(ky)
       ! independent sets of connected kx values
       if (.not.associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

       if(proc0.and.debug) then
         call time_message(.false., time_response_matrix_dgdphi, message_dgdphi)
       end if
     
       ! loop over the sets of connected kx values
       do ie = 1, neigen(iky)
        
          ! number of zeds x number of segments
          nz_ext = nsegments(ie,iky)*nzed_segment+1
        
          ! treat zonal mode specially to avoid double counting
          ! as it is periodic
          if (periodic(iky)) then
             nresponse = nz_ext-1
          else
             nresponse = nz_ext
          end if
        
          if (proc0.and.mat_gen) then
             write(unit=mat_unit) ie, nresponse
          end if

#if !defined ISO_C_BINDING || !defined MPI
          ! for each ky and set of connected kx values,
          ! must have a response matrix that is N x N
          ! with N = number of zeds per 2pi segment x number of 2pi segments
          if (.not.associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse,nresponse))

          ! response_matrix%idx is needed to keep track of permutations
          ! to the response matrix made during LU decomposition
          ! it will be input to LU back substitution during linear solve
          if (.not.associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
#else
          !exploit MPIs shared memory framework to reduce memory consumption of 
          !response matrices

          if (.not.associated(response_matrix(iky)%eigen(ie)%zloc)) then
            cptr = transfer(cur_pos,cptr)
            call c_f_pointer(cptr,response_matrix(iky)%eigen(ie)%zloc,(/nresponse,nresponse/))
            cur_pos = cur_pos + nresponse**2*2*nbytes_real
          endif

          if (.not.associated(response_matrix(iky)%eigen(ie)%idx)) then
            cptr = transfer(cur_pos,cptr)
            call c_f_pointer(cptr,response_matrix(iky)%eigen(ie)%idx,(/nresponse/))
            cur_pos = cur_pos + nresponse*4
          endif
#endif

          allocate (gext(nz_ext,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          allocate (phiext(nz_ext))
          ! idx is the index in the extended zed domain
          ! that we are giving a unit impulse
          idx = 0

          ! loop over segments, starting with 1
          ! first segment is special because it has
          ! one more unique zed value than all others
          ! since domain is [z0-pi:z0+pi], including both endpoints
          ! i.e., one endpoint is shared with the previous segment
          iseg = 1
          ! ikxmod gives the kx corresponding to iseg,ie,iky
          ikx = ikxmod(iseg,ie,iky)
          izl_offset = 0
          ! avoid double-counting of periodic points for zonal mode (and other periodic modes)
          if (periodic(iky)) then
             izup = iz_up(iseg)-1
          else
             izup = iz_up(iseg)
          end if
          ! no need to obtain response to impulses at negative kx values
          do iz = iz_low(iseg), izup
             idx = idx + 1
             call get_dgdphi_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
          end do
          ! once we have used one segments, remaining segments
          ! have one fewer unique zed point
          izl_offset = 1
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                ikx = ikxmod(iseg,ie,iky)
                do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
                   idx = idx + 1
                   call get_dgdphi_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
                end do
                if (izl_offset == 0) izl_offset = 1
             end do
          end if
          deallocate (gext,phiext)
       enddo
       !DSO - This ends parallelization over velocity space.
       !      At this point every processor has int dv dgdphi for a given ky
       !      and so the quasineutrality solve and LU decomposition can be
       !      parallelized locally if need be.
       !      This is preferable to parallelization over ky as the LU
       !      decomposition (and perhaps QN) will be dominated by the
       !      ky with the most connections

       if(proc0.and.debug) then
         call time_message(.true. , time_response_matrix_dgdphi, message_dgdphi)
         call time_message(.false., time_response_matrix_QN, message_QN)
       end if

#ifdef ISO_C_BINDING 
       call mpi_win_fence(0,window,ierr)
#endif

       ! solve quasineutrality
       ! for local stella, this is a diagonal process, but global stella
       ! may require something more sophisticated

       ! loop over the sets of connected kx values
       do ie = 1, neigen(iky)
#if defined ISO_C_BINDING && defined MPI
         if(sgproc0) then
#endif
           ! number of zeds x number of segments
           nz_ext = nsegments(ie,iky)*nzed_segment+1
        
           ! treat zonal mode specially to avoid double counting
           ! as it is periodic
           if (periodic(iky)) then
             nresponse = nz_ext-1
           else
             nresponse = nz_ext
           end if

           allocate (phiext(nz_ext))

           do idx = 1, nresponse
             phiext(nz_ext) = 0.0
             phiext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:,idx)
             call get_fields_for_response_matrix (phiext, iky, ie)

             ! next need to create column in response matrix from phiext
             ! negative sign because matrix to be inverted in streaming equation
             ! is identity matrix - response matrix
             ! add in contribution from identity matrix
             phiext(idx) = phiext(idx)-1.0
             response_matrix(iky)%eigen(ie)%zloc(:,idx) = -phiext(:nresponse)

           end do
           deallocate (phiext)
#if defined ISO_C_BINDING && defined MPI
         endif
#endif
       enddo

#ifdef ISO_C_BINDING       
       call mpi_win_fence(0,window,ierr)
#endif

       if(proc0.and.debug) then
         call time_message(.true. , time_response_matrix_QN, message_QN)
         call time_message(.false., time_response_matrix_lu, message_lu)
       end if

       !now we have the full response matrix. Finally, perform its LU decomposition
#ifdef MPI
       select case (lu_option_switch)
       case (lu_option_global)
        call parallel_LU_decomposition_global(iky)
       case (lu_option_local)
#ifdef ISO_C_BINDING       
         call parallel_LU_decomposition_local(iky)
#else
         call mp_abort('Stella must be built with HAS_ISO_BINDING in order to use local parallel LU decomposition.')
#endif
       case default
#endif
         do ie = 1, neigen(iky)
#ifdef ISO_C_BINDING       
           if(sgproc0) then
#endif
             ! now that we have the reponse matrix for this ky and set of connected kx values
             !get the LU decomposition so we are ready to solve the linear system
             call lu_decomposition (response_matrix(iky)%eigen(ie)%zloc, &
                                    response_matrix(iky)%eigen(ie)%idx,dum)
        
#ifdef ISO_C_BINDING       
           endif
#endif
        enddo
#ifdef MPI
       end select
#endif

     
       if(proc0.and.debug) then
         call time_message(.true., time_response_matrix_lu, message_lu)
       end if

       time_response_matrix_dgdphi = 0
       time_response_matrix_QN     = 0
       time_response_matrix_lu     = 0
     
       do ie = 1, neigen(iky)
          if (proc0.and.mat_gen) then
             write(unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
             write(unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
          end if
       end do

       !if(proc0)  write (*,*) 'job', iky, iproc, response_matrix(iky)%eigen(1)%zloc(5,:)
    end do

#ifdef ISO_C_BINDING       
    call mpi_win_fence(0,window,ierr)
#endif

    if (proc0.and.mat_gen) then
       close(unit=mat_unit)
    end if

    if (proc0.and.debug) then
      write (*,'(A)') "    ############################################################"
      write (*,'(A)') " "
    end if

  end subroutine init_response_matrix

  subroutine read_response_matrix
  
    use fields_arrays, only: response_matrix
    use common_types, only: response_matrix_type
    use kt_grids, only: naky
    use extended_zgrid, only: neigen
    use extended_zgrid, only: nsegments
    use extended_zgrid, only: nzed_segment
    use extended_zgrid, only: periodic
    use mp, only: proc0, job, broadcast, mp_abort

    implicit none

    integer :: iky, ie, nz_ext
    integer :: iky_dump, neigen_dump, naky_dump, nresponse_dump
    integer :: nresponse
    character(len=15) :: fmt, job_str
    character(len=100) :: file_name
    integer :: istatus, ie_dump, istat
    logical, parameter :: debug=.false.
    istatus = 0

!   All matrices handled for the job i_job are read
!   from a single file named: responst_mat.ijob by that
!   jobs root process

    if(proc0) then
      fmt = '(I5.5)'
      write (job_str, '(I1.1)') job
      file_name = './mat/response_mat.'//trim(job_str)

      open(unit=mat_unit, status='old', file=file_name, &
         action='read', form='unformatted', iostat=istat)
      if (istat /= 0) then
         print *, 'Error opening response_matrix by root processor for job ', job_str
      end if

      read(unit=mat_unit) naky_dump
      if(naky.ne.naky_dump) call mp_abort('mismatch in naky and naky_dump')
    endif

    if (.not.allocated(response_matrix)) allocate (response_matrix(naky))

    do iky = 1, naky
       if(proc0) then 
         read(unit=mat_unit) iky_dump, neigen_dump
         if(iky_dump.ne.iky.or.neigen_dump.ne.neigen(iky)) &
           call mp_abort('mismatch in iky_dump/neigen_dump')
       endif

       if (.not.associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

       ! loop over the sets of connected kx values
       do ie = 1, neigen(iky)
          ! number of zeds x number of segments
          nz_ext = nsegments(ie,iky)*nzed_segment+1
        
          ! treat zonal mode specially to avoid double counting
          ! as it is periodic
          if (periodic(iky)) then
             nresponse = nz_ext-1
          else
             nresponse = nz_ext
          end if

          if(proc0) then 
            read(unit=mat_unit) ie_dump, nresponse_dump
            if(ie_dump.ne.ie.or.nresponse.ne.nresponse_dump) &
              call mp_abort('mismatch in ie/nresponse_dump')
          endif

          ! for each ky and set of connected kx values,
          ! must have a response matrix that is N x N
          ! with N = number of zeds per 2pi segment x number of 2pi segments
          if (.not.associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse,nresponse))

          ! response_matrix%idx is needed to keep track of permutations
          ! to the response matrix made during LU decomposition
          ! it will be input to LU back substitution during linear solve
          if (.not.associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
          if(proc0) then
            read(unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
            read(unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
          endif

          call broadcast(response_matrix(iky)%eigen(ie)%idx)
          call broadcast(response_matrix(iky)%eigen(ie)%zloc)

       enddo
    enddo

    if (proc0) close (mat_unit)
  
    if (debug) then
       print *, 'File', file_name, ' successfully read by root proc for job: ', job_str
    end if
  end subroutine read_response_matrix

  subroutine get_dgdphi_matrix_column (iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use zgrid, only: delzed, nzgrid
    use extended_zgrid, only: periodic
    use species, only: spec
    use stella_geometry, only: gradpar, dbdzed
    use vpamu_grids, only: vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use fields_arrays, only: response_matrix
    use gyro_averages, only: aj0x
    use run_parameters, only: driftkinetic_implicit
    use run_parameters, only: maxwellian_inside_zed_derivative
    use parallel_streaming, only: stream_tridiagonal_solve
    use parallel_streaming, only: stream_sign
    use run_parameters, only: zed_upwind, time_upwind
#if defined ISO_C_BINDING && defined MPI
    use mp, only: sgproc0
#endif

    implicit none

    integer, intent (in) :: iky, ikx, iz, ie, idx, nz_ext, nresponse
    complex, dimension (:), intent (in out) :: phiext
    complex, dimension (:,vmu_lo%llim_proc:), intent (in out) :: gext

    integer :: ivmu, iv, imu, is, ia
    integer :: izp, izm
    real :: mu_dbdzed_p, mu_dbdzed_m
    real :: fac, fac0, fac1, gyro_fac

    ia = 1

    if (.not.maxwellian_inside_zed_derivative) then
       ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          ! initialize g to zero everywhere along extended zed domain
          gext(:,ivmu) = 0.0
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
        
          ! give unit impulse to phi at this zed location
          ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)
        
          ! NB:  assuming equal spacing in zed below
          ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*J0/dz
          ! b.gradz left out because needs to be centred in zed
          if (driftkinetic_implicit) then
             gyro_fac = 1.0
          else
             gyro_fac = aj0x(iky,ikx,iz,ivmu)
          end if
        
          ! 0.125 to account for two linear interpolations
          fac = -0.125*(1.+time_upwind)*code_dt*vpa(iv)*spec(is)%stm_psi0 &
               *gyro_fac*spec(is)%zt/delzed(0)*maxwell_vpa(iv,is)*maxwell_fac(is)
          
          ! In the following, gradpar and maxwell_mu are interpolated separately
          ! to ensure consistency to what is done in parallel_streaming.f90

          ! stream_sign < 0 corresponds to positive advection speed
          if (stream_sign(iv)<0) then
             if (iz > -nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                            + (1.-zed_upwind)*gradpar(iz-1)) &
                          *((1.+zed_upwind)*maxwell_mu(ia,iz,imu,is) &
                            + (1.-zed_upwind)*maxwell_mu(ia,iz-1,imu,is))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                if (iz < nzgrid) then
                   fac1 = fac*((1.+zed_upwind)*gradpar(iz+1) &
                               + (1.-zed_upwind)*gradpar(iz)) &
                             *((1.+zed_upwind)*maxwell_mu(ia,iz+1,imu,is) &
                               + (1.-zed_upwind)*maxwell_mu(ia,iz,imu,is))
                else
                   fac1 = fac*((1.+zed_upwind)*gradpar(-nzgrid+1) &
                               + (1.-zed_upwind)*gradpar(nzgrid)) &
                             *((1.+zed_upwind)*maxwell_mu(ia,-nzgrid+1,imu,is) &
                               + (1.-zed_upwind)*maxwell_mu(ia,nzgrid,imu,is))
                end if
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                            + (1.-zed_upwind)*gradpar(nzgrid-1)) &
                          *((1.+zed_upwind)*maxwell_mu(ia,iz,imu,is) &
                            + (1.-zed_upwind)*maxwell_mu(ia,nzgrid-1,imu,is))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar(iz+1)  &
                            + (1.-zed_upwind)*gradpar(iz)) &
                          *((1.+zed_upwind)*maxwell_mu(ia,iz+1,imu,is) &
                            + (1.-zed_upwind)*maxwell_mu(ia,iz,imu,is))
             end if
             gext(idx,ivmu) = fac0
             if (idx < nz_ext) gext(idx+1,ivmu) = -fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (periodic(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = fac0
                else if (idx == nz_ext-1) then
                   gext(1,ivmu) = -fac1
                end if
             end if
          else
             if (iz < nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                            + (1.-zed_upwind)*gradpar(iz+1)) &
                          *((1.+zed_upwind)*maxwell_mu(ia,iz,imu,is) &
                            + (1.-zed_upwind)*maxwell_mu(ia,iz+1,imu,is))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                if (iz > -nzgrid) then
                   fac1 = fac*((1.+zed_upwind)*gradpar(iz-1) &
                               + (1.-zed_upwind)*gradpar(iz)) &
                             *((1.+zed_upwind)*maxwell_mu(ia,iz-1,imu,is) &
                               + (1.-zed_upwind)*maxwell_mu(ia,iz,imu,is))
                else
                   fac1 = fac*((1.+zed_upwind)*gradpar(nzgrid-1) &
                               + (1.-zed_upwind)*gradpar(iz)) &
                             *((1.+zed_upwind)*maxwell_mu(ia,nzgrid-1,imu,is) &
                               + (1.-zed_upwind)*maxwell_mu(ia,iz,imu,is))
                end if
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                            + (1.-zed_upwind)*gradpar(-nzgrid+1)) &
                          *((1.+zed_upwind)*maxwell_mu(ia,iz,imu,is) &
                            + (1.-zed_upwind)*maxwell_mu(ia,-nzgrid+1,imu,is))
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar(iz-1) &
                            + (1.-zed_upwind)*gradpar(iz)) &
                          *((1.+zed_upwind)*maxwell_mu(ia,iz-1,imu,is) &
                            + (1.-zed_upwind)*maxwell_mu(ia,iz,imu,is))
             end if
             gext(idx,ivmu) = -fac0
             if (idx > 1) gext(idx-1,ivmu) = fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (periodic(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = -fac0
                   gext(nz_ext-1,ivmu) = fac1
                else if (idx == 2) then
                   gext(nz_ext,ivmu) = fac1
                end if
             end if
          end if

          ! hack for now (duplicates much of the effort from sweep_zed_zonal)
          if (periodic(iky)) then
             call sweep_zed_zonal_response (iv, is, stream_sign(iv), gext(:,ivmu))
          else
             ! invert parallel streaming equation to get g^{n+1} on extended zed grid
             ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
             call stream_tridiagonal_solve (iky, ie, iv, is, gext(:,ivmu))
          end if

       end do
    else
       ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          ! initialize g to zero everywhere along extended zed domain
          gext(:,ivmu) = 0.0
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
        
          ! give unit impulse to phi at this zed location
          ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)
        
          ! NB:  assuming equal spacing in zed below
          ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*J0/dz
          ! b.gradz left out because needs to be centred in zed
          if (driftkinetic_implicit) then
             gyro_fac = 1.0
          else
             gyro_fac = aj0x(iky,ikx,iz,ivmu)
          end if
        
          fac = -0.25*(1.+time_upwind)*code_dt*vpa(iv)*spec(is)%stm_psi0 &
               *gyro_fac*spec(is)%zt*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)

          mu_dbdzed_p = 1./delzed(0)+mu(imu)*dbdzed(ia,iz)*(1.+zed_upwind)
          mu_dbdzed_m = 1./delzed(0)+mu(imu)*dbdzed(ia,iz)*(1.-zed_upwind)
        
          ! stream_sign < 0 corresponds to positive advection speed
          if (stream_sign(iv)<0) then
             if (iz > -nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(iz-1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                if (iz < nzgrid) then
                   izp = iz+1
                else
                   izp = -nzgrid+1
                end if
                fac1 = fac*((1.+zed_upwind)*gradpar(izp) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(nzgrid-1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the right of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar(iz+1) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             end if
             gext(idx,ivmu) = fac0
             if (idx < nz_ext) gext(idx+1,ivmu) = -fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (periodic(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = fac0
                else if (idx == nz_ext-1) then
                   gext(1,ivmu) = -fac1
                end if
             end if
          else
             if (iz < nzgrid) then
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(iz+1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                if (iz > -nzgrid) then
                   izm = iz-1
                else
                   izm = nzgrid-1
                end if
                fac1 = fac*((1.+zed_upwind)*gradpar(izm) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             else
                ! fac0 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at this zed index
                fac0 = fac*((1.+zed_upwind)*gradpar(iz) &
                     + (1.-zed_upwind)*gradpar(-nzgrid+1))*mu_dbdzed_p
                ! fac1 is the factor multiplying delphi on the RHS
                ! of the homogeneous GKE at the zed index to the left of
                ! this one
                fac1 = fac*((1.+zed_upwind)*gradpar(iz-1) &
                     + (1.-zed_upwind)*gradpar(iz))*mu_dbdzed_m
             end if
             gext(idx,ivmu) = -fac0
             if (idx > 1) gext(idx-1,ivmu) = fac1
             ! zonal mode BC is periodic instead of zero, so must
             ! treat specially
             if (periodic(iky)) then
                if (idx == 1) then
                   gext(nz_ext,ivmu) = -fac0
                   gext(nz_ext-1,ivmu) = fac1
                else if (idx == 2) then
                   gext(nz_ext,ivmu) = fac1
                end if
             end if
          end if

          ! hack for now (duplicates much of the effort from sweep_zed_zonal)
          if (periodic(iky)) then
             call sweep_zed_zonal_response (iv, is, stream_sign(iv), gext(:,ivmu))
          else
             ! invert parallel streaming equation to get g^{n+1} on extended zed grid
             ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
             call stream_tridiagonal_solve (iky, ie, iv, is, gext(:,ivmu))
          end if
        
       end do
    end if

    ! we now have g on the extended zed domain at this ky and set of connected kx values
    ! corresponding to a unit impulse in phi at this location
    ! now integrate over velocities to get a square response matrix
    ! (this ends the parallelization over velocity space, so every core should have a
    !  copy of phiext)
    call integrate_over_velocity (gext, phiext, iky, ie)

#if !defined ISO_C_BINDING || !defined MPI
    response_matrix(iky)%eigen(ie)%zloc(:,idx) = phiext(:nresponse)
#else
    if(sgproc0) response_matrix(iky)%eigen(ie)%zloc(:,idx) = phiext(:nresponse)
#endif

  end subroutine get_dgdphi_matrix_column

! subroutine get_phi_matrix
! end subroutine get_phi_matrix

  subroutine sweep_zed_zonal_response (iv, is, sgn, g)

    use zgrid, only: nzgrid, delzed, nztot
    use run_parameters, only: zed_upwind, time_upwind
    use parallel_streaming, only: stream_c

    implicit none

    integer, intent (in) :: iv, is, sgn
    complex, dimension (:), intent (in out) :: g

    integer :: iz, iz1, iz2
    real :: fac1, fac2
    complex, dimension (:), allocatable :: gcf, gpi

    allocate (gpi(-nzgrid:nzgrid))
    allocate (gcf(-nzgrid:nzgrid))
    ! ky=0 is 2pi periodic (no extended zgrid)
    ! decompose into complementary function + particular integral
    ! zero BC for particular integral
    ! unit BC for complementary function (no source)
    if (sgn < 0) then
       iz1 = -nzgrid ; iz2 = nzgrid
    else
       iz1 = nzgrid ; iz2 = -nzgrid
    end if
    gpi(iz1) = 0. ; gcf(iz1) = 1.
    do iz = iz1-sgn, iz2, -sgn
       fac1 = 1.0+zed_upwind+sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       fac2 = 1.0-zed_upwind-sgn*(1.0+time_upwind)*stream_c(iz,iv,is)/delzed(0)
       gpi(iz) = (-gpi(iz+sgn)*fac2 + 2.0*g(iz+nzgrid+1))/fac1
       gcf(iz) = -gcf(iz+sgn)*fac2/fac1
    end do
    ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
    g = gpi + (spread(gpi(iz2),1,nztot)/(1.-gcf(iz2)))*gcf
    deallocate (gpi, gcf)

  end subroutine sweep_zed_zonal_response

  subroutine integrate_over_velocity(g,phi,iky,ie)

    use stella_layouts, only: vmu_lo
    use species, only: nspec, spec
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: nsegments
    use kt_grids, only: zonal_mode, akx
    use vpamu_grids, only: integrate_species
    use gyro_averages, only: gyro_average
    use mp, only: sum_allreduce

    implicit none

    complex, dimension (:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:), intent (out) :: phi
    integer, intent (in) :: iky, ie

    integer :: idx, iseg, ikx, iz, ia
    integer :: izl_offset
    real, dimension (nspec) :: wgt
    complex, dimension (:), allocatable :: g0

    ia = 1

    allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    wgt = spec%z*spec%dens_psi0
    phi = 0.

    idx = 0 ; izl_offset = 0
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    if(zonal_mode(iky).and.abs(akx(ikx)) < epsilon(0.)) then
      phi(:) = 0.0
      return
    endif
    do iz = iz_low(iseg), iz_up(iseg)
       idx = idx + 1
        call gyro_average (g(idx,:), iky, ikx, iz, g0)
        call integrate_species (g0, iz, wgt, phi(idx),reduce_in=.false.)
    end do
    izl_offset = 1
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
             idx = idx + 1
             call gyro_average (g(idx,:), iky, ikx, iz, g0)
             call integrate_species (g0, iz, wgt, phi(idx),reduce_in=.false.)
          end do
          if (izl_offset == 0) izl_offset = 1
       end do
    end if

    call sum_allreduce(phi)

  end subroutine integrate_over_velocity

  subroutine get_fields_for_response_matrix (phi, iky, ie)

    use stella_layouts, only: vmu_lo
    use species, only: spec
    use species, only: has_electron_species
    use stella_geometry, only: dl_over_b
    use extended_zgrid, only: iz_low, iz_up
    use extended_zgrid, only: ikxmod
    use extended_zgrid, only: nsegments
    use kt_grids, only: zonal_mode, akx
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use fields, only: gamtot3
    use fields_arrays, only: gamtot

    implicit none

    complex, dimension (:), intent (inout) :: phi
    integer, intent (in) :: iky, ie
  
    integer :: idx, iseg, ikx, iz, ia
    integer :: izl_offset
    complex, dimension (:), allocatable :: g0
    complex :: tmp

    ia = 1

    allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    idx = 0 ; izl_offset = 0
    iseg = 1
    ikx = ikxmod(iseg,ie,iky)
    if(zonal_mode(iky).and.abs(akx(ikx)) < epsilon(0.)) then
      phi(:) = 0.0
      return
    endif
    do iz = iz_low(iseg), iz_up(iseg)
       idx = idx + 1
       phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
    end do
    izl_offset = 1
    if (nsegments(ie,iky) > 1) then
       do iseg = 2, nsegments(ie,iky)
          ikx = ikxmod(iseg,ie,iky)
          do iz = iz_low(iseg)+izl_offset, iz_up(iseg)
             idx = idx + 1
             phi(idx) = phi(idx)/gamtot(iky,ikx,iz)
          end do
          if (izl_offset == 0) izl_offset = 1
       end do
    end if

    if (.not.has_electron_species(spec) .and. &
         adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       if (zonal_mode(iky)) then
          ! no connections for ky = 0
          iseg = 1
          tmp = sum(dl_over_b(ia,:)*phi)
          phi = phi + tmp*gamtot3(ikxmod(1,ie,iky),:)
       end if
    end if

    deallocate (g0)

  end subroutine get_fields_for_response_matrix

  subroutine finish_response_matrix

    use fields_arrays, only: response_matrix
#if !defined ISO_C_BINDING 

    implicit none
  
#else
    use mpi

    implicit none

    integer :: ierr

    if (window.ne.MPI_WIN_NULL) call mpi_win_free(window,ierr)
#endif 

    if (allocated(response_matrix)) deallocate (response_matrix)
    response_matrix_initialized = .false.

  end subroutine finish_response_matrix


#ifdef MPI
!----------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!PARALLEL LU DECOMPOSITIONS!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
!----------------------------------------------------------!

#ifdef ISO_C_BINDING

  !this subroutine parallelizes the LU decomposition on a single
  !node using MPIs shared memory interface
  !It also splits up jtwist the independent matrices across nodes
  !Ideal speed up: cores_per_node*min(jtwist,ncores)
  subroutine parallel_LU_decomposition_local (iky)


    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use fields_arrays, only: response_matrix
    use mp, only: barrier, broadcast, sum_allreduce
    use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
    use mp, only: scrossdomprocs, sgproc0, mp_abort, real_size
    use mp, only: job, iproc, proc0, nproc, numnodes, inode
    use mp_lu_decomposition, only: lu_decomposition_local
    use job_manage, only: njobs
    use extended_zgrid, only: neigen
    use mpi

    implicit none

    integer, intent (in) :: iky

    integer, dimension (:,:), allocatable :: eig_limits
    integer, dimension (:), allocatable :: job_list
    integer, dimension (:), allocatable :: row_limits
    logical, dimension (:,:), allocatable :: node_jobs

    complex, dimension (:,:), pointer :: lu

    type(c_ptr) :: bptr

    logical :: needs_send = .false.

    integer :: prior_focus, nodes_on_job
    integer :: ijob, j, ie, n, ediv, emod
    integer :: jroot, neig, ierr, win, nroot
    integer (kind=MPI_ADDRESS_KIND) :: win_size
    integer :: disp_unit = 1
    real :: dmax

    prior_focus = curr_focus

    call scope (sharedprocs)

    allocate (node_jobs(numnodes,njobs)); node_jobs = .false.
    allocate (job_list(nproc)); job_list = 0
    allocate (row_limits(0:nproc))
    allocate (eig_limits(0:numnodes,njobs)); eig_limits = 0

    job_list(iproc+1) = job
    call sum_allreduce(job_list)

    if(proc0) then
      do j = 1, nproc
        node_jobs(inode+1,job_list(j)+1) = .true. !create a map of which nodes have which jobs
      enddo
    endif

    !make sure all processors have this map
    call scope (allprocs)
    call mpi_allreduce &
         (MPI_IN_PLACE, node_jobs,size(node_jobs), MPI_LOGICAL, MPI_LOR, mp_comm, ierr)
    call scope (sharedprocs)

    do ijob = 0, njobs-1
      jroot=-1
      do j = 1, nproc
        if(job_list(j).eq.ijob) then
          jroot = j-1 !the first processor on this job will be the root process
          exit
        endif
      enddo

      if(jroot.eq.-1) cycle !no processors on this node are on this job

      if(iproc.eq.jroot) neig = neigen(iky)

      ! broadcast number of matrices
      call broadcast(neig,jroot)

      ! split up neig across nodes that have the current job
      nodes_on_job = count(node_jobs(:,ijob+1))
      ediv = neig/nodes_on_job
      emod = mod(neig,nodes_on_job)
      
      eig_limits(0,ijob+1) = 1
      do j = 1, numnodes
        if(node_jobs(j,ijob+1)) then
          eig_limits(j,ijob+1) = eig_limits(j-1,ijob+1) + ediv
          if(emod.gt.0) then
            eig_limits(j,ijob+1) = eig_limits(j,ijob+1) + 1
            emod = emod - 1
          endif
        else
          eig_limits(j,ijob+1) = eig_limits(j-1,ijob+1)
        endif
      enddo

      do ie = eig_limits(inode,ijob+1), eig_limits(inode+1,ijob+1)-1
        win_size = 0
        if(iproc.eq.jroot) then
          needs_send = .true.
          n = size(response_matrix(iky)%eigen(ie)%idx)
          win_size = int(n*n,MPI_ADDRESS_KIND)*2*real_size !complex size
        endif

        !broadcast size of matrix
        call broadcast(n, jroot)

        !allocate the window
        call mpi_win_allocate_shared(win_size,disp_unit,MPI_INFO_NULL,mp_comm,bptr,win,ierr)

        if(iproc.ne.jroot) then
          !make sure all the procs have the right memory address
          call mpi_win_shared_query(win, jroot, win_size, disp_unit, bptr, ierr)
        endif

        ! bind this c_ptr to our fortran matrix
        call c_f_pointer(bptr,lu,(/n,n/))

        !load the matrix
        if(iproc.eq.jroot) lu = response_matrix(iky)%eigen(ie)%zloc

        !syncronize the processors
        call mpi_win_fence(0,win,ierr)

        ! All the processors have the matrix.
        ! Now perform LU decomposition
        call lu_decomposition_local (mp_comm, jroot, win, lu, &
                                     response_matrix(iky)%eigen(ie)%idx, dmax)

        !copy the decomposed matrix over
        if(iproc.eq.jroot) response_matrix(iky)%eigen(ie)%zloc = lu

        call mpi_win_free(win,ierr)
      enddo
    enddo

    call scope(scrossdomprocs)

    !copy all the matrices across all nodes
    if(sgproc0) then
      do ie = 1, neigen(iky)
        nroot = 0
        if(needs_send.and. &
          (ie.ge.eig_limits(inode,job+1).and.ie.lt.eig_limits(inode+1,job+1))) nroot = iproc
        !first let processors know who is sending the data
        call sum_allreduce(nroot)
        !now send the data
        call broadcast(response_matrix(iky)%eigen(ie)%zloc,nroot)
        call broadcast(response_matrix(iky)%eigen(ie)%idx, nroot)
      enddo
    endif

    call scope(prior_focus)

    deallocate (node_jobs,job_list,row_limits,eig_limits)
  end subroutine parallel_LU_decomposition_local

#endif /* ISO_C_BINDING */

  !this subroutine parallelizes the LU decomposition across
  !all cores. Ideal speed up: ncores
  subroutine parallel_LU_decomposition_global (iky)

    use fields_arrays, only: response_matrix
    use mp, only: barrier, broadcast, sum_allreduce
    use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
    use mp, only: job, iproc, proc0, nproc, mpicmplx
#ifdef ISO_C_BINDING
    use mp, only: sgproc0, scrossdomprocs
#endif
    use job_manage, only: njobs
    use extended_zgrid, only: neigen
    use mpi
    use linear_solve, only: imaxloc
  
    implicit none

    integer, intent (in) :: iky

    integer, dimension (:), allocatable :: job_roots,  eig_roots
    integer, dimension (:), allocatable :: row_limits, eig_limits
    integer, dimension (MPI_STATUS_SIZE) :: status

    real, parameter :: zero = 1.0e-20
    integer, dimension (:), allocatable :: idx
    complex, dimension (:,:), allocatable :: lu
    real, dimension (:), allocatable :: vv
    complex, dimension (:), allocatable :: dum

    integer :: sproc
    logical :: sproc0
  
    integer :: eig_comm, ceig_comm !c for 'cross'
    integer :: ieig_core, ceig_core, eig_cores
    integer :: ncomm
    integer :: prior_focus
    integer :: ie, ie_hi, r_lo, r_hi
    integer :: ijob, i,j,k,n, n_send, rsize
    integer :: imax, neig, ierr
    integer :: istage, nstage
    integer :: rdiv, rmod
    integer :: ediv, emod

    real :: dmax, tmp

    prior_focus = curr_focus

    sproc  = iproc
    sproc0 = proc0

    call scope (allprocs)

    allocate (job_roots(0:njobs-1)); job_roots = 0

    if(sproc0) job_roots(job) = iproc

    call sum_allreduce(job_roots)
    
    do ijob = 0, njobs-1
      if(job.eq.ijob.and.sproc0) then
        neig = neigen(iky)
      endif
    
      ! broadcast number of matrices for this job
      call broadcast(neig,job_roots(ijob))

      !set up communicator for cores working on a single matrix
      call mpi_comm_split(mp_comm,mod(iproc,neig),iproc,eig_comm,ierr)
      call mpi_comm_size(eig_comm,eig_cores,ierr)
      call mpi_comm_rank(eig_comm,ieig_core,ierr)

      !set up a communicator that crosses the previous one
      call mpi_comm_split(mp_comm,ieig_core,iproc,ceig_comm,ierr)
      call mpi_comm_rank(ceig_comm,ceig_core,ierr)

      call mpi_bcast(ceig_core,1,MPI_INT,0,eig_comm,ierr)

      ncomm = min(neig,nproc) !number of communicators

      allocate (eig_roots(0:ncomm-1)); eig_roots = 0
      allocate (eig_limits(0:ncomm))
      allocate (row_limits(0:eig_cores))

      if(ieig_core.eq.0) eig_roots(ceig_core) = iproc

      call sum_allreduce(eig_roots)

      ! split up neigen across cores
      ediv = neig/ncomm
      emod = mod(neig,ncomm)

      !how many stages will the LU decomposition take?
      nstage = ediv
      if(emod.gt.0) nstage = nstage + 1


      !determine which parts of neigen this communicator processes
      eig_limits(0) = 1
      do j = 1, ncomm
        eig_limits(j) = eig_limits(j-1) + ediv
        if(j.le.emod) then
          eig_limits(j) = eig_limits(j) + 1
        endif
      enddo

      do istage = 0, nstage-1
        !transfer the data from job root to root of subcommunicator
        do j = 0,ncomm-1
          ie    = eig_limits(j) + istage
          ie_hi = eig_limits(j+1) - 1
          if(ie.gt.ie_hi) cycle

          if(iproc.eq.job_roots(ijob).and.iproc.eq.eig_roots(j)) then !no need for data transfer
            n=size(response_matrix(iky)%eigen(ie)%idx)
            allocate(lu(n,n))
            lu=response_matrix(iky)%eigen(ie)%zloc
          else if(iproc.eq.job_roots(ijob)) then !send data to subroots
            !send size of matrix
            n_send=size(response_matrix(iky)%eigen(ie)%idx)
            call mpi_send(n_send,1,MPI_INT,eig_roots(j),j,mp_comm,ierr)
            !send matrix
            call mpi_send(response_matrix(iky)%eigen(ie)%zloc, &
                          n_send*n_send,mpicmplx,eig_roots(j),nproc+j,mp_comm,ierr)
          else if(iproc.eq.eig_roots(j)) then !subroot gets the data
            !receive size of matrix
            call mpi_recv(n,1,MPI_INT,job_roots(ijob),j,mp_comm,status, ierr)
            allocate(lu(n,n))
            !receive matrix
            call mpi_recv(lu,n*n,mpicmplx,job_roots(ijob),nproc+j,mp_comm,status,ierr)
          endif
        enddo

        if(istage.ge.(eig_limits(ceig_core+1)-eig_limits(ceig_core))) cycle !nothing for this communicator to do

        !broadcast matrix and its size across the communicator
        call mpi_bcast(n,1,MPI_INT,0,eig_comm,ierr)
        if(.not.allocated(lu)) allocate (lu(n,n))
        if(.not.allocated(vv)) allocate (vv(n))

        call mpi_bcast(lu,n*n,mpicmplx,0,eig_comm,ierr)

        allocate (dum(n))
        allocate (idx(n))

        ! All the processors have the matrix.
        ! Now perform LU decomposition
        vv = maxval(cabs(lu),dim=2)
        if (any(vv==0.0)) &
          write (*,*) 'singular matrix in lu_decomposition on job ', job, ', process ', iproc
        vv = 1.0/vv
        do j = 1, n
          !divide up the work using row_limits
          rdiv = (n-j)/eig_cores
          rmod = mod(n-j,eig_cores)
          row_limits(0) = j+1
          if(rdiv.eq.0) then
            row_limits(rmod+1:) = -1
            do k=1,rmod
              row_limits(k)  = row_limits(k-1) + 1
            enddo
          else
            do k=1,eig_cores
              row_limits(k) = row_limits(k-1) + rdiv
              if(k.le.rmod) row_limits(k) = row_limits(k) + 1
            enddo
          endif

          !pivot if needed
          dmax = -1.0
          do k = j, n
            tmp = vv(k)*abs(lu(k,j))
            if(tmp.gt.dmax) then 
              dmax = tmp
              imax = k
            endif
          enddo
!         imax = (j-1) + imaxloc(vv(j:n)*cabs(lu(j:n,j)))
          if (j /= imax) then
            dum = lu(imax,:)
            lu(imax,:) = lu(j,:)
            lu(j,:) = dum
            vv(imax) = vv(j)
          end if
          if(ieig_core.eq.0) idx(j) = imax

           !get the lead multiplier
           if (lu(j,j)==0.0) lu(j,j) = zero
           do i = j+1, n
             lu(i,j) = lu(i,j)/lu(j,j)
           enddo

           r_lo = row_limits(ieig_core)
           r_hi = row_limits(ieig_core+1)-1

           do k=r_lo, r_hi
             do i = j+1,n
                lu(i,k) = lu(i,k) - lu(i,j)*lu(j,k)
             enddo
           enddo

           do i=0,eig_cores-1
             r_lo = row_limits(i)
             r_hi = row_limits(i+1)-1
             rsize = (r_hi-r_lo+1)*(n-j)
             if(r_lo.gt.r_hi) cycle
             !call mpi_bcast(lu(j+1:n,r_lo:r_hi),rsize,mpicmplx,i,eig_comm,ierr)
             do k=r_lo,r_hi
               call mpi_bcast(lu(j+1:n,k),n-j,mpicmplx,i,eig_comm,ierr)
             enddo
           enddo
        enddo
        !LU decomposition ends here

        !copy the decomposed matrix over
        do j = 0,ncomm-1

          ie = eig_limits(j) + istage
          ie_hi = eig_limits(j+1) - 1
          if(ie.gt.ie_hi) cycle

          if(iproc.eq.job_roots(ijob).and.iproc.eq.eig_roots(j)) then !no need for data transfer
            response_matrix(iky)%eigen(ie)%zloc = lu
            response_matrix(iky)%eigen(ie)%idx  = idx
          else if(iproc.eq.eig_roots(j)) then !subroot sends the data
            !send indices
            call mpi_send(idx,n,MPI_INT,job_roots(ijob),j,mp_comm,ierr)
            !send matrix
            call mpi_send(lu,n*n,mpicmplx,job_roots(ijob),nproc+j,mp_comm,ierr)
          else if(iproc.eq.job_roots(ijob)) then !receive data from subroot
            !receive indices
            call mpi_recv(response_matrix(iky)%eigen(ie)%idx, &
                          n,MPI_INT,eig_roots(j),j,mp_comm,status,ierr)
            !receive matrix
            call mpi_recv(response_matrix(iky)%eigen(ie)%zloc, &
                          n*n,mpicmplx,eig_roots(j),nproc+j,mp_comm,status,ierr)
          endif
        enddo
        deallocate (vv, lu, idx, dum)
      enddo
      deallocate (eig_roots, eig_limits, row_limits)
    enddo

#ifdef ISO_C_BINDING
    if(sgproc0) then
      call scope(scrossdomprocs)
      !copy all the matrices across all nodes
      do ie = 1, neigen(iky)
        call broadcast(response_matrix(iky)%eigen(ie)%zloc)
        call broadcast(response_matrix(iky)%eigen(ie)%idx)
      enddo
    endif

    call scope(prior_focus)
#else
    call scope(prior_focus)
    
    !copy all the matrices across all nodes
    do ie = 1, neigen(iky)
      call broadcast(response_matrix(iky)%eigen(ie)%zloc)
      call broadcast(response_matrix(iky)%eigen(ie)%idx)
    enddo
#endif
    deallocate (job_roots)
  end subroutine parallel_LU_decomposition_global

#endif /* MPI */

end module response_matrix
