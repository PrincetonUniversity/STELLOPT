program stella

  use mp, only: proc0
  use redistribute, only: scatter
  use job_manage, only: time_message, checkstop, job_fork
  use run_parameters, only: nstep, fphi, fapar
  use stella_time, only: update_time, code_time, code_dt
  use dist_redistribute, only: kxkyz2vmu
  use time_advance, only: advance_stella
  use stella_diagnostics, only: diagnose_stella, nsave
  use stella_save, only: stella_save_for_restart
  use dist_fn_arrays, only: gnew, gvmu
  use file_utils, only: error_unit, flush_output_file
  use stella_stellopt_mod

  implicit none

  ! Add the version number and date of last change when uploading to github
  character(len=4), parameter :: VERNUM = '0.3'
  character(len=10), parameter :: VERDATE = '2021.03.26'

  logical :: debug = .false.
  logical :: stop_stella = .false.
  logical :: mpi_initialized = .false.

  integer :: istep0, istep, ierr
  integer :: istatus
  real, dimension (2) :: time_init = 0.
  real, dimension (2) :: time_diagnostics = 0.
  real, dimension (2) :: time_total = 0.

  ! Initiate stella
  call init_stella(istep0, VERNUM, VERDATE)

  ! Add a header to the output file
  if (proc0) then
    write (*,'(A)') "############################################################"
    write (*,'(A)') "                OVERVIEW OF THE SIMULATION"
    write (*,'(A)') "############################################################"
    write (*,'(A)') " "
    write (*,'(A)') "    istep       time           dt         |phi|^2"
    write (*,'(A)') "------------------------------------------------------------"
  end if

  ! Diagnose stella
  if (debug) write(*,*) 'stella::diagnose_stella'
  if (istep0.eq.0) call diagnose_stella (istep0)

  ! Advance stella until istep=nstep
  if (debug) write(*,*) 'stella::advance_stella'
  do istep = (istep0+1), nstep
     if (debug) write(*,*) 'istep = ', istep
     if (mod(istep,10)==0) call checkstop (stop_stella)
     if (stop_stella) exit
     call advance_stella(istep)
     call update_time
     if (nsave > 0 .and. mod(istep,nsave)==0) then
        call scatter (kxkyz2vmu, gnew, gvmu)
        call stella_save_for_restart (gvmu, istep, code_time, code_dt, istatus)
     end if
     call time_message(.false.,time_diagnostics,' diagnostics')
     call diagnose_stella (istep)
     call time_message(.false.,time_diagnostics,' diagnostics')
     ierr = error_unit()
     call flush_output_file (ierr)
  end do

  ! Finish stella
  if (debug) write(*,*) 'stella::finish_stella'
  call finish_stella(last_call=.true.)

contains

  !> Initialise stella
  !>
  !> Calls the initialisation routines for all the geometry, physics, and
  !> diagnostic modules
  subroutine init_stella(istep0, VERNUM, VERDATE)
    use mp, only: init_mp, broadcast, sum_allreduce
    use mp, only: proc0,job, scope, subprocs, crossdomprocs
    use file_utils, only: init_file_utils
    use file_utils, only: runtype_option_switch, runtype_multibox
    use file_utils, only: run_name, init_job_name
    use file_utils, only: flush_output_file, error_unit
    use job_manage, only: checktime, time_message, njobs
    use physics_parameters, only: init_physics_parameters, g_exb, g_exbfac
    use physics_flags, only: init_physics_flags
    use physics_flags, only: nonlinear, include_parallel_nonlinearity
    use physics_flags, only: full_flux_surface, radial_variation
    use physics_flags, only: hammett_flow_shear
    use run_parameters, only: init_run_parameters
    use run_parameters, only: avail_cpu_time, nstep, rng_seed, delt
    use run_parameters, only: stream_implicit, driftkinetic_implicit
    use run_parameters, only: delt_option_switch, delt_option_auto
    use run_parameters, only: mat_gen, mat_read
    use species, only: init_species, read_species_knobs
    use species, only: nspec, communicate_species_multibox
    use zgrid, only: init_zgrid
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: init_geometry, communicate_geo_multibox
    use stella_geometry, only: finish_init_geometry
    use stella_layouts, only: init_stella_layouts, init_dist_fn_layouts
    use response_matrix, only: init_response_matrix, read_response_matrix
    use init_g, only: ginit, init_init_g, phiinit, scale_to_phiinit
    use init_g, only: tstart
    use fields, only: init_fields, advance_fields, get_radial_correction, fields_updated
    use stella_time, only: init_tstart, init_delt
    use stella_diagnostics, only: init_stella_diagnostics
    use fields_arrays, only: phi, apar
    use dist_fn_arrays, only: gnew, gvmu
    use dist_fn, only: init_gxyz, init_dist_fn
    use dist_redistribute, only: init_redistribute
    use time_advance, only: init_time_advance
    use extended_zgrid, only: init_extended_zgrid
    use kt_grids, only: init_kt_grids, read_kt_grids_parameters, communicate_ktgrids_multibox
    use kt_grids, only: naky, nakx, ny, nx, nalpha
    use vpamu_grids, only: init_vpamu_grids, read_vpamu_grids_parameters
    use vpamu_grids, only: nvgrid, nmu
    use stella_transforms, only: init_transforms
    use stella_save, only: init_dt
    use multibox, only: read_multibox_parameters, init_multibox, rhoL, rhoR
    use multibox, only: communicate_multibox_parameters, multibox_communicate
    use ran, only: get_rnd_seed_length, init_ranf
    use dissipation, only: init_dissipation
    use sources, only: init_sources
    use volume_averages, only: init_volume_averages, volume_average

    implicit none

    !> Starting timestep: zero unless the simulation has been restarted
    integer, intent (out) :: istep0
    !> stella version number
    character(len=4), intent (in) :: VERNUM
    !> Release date
    character(len=10), intent (in) :: VERDATE
    logical :: exit, list, restarted, needs_transforms
    character (500), target :: cbuff
    integer, dimension (:), allocatable  :: seed
    integer :: i, n, ierr
    real :: delt_saved, phi2, rescale

    ! initialize mpi message passing
    if (.not.mpi_initialized) call init_mp
    mpi_initialized = .true.
    debug = debug .and. proc0

    ! initialize timer
    if (debug) write (*,*) 'stella::init_stella::check_time'
    call checktime(avail_cpu_time,exit)

    if (proc0) then
       ! write message to screen with useful info regarding start of simulation
       if (debug) write (*,*) 'stella::init_stella::write_start_message'
       call write_start_message(VERNUM, VERDATE)
       ! initialize file i/o
       if (debug) write (*,*) 'stella::init_stella::init_file_utils'
       call init_file_utils (list)
       call time_message(.false.,time_total,' Total')
       call time_message(.false.,time_init,' Initialization')
    end if

    call broadcast (list)
    call broadcast (runtype_option_switch)
    if(list) call job_fork

    !proc0 may have changed
    debug = debug .and. proc0

    if (proc0) cbuff = trim(run_name)
    call broadcast (cbuff)
    if (.not. proc0) call init_job_name(cbuff)


    if (debug) write(6,*) "stella::init_stella::init_physics_flags"
    call init_physics_flags
    if (debug) write(6,*) "stella::init_stella::init_physics_parameters"
    call init_physics_parameters
    if (debug) write(6,*) "stella::init_stella::init_zgrid"
    call init_zgrid
    if (debug) write (6,*) "stella::init_stella::read_species_knobs"
    call read_species_knobs
    if (debug) write (6,*) "stella::init_stella::read_multibox_parameters"
    call read_kt_grids_parameters
    if (debug) write (6,*) "stella::init_stella::read_vpamu_grids_parameters"
    call read_multibox_parameters
    if (debug) write (6,*) "stella::init_stella::read_kt_grids_parameters"
    call read_vpamu_grids_parameters
    if (debug) write (6,*) "stella::init_stella::init_dist_fn_layouts"
    call init_dist_fn_layouts (nzgrid, ntubes, naky, nakx, nvgrid, nmu, nspec, ny, nx, nalpha)
    if (debug) write(6,*) "stella::init_stella::init_geometry"
    call init_geometry (nalpha)
    if (debug) write (6,*) 'stella::init_stella::init_species'
    call init_species
    if (debug) write(6,*) "stella::init_stella::init_init_g"
    call init_init_g
    if (debug) write(6,*) "stella::init_stella::init_run_parameters"
    call init_run_parameters

    if (debug) write(6,*) "stella::init_stella::init_ranf"
    n=get_rnd_seed_length()
    allocate(seed(n))
    if(rng_seed .lt. 0) then
      call init_ranf(.true.,seed,job+2)
    else
      seed = rng_seed + 37 * (/ ( i - 1, i = 1, n) /)
      call init_ranf(.false.,seed,job+2)
    endif
    deallocate(seed)

    if (debug) write (6,*) 'stella::init_stella::init_stella_layouts'
    call init_stella_layouts
    if (debug) write (6,*) 'stella::init_stella::init_kt_grids'
    call init_kt_grids
    !if (nonlinear .or. full_flux_surface .or. include_parallel_nonlinearity & 
    !    .or. radial_variation .or. (g_exb*g_exb).gt.epsilon(0.0).or. &
    !    runtype_option_switch.eq.runtype_multibox) then
    needs_transforms = .false.
    if(nonlinear.or.include_parallel_nonlinearity) needs_transforms = .true.
    if(radial_variation.or.full_flux_surface)      needs_transforms = .true.
    if(runtype_option_switch.eq.runtype_multibox)  needs_transforms = .true.
    if(abs(g_exb*g_exbfac).gt.epsilon(0.).and..not.hammett_flow_shear) & 
      needs_transforms = .true.
    if (needs_transforms) then
       if (debug) write (*,*) "stella::init_stella::init_transforms"
       call init_transforms
    end if
    if (debug) write (6,*) 'stella::init_stella::init_multibox'
    call init_multibox
    if (proc0.and.runtype_option_switch.eq.runtype_multibox &
             .and.(job.eq.1).and.radial_variation) then
      if (debug) write (6,*) 'stella::init_stella::init_multibox_geo'
      call communicate_geo_multibox(rhoL,rhoR)
      if (debug) write (6,*) 'stella::init_stella::init_multibox_spec'
      call communicate_species_multibox(rhoL,rhoR)
    endif
    if (runtype_option_switch.eq.runtype_multibox.and.(job.eq.1)) then
      call communicate_multibox_parameters
    endif
    if (runtype_option_switch.eq.runtype_multibox.and.radial_variation) then
      if (debug) write (6,*) 'stella::init_stella::init_multibox_ktgrid'
      call communicate_ktgrids_multibox
    endif
    if (debug) write (6,*) 'stella::init_stella::finish_init_geometry'
    call finish_init_geometry
    if (debug) write (6,*) 'stella::init_stella::init_vpamu_grids'
    call init_vpamu_grids
    if (debug) write (6,*) 'stella::init_stella::init_extended_zgrid'
    call init_extended_zgrid
    if (debug) write (6,*) 'stella::init_stella::init_volume_averages'
    call init_volume_averages
    if (debug) write(6,*) "stella::init_stella::init_dist_fn"
    call init_dist_fn
    if (debug) write(6,*) "stella::init_stella::init_redistribute"
    call init_redistribute
    if (debug) write (6,*) 'stella::init_stella::init_dissipation'
    call init_dissipation
    if (debug) write (6,*) 'stella::init_stella::init_sources'
    call init_sources
    if (debug) write (6,*) 'stella::init_stella::init_fields'
    call init_fields
    if (debug) write(6,*) "stella::init_stella::ginit"
    call ginit (restarted,istep0)
    if (debug) write(6,*) "stella::init_stella::init_gxyz"
    call init_gxyz (restarted)

    if(restarted.and.delt_option_switch == delt_option_auto) then
      delt_saved = delt
      if (debug) write(6,*) "stella::init_stella::init_dt"
      call init_dt(delt_saved, istatus)
      if(istatus == 0) delt = delt_saved
    endif
    if (debug) write(6,*) "stella::init_stella::init_delt"
    call init_delt(delt)
    if (debug) write (6,*) 'stella::init_stella::init_time_advance'
    call init_time_advance
    if (stream_implicit .or. driftkinetic_implicit) then
       if (mat_read) then
          if (debug) write (6,*) "stella::init_stella::read_response_matrix"
          call read_response_matrix
       else
          if (debug) write (6,*) "stella::init_stella::init_response_matrix"
          call init_response_matrix
       end if
    end if

    if (debug) write (6,*) 'stella::init_stella::get_fields'
    ! get initial field from initial distribution function
    call advance_fields (gnew, phi, apar, dist='gbar')
    if(radial_variation) then
      if (debug) write (6,*) 'stella::init_stella::get_radial_correction'
      call get_radial_correction(gnew,phi,dist='gbar')
    endif

    if(runtype_option_switch.eq.runtype_multibox) then
      if (debug) write (6,*) 'stella::init_stella:multibox_communicate'
      call multibox_communicate (gnew)
      if(job.eq.1) then
        fields_updated=.false.
        call advance_fields (gnew, phi, apar, dist='gbar')
      endif
    endif

    ! FLAG - the following code should probably go elsewhere
    if(.not.restarted.and.scale_to_phiinit) then
      call volume_average(phi,phi2)
      if(runtype_option_switch.eq.runtype_multibox) then
        call scope(crossdomprocs)
        call sum_allreduce(phi2)
        call scope(subprocs)
        phi2=phi2/njobs
      endif
      rescale=phiinit/sqrt(phi2)
      phi  = rescale*phi
      gnew = rescale*gnew
      gvmu = rescale*gvmu
    endif

    if (debug) write (6,*) 'stella::init_stella::init_stella_diagnostics'
    call init_stella_diagnostics (restarted,tstart)
    if (debug) write (6,*) 'stella::init_stella::init_tstart'
    call init_tstart (tstart)

    ierr = error_unit()
    if (proc0) call flush_output_file (ierr)

    if (proc0) call time_message(.false.,time_init,' Initialization')

  end subroutine init_stella

  !> Write the start message to screen
  subroutine write_start_message(VERNUM, VERDATE)
    use mp, only: proc0, nproc

    implicit none

    !> stella version number
    character(len=4), intent (in) :: VERNUM
    !> Release date
    character(len=10), intent (in) :: VERDATE
    character(len=23) :: str

    if (proc0) then
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) ''//achar(27)//'[32m'
      write (*,*) "            I8            ,dPYb, ,dPYb,            "
      write (*,*) "            I8            IP'`Yb IP'`Yb            "
      write (*,*) "         88888888         I8  8I I8  8I            "
      write (*,*) "            I8            I8  8' I8  8'            "
      write (*,*) "   ,g,      I8    ,ggg,   I8 dP  I8 dP    ,gggg,gg "
      write (*,*) "  ,8'8,     I8   i8' '8i  I8dP   I8dP    dP'  'Y8I "
      write (*,*) " ,8'  Yb   ,I8,  I8, ,8I  I8P    I8P    i8'    ,8I "
      write (*,*) ",8'_   8) ,d88b, `YbadP' ,d8b,_ ,d8b,_ ,d8,   ,d8b,"
      write (*,*) 'P` "YY8P8P8P""Y8888P"Y8888P`"Y888P`"Y88P"Y8888P"`Y8'
      write (*,*) ''//achar(27)//'[0m'
      write (*,*) ' '
      write (*,*) ' '
      write (*,*) '                       Version ', VERNUM
      write (*,*) '                        ', VERDATE
      write (*,*) ' '
      write (*,*) '                   Add author names.,'
      write (*,*) '                  More author names...,'
      write (*,*) ' '
      write (*,*) '                   Add institutions...'
      write (*,*) ' '
      write (*,*) ' '
      write (*,'(A)') "############################################################"
      write (*,'(A)') "                     PARALLEL COMPUTING"
      write (*,'(A)') "############################################################"
      if (nproc==1) then
         write (str,'(I10, A)') nproc, " processor."
         write (*,'(A,A,A)') " Running on ", adjustl(trim(str))
      else
         write (str,'(I10, A)') nproc, " processors."
         write (*,'(A,A,A)') " Running on ", adjustl(trim(str))
      end if
      write(*,*)
    end if

  end subroutine write_start_message

  !> Finish a simulation, call the finialisation routines of all modules
  subroutine finish_stella (last_call)

    use mp, only: finish_mp
    use mp, only: proc0
    use file_utils, only: finish_file_utils
    use job_manage, only: time_message
    use physics_parameters, only: finish_physics_parameters
    use physics_flags, only: finish_physics_flags
    use run_parameters, only: finish_run_parameters
    use zgrid, only: finish_zgrid
    use species, only: finish_species
    use time_advance, only: time_gke, time_parallel_nl
    use time_advance, only: finish_time_advance
    use parallel_streaming, only: time_parallel_streaming
    use mirror_terms, only: time_mirror
    use dissipation, only: time_collisions
    use sources, only: finish_sources
    use init_g, only: finish_init_g
    use dist_fn, only: finish_dist_fn
    use dist_redistribute, only: finish_redistribute
    use fields, only: finish_fields
    use fields, only: time_field_solve
    use stella_diagnostics, only: finish_stella_diagnostics
    use response_matrix, only: finish_response_matrix
    use stella_geometry, only: finish_geometry
    use extended_zgrid, only: finish_extended_zgrid
    use vpamu_grids, only: finish_vpamu_grids
    use kt_grids, only: finish_kt_grids
    use volume_averages, only: finish_volume_averages

    implicit none

    logical, intent (in), optional :: last_call

    if (debug) write (*,*) 'stella::finish_stella::finish_stella_diagnostics'
    call finish_stella_diagnostics(istep)
    if (debug) write (*,*) 'stella::finish_stella::finish_response_matrix'
    call finish_response_matrix
    if (debug) write (*,*) 'stella::finish_stella::finish_fields'
    call finish_fields
    if (debug) write (*,*) 'stella::finish_stella::finish_time_advance'
    call finish_time_advance
    if (debug) write (*,*) 'stella::finish_stella::finish_sources'
    call finish_sources
    if (debug) write (*,*) 'stella::finish_stella::finish_volume_averages'
    call finish_volume_averages
    if (debug) write (*,*) 'stella::finish_stella::finish_extended_zgrid'
    call finish_extended_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_dist_fn'
    call finish_dist_fn
    if (debug) write (*,*) 'stella::finish_stella::finish_redistribute'
    call finish_redistribute
    if (debug) write (*,*) 'stella::finish_stella::finish_init_g'
    call finish_init_g
    if (debug) write (*,*) 'stella::finish_stella::finish_vpamu_grids'
    call finish_vpamu_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_kt_grids'
    call finish_kt_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_run_parameters'
    call finish_run_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_species'
    call finish_species
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_flags'
    call finish_physics_flags
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_parameters'
    call finish_physics_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_geometry'
    call finish_geometry
    if (debug) write (*,*) 'stella::finish_stella::finish_zgrid'
    call finish_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_file_utils'
    if (proc0) then
       call finish_file_utils
       call time_message(.false.,time_total,' Total')
       write (*,*)
       write (*,'(A)') "############################################################"
       write (*,'(A)') "                        ELAPSED TIME"
       write (*,'(A)') "############################################################"
       write (*,fmt=101) 'initialization:', time_init(1)/60., 'min'
       write (*,fmt=101) 'diagnostics:', time_diagnostics(1)/60., 'min'
       write (*,fmt=101) 'fields:', time_field_solve(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_field_solve(1,2)/60., 'min'
       write (*,fmt=101) '(int_dv_g):', time_field_solve(1,3)/60., 'min'
       write (*,fmt=101) '(get_phi):', time_field_solve(1,4)/60., 'min'
       write (*,fmt=101) '(phi_adia_elec):', time_field_solve(1,5)/60., 'min'
       write (*,fmt=101) 'mirror:', time_mirror(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_mirror(1,2)/60., 'min'
       write (*,fmt=101) 'stream:', time_parallel_streaming(1)/60., 'min'
       write (*,fmt=101) 'dgdx:', time_gke(1,5)/60., 'min'
       write (*,fmt=101) 'dgdy:', time_gke(1,4)/60., 'min'
       write (*,fmt=101) 'wstar:', time_gke(1,6)/60., 'min'
       write (*,fmt=101) 'collisions:', time_collisions(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_collisions(1,2)/60., 'min'
       write (*,fmt=101) 'ExB nonlin:', time_gke(1,7)/60., 'min'
       write (*,fmt=101) 'parallel nonlin:', time_parallel_nl(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_parallel_nl(1,2)/60., 'min'
       write (*,fmt=101) 'radial var:', time_gke(1,10)/60., 'min'
       write (*,fmt=101) 'total implicit: ', time_gke(1,9)/60., 'min'
       write (*,fmt=101) 'total explicit: ', time_gke(1,8)/60., 'min'
       write (*,fmt=101) 'total:', time_total(1)/60., 'min'
       write (*,*)
    end if
101 format (a17,0pf8.2,a4)

    if (debug) write (*,*) 'stella::finish_stella::finish_mp'
    ! finish (clean up) mpi message passing
    if (present(last_call)) then
       call finish_mp
       mpi_initialized = .false.
    end if

  end subroutine finish_stella

end program stella
