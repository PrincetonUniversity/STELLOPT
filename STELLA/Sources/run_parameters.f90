!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters. In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.
 
module run_parameters

  implicit none

  public :: init_run_parameters, finish_run_parameters
  public :: fphi, fapar, fbpar
  public :: code_delt_max
  public :: nstep, delt
  public :: cfl_cushion, delt_adjust
  public :: avail_cpu_time
  public :: stream_implicit, mirror_implicit
  public :: drifts_implicit
  public :: driftkinetic_implicit
  public :: fully_explicit
  public :: ky_solve_radial, ky_solve_real
  public :: maxwellian_inside_zed_derivative
  public :: stream_matrix_inversion
  public :: mirror_semi_lagrange, mirror_linear_interp
  public :: zed_upwind, vpa_upwind, time_upwind
  public :: fields_kxkyz, mat_gen, mat_read
  public :: rng_seed
  
  private

  real :: cfl_cushion, delt_adjust
  real :: fphi, fapar, fbpar
  real :: delt, code_delt_max
  real :: zed_upwind, vpa_upwind, time_upwind
  logical :: stream_implicit, mirror_implicit
  logical :: driftkinetic_implicit
  logical :: fully_explicit, drifts_implicit
  logical :: maxwellian_inside_zed_derivative
  logical :: stream_matrix_inversion
  logical :: mirror_semi_lagrange, mirror_linear_interp
  logical :: fields_kxkyz, mat_gen, mat_read
  logical :: ky_solve_real
  real :: avail_cpu_time
  integer :: nstep, ky_solve_radial
  integer :: rng_seed
  integer, public :: delt_option_switch, lu_option_switch
  integer, public, parameter :: delt_option_hand = 1, delt_option_auto = 2
  integer, public, parameter :: lu_option_none   = 1, &
                                lu_option_local  = 2, &
                                lu_option_global = 3
  logical :: initialized = .false.
  logical :: knexist

contains

  subroutine init_run_parameters

    implicit none

    if (initialized) return
    initialized = .true.

    call read_parameters

  end subroutine init_run_parameters

  subroutine read_parameters

    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use text_options, only: text_option, get_option_value
    use physics_flags, only: include_mirror, full_flux_surface

    implicit none

    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_auto), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    type (text_option), dimension (4), parameter :: lu_opts = &
         (/ text_option('default', lu_option_none), &
            text_option('none',    lu_option_none), &
            text_option('local',   lu_option_local), &
            text_option('global',  lu_option_global) /)

    character(20) :: delt_option, lu_option

    integer :: ierr, in_file

    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, &
         delt_option, lu_option, &
         avail_cpu_time, cfl_cushion, delt_adjust, &
         stream_implicit, mirror_implicit, driftkinetic_implicit, &
         drifts_implicit, &
         stream_matrix_inversion, maxwellian_inside_zed_derivative, &
         mirror_semi_lagrange, mirror_linear_interp, &
         zed_upwind, vpa_upwind, time_upwind, &
         fields_kxkyz, mat_gen, mat_read, rng_seed, &
         ky_solve_radial, ky_solve_real

    if (proc0) then
       fphi = 1.0
       fapar = 1.0
       fbpar = -1.0
       fields_kxkyz = .false.
       stream_implicit = .true.
       mirror_implicit = .true.
       drifts_implicit = .false.
       driftkinetic_implicit = .false.
       maxwellian_inside_zed_derivative = .false.
       mirror_semi_lagrange = .true.
       mirror_linear_interp = .false.
       stream_matrix_inversion = .false.
       delt_option = 'default'
       lu_option = 'default'
       zed_upwind = 0.02
       vpa_upwind = 0.02
       time_upwind = 0.02
       avail_cpu_time = 1.e10
       cfl_cushion = 0.5
       delt_adjust = 2.0
       rng_seed = -1 !negative values use current time as seed
       ky_solve_radial = 0
       ky_solve_real   = .false.
       mat_gen = .true.
       mat_read = .false.

       in_file = input_unit_exist("knobs", knexist)
       if (knexist) read (unit=in_file, nml=knobs)

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs")

       call get_option_value &
            (lu_option, lu_opts, lu_option_switch, ierr, &
            "lu_option in knobs")

    end if

    call broadcast (fields_kxkyz)
    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (lu_option_switch)
    call broadcast (cfl_cushion)
    call broadcast (delt_adjust)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (fbpar)
    call broadcast (stream_implicit)
    call broadcast (mirror_implicit)
    call broadcast (drifts_implicit)
    call broadcast (driftkinetic_implicit)
    call broadcast (maxwellian_inside_zed_derivative)
    call broadcast (mirror_semi_lagrange)
    call broadcast (mirror_linear_interp)
    call broadcast (stream_matrix_inversion)
    call broadcast (zed_upwind)
    call broadcast (vpa_upwind)
    call broadcast (time_upwind)
    call broadcast (nstep)
    call broadcast (avail_cpu_time)
    call broadcast (rng_seed)
    call broadcast (ky_solve_radial)
    call broadcast (ky_solve_real)
    call broadcast (mat_gen)
    call broadcast (mat_read)
    
    if (.not.include_mirror) mirror_implicit = .false.

    code_delt_max = delt

    if (driftkinetic_implicit) then
       stream_implicit = .false.
    else if (stream_implicit .and. full_flux_surface) then
       stream_implicit = .false.
    end if

    if (mirror_implicit.or.stream_implicit.or.driftkinetic_implicit.or.drifts_implicit) then
       fully_explicit = .false.
    else
       fully_explicit = .true.
    end if

    if (fields_kxkyz .and. full_flux_surface) then
       write (*,*)
       write (*,*) 'WARNING!!!'
       write (*,*) 'The option fields_kxkyz=T is not currently supported for full_flux_surface=T.'
       write (*,*) 'Forcing fields_kxkyz=F.'
       write (*,*) 'WARNING!!!'
       write (*,*)
       fields_kxkyz = .false.
    end if

  end subroutine read_parameters

  subroutine finish_run_parameters

    implicit none

    initialized = .false.

  end subroutine finish_run_parameters

end module run_parameters
