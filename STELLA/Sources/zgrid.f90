module zgrid

  implicit none

  public :: init_zgrid, finish_zgrid
  public :: nzed, nzgrid, nperiod, ntubes
  public :: nztot, nz2pi
  public :: zed
  public :: delzed
  public :: zed_equal_arc
  public :: get_total_arc_length
  public :: get_arc_length_grid
  public :: shat_zero
  public :: boundary_option_switch
  public :: boundary_option_zero
  public :: boundary_option_self_periodic
  public :: boundary_option_linked

  private

  integer :: nzed, nzgrid, nztot, nz2pi
  integer :: nperiod, ntubes
  logical :: zed_equal_arc
  real :: shat_zero
  real, dimension (:), allocatable :: zed, delzed

  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_linked = 3

  logical :: zgridinit = .false.

contains

  subroutine init_zgrid

    use mp, only: proc0
    use constants, only: pi

    implicit none

    integer :: i

    if (zgridinit) return
    zgridinit = .true.

    if (proc0) then
       call read_parameters
    end if
    call broadcast_parameters

    if (.not.allocated(zed)) allocate (zed(-nzgrid:nzgrid))
    if (.not.allocated(delzed)) allocate (delzed(-nzgrid:nzgrid))

    zed = (/ (i*pi/real(nzed/2), i=-nzgrid, nzgrid ) /)
    delzed(:nzgrid-1) = zed(-nzgrid+1:) - zed(:nzgrid-1)
    delzed(nzgrid) = delzed(-nzgrid)

    nztot = 2*nzgrid+1
    ! number of zed in a 2*pi segment, including points at +/- pi
    nz2pi = 2*(nzed/2)+1

  end subroutine init_zgrid

  subroutine read_parameters

    use file_utils, only: input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value
    use physics_flags, only: full_flux_surface

    implicit none

    integer :: in_file, ierr
    logical :: exist

    type (text_option), dimension (6), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked) /)
    character(20) :: boundary_option

    namelist /zgrid_parameters/ nzed, nperiod, ntubes, &
         shat_zero, boundary_option, zed_equal_arc

    nzed = 24
    nperiod = 1
    ntubes = 1
    boundary_option = 'default'
    ! if zed_equal_arc = T, then zed is chosen to be arc length
    ! if zed_equal_arc = F, then zed is poloidal (axisymmetric)
    ! or zeta (toroidal) angle
    zed_equal_arc = .false.
    ! set minimum shat value below which we assume
    ! periodic BC
    shat_zero = 1.e-5

    in_file = input_unit_exist("zgrid_parameters", exist)
    if (exist) read (unit=in_file, nml=zgrid_parameters)

    ierr = error_unit()
    call get_option_value &
         (boundary_option, boundaryopts, boundary_option_switch, &
         ierr, "boundary_option in dist_fn_knobs")

    ! note that boundary_option may be changed to self-periodic later
    ! if magnetic shear is smaller than shat_zero
    ! cannot do this here as magnetic shear has yet to be input

    nzgrid = nzed/2 + (nperiod-1)*nzed

    ! force use of equal arc grid to ensure gradpar alpha-independent
    ! necessary to obtain efficient numerical solution of parallel streaming
    if (full_flux_surface) zed_equal_arc = .true.

  end subroutine read_parameters

  subroutine broadcast_parameters

    use mp, only: broadcast

    implicit none

    call broadcast (nzed)
    call broadcast (nzgrid)
    call broadcast (nperiod)
    call broadcast (ntubes)
    call broadcast (zed_equal_arc)
    call broadcast (shat_zero)
    call broadcast (boundary_option_switch)

  end subroutine broadcast_parameters

  subroutine finish_zgrid

    implicit none

    if (allocated(zed)) deallocate (zed)
    if (allocated(delzed)) deallocate (delzed)

    zgridinit = .false.

  end subroutine finish_zgrid

  subroutine get_total_arc_length (nz, gp, dz, length)

    implicit none

    integer, intent (in) :: nz
    real, dimension (-nz:), intent (in) :: gp
    real, intent (in) :: dz
    real, intent (out) :: length

    call integrate_zed (nz, dz, 1./gp, length)

  end subroutine get_total_arc_length

  subroutine get_arc_length_grid (nz_max, nzext_max, zboundary, gp, dz, zarc)

    implicit none

    integer, intent (in) :: nz_max, nzext_max
    real, intent (in) :: zboundary
    real, dimension (-nzext_max:), intent (in) :: gp
    real, intent (in) :: dz
    real, dimension (-nzext_max:), intent (out) :: zarc

    integer :: iz

    zarc(-nz_max) = zboundary
    if (nz_max /= nzext_max) then
       do iz = -nzext_max, -nz_max-1
          call integrate_zed (nzext_max, dz, 1./gp(iz:-nz_max), zarc(iz))
          zarc(iz) = zarc(-nz_max) - zarc(iz)
       end do
    end if
    do iz = -nz_max+1, nzext_max
       call integrate_zed (nz_max, dz, 1./gp(-nz_max:iz), zarc(iz))
       zarc(iz) = zarc(-nz_max) + zarc(iz)
    end do
    
  end subroutine get_arc_length_grid

  ! trapezoidal rule to integrate in zed
  subroutine integrate_zed (nz, dz, f, intf)
    
    implicit none
    
    integer, intent (in) :: nz
    real, intent (in) :: dz
    real, dimension (-nz:), intent (in) :: f
    real, intent (out) :: intf
    
    integer :: iz, iz_max
    
    iz_max = -nz + size(f) - 1
    intf = 0.
    do iz = -nz+1, iz_max
       intf = intf + dz*(f(iz-1)+f(iz))
    end do
    intf = 0.5*intf
    
  end subroutine integrate_zed

end module zgrid
