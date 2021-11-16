! Set up the perpendicular wavenumbers by calling the appropriate sub-modules. 
module kt_grids

  implicit none

  public :: init_kt_grids, finish_kt_grids
  public :: read_kt_grids_parameters, box
  public :: aky, theta0, akx, zed0
  public :: naky, nakx, nx, ny, reality
  public :: dx,dy,dkx, dky, dx_d
  public :: jtwist, jtwistfac, ikx_twist_shift, x0, y0
  public :: x, x_d
  public :: rho, rho_d, rho_clamped, rho_d_clamped
  public :: nalpha
  public :: ikx_max, naky_all
  public :: phase_shift_fac
  public :: zonal_mode
  public :: swap_kxky, swap_kxky_back
  public :: swap_kxky_ordered, swap_kxky_back_ordered
  public :: multiply_by_rho, centered_in_rho
  public :: periodic_variation
  public :: communicate_ktgrids_multibox

  private

  interface swap_kxky
     module procedure swap_kxky_real
     module procedure swap_kxky_complex
  end interface

  real, dimension (:,:), allocatable :: theta0, zed0
  real, dimension (:), allocatable :: aky, akx
  real, dimension (:), allocatable :: x, x_d
  real, dimension (:), allocatable :: rho, rho_d, rho_clamped, rho_d_clamped
  complex, dimension (:,:), allocatable :: g0x
  real :: dx, dy, dkx, dky, dx_d
  real :: jtwistfac, phase_shift_fac
  integer :: naky, nakx, nx, ny, nalpha
  integer :: jtwist, ikx_twist_shift
  integer :: ikx_max, naky_all
  logical :: reality = .false.
  logical :: centered_in_rho, periodic_variation, randomize_phase_shift
  character(20) :: grid_option
  logical, dimension (:), allocatable :: zonal_mode

  namelist /kt_grids_knobs/ grid_option

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_range = 1, gridopt_box = 2

  real :: aky_min, aky_max
  real :: akx_min, akx_max
  real :: theta0_min, theta0_max
  real :: x0, y0
  logical :: read_kt_grids_initialized = .false.
  logical :: init_kt_grids_initialized = .false.
  logical :: box

contains
  
  subroutine read_kt_grids_parameters

    use mp, only: proc0

    implicit none

    if (read_kt_grids_initialized) return
    read_kt_grids_initialized = .true.

    if (proc0) then
       call read_grid_option
       select case (gridopt_switch)
       case (gridopt_range)
          call read_kt_grids_range
       case (gridopt_box)
          call read_kt_grids_box
       end select
    end if

    call broadcast_input

    call allocate_arrays

  end subroutine read_kt_grids_parameters

  subroutine read_grid_option

    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value

    implicit none

    type (text_option), dimension (5), parameter :: gridopts = &
         (/ text_option('default', gridopt_range), &
         text_option('range', gridopt_range), &
         text_option('box', gridopt_box), &
         text_option('annulus', gridopt_box), &
         text_option('nonlinear', gridopt_box) /)
    
    integer :: ierr, in_file
    logical :: nml_exist
    
    grid_option = 'default'
    
    in_file = input_unit_exist ("kt_grids_knobs", nml_exist)
    if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)
    
    ierr = error_unit()
    call get_option_value (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")
    
  end subroutine read_grid_option

  subroutine read_kt_grids_box

    use file_utils, only: input_unit_exist
    use physics_flags, only: full_flux_surface

    implicit none

    integer :: in_file
    logical :: exist

    namelist /kt_grids_box_parameters/ nx, ny, jtwist, jtwistfac, y0, &
                                       centered_in_rho, periodic_variation, &
                                       randomize_phase_shift, phase_shift_fac

    ! note that jtwist and y0 will possibly be modified
    ! later in init_kt_grids_box if they make it out
    ! of this subroutine with negative values
    ! it is necessary to wait until later to do this check
    ! because the values to which they may be set will
    ! depend on information from the geometry module,
    ! which itself may rely on ny from here (number of alphas)
    nx = 1
    ny = 1
    jtwist = -1
    jtwistfac = 1.
    phase_shift_fac = 0.
    y0 = -1.0
    nalpha = 1
    centered_in_rho = .true.
    randomize_phase_shift = .false.
    periodic_variation = .false.

    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_box_parameters)

    ! get the number of de-aliased modes in y and x
    naky = (ny-1)/3 + 1
    nakx = 2*((nx-1)/3) +  1

    reality = .true.

    if (full_flux_surface) nalpha = ny

  end subroutine read_kt_grids_box

  subroutine read_kt_grids_range

    use file_utils, only: input_unit, input_unit_exist

    implicit none

    integer :: in_file
    logical :: exist

    namelist /kt_grids_range_parameters/ naky, nakx,  &
         aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max

    nalpha = 1
    naky = 1
    nakx = 1
    aky_min = 0.0
    aky_max = 0.0
    ! set these to be nonsense values
    ! so can check later if they've been set
    akx_min = 0.0
    akx_max = -1.0
    theta0_min = 0.0
    theta0_max = -1.0    

    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_range_parameters)

  end subroutine read_kt_grids_range

  subroutine init_kt_grids 

    use common_types, only: flux_surface_type
    use zgrid, only: init_zgrid

    implicit none

    if (init_kt_grids_initialized) return
    init_kt_grids_initialized = .true.

    call init_zgrid

    select case (gridopt_switch)
    case (gridopt_range)
       call init_kt_grids_range 
    case (gridopt_box)
       call init_kt_grids_box 
    end select

    ! determine if iky corresponds to zonal mode
    if (.not.allocated(zonal_mode)) allocate (zonal_mode(naky))
    zonal_mode = .false.
    if (abs(aky(1)) < epsilon(0.)) zonal_mode(1) = .true.


  end subroutine init_kt_grids

  subroutine init_kt_grids_box

    use mp, only: mp_abort, proc0, broadcast
    use common_types, only: flux_surface_type
    use constants, only: pi, zi
    use stella_geometry, only: geo_surf, twist_and_shift_geo_fac, dydalpha
    use stella_geometry, only: q_as_x, get_x_to_rho, dxdXcoord, drhodpsi
    use physics_parameters, only: rhostar
    use physics_flags, only: full_flux_surface, radial_variation
    use file_utils, only: runtype_option_switch, runtype_multibox
    use zgrid, only: shat_zero, nperiod
    use ran, only: ranf

    implicit none
    
    integer :: ikx, iky
    real :: x_shift, dqdrho, pfac, norm

    box = .true.

    ! set jtwist and y0 for cases where they have not been specified
    ! and for which it makes sense to set them automatically
    if (jtwist < 1) then 
      jtwist = max(1,int(abs(twist_and_shift_geo_fac)+0.5))
      jtwist = max(1,int(jtwistfac*jtwist+0.5))
    endif
    ! signed version of jtwist, with sign determined by, e.g., magnetic shear
    ikx_twist_shift = -jtwist*int(sign(1.0,twist_and_shift_geo_fac))

    if (y0 < 0.) then
       if (full_flux_surface) then
          ! if simulating a flux annulus, then
          ! y0 determined by the physical
          ! extent of the device
          if (rhostar > 0.) then
             y0 = 1./(rhostar*geo_surf%rhotor)
          else
             call mp_abort ('must set rhostar if simulating a full flux surface. aborting.')
          end if
       else
          ! if simulating a flux tube
          ! makes no sense to have y0 < 0.0
          ! so abort
          call mp_abort ('y0 negative only makes sense when simulating a flux annulus.  aborting.')
       end if
    end if

    ! get the grid spacing in ky and then in kx using twist-and-shift BC
    dky = 1./y0
    ! non-quantized b/c assumed to be periodic instead 
    ! of linked boundary conditions if zero magnetic shear
    if (abs(geo_surf%shat) <= shat_zero) then
       dkx = dky / real(jtwist)
    else
       dkx = (2*nperiod - 1) * dky * abs(twist_and_shift_geo_fac) / real(jtwist)
    end if

    x0 = 1./dkx

    ! ky goes from zero to ky_max
    do iky = 1, naky
       aky(iky) = real(iky-1)*dky
    end do

    ! get the ikx index corresponding to kx_max
    ikx_max = nakx/2+1

    ! get the total number of ky values, including negative ky
    naky_all = 2*naky-1

    ! kx goes from zero to kx_max down to zero...
    do ikx = 1, ikx_max
       akx(ikx) = real(ikx-1)*dkx
    end do
    ! and then from -kx_max to -|kx_min|
    do ikx = ikx_max+1, nakx
       akx(ikx) = real(ikx-nakx-1)*dkx
    end do

    ! set theta0=0 for ky=0
    theta0(1,:) = 0.0
    if (q_as_x) then
       do ikx = 1, nakx
          ! theta0 = kx/ky
          theta0(2:,ikx) = akx(ikx)/aky(2:)
       end do
    else if (abs(geo_surf%shat) > shat_zero) then
       do ikx = 1, nakx
          ! theta0 = kx/ky/shat
          theta0(2:,ikx) = akx(ikx)/(aky(2:)*geo_surf%shat)
       end do
    else
       do ikx = 1, nakx
          ! if shat=0, theta0 is meaningless, so be careful
          theta0(2:,ikx) = - akx(ikx)/aky(2:)
       end do
    end if

    norm = 1.
    if (nakx.gt.1) norm = aky(2)
    if (rhostar.gt.0.) then
      phase_shift_fac =-2.*pi*(2*nperiod-1)*geo_surf%qinp_psi0*dydalpha/rhostar
    else if (randomize_phase_shift) then
      if (proc0) phase_shift_fac = 2.*pi*ranf()/norm
      call broadcast (phase_shift_fac)
    else
      phase_shift_fac = phase_shift_fac/norm
    endif


    ! for radial variation
    if(.not.allocated(x)) allocate (x(nx))
    if(.not.allocated(x_d)) allocate (x_d(nakx))
    if(.not.allocated(rho)) allocate (rho(nx))
    if(.not.allocated(rho_d)) allocate (rho_d(nakx))
    if(.not.allocated(rho_clamped)) allocate (rho_clamped(nx))
    if(.not.allocated(rho_d_clamped)) allocate (rho_d_clamped(nakx))

    dx = (2*pi*x0)/nx
    dy = (2*pi*y0)/ny

    x_shift = pi*x0
    pfac = 1.0
    if (periodic_variation) pfac = 0.5
    if(centered_in_rho) then
      if(q_as_x) then
        dqdrho = geo_surf%shat*geo_surf%qinp/geo_surf%rhoc
        x_shift = pi*x0*(1.0 &
                - 0.5*pfac*rhostar*pi*x0*geo_surf%d2qdr2/(dqdrho**2*dxdXcoord))
      else
        x_shift = pi*x0*(1.0 &
                - 0.5*pfac*rhostar*pi*x0*geo_surf%d2psidr2*drhodpsi**2/dxdXcoord)
      endif
    endif

    do ikx = 1, nx
      if (radial_variation.or.runtype_option_switch.eq.runtype_multibox) then
        if(periodic_variation) then
          if(ikx.le.nx/2) then
            x(ikx) = (ikx-1)*dx - 0.5*x_shift
          else
            x(ikx) = x(nx-ikx+1)
          endif
        else
          x(ikx) = (ikx-0.5)*dx - x_shift
        endif
      else
        x(ikx) = (ikx-1)*dx
      endif
    enddo

    dx_d = (2*pi*x0)/nakx
    do ikx = 1, nakx
      if (radial_variation.or.runtype_option_switch.eq.runtype_multibox) then
        if(periodic_variation) then
          if(ikx.le.(nakx+1)/2) then
            x_d(ikx) = (ikx-1)*dx_d - 0.5*x_shift
          else
            x_d(ikx) = x_d(nakx-ikx+1)
          endif
        else
          x_d(ikx) = (ikx-0.5)*dx_d - x_shift
        endif
      else
        x_d(ikx) = (ikx-1)*dx_d
      endif
    enddo

    call get_x_to_rho(1,x,rho)
    call get_x_to_rho(1,x_d,rho_d)

    if(.not.allocated(rho_clamped)) allocate(rho_clamped(nx)); rho_clamped = rho
    if(.not.allocated(rho_d_clamped)) allocate(rho_d_clamped(nakx)); rho_d_clamped = rho_d

    zed0 = theta0*geo_surf%zed0_fac

    if(radial_variation) call dump_radial_grid

    if(radial_variation.and.(any((rho+geo_surf%rhoc).lt.0.0) & 
                             .or.any((rho+geo_surf%rhoc).gt.1.0))) then
      call mp_abort ('rho(x) is beyond range [0,1]. Try changing rhostar or q/psi profiles')
    endif

  end subroutine init_kt_grids_box

  subroutine init_kt_grids_range

    use mp, only: mp_abort
    use common_types, only: flux_surface_type
    use stella_geometry, only: geo_surf, q_as_x
    use zgrid, only: shat_zero

    implicit none

    integer :: i, j
    real :: dkx, dky, dtheta0, tfac
    real :: zero

    box = .false.

    ! NB: we are assuming here that all ky are positive
    ! when running in range mode
    dky = 0.0
    if (naky > 1) dky = (aky_max - aky_min)/real(naky - 1)
    aky = (/ (aky_min + dky*real(i), i = 0,naky-1) /)

    ! set default akx and theta0 to 0
    akx = 0.0 ; theta0 = 0.0

    if (q_as_x) then
      tfac = 1.0
    else
      tfac = geo_surf%shat
    endif

    zero = 100.*epsilon(0.)

    ! if theta0_min and theta0_max have been specified,
    ! use them to determine akx_min and akx_max
    if (theta0_max > theta0_min-zero) then
       if (geo_surf%shat > epsilon(0.)) then
          akx_min = theta0_min * tfac * aky(1)
          akx_max = theta0_max * tfac * aky(1)
       else
          akx_min = theta0_max * tfac * aky(1)
          akx_max = theta0_min * tfac * aky(1)
       end if
    end if

    ! shat_zero is minimum shat value below which periodic BC is enforced
    if (abs(geo_surf%shat) > shat_zero) then  ! ie assumes boundary_option .eq. 'linked'
       ! if akx_min and akx_max specified in input
       ! instead of theta0_min and theta0_max,
       ! use them to get theta0_min and theta0_max
       if (theta0_min > theta0_max+zero .and. abs(aky(1)) > zero) then
          theta0_min = akx_min/(tfac*aky(1))
          theta0_max = akx_max/(tfac*aky(1))
          dtheta0 = 0.0
          if (nakx > 1) dtheta0 = (theta0_max - theta0_min)/real(nakx - 1)
          
          do j = 1, naky
             theta0(j,:) &
                  = (/ (theta0_min + dtheta0*real(i), i=0,nakx-1) /)
          end do
          akx = theta0(1,:) * tfac * aky(1)
       else if (akx_max > akx_min-zero .or. nakx.eq.1) then
          dkx = 0.0
          if (nakx > 1) dkx = (akx_max - akx_min)/real(nakx - 1)
          akx = (/ (akx_min + dkx*real(i), i = 0,nakx-1) /)

          dtheta0 = 0.0
          if (nakx > 1) dtheta0 = (theta0_max - theta0_min)/real(nakx - 1)

          if (geo_surf%shat > epsilon(0.)) then
             do j = 1, naky
                theta0(j,:) &
                     = (/ (theta0_min + dtheta0*real(i), i=0,nakx-1) /)
             end do
          else
             do j = 1, naky
                theta0(j,:) &
                     = (/ (theta0_min + dtheta0*real(i), i=nakx-1,0,-1) /)
             end do
          end if
       else
          call mp_abort ('ky=0 is inconsistent with kx_min different from kx_max. aborting.')
       end if
       
    else
       ! here assume boundary_option .eq. 'periodic'
       ! used for periodic finite kx ballooning space runs with shat=0
       dkx = 0.0
       if (nakx > 1) dkx = (akx_max - akx_min)/real(nakx - 1)
       akx = (/ (akx_min + dkx*real(i), i = 0,nakx-1) /)
    endif

    ikx_max = nakx
    naky_all = naky

  end subroutine init_kt_grids_range

  subroutine broadcast_input

    use mp, only: broadcast

    implicit none

    call broadcast (gridopt_switch)
    call broadcast (centered_in_rho)
    call broadcast (periodic_variation)
    call broadcast (naky)
    call broadcast (nakx)
    call broadcast (ny)
    call broadcast (nx)
    call broadcast (nalpha)
    call broadcast (reality)
    call broadcast (jtwist)
    call broadcast (jtwistfac)
    call broadcast (y0)
    call broadcast (aky_min)
    call broadcast (aky_max)
    call broadcast (akx_min)
    call broadcast (akx_max)
    call broadcast (theta0_min)
    call broadcast (theta0_max)
    call broadcast (randomize_phase_shift)
    call broadcast (phase_shift_fac)

  end subroutine broadcast_input

  subroutine dump_radial_grid

    use file_utils, only: run_name
    use physics_parameters, only: rhostar
    use stella_geometry, only: q_as_x, geo_surf, dxdXcoord, drhodpsi

    implicit none

    integer :: ix
    character (300) :: filename

    filename = trim(trim(run_name)//'.radial_grid')
    open (1047,file=filename,status='unknown')
    if(q_as_x) then
      write (1047,'(1a12,1e12.4,1a12,1e12.4,1a12,1e12.4,1a12,1e12.4)') & 
            '#dxdXcoord = ', dxdXcoord, &
            ' q    = '  , geo_surf%qinp, &
            ' dqdr = '  , geo_surf%shat*geo_surf%qinp/geo_surf%rhoc, &
            ' d2qdr2 = ', geo_surf%d2qdr2
      write (1047,'(3a12)') '#1.x', '2.q', '3.rho'
      do ix = 1, nx
        write (1047,'(3e12.4,i9)') &
          x(ix), &
          rhostar*x(ix)/dxdXcoord + geo_surf%qinp, &
          rho(ix) + geo_surf%rhoc
      enddo
    else
      write (1047,'(1a12,1e12.4,1a12,1e12.4,1a12,1e12.4,1a12,1e12.4)') & 
            '#dxdXcoord = ', dxdXcoord, &
            ' dpsidr    = ' , 1.0/drhodpsi, &
            ' d2psidr2 = '  , geo_surf%d2psidr2
      write (1047,'(3a12,a9)') '#1.x', '2.psi-psi0', '3.rho'
      do ix = 1, nx
        write (1047,'(3e12.4,i9)') &
            x(ix), &
            rhostar*x(ix)/dxdXcoord, &
            rho(ix) + geo_surf%rhoc
      enddo
    endif

    close (1047)

  end subroutine dump_radial_grid

  subroutine allocate_arrays

    implicit none

    allocate (akx(nakx))
    allocate (aky(naky))
    allocate (theta0(naky,nakx))
    allocate (zed0(naky,nakx))

  end subroutine allocate_arrays

  ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
  subroutine swap_kxky_complex (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:naky,:) = gin(:,:ikx_max)
    ! next fill in ky < 0, kx >= 0 elements of array using reality
    ikx = 1
    ikxneg = ikx
    do iky = naky+1, naky_all
       ! this is the ky index corresponding to +ky in original array
       ikyneg = naky_all-iky+2
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = naky+1, naky_all
          ! this is the ky index corresponding to +ky in original array
          ikyneg = naky_all-iky+2
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_complex

  ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
  subroutine swap_kxky_real (gin, gout)

    implicit none

    real, dimension (:,:), intent (in) :: gin
    real, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:naky,:) = gin(:,:ikx_max)
    ! next fill in ky < 0, kx >= 0 elements of array using reality
    ikx = 1
    ikxneg = ikx
    do iky = naky+1, naky_all
       ! this is the ky index corresponding to +ky in original array
       ikyneg = naky_all-iky+2
       gout(iky,ikx) = gin(ikyneg,ikxneg)
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = naky+1, naky_all
          ! this is the ky index corresponding to +ky in original array
          ikyneg = naky_all-iky+2
          gout(iky,ikx) = gin(ikyneg,ikxneg)
       end do
    end do

  end subroutine swap_kxky_real

  ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
  subroutine swap_kxky_ordered (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(naky:,:) = gin(:,:ikx_max)
    ! next fill in ky < 0, kx >= 0 elements of array using reality
    ikx = 1
    ikxneg = ikx
    do iky = 1, naky-1
       ! this is the ky index corresponding to +ky in original array
       ikyneg = naky-iky+1
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = 1, naky-1
          ! this is the ky index corresponding to +ky in original array
          ikyneg = naky-iky+1
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_ordered

  ! take an array with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
  ! and uses reality condition to return array
  ! with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  subroutine swap_kxky_back (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:,:ikx_max) = gin(:naky,:)
    ! next fill in kx < 0, ky >= 0 elements of array using reality
    do ikx = ikx_max+1, nakx
       ikxneg = nakx-ikx+2
       iky = 1
       ikyneg = iky
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       do iky = 2, naky
          ikyneg = naky_all-iky+2
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_back

  ! take an array with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
  ! and uses reality condition to return array
  ! with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  subroutine swap_kxky_back_ordered (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:,:ikx_max) = gin(naky:,:)
    ! next fill in kx < 0, ky >= 0 elements of array using reality
    do ikx = ikx_max+1, nakx
       ikxneg = nakx-ikx+2
       do iky = 1, naky
          ikyneg = naky-iky+1
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_back_ordered
  
  subroutine communicate_ktgrids_multibox
      use job_manage, only: njobs
      use mp, only: job, scope, &
                  crossdomprocs, subprocs,  &
                  send, receive

    implicit none

    call scope(crossdomprocs)

    if(job==1) then
      call send(phase_shift_fac, 0, 120)
      call send(phase_shift_fac, njobs-1, 130)
    elseif(job == 0) then
      call receive(phase_shift_fac, 1, 120)
    elseif(job == njobs-1) then
      call receive(phase_shift_fac, 1, 130)
    endif

    call scope(subprocs)

  end subroutine communicate_ktgrids_multibox

  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (aky)
    if (allocated(akx)) deallocate (akx)
    if (allocated(theta0)) deallocate (theta0)
    if (allocated(zed0)) deallocate (zed0)

    if (allocated(x)) deallocate (x)
    if (allocated(x_d)) deallocate (x_d)
    if (allocated(rho)) deallocate (rho)
    if (allocated(rho_d)) deallocate (rho_d)
    if (allocated(rho_clamped)) deallocate (rho_clamped)
    if (allocated(rho_d_clamped)) deallocate (rho_d_clamped)

    if (allocated(g0x)) deallocate (g0x)

    reality = .false.
    read_kt_grids_initialized = .false.
    init_kt_grids_initialized = .false.

  end subroutine finish_kt_grids

  subroutine multiply_by_rho (gin)

    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
!   use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst

    implicit none

    complex, dimension (:,:), intent (inout) :: gin

    if(.not.allocated(g0x)) allocate(g0x(naky,nakx))

    call transform_kx2x_unpadded(gin,g0x)
    g0x = spread(rho_d_clamped,1,naky)*g0x
    if(zonal_mode(1)) g0x(1,:) = real(g0x(1,:))
    call transform_x2kx_unpadded(g0x,gin)

!   if(.not.allocated(g0x)) allocate(g0x(naky,nx))

!   call transform_kx2x_xfirst(gin,g0x)
!   g0x = spread(rho_clamped,1,naky)*g0x
!   if(zonal_mode(1)) g0x(1,:) = real(g0x(1,:))
!   call transform_x2kx_xfirst(g0x,gin)

  end subroutine multiply_by_rho

end module kt_grids

