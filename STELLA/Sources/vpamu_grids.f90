module vpamu_grids

  implicit none

  public :: init_vpamu_grids, finish_vpamu_grids
  public :: read_vpamu_grids_parameters
  public :: calculate_velocity_integrals
  public :: integrate_vmu, integrate_species
  public :: integrate_mu
  public :: vpa, nvgrid, nvpa
  public :: wgts_vpa, dvpa
  public :: mu, nmu, wgts_mu, wgts_mu_bare, dmu
  public :: dmu_ghost, dmu_cell, mu_cell
  public :: maxwell_vpa, maxwell_mu, ztmax
  public :: maxwell_fac
  public :: int_unit, int_vpa2, int_vperp2, int_vfrth
  public :: vperp2
  public :: equally_spaced_mu_grid
  public :: set_vpa_weights

  logical :: vpamu_initialized = .false.
  
  integer :: nvgrid, nvpa
  integer :: nmu
  real :: vpa_max, vperp_max

  ! arrays that are filled in vpamu_grids
  real, dimension (:), allocatable :: vpa, wgts_vpa, wgts_vpa_default, wgts_mu_bare
  real, dimension (:), allocatable :: mu, maxwell_fac
  real, dimension (:,:), allocatable :: maxwell_vpa
  real, dimension (:,:,:), allocatable :: int_unit, int_vpa2, int_vperp2, int_vfrth
  real, dimension (:,:,:), allocatable :: wgts_mu
  real, dimension (:,:,:,:), allocatable :: maxwell_mu
  real, dimension (:,:), allocatable :: ztmax
  real :: dvpa
  real, dimension (:), allocatable :: dmu
  real, dimension (:), allocatable :: dmu_ghost, dmu_cell, mu_cell
  complex, dimension (:), allocatable :: rbuffer
  logical :: equally_spaced_mu_grid, conservative_wgts_vpa

  ! vpa-mu related arrays that are declared here
  ! but allocated and filled elsewhere because they depend on z, etc.
  real, dimension (:,:,:), allocatable :: vperp2

  interface integrate_species
!     module procedure integrate_species_vmu
     module procedure integrate_species_vmu_single
     module procedure integrate_species_vmu_single_real
     module procedure integrate_species_vmu_block_complex
     module procedure integrate_species_vmu_block_real
!     module procedure integrate_species_local_complex
!     module procedure integrate_species_local_real
  end interface

  interface integrate_vmu
     module procedure integrate_vmu_local_real
     module procedure integrate_vmu_local_complex
     module procedure integrate_vmu_vmulo_complex
     module procedure integrate_vmu_vmulo_ivmu_only_real
  end interface

  interface integrate_mu
     module procedure integrate_mu_local
     module procedure integrate_mu_nonlocal
  end interface

contains

  subroutine read_vpamu_grids_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    namelist /vpamu_grids_parameters/ nvgrid, nmu, vpa_max, vperp_max, &
         equally_spaced_mu_grid, conservative_wgts_vpa

    integer :: in_file
    logical :: exist

    if (proc0) then

       nvgrid = 24
       vpa_max = 3.0
       nmu = 12
       vperp_max = 3.0
       equally_spaced_mu_grid = .false.
       conservative_wgts_vpa = .false.


       in_file = input_unit_exist("vpamu_grids_parameters", exist)
       if (exist) read (unit=in_file, nml=vpamu_grids_parameters)

    end if

    call broadcast (nvgrid)
    call broadcast (vpa_max)
    call broadcast (nmu)
    call broadcast (vperp_max)
    call broadcast (equally_spaced_mu_grid)
    call broadcast (conservative_wgts_vpa)

    nvpa = 2*nvgrid

  end subroutine read_vpamu_grids_parameters

  subroutine init_vpamu_grids

    use species, only: spec, nspec

    implicit none

    if (vpamu_initialized) return
    vpamu_initialized = .true.

    call init_vpa_grid
    call init_mu_grid

    if(.not.allocated(maxwell_fac)) then
      allocate(maxwell_fac(nspec)) ; maxwell_fac = 1.0
    endif

    maxwell_fac = spec%dens/spec%dens_psi0*(spec%temp_psi0/spec%temp)**1.5

!   do is = 1, nspec
!     norm = sum(maxwell_vpa(:,is)*wgts_vpa)
!     maxwell_vpa(:,is) = maxwell_vpa(:,is)/norm

!   
!     do iz = -nzgrid, nzgrid
!       norm = sum(maxwell_mu(1,iz,:,is)*wgts_mu(1,iz,:))
!       maxwell_mu(1,iz,:,is) = maxwell_mu(1,iz,:,is) / norm
!     enddo
!   enddo

  end subroutine init_vpamu_grids

  subroutine init_vpa_grid

    use mp, only: mp_abort
    use constants, only: pi
    use species, only: spec, nspec

    implicit none

    integer :: iv, idx, iseg, nvpa_seg
    real :: del

    if (.not. allocated(vpa)) then
       ! vpa is the parallel velocity at grid points
       allocate (vpa(nvpa)) ; vpa = 0.0
       ! wgts_vpa are the integration weights assigned
       ! to the parallel velocity grid points
       allocate (wgts_vpa(nvpa)) ; wgts_vpa = 0.0
       allocate (wgts_vpa_default(nvpa)) ; wgts_vpa_default = 0.0
       ! this is the Maxwellian in vpa
       allocate (maxwell_vpa(nvpa,nspec)) ; maxwell_vpa = 0.0
       allocate (ztmax(nvpa,nspec)) ; ztmax = 0.0
    end if

    ! velocity grid goes from -vpa_max to vpa_max
    ! with a point at vpa = 0

    ! equal grid spacing in vpa
    dvpa = 2.*vpa_max/(nvpa-1)

    ! obtain vpa grid for vpa > 0
    do iv = nvgrid+1, nvpa
       vpa(iv) = real(iv-nvgrid-0.5)*dvpa
    end do
    ! fill in vpa grid for vpa < 0
    vpa(:nvgrid) = -vpa(nvpa:nvgrid+1:-1)

    ! this is the equilibrium Maxwellian in vpa
    maxwell_vpa = exp(-spread(vpa*vpa,2,nspec)*spread(spec%temp_psi0/spec%temp,1,nvpa))
    ztmax = spread(spec%zt,1,nvpa)*maxwell_vpa

    ! get integration weights corresponding to vpa grid points
    ! for now use Simpson's rule; 
    ! i.e. subdivide grid into 3-point segments, with each segment spanning vpa_low to vpa_up
    ! then the contribution of each segment to the integral is
    ! (vpa_up - vpa_low) * (f1 + 4*f2 + f3) / 6
    ! inner boundary points are used in two segments, so they get double the weight

    if (nvpa < 6) &
         call mp_abort ('stella does not currently support nvgrid < 3.  aborting.')

    ! use simpson 3/8 rule at lower boundary and composite Simpson elsewhere
    del=0.375*dvpa
    wgts_vpa(1) = del
    wgts_vpa(2:3) = 3.*del
    wgts_vpa(4) = del
    ! composite simpson
    nvpa_seg = (nvpa-4)/2
    del = dvpa/3.
    do iseg = 1, nvpa_seg
       idx = 2*(iseg-1)+4
       wgts_vpa(idx) = wgts_vpa(idx) + del
       wgts_vpa(idx+1) = wgts_vpa(idx+1) +4.*del
       wgts_vpa(idx+2) = wgts_vpa(idx+2) + del
    end do

    ! for the sake of symmetry, do the same thing with 3/8 rule at upper boundary
    ! and composite elsewhere.
    del = 0.375*dvpa
    wgts_vpa(nvpa-3) = wgts_vpa(nvpa-3) + del
    wgts_vpa(nvpa-2:nvpa-1) = wgts_vpa(nvpa-2:nvpa-1) + 3.*del
    wgts_vpa(nvpa) = wgts_vpa(nvpa) + del
    nvpa_seg = (nvpa-4)/2
    del = dvpa/3.
    do iseg = 1, nvpa_seg
       idx = 2*(iseg-1)+1
       wgts_vpa(idx) = wgts_vpa(idx) + del
       wgts_vpa(idx+1) = wgts_vpa(idx+1) +4.*del
       wgts_vpa(idx+2) = wgts_vpa(idx+2) + del
    end do

    ! divide by 2 to account for double-counting
    wgts_vpa = 0.5*wgts_vpa/sqrt(pi)

    wgts_vpa_default = wgts_vpa

  end subroutine init_vpa_grid

  subroutine set_vpa_weights (conservative)

    use constants, only: pi

    implicit none

    logical, intent (in) :: conservative

    if (conservative) then
       wgts_vpa = dvpa/sqrt(pi)
    else if (conservative_wgts_vpa) then ! AVB: added option for density conserving form of collision operator
       wgts_vpa = dvpa/sqrt(pi)
    else if ((.not.conservative_wgts_vpa).and.(.not.conservative)) then
       wgts_vpa = wgts_vpa_default
    end if

  end subroutine set_vpa_weights

  subroutine integrate_mu_local (iz, g, total)

    use species, only: nspec

    implicit none

    integer, intent (in) :: iz
    real, dimension (:,:), intent (in) :: g
    real, dimension (:), intent (out) :: total

    integer :: is, imu, ia

    total = 0.

    ia = 1
    do is = 1, nspec
       ! sum over mu
       do imu = 1, nmu
          total(is) = total(is) + wgts_mu(ia,iz,imu)*g(imu,is)
       end do
    end do

  end subroutine integrate_mu_local

  subroutine integrate_mu_nonlocal (iz, g, total)

    use mp, only: nproc, sum_reduce
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, imu_idx, iv_idx

    implicit none

    integer, intent (in) :: iz
    real, dimension (vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,:), intent (out) :: total

    integer :: is, imu, iv, ivmu, ia

    total = 0.
    
    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       total(iv,is) = total(iv,is) + wgts_mu(ia,iz,imu)*g(ivmu)
    end do

    if (nproc > 1) call sum_reduce (total,0)

  end subroutine integrate_mu_nonlocal

  subroutine integrate_vmu_local_real (g, iz, total)

    implicit none

    real, dimension (:,:), intent (in) :: g
    integer, intent (in) :: iz
    real, intent (out) :: total

    integer :: iv, imu, ia

    total = 0.
    
    ia = 1
    do imu = 1, nmu
       do iv = 1, nvpa
          total = total + wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(iv,imu)
       end do
    end do

  end subroutine integrate_vmu_local_real

  subroutine integrate_vmu_local_complex (g, iz, total)

    implicit none

    complex, dimension (:,:), intent (in) :: g
    integer, intent (in) :: iz
    complex, intent (out) :: total

    integer :: iv, imu, ia

    total = 0.
    
    ia = 1
    do imu = 1, nmu
       do iv = 1, nvpa
          total = total + wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(iv,imu)
       end do
    end do

  end subroutine integrate_vmu_local_complex

  ! integrave over v-space in vmu_lo
  subroutine integrate_vmu_vmulo_complex (g, weights, total)

    use mp, only: sum_allreduce
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
    use zgrid, only: nzgrid

    implicit none

    integer :: ivmu, iv, iz, is, imu, ia

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:), intent (in) :: weights
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: total

    total = 0.

    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do iz = -nzgrid, nzgrid
          total(:,:,iz,:,is) = total(:,:,iz,:,is) + &
               wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(:,:,iz,:,ivmu)*weights(is)
       end do
    end do

    call sum_allreduce (total)

  end subroutine integrate_vmu_vmulo_complex

  ! integrave over v-space in vmu_lo
  subroutine integrate_vmu_vmulo_ivmu_only_real (g, ia, iz, total)

    use mp, only: sum_allreduce
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    implicit none

    integer :: ivmu, iv, is, imu

    real, dimension (vmu_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: ia, iz
    real, dimension (:), intent (out) :: total

    total = 0.

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       total(is) = total(is) + &
               wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(ivmu)
    end do

    call sum_allreduce (total)

  end subroutine integrate_vmu_vmulo_ivmu_only_real

!   subroutine integrate_species_local_real (g, weights, iz, total)

!     use species, only: nspec
!     use stella_geometry, only: bmag

!     implicit none

!     real, dimension (:,:,:), intent (in) :: g
!     real, dimension (:), intent (in) :: weights
!     integer, intent (in) :: iz
!     real, intent (out) :: total

!     integer :: iv, imu, is

!     total = 0.

!     do is = 1, nspec
!        do imu = 1, nmu
!           do iv = 1, nvpa
!              total = total + wgts_mu(imu)*wgts_vpa(iv)*bmag(1,iz)*g(iv,imu,is)*weights(is)
!           end do
!        end do
!     end do

!   end subroutine integrate_species_local_real

!   subroutine integrate_species_local_complex (g, weights, iz, total)

!     use species, only: nspec
!     use stella_geometry, only: bmag

!     implicit none

!     complex, dimension (:,:,:), intent (in) :: g
!     real, dimension (:), intent (in) :: weights
!     integer, intent (in) :: iz
!     complex, intent (out) :: total

!     integer :: iv, imu, is

!     total = 0.

!     do is = 1, nspec
!        do imu = 1, nmu
!           do iv = 1, nvpa
!              total = total + wgts_mu(imu)*wgts_vpa(iv)*bmag(1,iz)*g(iv,imu,is)*weights(is)
!           end do
!        end do
!     end do

!   end subroutine integrate_species_local_complex

!   ! integrave over v-space and sum over species
!   subroutine integrate_species_vmu (g, weights, total)

!     use mp, only: sum_allreduce
!     use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
!     use zgrid, only: nzgrid
!     use stella_geometry, only: bmag

!     implicit none

!     integer :: ivmu, iv, iz, is, imu

!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: g
!     real, dimension (:), intent (in) :: weights
!     complex, dimension (:,:,-nzgrid:), intent (out) :: total

!     total = 0.

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        imu = imu_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)
!        do iz = -nzgrid, nzgrid
!           total(:,:,iz) = total(:,:,iz) + &
!                wgts_mu(imu)*wgts_vpa(iv)*bmag(1,iz)*g(:,:,iz,ivmu)*weights(is)
!        end do
!     end do

!     call sum_allreduce (total)

!   end subroutine integrate_species_vmu

  ! integrave over v-space and sum over species for given (ky,kx,z) point
  subroutine integrate_species_vmu_single (g, iz, weights, total, ia_in, reduce_in)

    use mp, only: sum_allreduce
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    implicit none

    integer :: ivmu, iv, is, imu, ia
    logical :: reduce

    complex, dimension (vmu_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: iz
    real, dimension (:), intent (in) :: weights
    complex, intent (out) :: total
    integer, intent (in), optional :: ia_in
    logical, intent (in), optional :: reduce_in

    total = 0.
    
    if (present(ia_in)) then
       ia = ia_in
    else
       ia = 1
    end if
    if (present(reduce_in)) then
       reduce = reduce_in
    else
       reduce = .true.
    end if

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       total = total + &
            wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(ivmu)*weights(is)
    end do

    if (reduce) call sum_allreduce (total)

  end subroutine integrate_species_vmu_single

  ! integrave over v-space and sum over species for given (ky,kx,z) point
  subroutine integrate_species_vmu_single_real (g, iz, weights, total, ia_in,reduce_in)

    use mp, only: sum_allreduce
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    implicit none

    integer :: ivmu, iv, is, imu, ia
    logical :: reduce

    real, dimension (vmu_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: iz
    real, dimension (:), intent (in) :: weights
    real, intent (out) :: total
    integer, intent (in), optional :: ia_in
    logical, intent (in), optional :: reduce_in

    total = 0.
    
    if (present(ia_in)) then
       ia = ia_in
    else
       ia = 1
    end if
    if (present(reduce_in)) then
       reduce = reduce_in
    else
       reduce = .true.
    end if

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       total = total + &
            wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(ivmu)*weights(is)
    end do

    if (reduce) call sum_allreduce (total)

  end subroutine integrate_species_vmu_single_real

  subroutine integrate_species_vmu_block_complex (g, iz, weights, pout, ia_in, reduce_in)

    use mp, only: sum_allreduce
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    implicit none

    integer :: ivmu, iv, is, imu, ia
    logical :: reduce

    complex, dimension (:,:,vmu_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: iz
    integer, intent (in), optional :: ia_in
    logical, intent (in), optional :: reduce_in
    real, dimension (:), intent (in) :: weights
    complex, dimension (:,:), intent (out) :: pout

    pout =0.

    if (present(ia_in)) then
       ia = ia_in
    else
       ia = 1
    end if
    if (present(reduce_in)) then
       reduce = reduce_in
    else
       reduce = .true.
    end if

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       pout = pout + wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(:,:,ivmu)*weights(is)
    end do

    if (reduce) call sum_allreduce (pout)

  end subroutine integrate_species_vmu_block_complex

  subroutine integrate_species_vmu_block_real (g, iz, weights, pout, ia_in,reduce_in)

    use mp, only: sum_allreduce
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

    implicit none

    integer :: ivmu, iv, is, imu, ia
    logical :: reduce

    real, dimension (:,:,vmu_lo%llim_proc:), intent (in) :: g
    integer, intent (in) :: iz
    integer, intent (in), optional :: ia_in
    logical, intent (in), optional :: reduce_in
    real, dimension (:), intent (in) :: weights
    real, dimension (:,:), intent (out) :: pout

    pout =0.

    if (present(ia_in)) then
       ia = ia_in
    else
       ia = 1
    end if
    if (present(reduce_in)) then
       reduce = reduce_in
    else
       reduce = .true.
    end if

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       pout = pout + wgts_mu(ia,iz,imu)*wgts_vpa(iv)*g(:,:,ivmu)*weights(is)
    end do

    if (reduce) call sum_allreduce (pout)

  end subroutine integrate_species_vmu_block_real

  subroutine finish_vpa_grid

    implicit none

    if (allocated(vpa)) deallocate (vpa)
    if (allocated(wgts_vpa)) deallocate (wgts_vpa)
    if (allocated(wgts_vpa_default)) deallocate (wgts_vpa_default)
    if (allocated(maxwell_vpa)) deallocate (maxwell_vpa)
    if (allocated(ztmax)) deallocate (ztmax)

  end subroutine finish_vpa_grid

  subroutine init_mu_grid

    use gauss_quad, only: get_laguerre_grids
    use zgrid, only: nzgrid, nztot
    use kt_grids, only: nalpha
    use species, only: spec, nspec
    use stella_geometry, only: bmag, bmag_psi0
    
    implicit none

    integer :: imu
    real :: mu_max

    ! allocate arrays and initialize to zero
    if (.not. allocated(mu)) then
       allocate (mu(nmu)) ; mu = 0.0
       allocate (wgts_mu(nalpha,-nzgrid:nzgrid,nmu)) ; wgts_mu = 0.0
       allocate (wgts_mu_bare(nmu)) ; wgts_mu_bare = 0.0
       allocate (maxwell_mu(nalpha,-nzgrid:nzgrid,nmu,nspec)) ; maxwell_mu = 0.0
       allocate (dmu(nmu-1))
       allocate (dmu_ghost(nmu))
       allocate (mu_cell(nmu))
       allocate (dmu_cell(nmu))
    end if

    ! dvpe * vpe = d(2*mu*B0) * B/2B0
    if (equally_spaced_mu_grid) then
       ! first get equally spaced grid in mu with max value
       ! mu_max = vperp_max**2/(2*max(bmag))
       mu_max = vperp_max**2/(2.*maxval(bmag_psi0))
       ! want first grid point at dmu/2 to avoid mu=0 special point
       ! dmu/2 + (nmu-1)*dmu = mu_max
       ! so dmu = mu_max/(nmu-1/2)
       dmu = mu_max/(nmu-0.5)
       mu(1) = 0.5*dmu(1)
       do imu = 2, nmu
          mu(imu) = mu(1)+(imu-1)*dmu(1)
       end do
       ! do simplest thing to start
       wgts_mu_bare = dmu(1)
    else
       !    ! use Gauss-Laguerre quadrature in 2*mu*bmag(z=0)
       ! use Gauss-Laguerre quadrature in 2*mu*min(bmag)*max(
       call get_laguerre_grids (mu, wgts_mu_bare)
       if(vperp_max.lt.0) vperp_max = sqrt(mu(nmu))
       wgts_mu_bare = wgts_mu_bare*exp(mu)/(2.*minval(bmag_psi0)*mu(nmu)/vperp_max**2)
    
       !    mu = mu/(2.*bmag(1,0))
       mu = mu/(2.*minval(bmag_psi0)*mu(nmu)/vperp_max**2)

       dmu(:nmu-1) = mu(2:)-mu(:nmu-1)
       ! leave dmu(nmu) uninitialized. should never be used, so want 
       ! valgrind or similar to return error if it is
    end if

    ! this is the mu part of the v-space Maxwellian
    maxwell_mu = exp(-2.*spread(spread(spread(mu,1,nalpha),2,nztot)*spread(bmag,3,nmu),4,nspec) &
                       *spread(spread(spread(spec%temp_psi0/spec%temp,1,nalpha),2,nztot),3,nmu))
       
    ! factor of 2. necessary to account for 2pi from 
    ! integration over gyro-angle and 1/pi^(3/2) normalization
    ! of velocity space Jacobian
    wgts_mu = 2.*spread(spread(wgts_mu_bare,1,nalpha),2,nztot)*spread(bmag,3,nmu)


    ! add ghost cell at mu=0 and beyond mu_max for purposes of differentiation
    ! note assuming here that grid spacing for ghost cell is equal to
    ! grid spacing for last non-ghost cell
    dmu_ghost(:nmu-1) = dmu ; dmu_ghost(nmu) = dmu(nmu-1)
    ! this is mu at cell centres (including to left and right of mu grid boundary points)
    mu_cell(:nmu-1) = 0.5*(mu(:nmu-1)+mu(2:))
    mu_cell(nmu) = mu(nmu)+0.5*dmu(nmu-1)
    ! this is mu_{j+1/2} - mu_{j-1/2}
    dmu_cell(1) = mu_cell(1)
    dmu_cell(2:) = mu_cell(2:)-mu_cell(:nmu-1)

  end subroutine init_mu_grid

  subroutine finish_mu_grid

    implicit none

    if (allocated(mu)) deallocate (mu)
    if (allocated(mu_cell)) deallocate (mu_cell)
    if (allocated(wgts_mu)) deallocate (wgts_mu)
    if (allocated(wgts_mu_bare)) deallocate (wgts_mu_bare)
    if (allocated(maxwell_mu)) deallocate (maxwell_mu)
    if (allocated(dmu)) deallocate (dmu)
    if (allocated(dmu_cell)) deallocate (dmu_cell)
    if (allocated(dmu_ghost)) deallocate (dmu_ghost)
    if (allocated(rbuffer)) deallocate (rbuffer)

  end subroutine finish_mu_grid

  subroutine calculate_velocity_integrals

    use zgrid, only: nzgrid
    use species, only: nspec

    implicit none

    real, dimension (nvpa,nmu) :: moment
    integer :: ia, is, iz

    ia = 1
    is = 1
    
    if(.not.allocated(int_unit))   allocate (int_unit(1,-nzgrid:nzgrid,nspec))
    if(.not.allocated(int_vpa2))   allocate (int_vpa2(1,-nzgrid:nzgrid,nspec))
    if(.not.allocated(int_vperp2)) allocate (int_vperp2(1,-nzgrid:nzgrid,nspec))
    if(.not.allocated(int_vfrth))  allocate (int_vfrth(1,-nzgrid:nzgrid,nspec))


    do is = 1, nspec
      do iz = -nzgrid, nzgrid
        moment  = spread(maxwell_mu(ia,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)
        call integrate_vmu (moment,iz,int_unit(ia,iz,is))

        moment  = spread(maxwell_mu(ia,iz,:,is),1,nvpa)*spread(vpa(:)**2*maxwell_vpa(:,is),2,nmu)
        call integrate_vmu (moment,iz,int_vpa2(ia,iz,is))

        moment  = spread(vperp2(ia,iz,:)*maxwell_mu(ia,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)
        call integrate_vmu (moment,iz,int_vperp2(ia,iz,is))

        moment  = spread(maxwell_mu(ia,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)
        moment  = moment*(spread(vpa**2,2,nmu) + spread(vperp2(ia,iz,:),1,nvpa))**2
        call integrate_vmu (moment,iz,int_vfrth(ia,iz,is))
      enddo
    enddo

  end subroutine calculate_velocity_integrals

  subroutine finish_vpamu_grids

    implicit none
    
    call finish_vpa_grid
    call finish_mu_grid

    if (allocated(maxwell_fac)) deallocate (maxwell_fac)
    if (allocated(int_unit))    deallocate (int_unit)
    if (allocated(int_vpa2))    deallocate (int_vpa2)
    if (allocated(int_vperp2))  deallocate (int_vperp2)
    if (allocated(int_vfrth))   deallocate (int_vfrth)

    vpamu_initialized = .false.

  end subroutine finish_vpamu_grids

end module vpamu_grids
