module neoclassical_terms

  implicit none

  public :: init_neoclassical_terms
  public :: finish_neoclassical_terms
  public :: include_neoclassical_terms
  public :: dfneo_dzed, dfneo_dvpa, dfneo_drho, dfneo_dalpha
  public :: dphineo_dzed, dphineo_drho, dphineo_dalpha

  private

  logical :: include_neoclassical_terms
  integer :: nradii
  real :: drho

  integer :: neo_option_switch
  integer, parameter :: neo_option_sfincs = 1

  real, dimension (:,:,:), allocatable :: dfneo_dzed, dfneo_dvpa, dfneo_drho, dfneo_dalpha
  real, dimension (:,:), allocatable :: dphineo_dzed, dphineo_drho, dphineo_dalpha

  logical :: neoinit = .false.
  logical :: debug = .false.

contains

  subroutine init_neoclassical_terms

    use zgrid, only: nzgrid
    use kt_grids, only: nalpha
    use vpamu_grids, only: nvpa, nmu
    use species, only: nspec
    use stella_layouts, only: vmu_lo
    use sfincs_interface, only: get_neo_from_sfincs
    
    implicit none

    real, dimension (:,:,:,:,:,:), allocatable :: f_neoclassical
    real, dimension (:,:,:), allocatable :: phi_neoclassical
    real, dimension (:,:,:,:,:), allocatable :: dfneo_dalpha_local
    
    integer :: iz, ialpha

    if (neoinit) return
    neoinit = .true.

    call read_parameters
    if (include_neoclassical_terms) then
       allocate (f_neoclassical(nalpha,-nzgrid:nzgrid,nvpa,nmu,nspec,-nradii/2:nradii/2))
       allocate (phi_neoclassical(nalpha,-nzgrid:nzgrid,-nradii/2:nradii/2))
       allocate (dfneo_dalpha_local(nalpha,-nzgrid:nzgrid,nvpa,nmu,nspec))
       if (.not.allocated(dfneo_dvpa)) &
            allocate (dfneo_dvpa(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dfneo_drho)) &
            allocate (dfneo_drho(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dfneo_dzed)) &
            allocate (dfneo_dzed(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dfneo_dalpha)) &
            allocate (dfneo_dalpha(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(dphineo_dzed)) &
            allocate (dphineo_dzed(nalpha,-nzgrid:nzgrid))
       if (.not.allocated(dphineo_drho)) &
            allocate (dphineo_drho(nalpha,-nzgrid:nzgrid))
       if (.not.allocated(dphineo_dalpha)) &
            allocate (dphineo_dalpha(nalpha,-nzgrid:nzgrid))
       select case (neo_option_switch)
       case (neo_option_sfincs)
          call get_neo_from_sfincs (nradii, drho, f_neoclassical, phi_neoclassical, dfneo_dalpha_local, dphineo_dalpha)
       end select
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dzed'
       call get_dfneo_dzed (f_neoclassical(:,:,:,:,:,0), dfneo_dzed)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dvpa'
       call get_dfneo_dvpa (f_neoclassical(:,:,:,:,:,0), dfneo_dvpa)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_drho'
       call get_dfneo_drho (f_neoclassical, dfneo_drho)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dphineo_dzed'
       call get_dphineo_dzed (phi_neoclassical(:,:,0), dphineo_dzed)
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::get_dphineo_drho'
       call get_dphineo_drho (phi_neoclassical, dphineo_drho)
       do iz = -nzgrid, nzgrid
          do ialpha = 1, nalpha
             call distribute_vmus_over_procs (dfneo_dalpha_local(ialpha,iz,:,:,:), dfneo_dalpha(ialpha,iz,:))
          end do
       end do
       if (debug) write (6,*) 'neoclassical_terms::init_neoclassical_terms::write_neoclassical'
       call write_neoclassical (f_neoclassical, phi_neoclassical)
       deallocate (f_neoclassical, phi_neoclassical, dfneo_dalpha_local)
    end if

  end subroutine init_neoclassical_terms

  subroutine read_parameters

    use mp, only: proc0, broadcast
    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value

    implicit none

    type (text_option), dimension (2), parameter :: neoopts = (/ &
         text_option('default', neo_option_sfincs), &
         text_option('sfincs', neo_option_sfincs) /)
    character (10) :: neo_option

    namelist /neoclassical_input/ include_neoclassical_terms, &
         neo_option, nradii, drho

    logical :: exist
    integer :: ierr, in_file

    if (proc0) then
       ! set to .true. to include neoclassical terms in GK equation
       include_neoclassical_terms = .false.
       ! number of radial points used for radial derivatives
       ! of neoclassical quantities
       nradii = 5
       ! spacing in rhoc between points used for radial derivatives
       drho = 0.01
       ! option for obtaining neoclassical distribution function and potential
       neo_option = 'sfincs'

       in_file = input_unit_exist("neoclassical_input", exist)
       if (exist) read (unit=in_file, nml=neoclassical_input)

       ierr = error_unit()
       call get_option_value &
            (neo_option, neoopts, neo_option_switch, &
            ierr, "neo_option in neoclassical_input")

       if (nradii /= 3 .and. nradii /= 5) then
          write (*,*) 'WARNING: only nradii of 3 or 5 is currently supported in neoclassical_input namelist'
          write (*,*) 'WARNING: forcing nradii=5'
          nradii = 5
       end if
    end if

    call broadcast (include_neoclassical_terms)
    call broadcast (neo_option_switch)
    call broadcast (nradii)
    call broadcast (drho)

  end subroutine read_parameters

  subroutine distribute_vmus_over_procs (local, distributed)

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx

    implicit none

    real, dimension (:,:,:), intent (in) :: local
    real, dimension (vmu_lo%llim_proc:), intent (out) :: distributed

    integer :: ivmu, iv, imu, is
    
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       distributed(ivmu) = local(iv,imu,is)
    end do

  end subroutine distribute_vmus_over_procs

  subroutine get_dfneo_dvpa (fneo, dfneo)

    use finite_differences, only: fd5pt
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: dvpa
    use species, only: nspec
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:,:,:,:), intent (in) :: fneo
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dfneo

    integer :: ia, iz, imu, is
    real, dimension (:), allocatable :: tmp1, tmp2
    real, dimension (:,:,:,:,:), allocatable :: dfneo_local

    allocate (tmp1(nvpa), tmp2(nvpa))
    allocate (dfneo_local(nalpha,-nzgrid:nzgrid,nvpa,nmu,nspec))

    do is = 1, nspec
       do imu = 1, nmu
          do iz = -nzgrid, nzgrid
             do ia = 1, nalpha
                ! hack to avoid dealing with negative indices in fd5pt
                tmp1 = fneo(ia,iz,:,imu,is)
                call fd5pt (tmp1, tmp2, dvpa)
                dfneo_local(ia,iz,:,imu,is) = tmp2
             end do
          end do
       end do
    end do

    do iz = -nzgrid, nzgrid
       do ia = 1, nalpha
          call distribute_vmus_over_procs (dfneo_local(ia,iz,:,:,:), dfneo(ia,iz,:))
       end do
    end do

    deallocate (dfneo_local)
    deallocate (tmp1, tmp2)
    
  end subroutine get_dfneo_dvpa

  subroutine get_dfneo_dzed (fneo, dfneo)

    use finite_differences, only: fd5pt
    use zgrid, only: nztot, nzgrid, delzed
    use vpamu_grids, only: nvpa, nmu
    use species, only: nspec
    use stella_layouts, only: vmu_lo
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:,:,:,:), intent (in) :: fneo
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dfneo

    integer :: iv, imu, is, iz, ia
    real, dimension (:), allocatable :: tmp1, tmp2
    real, dimension (:), allocatable :: dfneo_local(:,:,:,:,:)

    allocate (tmp1(nztot), tmp2(nztot))
    allocate (dfneo_local(nalpha,-nzgrid:nzgrid,nvpa,nmu,nspec))

    do is = 1, nspec
       do imu = 1, nmu
          do iv = 1, nvpa
             do ia = 1, nalpha
                ! hack to avoid dealing with negative indices in fd5pt
                tmp1 = fneo(ia,:,iv,imu,is)
                call fd5pt (tmp1, tmp2,delzed(0))
                dfneo_local(ia,:,iv,imu,is) = tmp2
             end do
          end do
       end do
    end do

    do iz = -nzgrid, nzgrid
       do ia = 1, nalpha
          call distribute_vmus_over_procs (dfneo_local(ia,iz,:,:,:), dfneo(ia,iz,:))
       end do
    end do

    deallocate (dfneo_local)
    deallocate (tmp1, tmp2)

  end subroutine get_dfneo_dzed

  subroutine get_dfneo_drho (fneo, dfneo)

    use finite_differences, only: fd3pt, fd5pt
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nmu
    use species, only: nspec
    use stella_layouts, only: vmu_lo
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:,:,:,:,-nradii/2:), intent (in) :: fneo
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (out) :: dfneo

    integer :: ia, iz, iv, imu, is
    real, dimension (:), allocatable :: tmp1, tmp2
    real, dimension (:,:,:,:,:), allocatable :: dfneo_local

    allocate (tmp1(nradii), tmp2(nradii))
    allocate (dfneo_local(nalpha,-nzgrid:nzgrid,nvpa,nmu,nspec))

    do is = 1, nspec
       do imu = 1, nmu
          do iv = 1, nvpa
             do iz = -nzgrid, nzgrid
                do ia = 1, nalpha
                   ! hack to avoid dealing with negative indices in fd5pt
                   tmp1 = fneo(ia,iz,iv,imu,is,:)
                   if (nradii == 5) then
                      call fd5pt (tmp1, tmp2, drho)
                   else
                      call fd3pt (tmp1, tmp2, drho)
                   end if
                   dfneo_local(ia,iz,iv,imu,is) = tmp2(nradii/2+1)
                end do
             end do
          end do
       end do
    end do

    do iz = -nzgrid, nzgrid
       do ia = 1, nalpha
          call distribute_vmus_over_procs (dfneo_local(ia,iz,:,:,:), dfneo(ia,iz,:))
       end do
    end do

    deallocate (dfneo_local)
    deallocate (tmp1, tmp2)

  end subroutine get_dfneo_drho

  subroutine get_dphineo_dzed (phineo, dphineo)

    use finite_differences, only: fd5pt
    use zgrid, only: nztot, nzgrid, delzed
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:), intent (in) :: phineo
    real, dimension (:,-nzgrid:), intent (out) :: dphineo

    integer :: ia
    real, dimension (:), allocatable :: tmp1, tmp2

    allocate (tmp1(nztot), tmp2(nztot))

    do ia = 1, nalpha
       ! hack to avoid dealing with negative indices in fd5pt
       tmp1 = phineo(ia,:)
       call fd5pt (tmp1, tmp2, delzed(0))
       dphineo(ia,:) = tmp2
    end do

    deallocate (tmp1, tmp2)

  end subroutine get_dphineo_dzed

  subroutine get_dphineo_drho (phineo, dphineo)

    use finite_differences, only: fd3pt, fd5pt
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:,-nradii/2:), intent (in) :: phineo
    real, dimension (:,-nzgrid:), intent (out) :: dphineo

    integer :: iz, ia
    real, dimension (:), allocatable :: tmp1, tmp2

    allocate (tmp1(nradii), tmp2(nradii))

    do iz = -nzgrid, nzgrid
       do ia = 1, nalpha
          ! hack to avoid dealing with negative indices in fd5pt
          tmp1 = phineo(ia,iz,:)
          if (nradii == 5) then
             call fd5pt (tmp1, tmp2, drho)
          else
             call fd3pt (tmp1, tmp2, drho)
          end if
          dphineo(ia,iz) = tmp2(nradii/2+1)
       end do
    end do

    deallocate (tmp1, tmp2)

  end subroutine get_dphineo_drho

  subroutine write_neoclassical (fnc, phinc)

    use mp, only: proc0
    use mp, only: send, receive
    use file_utils, only: open_output_file, close_output_file
    use zgrid, only: nzgrid, zed
    use vpamu_grids, only: vpa, mu
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_layouts, only: idx_local, proc_id
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:,:,:,:,-nradii/2:), intent (in) :: fnc
    real, dimension (:,-nzgrid:,-nradii/2:), intent (in) :: phinc
    
    integer :: neo_unit
    integer :: irad, iz, ivmu, iv, imu, is, ia
    real, dimension (:,:), allocatable :: dfdv_local, dfdr_local, dfdz_local

    allocate (dfdv_local(nalpha,-nzgrid:nzgrid))
    allocate (dfdr_local(nalpha,-nzgrid:nzgrid))
    allocate (dfdz_local(nalpha,-nzgrid:nzgrid))

    if (proc0) then
       call open_output_file (neo_unit,'.neoclassical')
       write (neo_unit,'(3a8,10a13)') '#1.rad', '2.spec', '3.alpha', '4.zed', '5.mu', &
            '6.vpa', '7.f_neo', '8.dfdvpa_neo', '9.dfdrho_neo', '10.dfdzed_neo', '11.phi_neo', &
            '12.dphidrho', '13.dphidzed'
    end if
    do irad = -nradii/2, nradii/2
       do ivmu = vmu_lo%llim_world, vmu_lo%ulim_world
          iv = iv_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          is = is_idx(vmu_lo,ivmu)
          if (idx_local(vmu_lo, iv, imu, is)) then
             if (proc0) then
                dfdv_local = dfneo_dvpa(:,:,ivmu)
                dfdr_local = dfneo_drho(:,:,ivmu)
                dfdz_local = dfneo_dzed(:,:,ivmu)
             else
                call send (dfneo_dvpa(:,:,ivmu), 0)
                call send (dfneo_drho(:,:,ivmu), 0)
                call send (dfneo_dzed(:,:,ivmu), 0)
             end if
          else if (proc0) then
             call receive (dfdv_local, proc_id(vmu_lo,ivmu))
             call receive (dfdr_local, proc_id(vmu_lo,ivmu))
             call receive (dfdz_local, proc_id(vmu_lo,ivmu))
          end if
          if (proc0) then
             do iz = -nzgrid, nzgrid
                do ia = 1, nalpha
                   write (neo_unit,'(3i8,10e13.5)') irad, is, ia, zed(iz), mu(imu), vpa(iv), &
                        fnc(ia,iz,iv,imu,is,irad), &
                        dfdv_local(ia,iz), &
                        dfdr_local(ia,iz), &
                        dfdz_local(ia,iz), &
                        phinc(ia,iz,irad), dphineo_drho(ia,iz), dphineo_dzed(ia,iz)
                end do
             end do
          end if
       end do
       if (proc0) write (neo_unit,*)
    end do
    if (proc0) call close_output_file (neo_unit)

    deallocate (dfdv_local, dfdr_local, dfdz_local)

  end subroutine write_neoclassical

  subroutine finish_neoclassical_terms

    implicit none

    if (allocated(dfneo_dvpa)) deallocate (dfneo_dvpa)
    if (allocated(dfneo_drho)) deallocate (dfneo_drho)
    if (allocated(dfneo_dzed)) deallocate (dfneo_dzed)
    if (allocated(dfneo_dalpha)) deallocate (dfneo_dalpha)
    if (allocated(dphineo_dzed)) deallocate (dphineo_dzed)
    if (allocated(dphineo_drho)) deallocate (dphineo_drho)
    if (allocated(dphineo_dalpha)) deallocate (dphineo_dalpha)

    neoinit = .false.

  end subroutine finish_neoclassical_terms

end module neoclassical_terms
