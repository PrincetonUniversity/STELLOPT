module flow_shear

  implicit none

  public :: flow_shear_initialized
  public :: init_flow_shear, finish_flow_shear
  public :: prl_shear, prl_shear_p, prp_shear
  public :: advance_parallel_flow_shear, advance_perp_flow_shear
  public :: v_edge, v_shift
  public :: shift_times

  private

  logical :: flow_shear_initialized = .false.

  complex, dimension (:,:), allocatable :: upwind_advect
  real, dimension (:,:,:), allocatable :: prl_shear, prl_shear_p
  real, dimension (:), allocatable :: prp_shear, shift_times

  integer :: shift_sign, shift_start

  real :: v_edge, v_shift = 0.

contains

  subroutine init_flow_shear
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use constants, only: zi, pi
    use zgrid, only: nzgrid
    use kt_grids, only: x, x_d, nalpha, nx, nakx, naky, akx, aky, ikx_max, zonal_mode, box
    use fields_arrays, only: shift_state
    use stella_geometry, only: q_as_x, geo_surf, bmag, btor, rmajor, dBdrho, dIdrho
    use stella_geometry, only: dydalpha, drhodpsi
    use physics_parameters, only: g_exb, g_exbfac, omprimfac
    use vpamu_grids, only: vperp2, vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use physics_flags, only: radial_variation, prp_shear_enabled, hammett_flow_shear
    use file_utils, only: runtype_option_switch, runtype_multibox
    use job_manage, only: njobs
    use mp, only: job, send, receive, crossdomprocs, subprocs, scope

    implicit none

    integer :: is, imu, iv, ivmu, iz, ia
    real, dimension (:,:), allocatable :: energy

    if (flow_shear_initialized) return
    flow_shear_initialized = .true.

    if(abs(g_exb*g_exbfac) > epsilon(0.)) prp_shear_enabled = .true.
    if(runtype_option_switch .eq. runtype_multibox .and. job.eq.1) then 
      hammett_flow_shear = .false.
    endif


    if(runtype_option_switch .eq. runtype_multibox) then 
      call scope(crossdomprocs)
      if(job == 1) then
        call send(g_exbfac*g_exb*x(1) ,0,120)
        call send(g_exbfac*g_exb*x(nx),njobs-1,121)
        v_shift = 0.0
      elseif (job == 0) then
        call receive(v_edge, 1, 120)
        v_shift=v_edge-g_exbfac*g_exb*x(1)
      elseif (job == njobs-1) then
        call receive(v_edge, 1, 121)
        v_shift=v_edge-g_exbfac*g_exb*x(nx)
      endif
      call scope(subprocs)
    endif

    ia=1

    !parallel flow shear

    allocate (energy(nalpha,-nzgrid:nzgrid))

    if (.not.allocated(prl_shear)) then
      allocate (prl_shear(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      prl_shear = 0.0
    endif

    if (radial_variation.and..not.allocated(prl_shear_p)) &
      allocate (prl_shear_p(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))


    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      is = is_idx(vmu_lo,ivmu)
      iv = iv_idx(vmu_lo,ivmu)
      imu = imu_idx(vmu_lo,ivmu)
      do iz = -nzgrid,nzgrid
        prl_shear(ia,iz,ivmu) = -omprimfac*g_exb*code_dt*vpa(iv)*spec(is)%stm_psi0 &
               * dydalpha*drhodpsi  &
               *(geo_surf%qinp_psi0/geo_surf%rhoc_psi0)  & 
               *(btor(iz)*rmajor(iz)/bmag(ia,iz))*(spec(is)%mass/spec(is)%temp) &
               * maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
      enddo
      if(radial_variation) then
        energy = (vpa(iv)**2 + vperp2(:,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
        prl_shear_p(:,:,ivmu) = prl_shear(:,:,ivmu)*(dIdrho/spread(rmajor*btor,1,nalpha) &
                                - spread(dBdrho,1,nalpha)/bmag &
                                - spec(is)%fprim - spec(is)%tprim*(energy-2.5)  &
                                - 2.*mu(imu)*spread(dBdrho,1,nalpha))
      endif
    enddo

    if(q_as_x) prl_shear = prl_shear/geo_surf%shat_psi0

    deallocate(energy)

    !perpendicular flow shear

    if (.not.allocated(shift_times)) allocate (shift_times(naky)) 
    if (.not.allocated(upwind_advect)) allocate (upwind_advect(naky,nakx))

    if (.not.allocated(shift_state)) then
      allocate (shift_state(naky)) 
      shift_state = 0.
    endif

    if(nakx.gt.1.and.abs(g_exb*g_exbfac).gt.0) then
      shift_times = abs(akx(2)/(aky*g_exb*g_exbfac))
    endif
    if(zonal_mode(1)) shift_times(1) = huge(0.)

    if(g_exb*g_exbfac > 0.) then
      shift_sign = -1
      shift_start = ikx_max
    else
      shift_sign=1
      shift_start = ikx_max + 1
    endif

    if(box) upwind_advect = exp(-zi*g_exbfac*g_exb*code_dt*spread(aky,2,nakx)*spread(x_d,1,naky))

  end subroutine init_flow_shear

  subroutine advance_parallel_flow_shear (gout)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use fields, only: get_dchidy
    use fields_arrays, only: phi, apar

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    

    complex, dimension (:,:), allocatable :: g0k

    integer :: ivmu, iz, it, ia

    ia = 1

    allocate (g0k(naky,nakx))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
            
          !parallel flow shear
          gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + prl_shear(ia,iz,ivmu)*g0k

        enddo
      enddo
    enddo

    deallocate (g0k)

  end subroutine advance_parallel_flow_shear

  subroutine advance_perp_flow_shear (g)

    use stella_layouts, only: vmu_lo
    use constants, only: zi
    use physics_flags, only: prp_shear_enabled, hammett_flow_shear
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use zgrid, only: nzgrid, ntubes
    use fields_arrays, only: shift_state
    use kt_grids, only: aky, nakx, naky, ikx_max, zonal_mode
    use file_utils, only: runtype_option_switch, runtype_multibox
    use stella_time, only: code_dt

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:), allocatable :: g0k, g0x

    real :: shift_fac

    integer :: ivmu, iz, it, iky

    if(.not.prp_shear_enabled) return

    allocate (g0k(naky,nakx))
    allocate (g0x(naky,nakx))

    if(hammett_flow_shear) then
      !TODO (DSO) - This assumes the timestep is small enough so that a shift is never
      !             more than a single cell
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do iky=1,naky
              if(zonal_mode(iky)) cycle
              if(shift_state(iky) > 0.5*shift_times(iky)) then
                if(shift_sign < 0) then
                  !shift everything left by one
                  g(iky,(ikx_max+1):(nakx-1),iz,it,ivmu) = g(iky,ikx_max+2:,iz,it,ivmu)
                  g(iky,nakx,iz,it,ivmu) = g(iky,1,iz,it,ivmu)
                  g(iky,:ikx_max-1,iz,it,ivmu) = g(iky,2:ikx_max,iz,it,ivmu)
                else
                  !shift everything right by one
                  g(iky,2:ikx_max,iz,it,ivmu) = g(iky,1:(ikx_max-1),iz,it,ivmu)
                  g(iky,1,iz,it,ivmu) = g(iky,nakx,iz,it,ivmu)
                  g(iky,ikx_max+2:,iz,it,ivmu) = g(iky,(ikx_max+1):(nakx-1),iz,it,ivmu)
                endif
                g(iky,shift_start,iz,it,ivmu) = 0.
              endif
            enddo
          enddo
        enddo
      enddo
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid

            g0k = g(:,:,iz,it,ivmu)

            call transform_kx2x_unpadded (g0k, g0x)
            g0x = upwind_advect*g0x
         
            call transform_x2kx_unpadded (g0x, g0k)

            do iky = 1, naky
              if(zonal_mode(iky)) cycle
              if(shift_state(iky) > shift_times(iky)) then
                g0k(iky,shift_start)   = 0.0
              endif
            enddo

            g(:,:,iz,it,ivmu) = g0k

          enddo
        enddo
      enddo
    endif

    shift_fac = 1.0
    if(hammett_flow_shear) shift_fac = 0.5

    do iky=1,naky
      if(zonal_mode(iky)) cycle
      if(shift_state(iky) > shift_fac*shift_times(iky)) then
        shift_state(iky) = shift_state(iky) - shift_times(iky)
      endif
    enddo

    shift_state = shift_state + code_dt
    if(zonal_mode(1)) shift_state(1) = 0.

    if(runtype_option_switch .eq. runtype_multibox) then 
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g(:,:,iz,it,ivmu) = g(:,:,iz,it,ivmu)*exp(-code_dt*zi*spread(aky,2,nakx)*v_shift)
          enddo
        enddo
      enddo
    endif

    deallocate(g0k,g0x)

  end subroutine advance_perp_flow_shear

  subroutine finish_flow_shear
    use fields_arrays, only: shift_state

    implicit none

    if (allocated(prl_shear)) deallocate (prl_shear)
    if (allocated(prl_shear_p)) deallocate (prl_shear_p)
    if (allocated(shift_times)) deallocate (shift_times)
    if (allocated(shift_state)) deallocate (shift_state)
    if (allocated(upwind_advect)) deallocate (upwind_advect)

    flow_shear_initialized = .false.

  end subroutine finish_flow_shear

end module flow_shear
