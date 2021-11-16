module gyro_averages

  use common_types, only: coupled_alpha_type

  public :: aj0x, aj0v, aj1x, aj1v
  public :: init_bessel, finish_bessel
  public :: gyro_average
  public :: gyro_average_j1
  
  private

  interface gyro_average
     module procedure gyro_average_kxky_local
     module procedure gyro_average_kxkyz_local
     module procedure gyro_average_vmu_local
     module procedure gyro_average_vmus_nonlocal
  end interface

  interface gyro_average_j1
     module procedure gyro_average_j1_kxky_local
     module procedure gyro_average_j1_kxkyz_local
     module procedure gyro_average_j1_vmu_local
  end interface

  real, dimension (:,:,:,:), allocatable :: aj0x, aj1x
  ! (naky, nakx, nalpha, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:), allocatable :: aj0v, aj1v
  ! (nmu, -kxkyz-layout-)

!  integer, dimension (:,:,:,:), allocatable :: ia_max_aj0a
!  complex, dimension (:,:,:,:,:), allocatable :: aj0a

  type (coupled_alpha_type), dimension (:,:,:,:), allocatable :: aj0a

  integer, dimension (:,:,:), allocatable :: ia_max_gam0a
  complex, dimension (:,:,:,:), allocatable :: gam0a, lu_gam0a
  integer, dimension (:), allocatable :: lu_gam0a_idx

  logical :: bessinit = .false.

contains

  subroutine init_bessel

    use mp, only: sum_allreduce, proc0
    use dist_fn_arrays, only: kperp2
    use physics_flags, only: full_flux_surface
    use species, only: spec, nspec
    use stella_geometry, only: bmag
    use zgrid, only: nzgrid, nztot
    use vpamu_grids, only: vperp2, nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
!    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: integrate_species
    use vpamu_grids, only: mu
    use kt_grids, only: naky, nakx, nalpha
    use kt_grids, only: naky_all, ikx_max
    use kt_grids, only: swap_kxky
    use kt_grids, only: aky, akx
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, imu_idx, iv_idx
    use spfunc, only: j0, j1
    use stella_transforms, only: transform_alpha2kalpha

    implicit none

    integer :: iz, iky, ikx, imu, is, ia, iv
    integer :: ikxkyz, ivmu
    real :: arg, dum
    integer :: ia_max_aj0a_count, ia_max_gam0a_count
    real :: ia_max_aj0a_reduction_factor, ia_max_gam0a_reduction_factor

    real, dimension (:), allocatable :: wgts
    real, dimension (:), allocatable :: aj0_alpha
    complex, dimension (:), allocatable :: aj0_kalpha, gam0_kalpha
    real, dimension (:), allocatable :: gam0_alpha
    real, dimension (:,:,:), allocatable :: kperp2_swap

    if (bessinit) return
    bessinit = .true.

    if (.not.allocated(aj0v)) then
       allocate (aj0v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       aj0v = 0.
    end if
    if (.not.allocated(aj1v)) then
       allocate (aj1v(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       aj1v = 0.
    end if
    
    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
          aj0v(imu,ikxkyz) = j0(arg)
          ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
          aj1v(imu,ikxkyz) = j1(arg)
       end do
    end do

    if (full_flux_surface) then
       ! wgts are species-dependent factors apperaing in Gamma0 factor
       allocate (wgts(nspec))
       wgts = spec%dens*spec%z**2/spec%temp

       allocate (aj0_kalpha(naky))
       allocate (gam0_kalpha(naky))
       allocate (kperp2_swap(naky_all,ikx_max,nalpha))
!       if (.not.allocated(ia_max_aj0a)) allocate(ia_max_aj0a(naky_all,ikx_max,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       if (.not.allocated(ia_max_gam0a)) allocate(ia_max_gam0a(naky_all,ikx_max,-nzgrid:nzgrid))
       if (.not.allocated(aj0a)) then
!          allocate(aj0a(naky,naky_all,ikx_max,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          allocate(aj0a(naky_all,ikx_max,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!          aj0a = 0.
       end if
       if (.not.allocated(gam0a)) then
          allocate(gam0a(naky,naky_all,ikx_max,-nzgrid:nzgrid)) ; gam0a = 0.
          allocate(lu_gam0a(naky,naky_all,ikx_max,-nzgrid:nzgrid)) ; lu_gam0a = 0.
!          allocate(lu_gam0a_idx(naky,naky_all,ikx_max,-nzgrid:nzgrid)) ; lu_gam0a_idx = 0.
       end if
       
       ia_max_aj0a_count = 0 ; ia_max_gam0a_count = 0
       do iz = -nzgrid, nzgrid
          do ia = 1, nalpha
             call swap_kxky (kperp2(:,:,ia,iz), kperp2_swap(:,:,ia))
          end do
          allocate (aj0_alpha(nalpha))

          do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
             is = is_idx(vmu_lo,ivmu)
             imu = imu_idx(vmu_lo,ivmu)
             do ikx = 1, ikx_max
                do iky = 1, naky_all
                   do ia = 1, nalpha
                      arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2_swap(iky,ikx,ia))/bmag(ia,iz)
                      aj0_alpha(ia) = j0(arg)
                   end do
                   if (iz == 0 .and. ikx == 1 .and. iky == naky_all/2 .and. imu == nmu/2 .and. is == 1 .and. iv_idx(vmu_lo,ivmu)==1) then
                      write (*,*)
                      do ia = 1, nalpha
                         write (*,*) 'j0_alpha', ia, aky(iky), mu(imu), kperp2_swap(iky,ikx,ia), aj0_alpha(ia)
                      end do
                      write (*,*)
                   end if

                   ! fourier transform aj0_alpha
                   ! note that fourier coefficients aj0_kalpha have
                   ! been filter to avoid aliasing
                   call transform_alpha2kalpha (aj0_alpha, aj0_kalpha)
!                   call find_max_required_kalpha_index (aj0_kalpha, ia_max_aj0a(iky,ikx,iz,ivmu), imu, iz)
!                   ia_max_aj0a_count = ia_max_aj0a_count + ia_max_aj0a(iky,ikx,iz,ivmu)
                   call find_max_required_kalpha_index (aj0_kalpha, aj0a(iky,ikx,iz,ivmu)%max_idx, imu, iz, is)
                   ia_max_aj0a_count = ia_max_aj0a_count + aj0a(iky,ikx,iz,ivmu)%max_idx
                   if (.not.associated(aj0a(iky,ikx,iz,ivmu)%fourier)) &
                        allocate (aj0a(iky,ikx,iz,ivmu)%fourier(aj0a(iky,ikx,iz,ivmu)%max_idx))
!                   aj0a(:ia_max_aj0a(iky,ikx,iz,ivmu),iky,ikx,iz,ivmu) = aj0_kalpha(:ia_max_aj0a(iky,ikx,iz,ivmu))
                   aj0a(iky,ikx,iz,ivmu)%fourier = aj0_kalpha(:aj0a(iky,ikx,iz,ivmu)%max_idx)
                end do
             end do
          end do
          deallocate (aj0_alpha)

          allocate (aj0_alpha(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          allocate (gam0_alpha(nalpha))
          do ikx = 1, ikx_max
             do iky = 1, naky_all
                do ia = 1, nalpha
                   ! get J0 for all vpar, mu, spec values
                   do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                      is = is_idx(vmu_lo,ivmu)
                      imu = imu_idx(vmu_lo,ivmu)
                      iv = iv_idx(vmu_lo,ivmu)
                      arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2_swap(iky,ikx,ia))/bmag(ia,iz)
                      aj0_alpha(ivmu) = j0(arg)
                      ! form coefficient needed to calculate 1-Gamma_0
                      aj0_alpha(ivmu) = (1.0-aj0_alpha(ivmu)**2) &
                           * maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
                   end do

                   ! calculate gamma0 = int d3v (1-J0^2)*F_{Maxwellian}
                   call integrate_species (aj0_alpha, iz, wgts, gam0_alpha(ia), ia)
                end do
                if (iz == 0 .and. ikx == 1 .and. iky == naky_all/2) then
                   write (*,*)
                   do ia = 1, nalpha
                      write (*,*) 'gam0_alpha', ia, aky(iky), kperp2_swap(iky,ikx,ia), gam0_alpha(ia)
                   end do
                   write (*,*)
                end if
                ! fourier transform Gamma_0(alpha)
                call transform_alpha2kalpha (gam0_alpha, gam0_kalpha)
                call find_max_required_kalpha_index (gam0_kalpha, ia_max_gam0a(iky,ikx,iz))
                ia_max_gam0a_count = ia_max_gam0a_count + ia_max_gam0a(iky,ikx,iz)
                gam0a(:ia_max_gam0a(iky,ikx,iz),iky,ikx,iz) = gam0_kalpha(:ia_max_gam0a(iky,ikx,iz))
             end do
          end do
          deallocate (aj0_alpha, gam0_alpha)
       end do

!       lu_gam0a = gam0a
!       call lu_decomposition (lu_gam0a, lu_gam0a_idx, dum)

       ! calculate the reduction factor of Fourier modes
       ! used to represent J0
       call sum_allreduce (ia_max_aj0a_count)
       ia_max_aj0a_reduction_factor = real(ia_max_aj0a_count)/real(naky*nakx*nztot*nmu*nvpa*nspec*naky)
       call sum_allreduce (ia_max_gam0a_count)
       ia_max_gam0a_reduction_factor = real(ia_max_gam0a_count)/real(naky*nakx*nztot*naky)

       if (proc0) then
          write (*,*) 'average number of k-alphas needed to represent J0(kperp(alpha))=', ia_max_aj0a_reduction_factor*naky, 'out of ', naky
          write (*,*) 'average number of k-alphas needed to represent Gamma0(kperp(alpha))=', ia_max_gam0a_reduction_factor*naky, 'out of ', naky
          write (*,*)
       end if

       deallocate (wgts)
       deallocate (aj0_kalpha, gam0_kalpha)
       deallocate (kperp2_swap)
    else
       if (.not.allocated(aj0x)) then
          allocate (aj0x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          aj0x = 0.
       end if

       if (.not.allocated(aj1x)) then
          allocate (aj1x(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
          aj1x = 0.
       end if

       ia = 1
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          do iz = -nzgrid, nzgrid
             do ikx = 1, nakx
                do iky = 1, naky
                   arg = spec(is)%bess_fac*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                   aj0x(iky,ikx,iz,ivmu) = j0(arg)
                   ! note that j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)), 
                   aj1x(iky,ikx,iz,ivmu) = j1(arg)
                end do
             end do
          end do
       end do
    end if

  end subroutine init_bessel

  subroutine find_max_required_kalpha_index (ft, idx, imu, iz, is)

    use vpamu_grids, only: maxwell_mu

    implicit none

    complex, dimension (:), intent (in) :: ft
    integer, intent (out) :: idx
    integer, intent (in), optional :: imu, iz, is

    real, parameter :: tol_floor = 0.01
    integer :: i, n
    real :: subtotal, total
    real :: tol
    real, dimension (:), allocatable :: ftmod2

    n = size(ft)

    ! use conservative estimate
    ! when deciding number of modes to retain
    if (present(imu) .and. present(iz).and.present(is)) then
       tol = min(0.1,tol_floor/maxval(maxwell_mu(:,iz,imu,is)))
    else
       tol = tol_floor
    end if

    allocate (ftmod2(n))
    ! get spectral energy associated with each mode
    ftmod2 = real(ft*conjg(ft))
    ! get total spectral energy
    total = sum(ftmod2)
    subtotal = 0.

    ! find minimum spectral index for which
    ! desired percentage of spectral energy contained
    ! in modes with indices at or below it
    if (total > 0.) then
       i = 1
       do while (subtotal < total*(1.0-tol))
          idx = i
          subtotal = sum(ftmod2(:i))
          i = i + 1
       end do
    else
       idx = 1
    end if

    deallocate (ftmod2)

  end subroutine find_max_required_kalpha_index

  subroutine finish_bessel

    implicit none

    if (allocated(aj0v)) deallocate (aj0v)
    if (allocated(aj1v)) deallocate (aj1v)
    if (allocated(aj0x)) deallocate (aj0x)
    if (allocated(aj1x)) deallocate (aj1x)
    if (allocated(aj0a)) deallocate (aj0a)
    if (allocated(gam0a)) deallocate (gam0a)
    if (allocated(lu_gam0a)) deallocate (lu_gam0a)
!    if (allocated(ia_max_aj0a)) deallocate (ia_max_aj0a)
    if (allocated(ia_max_gam0a)) deallocate (ia_max_gam0a)

    bessinit = .false.

  end subroutine finish_bessel

  subroutine gyro_average_kxky_local (field, iz, ivmu, gyro_field)

    use physics_flags, only: full_flux_surface
    use kt_grids, only: naky, nakx
    use kt_grids, only: ikx_max
    use kt_grids, only: swap_kxky_ordered, swap_kxky_back_ordered

    implicit none

    complex, dimension (:,:), intent (in) :: field
    integer, intent (in) :: iz, ivmu
    complex, dimension (:,:), intent (out) :: gyro_field

    integer :: ia, iky, ikx
    integer :: idx
    integer :: naky_all
    complex, dimension (:,:), allocatable :: field_kyall, gyro_field_kyall

    if (full_flux_surface) then
!        naky_all = 2*naky-1
!        ! need to switch from ky>=0 and all kx
!        ! to kx>=0 and all ky (using reality condition)
!        allocate (field_kyall(naky_all,ikx_max))
!        allocate (gyro_field_kyall(naky_all,ikx_max)) ; gyro_field_kyall = 0.
!        call swap_kxky_ordered (field, field_kyall)
!        ! NB: J0(kx,ky) = J0(-kx,-ky)
!        do ikx = 1, ikx_max
!           do iky = 1, naky_all
!              ! account for contributions from less positive ky values (and zero)
!              do ia = 1, min(aj0a(iky,ikx,iz,ivmu)%max_idx,iky)
!                 idx = iky-ia+1
!                 gyro_field_kyall(iky,ikx) = gyro_field_kyall(iky,ikx) &
!                      + aj0a(idx,ikx,iz,ivmu)%fourier(ia)*field_kyall(idx,ikx)
!              end do
!              ! account for contributions from more positive ky values
!              if (aj0a(iky,ikx,iz,ivmu)%max_idx > 1 .and. iky /= naky_all) then
!                 do ia = 2, min(aj0a(iky,ikx,iz,ivmu)%max_idx,naky_all-iky+1)
!                    idx = iky+ia-1
!                    gyro_field_kyall(iky,ikx) = gyro_field_kyall(iky,ikx) &
!                         + aj0a(idx,ikx,iz,ivmu)%fourier(ia)*field_kyall(idx,ikx)
!                 end do
!              end if
!           end do
!        end do
!        call swap_kxky_back_ordered (gyro_field_kyall, gyro_field)
!        deallocate (field_kyall, gyro_field_kyall)
    else
       gyro_field = aj0x(:,:,iz,ivmu)*field
       ! INCLUSION OF BELOW ALLOCATE/DEALLOCATE STATEMENTS
       ! OBSERVED TO INCREASE RUN-TIME BY FACTOR OF 3 FOR
       ! MODERATE RESOLUTION LINEAR SIMULATION
       ! I DO NOT UNDERSTAND WHY
!       allocate (field_kyall(1,1)) ; deallocate (field_kyall)
!       allocate (gyro_field_kyall(1,1)) ; deallocate (gyro_field_kyall)
    end if

  end subroutine gyro_average_kxky_local

  subroutine gyro_average_kxkyz_local (field, ivmu, gyro_field)

    use zgrid, only: nzgrid, ntubes

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: gyro_field

    integer :: iz, it

    ! NEED TO FIGURE OUT WHY BELOW LOOP SLOWS DOWN CODE A BIT
    ! WILL HAVE TO USE IT WHEN DOING FULL FLUX SURFACE
!    do it = 1, ntubes
!       do iz = -nzgrid, nzgrid
!          call gyro_average (field(:,:,iz,it), iz, ivmu, gyro_field(:,:,iz,it))
!       end do
!    end do
    gyro_field = spread(aj0x(:,:,:,ivmu),4,ntubes)*field

  end subroutine gyro_average_kxkyz_local

  subroutine gyro_average_vmu_local (distfn, ikxkyz, gyro_distfn)

    use vpamu_grids, only: nvpa

    implicit none

    complex, dimension (:,:), intent (in) :: distfn
    integer, intent (in) :: ikxkyz
    complex, dimension (:,:), intent (out) :: gyro_distfn

    gyro_distfn = spread(aj0v(:,ikxkyz),1,nvpa)*distfn

  end subroutine gyro_average_vmu_local

  subroutine gyro_average_vmus_nonlocal (field, iky, ikx, iz, gyro_field)

    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (vmu_lo%llim_proc:), intent (in) :: field
    integer, intent (in) :: iky, ikx, iz
    complex, dimension (vmu_lo%llim_proc:), intent (out) :: gyro_field

    gyro_field = aj0x(iky,ikx,iz,:)*field

  end subroutine gyro_average_vmus_nonlocal

  subroutine gyro_average_j1_kxky_local (field, iz, ivmu, gyro_field)

    use physics_flags, only: full_flux_surface
    use kt_grids, only: naky, nakx
    use kt_grids, only: ikx_max
    use kt_grids, only: swap_kxky_ordered, swap_kxky_back_ordered

    implicit none

    complex, dimension (:,:), intent (in) :: field
    integer, intent (in) :: iz, ivmu
    complex, dimension (:,:), intent (out) :: gyro_field

    gyro_field = aj1x(:,:,iz,ivmu)*field

  end subroutine gyro_average_j1_kxky_local

  subroutine gyro_average_j1_kxkyz_local (field, ivmu, gyro_field)

    use zgrid, only: nzgrid, ntubes

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    integer, intent (in) :: ivmu
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: gyro_field

    integer :: iz, it

    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          call gyro_average_j1 (field(:,:,iz,it), iz, ivmu, gyro_field(:,:,iz,it))
       end do
    end do

  end subroutine gyro_average_j1_kxkyz_local

  subroutine gyro_average_j1_vmu_local (distfn, ikxkyz, gyro_distfn)

    use vpamu_grids, only: nvpa

    implicit none

    complex, dimension (:,:), intent (in) :: distfn
    integer, intent (in) :: ikxkyz
    complex, dimension (:,:), intent (out) :: gyro_distfn

    gyro_distfn = spread(aj1v(:,ikxkyz),1,nvpa)*distfn

  end subroutine gyro_average_j1_vmu_local

end module gyro_averages
