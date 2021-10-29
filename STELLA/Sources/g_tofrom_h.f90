module g_tofrom_h

!  public :: gbar_to_g
!  public :: gbar_to_h
!  public :: gstar_to_g
  public :: g_to_h

  private

!   interface gbar_to_g
!      module procedure gbar_to_g_kxkyz
!      module procedure gbar_to_g_vmu
!   end interface

!  interface gbar_to_h
!     module procedure gbar_to_h_kxkyz
!     module procedure gbar_to_h_vmu
!  end interface

  interface g_to_h
     module procedure g_to_h_kxkyz
     module procedure g_to_h_vmu
!     module procedure g_to_h_vmu_zext
  end interface

contains

!   subroutine gbar_to_h_vmu (g, phi, apar, facphi, facapar)

!     use species, only: spec
!     use zgrid, only: nzgrid
!     use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, imu_idx, is_idx
!     use kt_grids, only: naky, nakx
!     use gyro_averages, only: aj0x

!     implicit none
!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
!     real, intent (in) :: facphi, facapar

!     integer :: ivmu, iz, iky, ikx, is, imu, iv
!     complex :: adj

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        imu = imu_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)
!        do iz = -nzgrid, nzgrid
!           do ikx = 1, nakx
!              do iky = 1, naky
!                 adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
!                      * ( facphi*phi(iky,ikx,iz) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
!                 g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
!              end do
!           end do
!        end do
!     end do

!   end subroutine gbar_to_h_vmu

!   subroutine gbar_to_h_kxkyz (g, phi, apar, facphi, facapar)

!     use species, only: spec
!     use zgrid, only: nzgrid
!     use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
!     use vpamu_grids, only: nvpa, nmu
!     use stella_layouts, only: kxkyz_lo
!     use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
!     use gyro_averages, only: aj0v

!     implicit none
!     complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
!     real, intent (in) :: facphi, facapar

!     integer :: ikxkyz, iz, iky, ikx, is, imu, iv
!     complex :: adj

!     do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
!        iz = iz_idx(kxkyz_lo,ikxkyz)
!        ikx = ikx_idx(kxkyz_lo,ikxkyz)
!        iky = iky_idx(kxkyz_lo,ikxkyz)
!        is = is_idx(kxkyz_lo,ikxkyz)
!        do imu = 1, nmu
!           do iv = 1, nvpa
!              adj = aj0v(imu,ikxkyz)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
!                   * ( facphi*phi(iky,ikx,iz) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
!              g(iv,imu,ikxkyz) = g(iv,imu,ikxkyz) + adj
!           end do
!        end do
!     end do

!   end subroutine gbar_to_h_kxkyz

!   subroutine gbar_to_g_vmu (g, apar, facapar)

!     use species, only: spec
!     use zgrid, only: nzgrid
!     use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, imu_idx, is_idx
!     use kt_grids, only: naky, nakx
!     use gyro_averages, only: aj0x

!     implicit none
!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:), intent (in) :: apar
!     real, intent (in) :: facapar

!     integer :: ivmu, iz, iky, ikx, is, imu, iv
!     complex :: adj

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        imu = imu_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)
!        do iz = -nzgrid, nzgrid
!           do ikx = 1, nakx
!              do iky = 1, naky
!                 adj = -aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
!                      * ( facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
!                 g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
!              end do
!           end do
!        end do
!     end do

!   end subroutine gbar_to_g_vmu

!   subroutine gbar_to_g_kxkyz (g, apar, facapar)

!     use species, only: spec
!     use zgrid, only: nzgrid
!     use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa
!     use vpamu_grids, only: nvpa, nmu
!     use stella_layouts, only: kxkyz_lo
!     use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
!     use gyro_averages, only: aj0v

!     implicit none
!     complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:), intent (in) :: apar
!     real, intent (in) :: facapar

!     integer :: ikxkyz, iz, iky, ikx, is, imu, iv
!     complex :: adj

!     do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
!        iz = iz_idx(kxkyz_lo,ikxkyz)
!        ikx = ikx_idx(kxkyz_lo,ikxkyz)
!        iky = iky_idx(kxkyz_lo,ikxkyz)
!        is = is_idx(kxkyz_lo,ikxkyz)
!        do imu = 1, nmu
!           do iv = 1, nvpa
!              adj = -aj0v(imu,ikxkyz)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
!                   * ( facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
!              g(iv,imu,ikxkyz) = g(iv,imu,ikxkyz) + adj
!           end do
!        end do
!     end do

!   end subroutine gbar_to_g_kxkyz

  subroutine g_to_h_vmu (g, phi, facphi, phi_corr)

    use species, only: spec
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_geometry, only: bmag, dBdrho
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use kt_grids, only: naky, nakx, multiply_by_rho
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac, vperp2, mu, vpa
    use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
    use gyro_averages, only: gyro_average, aj0x, aj1x
    use physics_flags, only: radial_variation

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:), optional, intent (in) :: phi_corr
    real, intent (in) :: facphi

    integer :: ivmu, iz, it, is, imu, iv, ia
    complex, dimension (:,:), allocatable :: field, adjust, g0k

    allocate (field(naky,nakx))
    allocate (adjust(naky,nakx))
    if(radial_variation) then
      allocate (g0k(naky,nakx))
    endif

    ia = 1
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             field = spec(is)%zt*facphi*phi(:,:,iz,it) &
                  *maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
             if(radial_variation.and.present(phi_corr)) then
               g0k = field*( -spec(is)%tprim*(vpa(iv)**2+vperp2(ia,iz,imu) - 2.5) &
                             -spec(is)%fprim - 2.0*dBdrho(iz)*mu(imu)  &
                             -0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 &
                                 * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                                 * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)))

               call multiply_by_rho(g0k)

               field = field + g0k &
                       + phi_corr(:,:,iz,it)*spec(is)%zt*facphi &
                         *maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
             endif
             call gyro_average (field, iz, ivmu, adjust)
             g(:,:,iz,it,ivmu) = g(:,:,iz,it,ivmu) + adjust
          end do
       end do
    end do

    deallocate (field, adjust)
    if (allocated(g0k)) deallocate (g0k)

  end subroutine g_to_h_vmu

!   subroutine g_to_h_vmu_zext (gext, phiext, facphi, iky, ie)

!     use species, only: spec
!     use extended_zgrid, only: ikxmod
!     use extended_zgrid, only: iz_low, iz_up
!     use extended_zgrid, only: nsegments
!     use vpamu_grids, only: maxwell_vpa, maxwell_mu
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, imu_idx, is_idx
!     use gyro_averages, only: aj0x

!     implicit none
!     complex, dimension (:,vmu_lo%llim_proc:), intent (in out) :: gext
!     complex, dimension (:), intent (in) :: phiext
!     real, intent (in) :: facphi
!     integer, intent (in) :: iky, ie

!     integer :: ivmu, iseg, iz, ikx, is, imu, iv
!     integer :: idx
!     complex :: adj

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        imu = imu_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)

!        idx = 0
!        iseg = 1
!        ikx = ikxmod(iseg,ie,iky)
!        do iz = iz_low(iseg), iz_up(iseg)
!           idx = idx + 1
!           adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
!                * facphi*phiext(idx)
!           gext(idx,ivmu) = gext(idx,ivmu) + adj
!        end do
!        if (nsegments(ie,iky) > 1) then
!           do iseg = 2, nsegments(ie,iky)
!              do iz = iz_low(iseg)+1, iz_up(iseg)
!                 adj = aj0x(iky,ikx,iz,ivmu)*spec(is)%zt*maxwell_vpa(iv)*maxwell_mu(1,iz,imu) &
!                      * facphi*phiext(idx)
!                 gext(idx,ivmu) = gext(idx,ivmu) + adj
!                 idx = idx + 1
!              end do
!           end do
!        end if

!     end do

!   end subroutine g_to_h_vmu_zext

  subroutine g_to_h_kxkyz (g, phi, facphi)

    use species, only: spec
    use zgrid, only: nzgrid
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: gyro_average

    implicit none
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    real, intent (in) :: facphi

    integer :: ikxkyz, iz, it, iky, ikx, is, ia
    complex, dimension (:,:), allocatable :: field, adjust

    allocate (field(nvpa,nmu))
    allocate (adjust(nvpa,nmu))

    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       field = facphi*phi(iky,ikx,iz,it)*spec(is)%zt &
            *spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)
       call gyro_average (field, ikxkyz, adjust)
       g(:,:,ikxkyz) = g(:,:,ikxkyz) + adjust
    end do

    deallocate (field, adjust)

  end subroutine g_to_h_kxkyz

!   subroutine gstar_to_g (g, phi, apar, facphi, facapar)

!     use constants, only: zi
!     use species, only: spec
!     use zgrid, only: nzgrid
!     use vpamu_grids, only: vpa
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, is_idx
!     use kt_grids, only: naky, nakx
!     use kt_grids, only: aky
!     use gyro_averages, only: aj0x

!     implicit none
!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:), intent (in) :: phi, apar
!     real, intent (in) :: facphi, facapar

!     integer :: ivmu, iz, iky, ikx, is, iv
!     complex :: adj

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        iv = iv_idx(vmu_lo,ivmu)
!        is = is_idx(vmu_lo,ivmu)
!        do iz = -nzgrid, nzgrid
!           do ikx = 1, nakx
!              do iky = 1, naky
!                 ! BACKWARDS DIFFERENCE FLAG
! !                adj = zi*aj0x(iky,ikx,iz,ivmu)*aky(iky)*wstar(iz,ivmu) &
!                 adj = zi*aj0x(iky,ikx,iz,ivmu)*aky(iky)*2.0*wstar(1,iz,ivmu) &
!                      * ( facphi*phi(iky,ikx,iz) - facapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz) )
!                 g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) + adj
!              end do
!           end do
!        end do
!     end do

!   end subroutine gstar_to_g

end module g_tofrom_h
