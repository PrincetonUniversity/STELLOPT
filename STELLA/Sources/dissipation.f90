module dissipation

  implicit none

  public :: init_dissipation, finish_dissipation
  public :: init_collisions, collisions_initialized
  public :: include_collisions
  public :: advance_collisions_explicit, advance_collisions_implicit
  public :: time_collisions
  public :: hyper_dissipation
  public :: advance_hyper_dissipation
  public :: collisions_implicit
  public :: vpa_operator, mu_operator
  public :: cfl_dt_vpadiff, cfl_dt_mudiff
  public :: fieldpart

  private

  logical :: include_collisions, vpa_operator, mu_operator
  logical :: collisions_implicit
  logical :: momentum_conservation, energy_conservation
  logical :: hyper_dissipation
  logical :: use_physical_ksqr
  real :: D_hyper
  real :: cfl_dt_vpadiff, cfl_dt_mudiff
  logical :: density_conservation, density_conservation_field, density_conservation_tp, exact_conservation_tp, exact_conservation, spitzer_problem, no_j1l1, no_j1l2, no_j0l2
  logical :: fieldpart, testpart
  logical :: interspec, intraspec
  logical :: advfield_coll

  character(30) :: collision_model
  integer :: nresponse = 1
  integer :: nresponse_vpa = 1
  integer :: nresponse_mu = 1
  real :: cfac, cfac2
  real :: nuxfac
  real :: iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob
  logical :: eimassr_approx
  integer :: jmax = 1
  integer :: lmax = 1
  integer :: nvel_local

  real, dimension (:,:), allocatable :: aa_vpa, bb_vpa, cc_vpa
  real, dimension (:,:,:), allocatable :: aa_mu, cc_mu
  real, dimension (:,:), allocatable :: bb_mu
  complex, dimension (:,:,:), allocatable :: vpadiff_response
  integer, dimension (:,:), allocatable :: vpadiff_idx
  complex, dimension (:,:,:), allocatable :: mudiff_response
  integer, dimension (:,:), allocatable :: mudiff_idx
  complex, dimension (:,:,:), allocatable :: fp_response
  integer, dimension (:,:), allocatable :: diff_idx
  complex, dimension (:,:,:,:,:), allocatable :: response_vpamu

  complex, dimension (:,:,:), allocatable :: vpadiff_zf_response
  integer, dimension (:,:), allocatable :: vpadiff_zf_idx
  complex, dimension (:,:,:), allocatable :: mudiff_zf_response
  integer, dimension (:,:), allocatable :: mudiff_zf_idx

  complex, dimension (:,:,:,:,:), allocatable :: aa_blcs, cc_blcs
  complex, dimension (:,:,:,:,:), allocatable :: bb_blcs
  complex, dimension (:,:,:,:), allocatable :: dd_vecs
  complex, dimension (:,:,:,:,:,:), allocatable :: cdiffmat_band
  complex, dimension (:,:,:,:), allocatable :: blockmatrix
  complex, dimension (:,:,:), allocatable :: blockmatrix_sum
  integer, dimension (:,:,:,:,:), allocatable :: ipiv
  real, dimension (:,:,:,:,:), allocatable :: nus, nuD, nupa, nux
  real, dimension (:,:,:,:), allocatable :: mw, modmw
  real, dimension (:,:,:), allocatable :: velvpamu
  integer :: info

  real, dimension (:), allocatable :: wgts_v
  real, dimension (:), allocatable :: vel

  real, dimension (:,:), allocatable :: deltajs
  real, dimension (:,:), allocatable :: psijs
  real, dimension (:,:,:,:,:,:,:,:), allocatable :: deltaj, deltaj_tp
  real, dimension (:,:,:,:,:,:,:,:), allocatable :: deltajmod
  complex, dimension (:,:,:,:), allocatable :: deltajint
  real, dimension (:,:,:,:,:), allocatable :: psijnorm
  complex, dimension (:,:,:,:,:), allocatable :: deltajfield
  real, dimension (:,:,:,:,:), allocatable :: legendre_vpamu
  real, dimension (:,:,:,:,:,:), allocatable :: jm
  real, dimension (:,:,:,:,:), allocatable :: jm0
  real, dimension (:), allocatable :: mwnorm
  real, dimension (:), allocatable :: modmwnorm

  logical :: collisions_initialized = .false.
  real, dimension (2,2) :: time_collisions = 0.
  real :: i1fac, i2fac

contains

  subroutine init_dissipation

    use mp, only: proc0

    implicit none

    call read_parameters
    if (include_collisions) then
      if (collision_model == "dougherty") then
        write(*,*)
        write(*,*) 'Coll. model:     Dougherty'
        if (collisions_implicit) then
          write(*,*) 'Coll. algorithm: implicit'
        else
          write(*,*) 'Coll. algorithm: explicit'
        end if
      end if
      if (collision_model == "fokker-planck") then
        write(*,*) 'Coll. model:     Fokker-Planck'
        if (collisions_implicit) then
          write(*,*) 'Coll. algorithm: implicit'
        else
          write(*,*) 'Coll. algorithm: explicit'
        end if
      end if
      write(*,*)
    else
        if (proc0) then
           write (*,'(A)') "############################################################"
           write (*,'(A)') "                         COLLISIONS"
           write (*,'(A)') "############################################################"
           write (*,*) 'Coll. model:     None'
           write (*,*)
        end if
    end if

  end subroutine init_dissipation

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use physics_flags, only: full_flux_surface, radial_variation
    use mp, only: proc0, broadcast
    use run_parameters, only: fully_explicit

    implicit none

    namelist /dissipation/ collision_model, testpart, fieldpart, lmax, jmax, nvel_local, hyper_dissipation, D_hyper, &
         include_collisions, collisions_implicit, &
         interspec, intraspec, iiknob, ieknob, eeknob, eiknob, eiediffknob, eideflknob, deflknob, eimassr_approx, advfield_coll, spitzer_problem, &
         density_conservation, density_conservation_field, density_conservation_tp, exact_conservation, exact_conservation_tp, &
         momentum_conservation, energy_conservation, &
         momentum_conservation, energy_conservation, &
         vpa_operator, mu_operator, use_physical_ksqr, &
         cfac,  cfac2, nuxfac, i1fac, i2fac, no_j1l1, no_j1l2, no_j0l2

    integer :: in_file
    logical :: dexist

    if (proc0) then
       include_collisions = .false.
       collisions_implicit = .true.
       collision_model = "dougherty"        ! dougherty or fokker-planck
       !!! control parameters specific to the Fokker-Planck collision model
       testpart = .true.                    ! test particle component (TPO) of fokker-planck operator, must be True
       fieldpart = .false.                  ! enable the field particle component (FPO) of the fokker-planck operator
       intraspec = .true.                   ! intra-species collisions in the Fokker-Planck operator
       interspec = .true.                   ! inter-species
       iiknob = 1.                          ! control the ion-ion coll freq in Fokker-Planck operator
       ieknob = 1.                          ! ...ion-eon coll freq
       eeknob = 1.                          ! ...eon-eon coll freq
       eiknob = 1.                          ! ...eon-ion coll freq
       eiediffknob = 1.                     ! control the eon-ion energy diffusion in Fokker-Planck operator
       deflknob = 1.                        ! control pitch angle scattering in Fokker-Planck operator, must be 1 or 0
       eimassr_approx = .false.             ! use mass ratio approximation for test particle operator, beta
       advfield_coll = .true.               ! disable electrostatic potential terms in the field particle operator, beta
       density_conservation = .false.       ! if True and equally_spaced_mu_grid=True and conservative_wgts_vpa=True, then TPO conserves density to machine precision
       density_conservation_field = .false. ! if True and jmax, lmax < 2, then FPO conserves density to machine precision
       density_conservation_tp = .false.    ! if True add term to field particle operator to ensure density conservation, also on non-uniform grids
       exact_conservation = .false.         ! if True and fieldpart=True and lmax=jmax=1 then momentum and energy conserved to machine precision - in beta &
                                            ! & works only if nux = 0, need to correct the discretisation of nux terms in TPO
       exact_conservation_tp = .false.      ! if True and lmax=jmax=1 then momentum and energy conserved to machine precision, by using the test particle operator &
                                            ! to compute field particle terms; this is slower than exact_conservation
       momentum_conservation = .true.       ! momentum conservation for Dougherty operator
       energy_conservation = .true.         ! energy conservation for Dougherty operator
       spitzer_problem = .false.            ! to solve the Spitzer problem for tests of the collision operator
       cfac = 1                             ! scale gyrodiffusive term in test particle component of Fokker-Planck operator
       cfac2 = 1                            ! scale gyrodiffusive terms in field particle component of Fokker-Planck operator - in beta
       nuxfac = 1                           ! scale nux (mixed derivative) terms in test particle component of Fokker-Planck operator
       jmax = 1                             ! maximum j in Hirshman-Sigmar expansion of the field particle operator
       lmax = 1                             ! maximum l in spherical harmonic expansion of the field particle operator
       i1fac = 1                            ! for Spitzer problem
       i2fac = 0                            ! for Spitzer problem
       no_j1l1 = .true.                     ! disable j1l1 term in the field particle component of Fokker-Planck operator
       no_j1l2 = .false.                    ! disable j1l2 term
       no_j0l2 = .false.                    ! disable j0l2 term
       !!!
       vpa_operator = .true.                ! include vpa components in Dougherty or Fokker-Planck operator
       mu_operator = .true.                 ! include mu components in Dougherty or Fokker-Planck operator
       hyper_dissipation = .false.
       nvel_local = 512
       use_physical_ksqr = .not.(full_flux_surface.or.radial_variation)
       D_hyper = 0.05

       in_file = input_unit_exist("dissipation", dexist)
       if (dexist) read (unit=in_file, nml=dissipation)
    end if

    call broadcast (include_collisions)
    call broadcast (collisions_implicit)
    call broadcast (collision_model)
    call broadcast (fieldpart)
    call broadcast (testpart)
    call broadcast (interspec)
    call broadcast (intraspec)
    call broadcast (iiknob)
    call broadcast (ieknob)
    call broadcast (eeknob)
    call broadcast (eiknob)
    call broadcast (eiediffknob)
    call broadcast (deflknob)
    call broadcast (eimassr_approx)
    call broadcast (eideflknob)
    call broadcast (advfield_coll)
    call broadcast (density_conservation)
    call broadcast (density_conservation_field)
    call broadcast (density_conservation_tp)
    call broadcast (exact_conservation)
    call broadcast (exact_conservation_tp)
    call broadcast (momentum_conservation)
    call broadcast (energy_conservation)
    call broadcast (spitzer_problem)
    call broadcast (vpa_operator)
    call broadcast (mu_operator)
    call broadcast (hyper_dissipation)
    call broadcast (use_physical_ksqr)
    call broadcast (D_hyper)
    call broadcast (cfac)
    call broadcast (cfac2)
    call broadcast (nuxfac)
    call broadcast (jmax)
    call broadcast (lmax)
    call broadcast (nvel_local)
    call broadcast (i1fac)
    call broadcast (i2fac)
    call broadcast (no_j1l1)
    call broadcast (no_j1l2)
    call broadcast (no_j0l2)

    if (.not.include_collisions) collisions_implicit = .false.

  end subroutine read_parameters

  subroutine init_collisions

    use species, only: spec, nspec
    use vpamu_grids, only: dvpa, dmu, mu, nmu
!   use vpamu_grids, only: calculate_velocity_integrals
    use stella_geometry, only: bmag
    use stella_layouts
    use run_parameters, only: fully_explicit
    use common_types, only: spec_type

    implicit none

    integer :: is, is2
    integer, parameter :: ion_species = 1
    integer, parameter :: electron_species = 2
    integer, parameter :: impurity_species = 3 ! AVB: clear up difference between slowing down species ('3' in species.f90) and impurity species
    real :: vnew_max


    if (collisions_initialized) return
    collisions_initialized = .true.

!   call calculate_velocity_integrals

    if (collision_model == "dougherty") then
        if (collisions_implicit) then
           if (vpa_operator) then
              call init_vpadiff_matrix
              call init_vpadiff_conserve
           end if
           if (mu_operator) then
              call init_mudiff_matrix
              call init_mudiff_conserve
           end if
        else
           vnew_max = 0.0
           do is = 1, nspec
              vnew_max = max(vnew_max,maxval(spec(is)%vnew))
           end do
           cfl_dt_vpadiff = 2.0*dvpa**2/vnew_max
           cfl_dt_mudiff = minval(bmag)/(vnew_max*maxval(mu(2:)/dmu(:nmu-1)**2))
        end if
    end if

    if (collision_model == "fokker-planck") then
        ! disable inter-species collisions if interspec==false
        if (.not.interspec) then
            do is = 1, nspec
                do is2 = 1, nspec
                    if (is/=is2) then
                        spec(is)%vnew(is2) = 0.
                    end if
                end do
            end do
        end if

        ! disable intra-species collisions if intraspec==false
        if (.not.intraspec) then
            do is = 1, nspec
                do is2 = 1, nspec
                    if (is==is2) then
                        spec(is)%vnew(is2) = 0.
                    end if
                end do
            end do
        end if

        ! control inter-species collisions
        do is = 1, nspec
            do is2 = 1, nspec
                if ((spec(is)%type==ion_species).and.(spec(is2)%type==ion_species)) then
                    spec(is)%vnew(is2) = spec(is)%vnew(is2) * iiknob
                else if ((spec(is)%type==ion_species).and.(spec(is2)%type==electron_species)) then
                    spec(is)%vnew(is2) = spec(is)%vnew(is2) * ieknob
                else if ((spec(is)%type==electron_species).and.(spec(is2)%type==electron_species)) then
                    spec(is)%vnew(is2) = spec(is)%vnew(is2) * eeknob
                else if ((spec(is)%type==electron_species).and.(spec(is2)%type==ion_species)) then
                    spec(is)%vnew(is2) = spec(is)%vnew(is2) * eiknob
                else
                    spec(is)%vnew(is2) = spec(is)%vnew(is2)
                end if
                ! AVB: to do - add impurity collision control
            end do
        end do

        ! initialise speed dependent collision frequencies
        call init_nusDpa

        if (collisions_implicit) then
            write(*,*) 'Coll. algorithm: implicit'
            fully_explicit = .false.

            call init_legendre
            call init_vgrid
            call init_bessel_fn
            call init_fp_diffmatrix
            call init_deltaj_vmu
            call init_fp_conserve
        else
           vnew_max = 0.0
           do is = 1, nspec
              vnew_max = max(vnew_max,maxval(spec(is)%vnew))
           end do
           cfl_dt_vpadiff = 2.0*dvpa**2/vnew_max
           cfl_dt_mudiff = minval(bmag)/(vnew_max*maxval(mu(2:)/dmu(:nmu-1)**2))
        end if
    end if
  end subroutine init_collisions

  subroutine init_nusDpa

      ! AVB: compute the collision frequencies nuD, nus and nupa

      use zgrid, only: nzgrid
      use vpamu_grids, only: nmu, mu, dmu, vpa, nvpa, integrate_vmu
      use stella_geometry, only: bmag
      use species, only: spec, nspec
      use spfunc, only: erf => erf_ext
      use finite_differences, only: fd3pt
      use vpamu_grids, only: maxwell_mu, maxwell_vpa
      use constants, only: pi

      implicit none

      real, dimension (-nzgrid:nzgrid) :: v2mwint, v4mwint
      integer :: ia, imu, iv, iz, is, isb
      real :: x, Gf, massr

      if (.not.allocated(nus)) allocate (nus(nvpa,nmu,-nzgrid:nzgrid,nspec,nspec))
      if (.not.allocated(nuD)) allocate (nuD(nvpa,nmu,-nzgrid:nzgrid,nspec,nspec))
      if (.not.allocated(nupa)) allocate (nupa(nvpa,nmu,-nzgrid:nzgrid,nspec,nspec))
      if (.not.allocated(nux)) allocate (nux(nvpa,nmu,-nzgrid:nzgrid,nspec,nspec))
      if (.not.allocated(mw)) allocate (mw(nvpa,nmu,-nzgrid:nzgrid,nspec))
      if (.not.allocated(modmw)) allocate (modmw(nvpa,nmu,-nzgrid:nzgrid,nspec))
      if (.not.allocated(velvpamu)) allocate (velvpamu(nvpa,nmu,-nzgrid:nzgrid))

      ia = 1

      do is = 1, nspec
          do isb = 1, nspec

              massr = spec(is)%mass / spec(isb)%mass

              do iz = -nzgrid, nzgrid
                  do iv = 1, nvpa
                      do imu = 1, nmu
                          x = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))
                          Gf = (erf(x/sqrt(massr))-x/sqrt(massr)*(2/sqrt(pi))*exp(-x**2/massr)) / (2*x**2/massr)
                          nuD(iv,imu,iz,is,isb) = deflknob*spec(is)%vnew(isb)*(erf(x/sqrt(massr))-Gf)/x**3
                          nus(iv,imu,iz,is,isb) = spec(is)%vnew(isb)*2*(1+1./massr)*Gf/x ! nus is never used; note - have assumed T_a = T_b here
                          nupa(iv,imu,iz,is,isb)= spec(is)%vnew(isb)*2*Gf/x**3
                          velvpamu(iv,imu,iz)   = x
                          mw(iv,imu,iz,is)      = maxwell_vpa(iv,is)*maxwell_mu(1,iz,imu,is)
                      end do
                  end do
              end do

              ! electron-ion collisions
              ! approximation of Lorentz operator using mass-ratio expansion
              if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                  do iz = -nzgrid, nzgrid
                      do iv = 1, nvpa
                          do imu = 1, nmu
                              x = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))
                              nuD(iv,imu,iz,is,isb) = deflknob*spec(is)%vnew(isb)*1./x**3
                          end do
                      end do
                  end do
              end if

          end do
      end do

      ! get a function with vanishing energy moment, modmw
      do is = 1, nspec
          do iz = -nzgrid, nzgrid
              call integrate_vmu(velvpamu(:,:,iz)**2*mw(:,:,iz,is), iz, v2mwint(iz))
              call integrate_vmu(velvpamu(:,:,iz)**4*mw(:,:,iz,is), iz, v4mwint(iz))
              modmw(:,:,iz,is) = mw(:,:,iz,is) - velvpamu(:,:,iz)**2 * v2mwint(iz)/v4mwint(iz) * mw(:,:,iz,is)
          end do
      end do

      nux = nuxfac*(nupa - deflknob*nuD)

      if (nspec>1) then
          ! eiediffknob controls e-i energy diffusion, note that it is also used in blockmatrix
          nux(:,:,:,2,1) = nuxfac*(eiediffknob*nupa(:,:,:,2,1) - eideflknob*deflknob*nuD(:,:,:,2,1))
          nuD(:,:,:,2,1) = eideflknob*nuD(:,:,:,2,1)
      end if

  end subroutine init_nusDpa

  subroutine finish_nusDpa

      implicit none

      if (allocated(nus)) deallocate (nus)
      if (allocated(nuD)) deallocate (nuD)
      if (allocated(nupa)) deallocate (nupa)
      if (allocated(nux)) deallocate (nux)
      if (allocated(mw)) deallocate (mw)
      if (allocated(modmw)) deallocate (modmw)

      if (allocated(velvpamu)) deallocate (velvpamu)

  end subroutine finish_nusDpa

  subroutine init_fp_diffmatrix

    use stella_time, only: code_dt
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: dvpa, vpa, nvpa, mu, nmu, maxwell_mu, maxwell_vpa, dmu, equally_spaced_mu_grid
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use stella_geometry, only: bmag
    use dist_fn_arrays, only: kperp2
    use physics_parameters, only: zeff
    use constants, only: pi
    use common_types, only: spec_type
    use kt_grids, only: naky, nakx
    use spfunc, only: erf => erf_ext
    use file_utils, only: open_output_file, close_output_file

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is, isb, tmpunit
    integer :: imu, ia, iv, ivv, imm, imu2
    integer :: nc, nb, lldab, colbl, rowbl, bm_colind, bm_rowind
    real :: vpap, vpam, vfac, mum, mup
    real :: xpv, xmv, nupapv, nupamv, nuDpv, nuDmv, mwpv, mwmv, gam_mu, gam_mum, gam_mup
    real :: mwm, mwp, nuDm, nuDp, nupam, nupap, xm, xp
    real :: nuDfac, massr, eiediff, eidefl

    integer, parameter :: ion_species = 1
    integer, parameter :: electron_species = 2

    if (.not.allocated(aa_blcs)) allocate (aa_blcs(nvpa,nmu,nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nspec))
    if (.not.allocated(bb_blcs)) allocate (bb_blcs(nvpa,nmu,nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nspec))
    if (.not.allocated(cc_blcs)) allocate (cc_blcs(nvpa,nmu,nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nspec))
    if (.not.allocated(cdiffmat_band)) allocate (cdiffmat_band(3*(nmu+1)+1,nmu*nvpa,naky,nakx,-nzgrid:nzgrid,nspec))
    if (.not.allocated(ipiv)) allocate (ipiv(nvpa*nmu,naky,nakx,-nzgrid:nzgrid,nspec))

    ! AVB: calculate the discretisation matrix -\Delta t C_test^{ab}
    ! because of mixed vpa-mu derivatives in the test particle operator
    ! this matrix is block tri-diagonal, with dimension nmu*nvpa x nmu*nvpa
    ! store and operate with the matrix in band format

    ! aa_blcs stores subdiagonal blocks, bb_blcs diagonal blocks and cc_blcs superdiagonal blocks
    ! aa_blcs(1,:,:) and cc_blcs(nvpa,:,:) are never used
    ! mu-derivatives are contained within blocks, thus blocks have dimension nmu x nmu

    ia = 1
    vfac = 1 ! zero vpa-operator, in beta
    nuDfac = 1.
    aa_blcs = 0.
    bb_blcs = 0.
    cc_blcs = 0.

    do imu = 1, nmu
        do iv = 1, nvpa
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo,ikxkyz)
               ikx = ikx_idx(kxkyz_lo,ikxkyz)
               iz = iz_idx(kxkyz_lo,ikxkyz)
               is = is_idx(kxkyz_lo,ikxkyz)

               if (spitzer_problem) then
                  if (.not.(spec(is)%type==electron_species)) cycle ! add eon-eon and eon-ion collisions only for Spitzer problem
               end if

               do isb = 1, nspec
                    ! for Spitzer problem, disable e-i energy diffusion if eiediffknob = 0.
                    if (spitzer_problem) then
                        if ((is==2).and.(isb==1)) then
                            eiediff = eiediffknob
                            eidefl  = eideflknob
                        else
                            eiediff = 1.
                            eidefl  = 1.
                        end if
                    else
                        if ((is==2).and.(isb==1)) then
                            eiediff = eiediffknob
                            eidefl  = eideflknob
                        else
                            eiediff = 1.
                            eidefl  = 1.
                        end if
                    end if

                    massr = spec(is)%mass / spec(isb)%mass

                    if (iv == 1) then
                         vpap   = 0.5*(vpa(iv)+vpa(iv+1))
                         mwpv   = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
                         xpv    = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
                         nupapv = vfac*spec(is)%vnew(isb)*2*(erf(xpv/sqrt(massr))&
                            -xpv/sqrt(massr)*(2/sqrt(pi))*exp(-(xpv/sqrt(massr))**2)) / (2*(xpv/sqrt(massr))**2)/xpv**3

                         if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                             nuDpv  = vfac*spec(is)%vnew(isb)/xpv**3
                         else
                             nuDpv  = vfac*spec(is)%vnew(isb)*(erf(xpv/sqrt(massr))-(erf(xpv/sqrt(massr))&
                                -(xpv/sqrt(massr))*(2/sqrt(pi))*exp(-(xpv/sqrt(massr))**2)) / (2*(xpv/sqrt(massr))**2))/xpv**3
                         end if

                         if (imu == 1) then
                             ! one-sided difference for mu-derivative at imu=1:
                             !cc_blcs(iv,imu,imu+1,iz,is) = cc_blcs(iv,imu,imu+1,iz,is) &
                             !  - vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                             !cc_blcs(iv,imu,imu,iz,is)   = cc_blcs(iv,imu,imu,iz,is) &
                             !  - code_dt*0.5*(eiediff*nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                             !                                +vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu,iz,is) / (2*dvpa) / dmu(imu)
                             ! using ghost cell at mu=0:
                             cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                - 0.5*vfac*code_dt*vpa(iv+1)*0.5*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                /mw(iv+1,imu+1,iz,is) / (dvpa) / dmu(imu)
                             cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                + 0.5*vfac*code_dt*vpa(iv+1)*0.5*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                /mw(iv+1,imu  ,iz,is) / (dvpa) / dmu(imu)
                             bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                - 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                /mw(iv,imu+1,iz,is) / (dvpa) / dmu(imu)
                             bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                + 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                /mw(iv,imu  ,iz,is) / (dvpa) / dmu(imu)
                             !
                             bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv) / dvpa**2 / mw(iv,imu,iz,is)

                             ! mu operator
                             if (mu_operator) then
                                 ! use ghost cell at mu_{0} = 0
                                 mup = 0.5*(mu(imu)+mu(imu+1))
                                 mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                 xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb) / xp**3
                                 else
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xp/sqrt(massr))&
                                        -(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                 end if
                                 nupap = spec(is)%vnew(isb)*2*(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr) / xp**3
                                 gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mup = 2*(eiediff*nupap*mup**2+eidefl*deflknob*nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp
                                 bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*(gam_mu*-1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                    +(gam_mup*-1/dmu(imu) - gam_mu*-1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu,iz,is) / (dmu(imu)/2.+mu(imu))
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) - code_dt*(gam_mu * 1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                    +(gam_mup*1/dmu(imu) - gam_mu*1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu+1,iz,is) / (dmu(imu)/2.+mu(imu))
                                 ! mixed derivative:
                                 if (density_conservation) then
                                     ! to ensure density conservation, change discretisation of mixed derivative term at imu=1
                                     ! from explicit routine: Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah &
                                     !  + vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz,is)*Dvpah_p) / (2.*dmu(imu))
                                     cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )&
                                        *nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))*1/mw(iv+1,imu  ,iz,is)) / (2*dmu(imu)) / dvpa
                                     cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)+vpa(iv+1)*mu(imu+1)&
                                        *nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))*1/mw(iv+1,imu+1,iz,is)) / (2*dmu(imu)) / dvpa
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )&
                                        *nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))*1/mw(iv  ,imu  ,iz,is)) / (2*dmu(imu)) / dvpa
                                     bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)+vpa(iv+1)*mu(imu+1)&
                                        *nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))*1/mw(iv  ,imu+1,iz,is)) / (2*dmu(imu)) / dvpa
                                 else
                                     cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                        / (mu(imu)+dmu(imu)) / dvpa
                                     cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv+1,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                        / (mu(imu)+dmu(imu)) / dvpa
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                        / (mu(imu)+dmu(imu)) / dvpa
                                     bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                        / (mu(imu)+dmu(imu)) / dvpa
                                 end if
                             end if
                         else if (imu == nmu) then
                             ! AVB: one-sided difference for mu-derivative at imu=nmu:
                             !cc_blcs(iv,imu,imu-1,iz,is) = cc_blcs(iv,imu,imu-1,iz,is) &
                             !  + 1.000*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                             !cc_blcs(iv,imu,imu,iz,is)   = cc_blcs(iv,imu,imu,iz,is) &
                             !  - code_dt*0.5*(eiediff*nupapv*vpap**2 + 2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                             !  - 1.000*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu,iz,is) / (2*dvpa) / dmu(nmu-1)
                             ! AVB: second-order:
                             if (density_conservation) then
                                cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + 0.5*vfac*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                    /mw(iv+1,imu-1,iz,is) / (dvpa) / dmu(nmu-1)
                                cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    - 0.5*vfac*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                    /mw(iv+1,imu  ,iz,is) / (dvpa) / dmu(nmu-1)
                                bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                    /mw(iv,imu-1,iz,is) / (dvpa) / dmu(nmu-1)
                                bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                    /mw(iv,imu  ,iz,is) / (dvpa) / dmu(nmu-1)
                             else
                                cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + 1.0*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    - 1.0*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu,iz,is) / (2*dvpa) / dmu(nmu-1)
                             end if
                             bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv) / dvpa**2 / mw(iv,imu,iz,is)
                             ! mu operator
                             if (mu_operator) then
                                 mum = 0.5*(mu(imu)+mu(imu-1))
                                 mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                 xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb) / xm**3
                                 else
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))&
                                        -(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                                 end if
                                 nupam = spec(is)%vnew(isb)*2*(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr) / xm**3
                                 gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mum = 2*(eiediff*nupam*mum**2+eidefl*deflknob*nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                 if (density_conservation) then
                                     ! to ensure density conservation we assume that the argument of the outer derivative vanishes at nmu+1/2, ie
                                     ! d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2)
                                     ! where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)  = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                                         + (-gam_mu/dmu(imu-1))*dmu(imu-1)/dmu(imu-1)) / mw(iv,imu,iz,is) / (dmu(imu-1)/2.+dmu(imu-1)/2.)
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb)= bb_blcs(iv,imu,imu-1,ikxkyz,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                                         -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/dmu(imu-1)) / mw(iv,imu-1,iz,is) / (dmu(imu-1)/2.+dmu(imu-1)/2.)
                                 else
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)  = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                         + (-gam_mu/dmu(imu-1))*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu,iz,is) / (dmu(imu-1)/2.+dmu(imu-1))
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb)= bb_blcs(iv,imu,imu-1,ikxkyz,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                         -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu-1,iz,is) / (dmu(imu-1)/2.+dmu(imu-1))
                                 end if
                                 ! mixed derivative:
                                 !cc_blcs(iv,imu,imu,iz,is)   = cc_blcs(iv,imu,imu  ,iz,is)  &
                                 ! - 1.000*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 !cc_blcs(iv,imu,imu-1,iz,is) = cc_blcs(iv,imu,imu-1,iz,is)  &
                                 ! + 1.000*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 !bb_blcs(iv,imu,imu,ikxkyz)  = bb_blcs(iv,imu,imu  ,ikxkyz) &
                                 !+ 1.000*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 !bb_blcs(iv,imu,imu-1,ikxkyz)= bb_blcs(iv,imu,imu-1,ikxkyz) &
                                 !- 1.000*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 if (density_conservation) then
                                     cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                        *1/mw(iv+1,imu  ,iz,is)) / (2*dmu(imu-1)) / dvpa
                                     cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                        *1/mw(iv+1,imu-1,iz,is)) / (2*dmu(imu-1)) / dvpa
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                        *1/mw(iv  ,imu  ,iz,is)) / (2*dmu(imu-1)) / dvpa
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                        *1/mw(iv  ,imu-1,iz,is)) / (2*dmu(imu-1)) / dvpa
                                 else
                                     cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                     cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 end if

                             end if
                         else
                             ! interior mu points for iv = 1
                             ! use ghost cell for derivative
                             ! d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                             ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa
                             ! could be cleared up, by moving non-nux part of aa_blcs down to non-nux part of bb_bcls
                             cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                + 0.5*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is)*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                - 0.5*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is)*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu,ikxkyz,isb)   &
                                - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                - 0.5*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                + 0.5*vfac*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu-1,iz,is)*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                - 0.5*vfac*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu+1,iz,is)*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu,ikxkyz,isb)   &
                                - 0.5*vfac*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)

                             bb_blcs(iv,imu,imu,ikxkyz,isb)  = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv) / dvpa**2 / mw(iv,imu,iz,is)

                             ! mu operator
                             if (mu_operator) then
                                 ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                                 mup = 0.5*(mu(imu)+mu(imu+1))
                                 mum = 0.5*(mu(imu)+mu(imu-1))
                                 mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                 mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                 xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                 xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb) / xp**3
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb) / xm**3
                                 else
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xp/sqrt(massr))&
                                        -(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))&
                                        -(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                                 end if
                                 nupap = spec(is)%vnew(isb)*2*(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr) / xp**3
                                 nupam = spec(is)%vnew(isb)*2*(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr) / xm**3
                                 gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mum = 2*(eiediff*nupam*mum**2+eidefl*deflknob*nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                 gam_mup = 2*(eiediff*nupap*mup**2+eidefl*deflknob*nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                 bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                    +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu))&
                                    / mw(iv,imu,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                                 bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - code_dt*((gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) - gam_mum* -1/dmu(imu-1)) * dmu(imu)/dmu(imu-1) &
                                    -gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) * dmu(imu-1)/dmu(imu)) / mw(iv,imu-1,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - code_dt*( gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) * dmu(imu)/dmu(imu-1) &
                                    + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu))&
                                    / mw(iv,imu+1,iz,is) * 2/(dmu(imu-1)+dmu(imu))
                                 ! mixed derivative:
                                 cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                    *1/mw(iv+1,imu  ,iz,is)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                    *1/mw(iv+1,imu-1,iz,is)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)+vpa(iv+1)*mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                    *1/mw(iv+1,imu+1,iz,is)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                    *1/mw(iv  ,imu  ,iz,is)*(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                    *1/mw(iv  ,imu-1,iz,is)* dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)+vpa(iv+1)*mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                    *1/mw(iv  ,imu+1,iz,is)* dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             end if
                         end if

                     else if (iv == nvpa) then
                         vpam = 0.5*(vpa(iv)+vpa(iv-1))
                         mwmv = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
                         xmv = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
                         nupamv = vfac*spec(is)%vnew(isb)*2*(erf(xmv/sqrt(massr))-xmv/sqrt(massr)*(2/sqrt(pi))*exp(-xmv**2/massr)) / (2*xmv**2/massr)/xmv**3
                         if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                             nuDmv = eidefl*deflknob*vfac*spec(is)%vnew(isb) / xmv**3
                         else
                             nuDmv = eidefl*deflknob*vfac*spec(is)%vnew(isb)*(erf(xmv/sqrt(massr))&
                                -(erf(xmv/sqrt(massr))-xmv/sqrt(massr)*(2/sqrt(pi))*exp(-xmv**2/massr)) / (2*xmv**2/massr))/xmv**3
                         end if
                         if (imu == 1) then
                             ! one-sided difference for mu-derivative at imu=1:
                             ! ie d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                             ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa
                             ! could be cleared up, by moving non-nux part of aa_blcs down to non-nux part of bb_bcls
                             aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                + 0.5*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                /mw(iv-1,imu+1,iz,is) / (dvpa) / dmu(imu)
                             aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) - code_dt*0.5*(eiediff*nupamv*vpam**2 &
                                + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                - 0.5*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                /mw(iv-1,imu  ,iz,is) / (dvpa) / dmu(imu)
                             bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                + 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                /mw(iv,imu+1,iz,is) / (dvpa) / dmu(imu)
                             bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                - 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                /mw(iv,imu  ,iz,is) / (dvpa) / dmu(imu)

                             bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2  / mw(iv,imu,iz,is)

                             ! mu operator
                             if (mu_operator) then
                                 ! use ghost cell at mu_{0} = 0
                                 mup = 0.5*(mu(imu)+mu(imu+1))
                                 mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                 xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb) / xp**3
                                 else
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xp/sqrt(massr)) &
                                     -(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                 end if
                                 nupap = spec(is)%vnew(isb)*2*(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr) / xp**3
                                 gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mup = 2*(eiediff*nupap*mup**2+eidefl*deflknob*nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp

                                 bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*(gam_mu *-1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                    +(gam_mup*-1/dmu(imu) - gam_mu*-1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu,iz,is) / (dmu(imu)/2.+mu(imu))
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) - code_dt*(gam_mu * 1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                    +(gam_mup* 1/dmu(imu) - gam_mu*1/dmu(imu)) * mu(imu)/(dmu(imu)/2.))  / mw(iv,imu+1,iz,is) / (dmu(imu)/2.+mu(imu))
                                 ! mixed derivative:
                                 if (density_conservation) then
                                     ! to ensure density conservation, change discretisation of mixed derivative term at imu=1
                                     ! from explicit routine: Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah + vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz,is)*Dvpah_p) / (2.*dmu(imu))
                                     aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                        *1/mw(iv-1,imu  ,iz,is)) / (2*dmu(imu)) / (dvpa)
                                     aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                        *1/mw(iv-1,imu+1,iz,is)) / (2*dmu(imu)) / (dvpa)
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                        *1/mw(iv  ,imu  ,iz,is)) / (2*dmu(imu)) / (dvpa)
                                     bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                        *1/mw(iv  ,imu+1,iz,is)) / (2*dmu(imu)) / (dvpa)
                                 else
                                     aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                        / (mu(imu)+dmu(imu)) / (dvpa)
                                     aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv-1,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                        / (mu(imu)+dmu(imu)) / (dvpa)
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                        / (mu(imu)+dmu(imu)) / (dvpa)
                                     bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                        / (mu(imu)+dmu(imu)) / (dvpa)
                                 end if
                             end if
                         else if (imu == nmu) then
                             ! one-sided difference for mu-derivative at imu=nmu:
                             !aa_blcs(iv,imu,imu-1,iz,is) = aa_blcs(iv,imu,imu-1,iz,is) - 1.000*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                             !aa_blcs(iv,imu,imu,iz,is)   = aa_blcs(iv,imu,imu,iz,is) - code_dt*0.5*(eiediff*nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2  / mw(iv-1,imu,iz,is) &
                             !  + 1.000*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu,iz,is) / (2*dvpa) / dmu(nmu-1)
                             ! 24.02.21, second order:
                             if (density_conservation) then
                                 !aa_blcs(iv,imu,imu-1,iz,is) = aa_blcs(iv,imu,imu-1,iz,is) &
                                 !- 1.000*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))/mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                 !aa_blcs(iv,imu,imu,iz,is)   = aa_blcs(iv,imu,imu  ,iz,is) - code_dt*0.5*(eiediff*nupamv*vpam**2 + 2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                 ! + 1.000*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))/mw(iv-1,imu  ,iz,is) / (2*dvpa) / dmu(nmu-1)
                                 ! take vpa derivative using ghost cell at nvpa+0.5
                                 ! ie d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                                 ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa
                                 ! note this could be cleared up, by moving the non-nux part of aa_blcs down to the non-nux part of bb_bcls
                                 aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - 0.5*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                    /mw(iv-1,imu-1,iz,is) / (dvpa) / dmu(nmu-1)
                                 aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                    + 0.5*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                    /mw(iv-1,imu  ,iz,is) / (dvpa) / dmu(nmu-1)
                                 bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                    /mw(iv,imu-1,iz,is) / (dvpa) / dmu(nmu-1)
                                 bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    + 0.5*vfac*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                    /mw(iv,imu  ,iz,is) / (dvpa) / dmu(nmu-1)
                             else
                                 aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - 1.0*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                 aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu,ikxkyz,isb)   &
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2  / mw(iv-1,imu,iz,is) &
                                    + 1.0*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu,iz,is) / (2*dvpa) / dmu(nmu-1)
                             end if
                             bb_blcs(iv,imu,imu,ikxkyz,isb)  = bb_blcs(iv,imu,imu,ikxkyz,isb)&
                                + code_dt*(0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz,is)

                             ! mu operator
                             if (mu_operator) then
                                 mum = 0.5*(mu(imu)+mu(imu-1))
                                 mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                 xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb) / xm**3
                                 else
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))&
                                        -(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                                 end if
                                 nupam = spec(is)%vnew(isb)*2*(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr) / xm**3
                                 gam_mu  = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mum = 2*(eiediff*nupam*mum**2+eidefl*deflknob*nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                 if (density_conservation) then
                                     ! to ensure density conservation we assume that the argument of the outer derivative vanishes at nmu+1/2, ie
                                     ! d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2), where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
                                     bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                                        + (-gam_mu/dmu(imu-1))*dmu(imu-1)/dmu(imu-1)) / mw(iv,imu,iz,is) / (dmu(imu-1)/2.+dmu(imu-1)/2.)
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                                        -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/dmu(imu-1)) / mw(iv,imu-1,iz,is) / (dmu(imu-1)/2.+dmu(imu-1)/2.)
                                 else
                                     bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                        + (-gam_mu/dmu(imu-1))*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu,iz,is) / (dmu(imu-1)/2.+dmu(imu-1))
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                        -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu-1,iz,is) / (dmu(imu-1)/2.+dmu(imu-1))
                                 end if
                                 ! mixed derivative:
                                 ! 24.02.21, commenting out
                                 !aa_blcs(iv,imu,imu,iz,is)    = aa_blcs(iv,imu,imu  ,iz,is) &
                                 !  + 1.000*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 !aa_blcs(iv,imu,imu-1,iz,is)  = aa_blcs(iv,imu,imu-1,iz,is) &
                                 !  - 1.000*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 !bb_blcs(iv,imu,imu,ikxkyz)   = bb_blcs(iv,imu,imu  ,ikxkyz) &
                                 !  - 1.000*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 !bb_blcs(iv,imu,imu-1,ikxkyz) = bb_blcs(iv,imu,imu-1,ikxkyz) &
                                 !  + 1.000*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 if (density_conservation) then
                                     ! 04.03. removed an extra bracket at end of lines in following block
                                     aa_blcs(iv,imu,imu-1,ikxkyz,isb)= aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                        *1/mw(iv-1,imu-1,iz,is)) / (2*dmu(imu-1)) / dvpa
                                     aa_blcs(iv,imu,imu  ,ikxkyz,isb)= aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                        *1/mw(iv-1,imu  ,iz,is)) / (2*dmu(imu-1)) / dvpa
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb)= bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                        *1/mw(iv  ,imu-1,iz,is)) / (2*dmu(imu-1)) / dvpa
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)  = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                        *1/mw(iv  ,imu  ,iz,is)) / (2*dmu(imu-1)) / dvpa
                                 else
                                     aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                     aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                     bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                     bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        + 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                        / (dmu(imu-1)+dmu(imu-1)) / (dvpa)
                                 end if
                             end if
                         else
                             ! interior mu points for iv=nvpa
                             ! take vpa derivative using ghost cell at nvpa+0.5
                             ! ie d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                             ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa

                             ! note this could be cleared up, by moving the non-nux part of aa_blcs down to the non-nux part of bb_bcls
                             aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb)&
                                - 0.5*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb)&
                                + 0.5*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             aa_blcs(iv,imu,imu,ikxkyz,isb) = aa_blcs(iv,imu,imu  ,ikxkyz,isb)&
                                - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                + 0.5*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb)&
                                - 0.5*vfac*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb)&
                                + 0.5*vfac*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu  ,ikxkyz,isb)&
                                + 0.5*vfac*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)

                             bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb)&
                                + code_dt*(0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz,is)

                             ! mu operator
                             if (mu_operator) then
                                 ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                                 mup = 0.5*(mu(imu)+mu(imu+1))
                                 mum = 0.5*(mu(imu)+mu(imu-1))
                                 mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                 mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                 xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                 xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb) / xp**3
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb) / xm**3
                                 else
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xp/sqrt(massr))-(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))-(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                                 end if
                                 nupap = spec(is)%vnew(isb)*2*(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr) / xp**3
                                 nupam = spec(is)%vnew(isb)*2*(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr) / xm**3
                                 gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mum = 2*(eiediff*nupam*mum**2+eidefl*deflknob*nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                 gam_mup = 2*(eiediff*nupap*mup**2+eidefl*deflknob*nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp
                                 bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                    +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu))&
                                    / mw(iv,imu,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                                 bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - code_dt*((gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) - gam_mum* -1/dmu(imu-1)) * dmu(imu)/dmu(imu-1) &
                                    - gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) * dmu(imu-1)/dmu(imu)) / mw(iv,imu-1,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - code_dt*( gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) * dmu(imu)/dmu(imu-1) &
                                    + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu))&
                                    / mw(iv,imu+1,iz,is) * 2/(dmu(imu-1)+dmu(imu))
                                 ! mixed derivative, one-sided difference in vpa at iv = nvpa:
                                 ! use second order accurate treatment for vpa derivative at nvpa
                                 aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb)&
                                    + 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                    *1/mw(iv-1,imu  ,iz,is)*(dmu(imu  )/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb)&
                                    - 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                    *1/mw(iv-1,imu-1,iz,is)* dmu(imu  )/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb)&
                                    + 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                    *1/mw(iv-1,imu+1,iz,is)* dmu(imu-1)/dmu(imu  )) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu  ,ikxkyz,isb)&
                                    - 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                    *1/mw(iv  ,imu  ,iz,is)*(dmu(imu  )/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb)&
                                    + 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                    *1/mw(iv  ,imu-1,iz,is)* dmu(imu  )/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb)&
                                    - 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                    *1/mw(iv  ,imu+1,iz,is)* dmu(imu-1)/dmu(imu  )) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                             end if
                         end if

                     else ! interior vpa points
                         vpam = 0.5*(vpa(iv)+vpa(iv-1))
                         vpap = 0.5*(vpa(iv)+vpa(iv+1))
                         mwmv = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
                         mwpv = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
                         xpv = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
                         xmv = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
                         nupamv = vfac*spec(is)%vnew(isb)*2*(erf(xmv/sqrt(massr))-xmv/sqrt(massr)*(2/sqrt(pi))*exp(-xmv**2/massr)) / (2*xmv**2/massr)/xmv**3
                         nupapv = vfac*spec(is)%vnew(isb)*2*(erf(xpv/sqrt(massr))-xpv/sqrt(massr)*(2/sqrt(pi))*exp(-xpv**2/massr)) / (2*xpv**2/massr)/xpv**3
                         if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                             nuDmv = eidefl*deflknob*vfac*spec(is)%vnew(isb)/xmv**3
                             nuDpv = eidefl*deflknob*vfac*spec(is)%vnew(isb)/xpv**3
                         else
                             nuDmv = eidefl*deflknob*vfac*spec(is)%vnew(isb)*(erf(xmv/sqrt(massr))&
                                -(erf(xmv/sqrt(massr))-xmv/sqrt(massr)*(2/sqrt(pi))*exp(-xmv**2/massr)) / (2*xmv**2/massr))/xmv**3
                             nuDpv = eidefl*deflknob*vfac*spec(is)%vnew(isb)*(erf(xpv/sqrt(massr))&
                                -(erf(xpv/sqrt(massr))-xpv/sqrt(massr)*(2/sqrt(pi))*exp(-xpv**2/massr)) / (2*xpv**2/massr))/xpv**3
                         end if

                         if (imu == 1) then
                              ! one-sided difference for mu-derivative at imu=1:
                              if (.not.density_conservation) then
                                  aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb)&
                                    + vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                                  aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb)&
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                    - vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb) / (2*dvpa) / dmu(imu)
                                  cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb)&
                                    - vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                                  cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb)&
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    + vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb) / (2*dvpa) / dmu(imu)
                              else if (density_conservation) then
                                   ! to ensure density conservation, assume that nux*F0 vanishes at iv=1 and iv=nvpa
                                   aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                     + 0.5*vfac*code_dt*vpa(iv-1)*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                     /mw(iv-1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                                   aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                     - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                     - 0.5*vfac*code_dt*vpa(iv-1)*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                     /mw(iv-1,imu,iz,is) / (2*dvpa) / dmu(imu)
                                   cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                     - 0.5*vfac*code_dt*vpa(iv+1)*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                     /mw(iv+1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                                   cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                     - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                     + 0.5*vfac*code_dt*vpa(iv+1)*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                     /mw(iv+1,imu  ,iz,is) / (2*dvpa) / dmu(imu)
                              end if
                              bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv &
                                + 0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz,is)

                              ! mu operator
                              if (mu_operator) then
                                  ! use ghost cell at mu_{0} = 0, where term vanishes, so dmu(0) = mu(1).
                                  mup = 0.5*(mu(imu)+mu(imu+1))
                                  mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                  xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                  if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                      nuDp = eidefl*deflknob*spec(is)%vnew(isb) / xp**3
                                  else
                                      nuDp = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xp/sqrt(massr))&
                                        -(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                  end if
                                  nupap = spec(is)%vnew(isb)*2*(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr) / xp**3
                                  gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                  gam_mup = 2*(eiediff*nupap*mup**2+eidefl*deflknob*nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp
                                  bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) - code_dt*(gam_mu *-1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                     +(gam_mup*-1/dmu(imu) - gam_mu*-1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu,iz,is) / (dmu(imu)/2.+mu(imu))
                                  bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) - code_dt*(gam_mu * 1/dmu(imu) * dmu(imu)/2./mu(imu) &
                                     +(gam_mup* 1/dmu(imu) - gam_mu*1/dmu(imu)) * mu(imu)/(dmu(imu)/2.)) / mw(iv,imu+1,iz,is) / (dmu(imu)/2.+mu(imu))
                                  ! mixed derivative:
                                  if (density_conservation) then
                                      ! to ensure density conservation, we change discretisation of mixed derivative term at imu=1
                                      ! from explicit routine: Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah &
                                      !     + vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz,is)*Dvpah_p) / (2.*dmu(imu))
                                      aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)) / (2*dmu(imu)) / (2*dvpa)
                                      aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv-1,imu+1,iz,is)) / (2*dmu(imu)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)) / (2*dmu(imu)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv+1,imu+1,iz,is)) / (2*dmu(imu)) / (2*dvpa)
                                  else
                                      aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                        / (mu(imu)+dmu(imu)) / (2*dvpa)
                                      aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv-1,imu+1,iz,is)* mu(imu)/dmu(imu)) &
                                        / (mu(imu)+dmu(imu)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) &
                                        / (mu(imu)+dmu(imu)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv+1,imu+1,iz,is)* mu(imu)/dmu(imu)) &
                                        / (mu(imu)+dmu(imu)) / (2*dvpa)
                                  end if
                              end if

                         else if (imu == nmu) then
                              ! one-sided difference for mu-derivative at imu=nmu:
                              ! to be consistent with treatment of mixed mu operator; assume that nux(imu)=0.
                              if (density_conservation) then
                                  aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - 1.0*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                    /mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                  aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                    + 1.0*vfac*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                    /mw(iv-1,imu  ,iz,is) / (2*dvpa) / dmu(nmu-1)
                                  cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + 1.0*vfac*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                    /mw(iv+1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                  cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    - 1.0*vfac*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                    /mw(iv+1,imu  ,iz,is) / (2*dvpa) / dmu(nmu-1)
                              else
                                  aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - 1.0*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                  aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                    + 1.0*vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb) / (2*dvpa) / dmu(nmu-1)
                                  cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + 1.0*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                                  cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    - 1.0*vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb) / (2*dvpa) / dmu(nmu-1)
                              end if

                              bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv &
                                + 0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz,is)

                              ! mu operator
                              if (mu_operator) then
                                  mum = 0.5*(mu(imu)+mu(imu-1))
                                  mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                  xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                  if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                      nuDm = eidefl*deflknob*spec(is)%vnew(isb) / xm**3
                                  else
                                      nuDm = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))&
                                        -(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                                  end if
                                  nupam = spec(is)%vnew(isb)*2*(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr) / xm**3
                                  gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                  gam_mum = 2*(eiediff*nupam*mum**2+eidefl*deflknob*nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm

                                  if (density_conservation) then
                                      ! to ensure density conservation assume that argument of outer derivative vanishes at nmu+1/2, ie
                                      ! d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2)
                                      ! where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
                                      bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                        - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                                        + (-gam_mu/dmu(imu-1))*dmu(imu-1)/dmu(imu-1)) / mw(iv,imu,iz,is) / (dmu(imu-1)/2.+dmu(imu-1)/2.)
                                      bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                                        -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/dmu(imu-1)) / mw(iv,imu-1,iz,is) / (dmu(imu-1)/2.+dmu(imu-1)/2.)
                                  else
                                      bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                        - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                        + (-gam_mu/dmu(imu-1))*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu,iz,is) / (dmu(imu-1)/2.+dmu(imu-1))
                                      bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * -1/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                        -gam_mu * -1/dmu(imu-1)*dmu(imu-1)/2./dmu(imu-1)) / mw(iv,imu-1,iz,is) / (dmu(imu-1)/2.+dmu(imu-1))
                                  end if
                                  ! no distinction here between density_conserving and default
                                  if (density_conservation) then
                                      aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)) / (2*dmu(imu-1)) / (2*dvpa)
                                      aa_blcs(iv,imu,imu  ,ikxkyz,isb) = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)*1/mw(iv-1,imu  ,iz,is)) / (2*dmu(imu-1)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)) / (2*dmu(imu-1)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu  ,ikxkyz,isb) = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)*1/mw(iv+1,imu  ,iz,is)) / (2*dmu(imu-1)) / (2*dvpa)
                                  else
                                      aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) &
                                        / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                      aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) &
                                        / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                        - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) &
                                        / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                      cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                        + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) &
                                        / (dmu(imu-1)+dmu(imu-1)) / (2*dvpa)
                                  end if

                              end if
                         else ! interior mu points
                             if (density_conservation) then
                                 aa_blcs(iv,imu,imu,ikxkyz,isb) = aa_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                    +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu,ikxkyz,isb) = cc_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 ! vpa operator, mixed (interior treatment):
                                 aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    + vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                             else
                                 aa_blcs(iv,imu,imu,ikxkyz,isb) = aa_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupamv*vpam**2 + eidefl*deflknob*2*nuDmv*bmag(ia,iz)*mu(imu))*mwmv / dvpa**2 / mw(iv-1,imu,iz,is) &
                                    +vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu,ikxkyz,isb) = cc_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*0.5*(eiediff*nupapv*vpap**2 + eidefl*deflknob*2*nuDpv*bmag(ia,iz)*mu(imu))*mwpv / dvpa**2 / mw(iv+1,imu,iz,is) &
                                    -vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 ! vpa operator, mixed (interior treatment):
                                 aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    + vfac*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - vfac*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                             end if

                             bb_blcs(iv,imu,imu,ikxkyz,isb) = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                + code_dt*(0.5*(eiediff*nupapv*vpap**2 + 2*eidefl*deflknob*nuDpv*bmag(ia,iz)*mu(imu))*mwpv &
                                + 0.5*(eiediff*nupamv*vpam**2 + 2*eidefl*deflknob*nuDmv*bmag(ia,iz)*mu(imu))*mwmv) / dvpa**2 / mw(iv,imu,iz,is)

                             ! mu operator
                             if (mu_operator) then
                                 ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                                 mup = 0.5*(mu(imu)+mu(imu+1))
                                 mum = 0.5*(mu(imu)+mu(imu-1))
                                 mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
                                 mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
                                 xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
                                 xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
                                 if ((is==2).and.(isb==1).and.(eimassr_approx)) then
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb) / xp**3
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb) / xm**3
                                 else
                                     nuDp = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xp/sqrt(massr))&
                                        -(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                     nuDm = eidefl*deflknob*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))&
                                        -(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                                 end if
                                 nupap = spec(is)%vnew(isb)*2*(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr) / xp**3
                                 nupam = spec(is)%vnew(isb)*2*(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr) / xm**3
                                 gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*deflknob*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                                 gam_mum = 2*(eiediff*nupam*mum**2+eidefl*deflknob*nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm
                                 gam_mup = 2*(eiediff*nupap*mup**2+eidefl*deflknob*nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp
                                 ! mu_operator (interior treatment):
                                 bb_blcs(iv,imu,imu,ikxkyz,isb)   = bb_blcs(iv,imu,imu,ikxkyz,isb) &
                                    - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                    +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu)) / mw(iv,imu,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                                 bb_blcs(iv,imu,imu-1,ikxkyz,isb) = bb_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - code_dt*((gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) - gam_mum* -1/dmu(imu-1)) * dmu(imu)/dmu(imu-1) &
                                    -gam_mu * -1*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) * dmu(imu-1)/dmu(imu)) / mw(iv,imu-1,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                                 bb_blcs(iv,imu,imu+1,ikxkyz,isb) = bb_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - code_dt*( gam_mu * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) * dmu(imu)/dmu(imu-1) &
                                    + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu)) / mw(iv,imu+1,iz,is) * 2/(dmu(imu-1)+dmu(imu))
                                 ! mu operator, mixed (interior treatment):
                                 aa_blcs(iv,imu,imu,ikxkyz,isb)   = aa_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)*1/mw(iv-1,imu  ,iz,is)&
                                    *(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 aa_blcs(iv,imu,imu-1,ikxkyz,isb) = aa_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)&
                                    * dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 aa_blcs(iv,imu,imu+1,ikxkyz,isb) = aa_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv-1,imu+1,iz,is)&
                                    * dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu,ikxkyz,isb)   = cc_blcs(iv,imu,imu  ,ikxkyz,isb) &
                                    - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)*1/mw(iv+1,imu  ,iz,is)&
                                    *(dmu(imu)/dmu(imu-1)-dmu(imu-1)/dmu(imu))) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu-1,ikxkyz,isb) = cc_blcs(iv,imu,imu-1,ikxkyz,isb) &
                                    + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)&
                                    * dmu(imu)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                                 cc_blcs(iv,imu,imu+1,ikxkyz,isb) = cc_blcs(iv,imu,imu+1,ikxkyz,isb) &
                                    - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv+1,imu+1,iz,is)&
                                    * dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                            end if
                         end if
                     end if

                end do

            end do
        end do
    end do

    if (testpart .eqv. .false.) then
        aa_blcs = 0.
        cc_blcs = 0.
        bb_blcs = 0.
        do isb = 1, nspec
            do imu = 1, nmu
                do iv = 1, nvpa
                    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                        bb_blcs(iv,imu,imu,ikxkyz,isb) = 1. ! AVB: in beta
                    end do
               end do
           end do
       end do
    end if

    ! construct full block matrix
    if ((exact_conservation_tp).or.(density_conservation_tp)) then
        ! this is memory intensive, operating with blockmatrix is slow
        ! currently used for exact conservation scheme on non-uniform grids
        ! AVB: to do - replace this with band matrix operations

        if (.not.allocated(blockmatrix)) allocate (blockmatrix(nvpa*nmu,nvpa*nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc,nspec))
        if (.not.allocated(blockmatrix_sum)) allocate (blockmatrix_sum(nvpa*nmu,nvpa*nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

        blockmatrix = 0.
        do isb = 1, nspec
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                iz  = iz_idx(kxkyz_lo,ikxkyz)
                is  = is_idx(kxkyz_lo,ikxkyz)
                do iv = 1, nvpa
                    ! diagonal blocks:
                    blockmatrix(nmu*(iv-1)+1:nmu*iv, nmu*(iv-1)+1:nmu*iv, ikxkyz, isb) = bb_blcs(iv,:,:,ikxkyz,isb)
                    if (iv < nvpa) then
                        ! subdiagonal blocks:
                        blockmatrix(nmu*iv+1:nmu*(iv+1), nmu*(iv-1)+1:nmu*iv, ikxkyz, isb) = aa_blcs(iv+1,:,:,ikxkyz,isb)
                        ! superdiagonal blocks:
                        blockmatrix(nmu*(iv-1)+1:nmu*iv, nmu*iv+1:nmu*(iv+1), ikxkyz, isb) = cc_blcs(iv,:,:,ikxkyz,isb)
                    end if
                end do
            end do
        end do

    end if

    ! switch to band-storage for LAPACK banded solver routines:
    ! and sum the interspecies and intraspecies operators for each species
    ! a_ij is stored in aband(ku+1+i-j,j) for $\max(1,j-ku) \leq i \leq \min(m,j+kl)$
    cdiffmat_band = 0.

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        is  = is_idx(kxkyz_lo,ikxkyz)

        ! loop through aa_blcs, bb_blcs, cc_blcs
        ! find corresponding index of blockmatrix at each location within blcs
        ! calculate index of cdiffmat_band

        do iv = 1, nvpa

            ! bb_blcs
            do imu = 1, nmu
                bm_rowind = (iv-1)*nmu + imu
                do imu2 = 1, nmu
                    bm_colind = (iv-1)*nmu + imu2
                    if ((max(1, bm_colind-(nmu+1)) .le. bm_rowind) .and. (bm_rowind .le. min(nvpa*nmu, bm_colind+(nmu+1)))) then
                        cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) = bb_blcs(iv,imu,imu2,ikxkyz,is)
                    end if
                end do
            end do
            ! cc_blcs
            if (iv<nvpa) then ! aa_blcs and cc_blcs contain only (nvpa-1) blocks, since they are off-diagonal
                do imu = 1, nmu
                    bm_rowind = (iv-1)*nmu + imu
                    do imu2 = 1, nmu
                        bm_colind = nmu + (iv-1)*nmu + imu2 ! nvpa*nmu - nmu + imu2
                        if ((max(1, bm_colind-(nmu+1)) .le. bm_rowind) .and. (bm_rowind .le. min(nvpa*nmu, bm_colind+(nmu+1)))) then
                            cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) = cc_blcs(iv,imu,imu2,ikxkyz,is)
                        end if
                    end do
                end do
                ! aa_blcs
                do imu = 1, nmu
                    bm_rowind = nmu + (iv-1)*nmu + imu
                    do imu2 = 1, nmu
                        bm_colind = (iv-1)*nmu + imu2
                        if ((max(1, bm_colind-(nmu+1)) .le. bm_rowind) .and. (bm_rowind .le. min(nvpa*nmu, bm_colind+(nmu+1)))) then
                            cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) = aa_blcs(1+iv,imu,imu2,ikxkyz,is)
                        end if
                    end do
                end do
            end if

            ! inter-species test particle contributions:

            do isb = 1, nspec
                if (isb==is) cycle
                ! bb_blcs
                do imu = 1, nmu
                    bm_rowind = (iv-1)*nmu + imu
                    do imu2 = 1, nmu
                        bm_colind = (iv-1)*nmu + imu2
                        if ((max(1, bm_colind-(nmu+1)) .le. bm_rowind) .and. (bm_rowind .le. min(nvpa*nmu, bm_colind+(nmu+1)))) then
                            cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) = &
                                cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) + bb_blcs(iv,imu,imu2,ikxkyz,isb)
                        end if
                    end do
                end do
                ! cc_blcs
                if (iv<nvpa) then
                    do imu = 1, nmu
                        bm_rowind = (iv-1)*nmu + imu
                        do imu2 = 1, nmu
                            bm_colind = (iv-1)*nmu + nmu + imu2
                            if ((max(1, bm_colind-(nmu+1)) .le. bm_rowind) .and. (bm_rowind .le. min(nvpa*nmu, bm_colind+(nmu+1)))) then
                                cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) = &
                                    cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) + cc_blcs(iv,imu,imu2,ikxkyz,isb)
                            end if
                        end do
                    end do
                    ! aa_blcs
                    do imu = 1, nmu
                        bm_rowind = (iv-1)*nmu + nmu + imu
                        do imu2 = 1, nmu
                            bm_colind = (iv-1)*nmu + imu2
                            if ((max(1, bm_colind-(nmu+1)) .le. bm_rowind) .and. (bm_rowind .le. min(nvpa*nmu, bm_colind+(nmu+1)))) then
                                cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) = &
                                    cdiffmat_band(nmu+1 + nmu+1 + 1+bm_rowind-bm_colind,bm_colind,iky,ikx,iz,is) + aa_blcs(1+iv,imu,imu2,ikxkyz,isb)
                            end if
                        end do
                    end do
                end if
            end do

        end do

    end do

    ! add the gyro-diffusive term to cdiffmat
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        is  = is_idx(kxkyz_lo,ikxkyz)
        do iv = 1, nvpa
            do imu = 1, nmu
                ! diagonal indices in blockmatrix
                ivv = nmu*(iv-1)+imu
                imm = ivv
                if ((max(1, ivv-(nmu+1)) .le. imm) .and. (imm .le. min(nvpa*nmu, ivv+(nmu+1)))) then
                    ! intra-species:
                    cdiffmat_band(nmu+1 + nmu+1 + 1+imm-ivv,ivv,iky,ikx,iz,is) = cdiffmat_band(nmu+1 + nmu+1 + 1+imm-ivv,ivv,iky,ikx,iz,is) &
                            + code_dt*cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
                            + deflknob*nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))
                    ! inter-species:
                    do isb = 1, nspec
                        if (isb==is) cycle

                        if ((is==2).and.(isb==1)) then
                            eiediff = eiediffknob
                            eidefl = eideflknob
                        else
                            eiediff = 1
                            eidefl = 1
                        end if

                        cdiffmat_band(nmu+1 + nmu+1 + 1+imm-ivv,ivv,iky,ikx,iz,is) = cdiffmat_band(nmu+1 + nmu+1 + 1+imm-ivv,ivv,iky,ikx,iz,is) &
                            + code_dt*cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(eiediff*nupa(iv,imu,iz,is,isb)*bmag(ia,iz)*mu(imu) &
                            + eidefl*deflknob*nuD(iv,imu,iz,is,isb)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))
                    end do
                end if
            end do
        end do
    end do

    ! add 1 to the diagonal, since the matrix operator is 1 - C^{ab}[h_a]
    do is = 1, nspec
        do iv = 1, nmu*nvpa
            do imu = 1, nmu*nvpa
                if ((max(1, iv-(nmu+1)) .le. imu) .and. (imu .le. min(nvpa*nmu, iv+(nmu+1)))) then
                    if (iv==imu) then
                        cdiffmat_band(nmu+1 + nmu+1 + 1+imu-iv,iv,:,:,:,is) = cdiffmat_band(nmu+1 + nmu+1 + 1+imu-iv,iv,:,:,:,is) + 1.
                    end if
                end if
            end do
        end do
    end do

    ! to write matrix in band-storage, for debugging
    !do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
    !    iky = iky_idx(kxkyz_lo,ikxkyz)
    !    ikx = ikx_idx(kxkyz_lo,ikxkyz)
    !    iz = iz_idx(kxkyz_lo,ikxkyz)
    !    is  = is_idx(kxkyz_lo,ikxkyz)
    !    if (iz/=0) cycle
    !    if (iky/=1) cycle
    !    if (is==2) then
    !        call open_output_file (tmpunit,'.cdiffmatband')
    !        do iv = 1, nvpa*nmu
    !          write (tmpunit,'(9es15.4e3)') cdiffmat_band(iv,:,iky,ikx,iz,is)
    !        end do
    !        write (tmpunit,*)
    !        call close_output_file (tmpunit)
    !    end if
    !end do

    ! AVB: LU factorise cdiffmat, using LAPACK's zgbtrf routine for banded matrices
    nc = nvpa*nmu
    nb = nmu+1
    lldab = 3*(nmu+1)+1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        is  = is_idx(kxkyz_lo,ikxkyz)
        call zgbtrf(nc, nc, nb, nb, cdiffmat_band(:,:,iky,ikx,iz,is), lldab, ipiv(:,iky,ikx,iz,is), info)
    end do

end subroutine init_fp_diffmatrix

     elemental function associated_laguerre (n, alpha, x)

      integer, intent (in) :: n
      real, intent (in) :: x
      real, intent (in) :: alpha
      integer :: k
      real :: associated_laguerre, p, p1, p2

      p1 = dble(1.0)
      p2 = dble(1.0) + alpha - x

      if (n==0) then
         associated_laguerre = p1
         return
      else if (n==1) then
         associated_laguerre = p2
         return
      end if

      do k=2, n
         p = ((dble(2.0)*k-dble(1.0) + alpha - x) * p2 - (k-dble(1.0) + alpha) * p1) / k
         p1 = p2
         p2 = p
      end do

      associated_laguerre = p

    end function associated_laguerre

    elemental function associated_legendre (l, m, x)

      integer, intent (in) :: l, m
      double precision, intent (in) :: x
      integer :: k
      double precision :: associated_legendre, p, p1, p2, fac
      double precision :: pi

      pi = 3.14159265359

      ! to start the recursion, use that P_l^m = 0 for l < abs(m)
      ! and P_l^l = (-1)^l*(2l-1)!!(1-x^2)^(l/2)
      ! where (2l-1)!! = 2**l * Gamma(l+0.5) / sqrt(pi)
      p1 = 0.
      p2 = (-1)**abs(m) * 2**abs(m) * gamma(abs(m)+0.5) / sqrt(pi) * (1.-x**2)**(abs(m)/2.)

      if (abs(m)>l) then
          associated_legendre = 0.
          return
      end if

      if (l==0) then
         associated_legendre = 1.
         return
      end if

      if (l==m) then
         associated_legendre = p2
         return
      end if

      if (l==-m) then
         fac = (-1)**abs(m)*gamma(l-abs(m)+1.)/gamma(l+abs(m)+1.)
         associated_legendre = p2*fac
         return
     end if

      do k = abs(m)+1, l
         p = ((dble(2.0)*k-dble(1.0)) * x * p2 - (k-dble(1.0) + abs(m)) * p1) / (k-abs(m))
         p1 = p2
         p2 = p
      end do

      if (m < 0) then
          fac = (-1)**abs(m)*gamma(l-abs(m)+1.)/gamma(l+abs(m)+1.)
          p = p*fac
      end if

      associated_legendre = p

  end function associated_legendre

  subroutine init_legendre

      use vpamu_grids, only: mu, nmu, vpa, nvpa
      use zgrid, only: nzgrid
      use stella_geometry, only: bmag
      use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx
      use file_utils, only: open_output_file, close_output_file

      implicit none

      integer :: is, iv, imu, iz, ia, mm, ll, tmpunit
      double precision :: xi

      allocate (legendre_vpamu(0:lmax,-lmax:lmax,nvpa,nmu,-nzgrid:nzgrid))

      ia = 1

      ! note lmin = 0, lmax = nsph-1
      ! mmin = -lmax, mmax = lmax

      legendre_vpamu = 0.

      ia = 1

      do iv = 1, nvpa
            do imu = 1, nmu
                do iz = -nzgrid, nzgrid
                    do ll = 0, lmax
                        do mm = -lmax, lmax
                            xi = dble(vpa(iv)/sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu)))
                            legendre_vpamu(ll,mm,iv,imu,iz) = associated_legendre(ll,mm,xi)
                        end do
                    end do
                end do
            end do
      end do

  end subroutine init_legendre

    subroutine init_bessel_fn
        use zgrid, only: nzgrid
        use vpamu_grids, only: nvpa, nmu, vpa, mu, vperp2
        use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, it_idx, is_idx
        use gyro_averages, only: aj0v, aj1v
        use species, only: spec, nspec
        use stella_geometry, only: bmag
        use kt_grids, only: naky, nakx
        use dist_fn_arrays, only: kperp2
        use file_utils, only: open_output_file, close_output_file

        implicit none

        integer :: ikxkyz, iky, ikx, iz, is, ia, mm, tmpunit, imu
        real :: arg, aj1fac, aj1exp, aj0exp, iv

        allocate (jm(nmu,0:lmax,naky,nakx,-nzgrid:nzgrid,nspec))
        allocate (jm0(nmu,naky,nakx,-nzgrid:nzgrid,nspec))

        jm = 0
        ia = 1

        aj1fac = 1.0
        aj0exp = 1.0
        aj1exp = 1.0

        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo,ikxkyz)
            ikx = ikx_idx(kxkyz_lo,ikxkyz)
            iz  = iz_idx(kxkyz_lo,ikxkyz)
            is  = is_idx(kxkyz_lo,ikxkyz)
            jm0(:,iky,ikx,iz,is) = aj0v(:,ikxkyz)**aj0exp
            do mm = 0, lmax
                if (mm==0) then
                    jm(:,0,iky,ikx,iz,is) = aj0v(:,ikxkyz)**aj0exp
                else if (mm==1) then
                    !jm(:,1,iky,ikx,iz,is) = aj1fac*aj1v(:,ikxkyz)*spec(is)%smz_psi0*sqrt(vperp2(ia,iz,:)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                    !mu*spec(is)%smz*sqrt(2*bmag(ia,iz)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                    do imu = 1, nmu
                        arg = spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                        jm(imu,mm,iky,ikx,iz,is) = bessel_jn(mm,arg)**aj1exp*aj1fac
                    end do
                else
                    do imu = 1, nmu
                        arg = spec(is)%smz_psi0*sqrt(vperp2(ia,iz,imu)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                        !mu(imu)*spec(is)%smz*sqrt(2*bmag(ia,iz)*kperp2(iky,ikx,ia,iz))/bmag(ia,iz)
                        jm(imu,mm,iky,ikx,iz,is) = bessel_jn(mm,arg)
                    end do
                end if
            end do
        end do

        if (cfac2==0.) then
            ! disable gyro-diffusive effects in the field particle operator
            jm(:,0,:,:,:,:) = 1.
            if (lmax > 0) then
                do mm = 1, lmax
                    jm(:,mm,:,:,:,:) = 0.
                end do
            end if
        end if

    end subroutine init_bessel_fn

    subroutine init_vgrid

        ! v grid used for writing various coll freqs to file for debugging

        use vpamu_grids, only: dmu, mu, nmu, vpa, dvpa, nvpa, vperp_max, vpa_max
        use stella_geometry, only: bmag

        integer :: ikxkyz, iky, ikx, iz, is
        integer :: ia, iv, idx, nv_seg, iseg, ll
        real :: del, delv, vmax, vmin

        allocate (wgts_v(nvel_local)); wgts_v = 0.0
        allocate (vel(nvel_local)); vel = 0.0
        ia =  1

        ! calculation of Delta_j[x^l L_j(x^2) exp(-x^2)] requires integrals over v, set these up below:
        ! velocity grid, equally spaced, from 0 to vmax:
        vmax = sqrt(maxval(vpa)**2 + 2*maxval(bmag(ia,:))*maxval(mu)) !sqrt(vperp_max**2 + vpa_max**2)
        !vmin = 1e-9
        delv = vmax/nvel_local ! (vmax-vmin)/(nvel_local-1)
        vmin = sqrt(minval(abs(vpa))**2 + 2*minval(bmag)*minval(mu)) !delv !1e-9
        do iv = 1, nvel_local
            vel(iv) = vmin+(iv-1)*delv
        end do

        ! Trapezoidal rule:
        wgts_v = 0.5*delv

        ! Composite Simpson's at interior nodes, average of Simpson's 3/8 and Composite at boundaries:
        ! Lower boundary, Simpson's 3/8:
        !del = 0.375*delv
        !wgts_v(1) = del
        !wgts_v(2:3) = 3.*del
        !wgts_v(4) = del
        ! Interior points, Composite:

        !nv_seg = (nvel_local-4)/2
        !del = delv/3.
        !do iseg = 1, nv_seg
        !   idx = 2*(iseg-1) + 4 ! for iseg = 1, idx = 4; for iseg = nv_seg, idx = nv-2.
        !   wgts_v(idx) = wgts_v(idx) + del
        !   wgts_v(idx+1) = wgts_v(idx+1) + 4.*del
        !   wgts_v(idx+2) = wgts_v(idx+2) + del
        !end do

        ! Upper boundary, Simpson's 3/8:
        !del = 0.375*delv
        !wgts_v(nvel_local-3) = wgts_v(nvel_local-3) + del
        !wgts_v(nvel_local-2:nvel_local-1) = wgts_v(nvel_local-2:nvel_local-1) + 3.*del
        !wgts_v(nvel_local) = wgts_v(nvel_local) + del
        ! Interior points, Composite:
        !nv_seg = (nvel_local-4)/2
        !del = delv/3.
        !do iseg = 1, nv_seg
        !   idx = 2*(iseg-1)+1 ! for iseg = 1, idx = 1; for iseg = nv_seg, idx = nv-5.
        !   wgts_v(idx) = wgts_v(idx) + del
        !   wgts_v(idx+1) = wgts_v(idx+1) + 4.*del
        !   wgts_v(idx+2) = wgts_v(idx+2) + del
        !end do

        ! AVB: the points idx = 4 and idx = nv-3 are counted three times
        ! divide by 2 to account for double-counting
        !wgts_v = 0.5*wgts_v

    end subroutine init_vgrid

    recursive subroutine gamlow (a, x, gl)

        ! recursive calculation of lower incomplete gamma function for half-integer a > 0

        use constants, only: pi
        use spfunc, only: erf => erf_ext

        implicit none
        real, intent (in) :: a
        real, intent (in) :: x
        real, intent (out) :: gl
        real :: glm1

        if (a==0.5) then
            gl = sqrt(pi)*erf(sqrt(x))
        else
            call gamlow(a-1., x, glm1)
            gl = (a-1.)*glm1 - x**(a-1.)*exp(-x)
        end if

    end subroutine gamlow

    recursive subroutine gamup (a, x, gu)

        ! recursive calculation of the upper incomplete gamma function for integer a > 0

        use constants, only: pi
        use spfunc, only: erf => erf_ext

        implicit none
        real, intent (in) :: a
        real, intent (in) :: x
        real, intent (out) :: gu
        real :: gum1

        if (a==1.0) then
            gu = exp(-x)
        else
            call gamup(a-1.,x,gum1)
            gu = (a-1.)*gum1 + x**(a-1.)*exp(-x)
        end if

    end subroutine gamup

    subroutine calc_delta0 (xa, jj, ll, isa, isb, delt0)

        ! calculate Delta0^{j,l,ab}(xa) (on xa grid)
        ! j and l denote the degree and index of the associated laguerre polynomial

        use species, only: nspec, spec
        use constants, only: pi

        implicit none

        integer, intent (in) :: jj, ll, isa, isb
        real, intent (in) :: xa
        real, intent (out) :: delt0

        real :: massr, ckjl, xb, gaml1, gaml2, gamu1, gamu2
        integer :: kk

        massr = spec(isa)%mass / spec(isb)%mass
        xb    = xa / sqrt(massr)

        delt0 = 0
        do kk = 0, jj
            ckjl = (-1)**kk*gamma(jj+ll+0.5+1)/(gamma(jj-kk+1.)*gamma(ll+kk+0.5+1)*gamma(kk+1.))
            call gamlow(1.5+ll+kk,xb**2,gaml1)
            call gamlow(2.5+ll+kk,xb**2,gaml2)
            call gamup(1.+kk,xb**2,gamu1)
            call gamup(2.+kk,xb**2,gamu2)
            delt0 = delt0 + ckjl*((2*ll+1.)*xb**(ll+2.*kk)*exp(-xb**2) &
                - xb*(1.-massr)*( -(ll+1.)/xb**(ll+2.)*gaml1 + ll*xb**(ll-1.)*gamu1 ) &
                - ( 1./xb**(ll+1.)*gaml1 + xb**ll*gamu1 ) &
                + massr*xb**2*( (ll+1.)*(ll+2.)/(2*ll+3.)*(xb**(-ll-3.)*gaml2 + xb**ll*gamu1) &
                                -ll*(ll-1.)/(2*ll-1.)*(xb**(-ll-1.)*gaml1 + xb**(ll-2.)*gamu2) ))
        end do
        delt0 = delt0*4*pi/(pi**1.5)*exp(-xa**2)*massr/(2*ll+1.)

    end subroutine calc_delta0

    recursive subroutine calc_deltaj_vmu (jj, nn, ll, isa, isb, deltj)

        ! calculate Delta_j^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)](xa) (on x_a grid)
        ! these are normalised, and calculated without the collision frequency
        ! in contrast to Hirshman & Sigmar 1976

        use species, only: nspec, spec
        use vpamu_grids, only: dmu, mu, nmu, vpa, dvpa, nvpa
        use zgrid, only: nzgrid
        use stella_geometry, only: bmag

        implicit none

        integer, intent (in) :: jj, nn, ll, isa, isb
        real, dimension (nvpa,nmu,-nzgrid:nzgrid), intent (out) :: deltj
        real, dimension (nvpa,nmu,-nzgrid:nzgrid) :: deltajm1_n, deltajm1_j
        integer :: iv, imu, iz, ia
        real :: v
        real, dimension (-nzgrid:nzgrid) :: psijm1_n

        ia = 1
        if (jj==0) then
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid
                        v = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))
                        call calc_delta0 (v, nn, ll, isa, isb, deltj(iv,imu,iz))
                    end do
                end do
            end do
        else
            ! get Delta_[j-1]^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)](xa)
            call calc_deltaj_vmu (jj-1, nn, ll, isa, isb, deltajm1_n)
            ! get Delta_[j-1]^{ab,l}[x_b^l L_[j-1]^{l+0.5}(x_b^2) exp(-x_b^2)](xa)
            call calc_deltaj_vmu (jj-1, jj-1, ll, isa, isb, deltajm1_j)
            ! get psi_[j-1]^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)](xa)
            call calc_psi_vmu (jj-1, nn, ll, isa, isb, psijm1_n)
            deltj = deltajm1_n - spread(spread(psijm1_n,1,nvpa),2,nmu)*deltajm1_j
        end if

    end subroutine calc_deltaj_vmu

    subroutine vLj_vmu (jj, ll, vLj)

        use species, only: nspec, spec
        use vpamu_grids, only: dmu, mu, nmu, vpa, dvpa, nvpa
        use zgrid, only: nzgrid
        use stella_geometry, only: bmag

        implicit none

        integer, intent (in) :: jj, ll
        real, dimension (nvpa,nmu,-nzgrid:nzgrid), intent (out) :: vLj
        integer :: iv, imu, iz, ia
        real :: v

        ia = 1
        do iv = 1, nvpa
            do imu = 1, nmu
                do iz = -nzgrid, nzgrid
                    v = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))
                    vLj(iv,imu,iz) = v**ll * associated_laguerre(jj,ll+1./2.,v**2)
                end do
            end do
        end do

    end subroutine vLj_vmu

    recursive subroutine calc_psi_vmu (jj, nn, ll, isa, isb, psij)

        ! calculate psi_j^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)]

        ! have defined deltaj without collision frequency
        ! and normalised everthing to species thermal speeds

        use species, only: nspec, spec
        use vpamu_grids, only: dmu, mu, nmu, vpa, dvpa, nvpa, integrate_vmu
        use zgrid, only: nzgrid

        implicit none

        integer, intent (in) :: jj, nn, ll, isa, isb
        real, dimension (nvpa,nmu,-nzgrid:nzgrid) :: deltj_j, vLj, vLn !deltj_n
        real, dimension (-nzgrid:nzgrid), intent (out) :: psij
        integer :: iv, imu, iz, ia
        real :: v, psijm1_n
        real, dimension (-nzgrid:nzgrid) :: num, den

        if (jj==0) then
            if ((ll==0).and.(nn==0)) then ! never used
                num = 1
                den = 1
            else
                ! get delta_j^{ba}(x_a^l L_j^{l+0.5}(x_a^2) exp(-x_a^2)), on x_b grid
                call calc_deltaj_vmu (jj, jj, ll, isb, isa, deltj_j)
                ! multiply by x_b^l L_n(x_b^2) and integrate
                call vLj_vmu(nn, ll, vLn)
                do iz = -nzgrid, nzgrid
                    call integrate_vmu(vLn(:,:,iz)*deltj_j(:,:,iz), iz, num(iz)) ! numerator in psijl
                end do

                ! if mb / ma < 1 use self-adjointness to avoid resolution problems
                if (spec(isb)%mass / spec(isa)%mass < 1.) then
                    call calc_deltaj_vmu (jj, jj, ll, isa, isb, deltj_j)
                    ! need to account for mass ratio normalisation
                    ! and collision frequency
                    ! \int x_a^l L_j^{l+1/2}(x_a^2) \Delta_j'^{ab}[\tilde{f}_b/F0b * F0b] x_a^2 dx_a
                    !   = ma^0.5/mb^0.5 * ma^3/mb^3 \int f_b/F0b \Delta_j^{ba}'[x_a^l L_j^{l+1/2}(x_a^2) \tilde{F0a}] x_b^2 dx_b
                    call vLj_vmu(jj, ll, vLj)
                    do iz = -nzgrid, nzgrid
                        call integrate_vmu(spec(isb)%mass**3.5/spec(isa)%mass**3.5*vLj(:,:,iz)*deltj_j(:,:,iz), iz, den(iz)) ! denominator in psijl
                    end do
                else
                    call vLj_vmu(jj, ll, vLj)
                    do iz = -nzgrid, nzgrid
                        call integrate_vmu(vLj(:,:,iz)*deltj_j(:,:,iz), iz, den(iz)) ! denominator in psijl
                    end do
                end if

            end if
        else
            ! get delta_j^{ba}(x_a^l L_j^{l+0.5}(x_a^2) exp(-x_a^2)), on x_b grid
            call calc_deltaj_vmu (jj, jj, ll, isb, isa, deltj_j)
            ! multiply by x_b^l L_n(x_b^2) and integrate
            call vLj_vmu(nn, ll, vLn)
            do iz = -nzgrid, nzgrid
                call integrate_vmu(vLn(:,:,iz)*deltj_j(:,:,iz), iz, num(iz)) ! numerator in psijl
            end do

            ! if mb / ma < 1 use self-adjointness to avoid resolution problems
            if (spec(isb)%mass / spec(isa)%mass < 1.) then
                call calc_deltaj_vmu (jj, jj, ll, isa, isb, deltj_j)
                ! need to account for mass ratio normalisation
                ! and collision frequency
                ! \int x_a^l L_j^{l+1/2}(x_a^2) \Delta_j'^{ab}[\tilde{f}_b/F0b * F0b] x_a^2 dx_a
                !   = ma^0.5/mb^0.5 * ma^3/mb^3 \int f_b/F0b \Delta_j^{ba}'[x_a^l L_j^{l+1/2}(x_a^2) \tilde{F0a}] x_b^2 dx_b
                call vLj_vmu(jj, ll, vLj)
                do iz = -nzgrid, nzgrid
                    call integrate_vmu(spec(isb)%mass**3.5/spec(isa)%mass**3.5*vLj(:,:,iz)*deltj_j(:,:,iz), iz, den(iz)) ! denominator in psijl
                end do
            else
                call vLj_vmu(jj, ll, vLj)
                do iz = -nzgrid, nzgrid
                    call integrate_vmu(vLj(:,:,iz)*deltj_j(:,:,iz), iz, den(iz)) ! denominator in psijl
                end do
            end if

        end if

        psij = num/den

    end subroutine calc_psi_vmu

    subroutine init_deltaj_vmu

        use species, only: nspec, spec
        use zgrid, only: nzgrid
        use stella_geometry, only: bmag
        use vpamu_grids, only: dmu, mu, nmu, vpa, dvpa, nvpa, vperp_max, vpa_max, integrate_vmu, set_vpa_weights, wgts_vpa
        use file_utils, only: open_output_file, close_output_file
        use constants, only: pi
        use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx, it_idx
        use stella_time, only: code_dt
        use kt_grids, only: naky

        implicit none

        integer :: ll, tmpunit, iv, jj, ia, imu, iz, is, isa, isb, ix, ikx, iky, ikxkyz
        real :: v, orth_test1, orth_test2, orth_test3, orth_test4, orth_test5
        real, dimension (-nzgrid:nzgrid) :: deltajint, deltajint_tp, orthtestz
        real, dimension (nvpa,nmu,-nzgrid:nzgrid) :: vLj
        real, dimension (0:lmax,0:jmax,nvpa,nmu,1,-nzgrid:nzgrid) :: vlaguerre_vmu, vlaguerre_vmu_deltaj
        real, dimension (0:lmax,0:jmax,nspec,nspec,nvpa,nmu,1,-nzgrid:nzgrid) :: vlaguerre_vmu_xb
        real, dimension (nvel_local) :: delt0test1, delt0test2, delt0test3, delt0test4, delt0test1_xb, delt0test2_xb
        real, dimension (nvpa*nmu,-nzgrid:nzgrid) :: vpaF0vec, v2F0vec
        real, dimension (nvpa*nmu,nvpa*nmu) :: ident
        real, dimension (nvpa,nmu,-nzgrid:nzgrid) :: vlLag

        logical :: conservative_wgts

        ia = 1

        allocate (deltaj(0:lmax,0:jmax,nspec,nspec,nvpa,nmu,ia,-nzgrid:nzgrid))
        allocate (deltaj_tp(0:lmax,0:jmax,nspec,nspec,nvpa,nmu,ia,-nzgrid:nzgrid))
        allocate (psijnorm(0:lmax,0:jmax,nspec,nspec,-nzgrid:nzgrid))

        allocate (mwnorm(-nzgrid:nzgrid))
        allocate (modmwnorm(-nzgrid:nzgrid))

        if (density_conservation) then
            conservative_wgts = .true.
            call set_vpa_weights (conservative_wgts)
        else if (exact_conservation_tp) then
            conservative_wgts = .false.
            call set_vpa_weights (conservative_wgts)
        else
            conservative_wgts = .false.
            call set_vpa_weights (conservative_wgts)
        end if

        ! AVB: to do - option for cons_wgts when exact_conservation is true

        ! get Delta_j^{l,ab}[x_b^l L_j^{l+0.5}(x_b^2)F_{0b}](x_a)
        ! and Delta_j^{l,ba}[x_a^l L_j^{l+0.5}(x_a^2)F_{0a}](x_b)

        deltaj = 0

        ! construct an identity matrix required below
        ident = 0
        forall (iv=1:nvpa*nmu) ident(iv,iv) = 1.

        do isa = 1, nspec
            do isb = 1, nspec
                do ll = 0, lmax
                    do jj = 0, jmax

                        call calc_deltaj_vmu(jj, jj, ll, isa, isb, deltaj(ll, jj, isa, isb, :, :, ia, :))
                        call vLj_vmu(jj, ll, vlaguerre_vmu(ll,jj,:,:,ia,:))

                        if (spitzer_problem) then
                             if ((exact_conservation).and.(ll==1).and.(jj==0)&
                                    .and..not.((isa==1).and.(isb==1))&
                                    .and..not.((isa==1).and.(isb==2))&
                                ) then

                                ! to ensure conservation of momentum to machine precision

                                !call open_output_file (tmpunit,'.deltaj_l1default')
                                !do iv = 1, nvpa
                                !    write (tmpunit,'(32es15.4e3)') ( deltaj(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0 ! electron-ion
                                !end do
                                !write (tmpunit,*)
                                !call close_output_file (tmpunit)

                                ! calculate delta_{j=0,l=1}^{ab} using the differential test particle operator C^{ab} acting on m_a v_\parallel F0a
                                do iv = 1, nvpa
                                    do imu = 1, nmu
                                        do iz = -nzgrid, nzgrid
                                            vpaF0vec(nmu*(iv-1)+imu,iz) = vpa(iv)*mw(iv,imu,iz,isa)
                                        end do
                                    end do
                                end do

                                do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                                   iky = iky_idx(kxkyz_lo,ikxkyz)
                                   ikx = ikx_idx(kxkyz_lo,ikxkyz)
                                   iz = iz_idx(kxkyz_lo,ikxkyz)
                                   is = is_idx(kxkyz_lo,ikxkyz)
                                   if (is/=isa) cycle
                                   vpaF0vec(:,iz) = matmul(blockmatrix(:,:,ikxkyz,isb)/code_dt/spec(is)%vnew(isb), (spec(isa)%mass/spec(isb)%mass)**2 * vpaF0vec(:,iz))
                                end do

                                do ix = 1, nvpa
                                    deltaj_tp(ll, jj, isa, isb, ix, :, ia, :) = vpaF0vec(nmu*(ix-1)+1 : nmu*(ix-1)+nmu,:)
                                end do

                                deltaj_tp(ll, jj, isa, isb, :, :, ia, :) = deltaj_tp(ll, jj, isa, isb, :, :, ia, :) * velvpamu / spread(spread(vpa,2,nmu),3,2*nzgrid+1)

                                !call open_output_file (tmpunit,'.deltaj_l1modified')
                                !do iv = 1, nvpa
                                !    write (tmpunit,'(32es15.4e3)') ( deltaj_tp(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0
                                !end do
                                !write (tmpunit,*)
                                !call close_output_file (tmpunit)

                            end if

                        else
                            if ((exact_conservation).and.(ll==1).and.(jj==0)) then

                               if (spec(isa)%vnew(isb)==0) cycle

                               !if ((isa==1).and.(isb==2)) then
                               !   call open_output_file (tmpunit,'.deltaj_l1default_ie')
                               !   do iv = 1, nvpa
                               !       write (tmpunit,'(32es15.4e3)') ( deltaj(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu )
                               !   end do
                               !   write (tmpunit,*)
                               !   call close_output_file (tmpunit)
                               !end if

                               ! calculate delta_{j=0,l=1}^{ab} using the differential test particle operator C^{ab} acting on m_a v_\parallel F0a
                               do iv = 1, nvpa
                                   do imu = 1, nmu
                                       do iz = -nzgrid, nzgrid
                                           vpaF0vec(nmu*(iv-1)+imu,iz) = vpa(iv)*mw(iv,imu,iz,isa)
                                       end do
                                   end do
                               end do

                               do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                                  iky = iky_idx(kxkyz_lo,ikxkyz)
                                  ikx = ikx_idx(kxkyz_lo,ikxkyz)
                                  iz = iz_idx(kxkyz_lo,ikxkyz)
                                  is = is_idx(kxkyz_lo,ikxkyz)
                                  if (iky/=naky) cycle
                                  if (is/=isa) cycle
                                  vpaF0vec(:,iz) = matmul(blockmatrix(:,:,ikxkyz,isb)/code_dt/spec(is)%vnew(isb), (spec(isa)%mass/spec(isb)%mass)**2 * vpaF0vec(:,iz))
                               end do

                               do ix = 1, nvpa
                                   deltaj_tp(ll, jj, isa, isb, ix, :, ia, :) = vpaF0vec(nmu*(ix-1)+1 : nmu*(ix-1)+nmu,:)
                               end do

                               deltaj_tp(ll, jj, isa, isb, :, :, ia, :) = deltaj_tp(ll, jj, isa, isb, :, :, ia, :) * velvpamu / spread(spread(vpa,2,nmu),3,2*nzgrid+1)

                               !if ((isa==1).and.(isb==2)) then
                               !call open_output_file (tmpunit,'.deltaj_l1modified_ie')
                               !   do iv = 1, nvpa
                               !       write (tmpunit,'(32es15.4e3)') ( deltaj_tp(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0
                               !  end do
                               !  write (tmpunit,*)
                               !  call close_output_file (tmpunit)
                               !end if

                           end if
                       end if

                        if (spitzer_problem) then
                            if ((exact_conservation).and.(ll==0).and.(jj==1)&
                                    .and..not.((isa==1).and.(isb==1))&
                                    .and..not.((isa==1).and.(isb==2))&
                                ) then

                                ! energy conservation to machine precision

                                !call open_output_file (tmpunit,'.deltaj_j1default')
                                !do iv = 1, nvpa
                                !    write (tmpunit,'(32es15.4e3)') ( deltaj(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0 ! electron-ion
                                !end do
                                !write (tmpunit,*)
                                !call close_output_file (tmpunit)

                                ! calculate delta_{j=1,l=0}^{ab} using the differential test particle operator C^{ab} acting on m_a v^2 F0a
                                do iv = 1, nvpa
                                    do imu = 1, nmu
                                        do iz = -nzgrid, nzgrid
                                            v2F0vec(nmu*(iv-1)+imu,iz) = (vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))*mw(iv,imu,iz,isa)
                                        end do
                                    end do
                                end do

                                do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                                   iky = iky_idx(kxkyz_lo,ikxkyz)
                                   ikx = ikx_idx(kxkyz_lo,ikxkyz)
                                   iz = iz_idx(kxkyz_lo,ikxkyz)
                                   is = is_idx(kxkyz_lo,ikxkyz)
                                   if (is/=isa) cycle
                                   v2F0vec(:,iz) = matmul(blockmatrix(:,:,ikxkyz,isb)/code_dt/spec(is)%vnew(isb), -(spec(isa)%mass/spec(isb)%mass)**1.5 * v2F0vec(:,iz))
                                end do

                                do ix = 1, nvpa
                                    deltaj_tp(ll, jj, isa, isb, ix, :, ia, :) = v2F0vec(nmu*(ix-1)+1 : nmu*(ix-1)+nmu,:)
                                end do

                                !call open_output_file (tmpunit,'.deltaj_j1modified')
                                !do iv = 1, nvpa
                                !    write (tmpunit,'(32es15.4e3)') ( deltaj_tp(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0
                                !end do
                                !write (tmpunit,*)
                                !call close_output_file (tmpunit)

                            end if

                        else
                            if ((exact_conservation).and.(ll==0).and.(jj==1)) then

                                if (spec(isa)%vnew(isb)==0) cycle

                                !if ((isa==1).and.(isb==2)) then
                                !    call open_output_file (tmpunit,'.deltaj_j1default_ie')
                                !    do iv = 1, nvpa
                                !        write (tmpunit,*) ( deltaj(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0 ! electron-ion
                                !    end do
                                !    write (tmpunit,*)
                                !    call close_output_file (tmpunit)
                                !end if

                                ! calculate delta_{j=1,l=0}^{ab} using the differential test particle operator C^{ab} acting on m_a v_\parallel^2 F0a
                                do iv = 1, nvpa
                                    do imu = 1, nmu
                                        do iz = -nzgrid, nzgrid
                                            v2F0vec(nmu*(iv-1)+imu,iz) = (vpa(iv)**2 + 2*bmag(ia,iz)*mu(imu))*mw(iv,imu,iz,isa)
                                        end do
                                    end do
                                end do

                                do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                                   iky = iky_idx(kxkyz_lo,ikxkyz)
                                   ikx = ikx_idx(kxkyz_lo,ikxkyz)
                                   iz = iz_idx(kxkyz_lo,ikxkyz)
                                   is = is_idx(kxkyz_lo,ikxkyz)
                                   if (is/=isa) cycle
                                   if (iky/=naky) cycle
                                   v2F0vec(:,iz) = matmul(blockmatrix(:,:,ikxkyz,isb)/code_dt/spec(is)%vnew(isb), -(spec(isa)%mass/spec(isb)%mass)**1.5 * v2F0vec(:,iz))
                                end do

                                do ix = 1, nvpa
                                    deltaj_tp(ll, jj, isa, isb, ix, :, ia, :) = v2F0vec(nmu*(ix-1)+1 : nmu*(ix-1)+nmu,:)
                                end do

                                !if ((isa==1).and.(isb==2)) then
                                !    call open_output_file (tmpunit,'.deltaj_j1modified_ie')
                                !    do iv = 1, nvpa
                                !        write (tmpunit,*) ( deltaj_tp(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0
                                !    end do
                                !    write (tmpunit,*)
                                !    call close_output_file (tmpunit)
                                !end if
                            end if
                        end if

                        !conservative_wgts = .false.
                        !call set_vpa_weights (conservative_wgts)

                        ! required integrals to ensure density conservation in field particle operator:
                        do iz = -nzgrid, nzgrid
                            call integrate_vmu(mw(:,:,iz,isa),iz,mwnorm(iz))
                        end do

                        do iz = -nzgrid, nzgrid
                            call integrate_vmu(modmw(:,:,iz,isa),iz,modmwnorm(iz))
                        end do

                        ! to ensure number conservation to machine precision, \Delta_{j=1} -> \Delta_{j=1} - F0/\int F0 dv * \int \Delta_{j=1} dv

                        ! AVB: to do - need to generalise this to higher-order terms in the field particle operator

                        if ((density_conservation_field).and.(ll==0).and.(jj==1)) then

                            do iz = -nzgrid, nzgrid
                                call integrate_vmu(deltaj(ll, jj, isa, isb, :, :, ia, iz), iz, deltajint(iz))
                            end do

                            do iz = -nzgrid, nzgrid
                                call integrate_vmu(deltaj_tp(ll, jj, isa, isb, :, :, ia, iz), iz, deltajint_tp(iz))
                            end do

                            ! check accuracy of deltaj:
                            !call open_output_file (tmpunit,'.deltaj_check')
                            !do iv = 1, nvpa
                            !    write (tmpunit,'(32es15.4e3)') ( deltaj(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0
                            !end do
                            !write (tmpunit,*)
                            !call close_output_file (tmpunit)

                            deltaj(ll, jj, isa, isb, :, :, ia, :) = deltaj(ll, jj, isa, isb, :, :, ia, :) &
                                - mw(:,:,:,isa) / spread(spread(mwnorm,1,nvpa),2,nmu) * spread(spread(deltajint,1,nvpa),2,nmu)

                            deltaj_tp(ll, jj, isa, isb, :, :, ia, :) = deltaj_tp(ll, jj, isa, isb, :, :, ia, :) &
                                - mw(:,:,:,isa) / spread(spread(mwnorm,1,nvpa),2,nmu) * spread(spread(deltajint_tp,1,nvpa),2,nmu)

                            !all open_output_file (tmpunit,'.deltaj_check_mod')
                            !do iv = 1, nvpa
                            !    write (tmpunit,'(32es15.4e3)') ( deltaj(ll, jj, isa, isb, iv, imu, ia, 0), imu=1, nmu ) ! at z = 0
                            !end do
                            !write (tmpunit,*)
                            !call close_output_file (tmpunit)

                            ! check density conservation to machine precision
                            !do iz = -nzgrid, nzgrid
                            !    call integrate_vmu(deltaj(ll, jj, isa, isb, :, :, ia, iz), iz, deltajint(iz))
                            !end do

                        end if

                    end do
                end do
            end do
        end do

        ! get normalisations for psijls
        do isa = 1, nspec
            do isb = 1, nspec
                do ll = 0, lmax
                    do iz = -nzgrid, nzgrid
                        do jj = 0, jmax
                            if ((ll==0).and.(jj==0)) then
                                psijnorm(ll,jj,isa,isb,iz) = 1 ! never used
                            else

                                ! note self-adjointness of deltaj after normalisation and without coll freqs (\Delta_j') is
                                ! \int x_a^l L_j^{l+1/2}(x_a^2) \Delta_j'^{ab}[\tilde{f}_b/F0b * F0b] x_a^2 dx_a
                                ! = ma^0.5/mb^0.5 * ma^3/mb^3 \int f_b/F0b \Delta_j^{ba}'[x_a^l L_j^{l+1/2}(x_a^2) \tilde{F0a}] x_b^2 dx_b

                                if ((exact_conservation).and.(ll==1).and.(jj==0)) then
                                    ! momentum conservation term for uniform mu grid
                                    call integrate_vmu( 3*legendre_vpamu(ll,0,:,:,iz)**2*vlaguerre_vmu(ll,jj,:,:,ia,iz)*deltaj(ll,jj,isa,isb,:,:,ia,iz) &
                                        * (spec(isa)%mass/spec(isb)%mass)**(-3) * (spec(isa)%mass/spec(isb)%mass)**-0.5, iz, psijnorm(ll,jj,isa,isb,iz) )
                                else if ((exact_conservation).and.(ll==0).and.(jj==1)) then
                                    ! energy conservation term for uniform mu grid
                                    call integrate_vmu( legendre_vpamu(ll,0,:,:,iz)**2*vlaguerre_vmu(ll,jj,:,:,ia,iz)*deltaj(ll,jj,isa,isb,:,:,ia,iz) &
                                        * (spec(isa)%mass/spec(isb)%mass)**(-3) * (spec(isa)%mass/spec(isb)%mass)**-0.5, iz, psijnorm(ll,jj,isa,isb,iz) )
                                else if ((exact_conservation_tp).and.(ll==0).and.(jj==1)) then
                                    ! energy conservation term for non-uniform mu grid
                                    call integrate_vmu( legendre_vpamu(ll,0,:,:,iz)**2*vlaguerre_vmu(ll,jj,:,:,ia,iz)*deltaj(ll,jj,isa,isb,:,:,ia,iz) &
                                        * (spec(isa)%mass/spec(isb)%mass)**(-3) * (spec(isa)%mass/spec(isb)%mass)**-0.5, iz, psijnorm(ll,jj,isa,isb,iz) )
                                else if ((exact_conservation_tp).and.(ll==1).and.(jj==0)) then
                                    ! momentum conservation term for non-uniform mu grid
                                    call integrate_vmu( 3*legendre_vpamu(ll,0,:,:,iz)**2*vlaguerre_vmu(ll,jj,:,:,ia,iz)*deltaj(ll,jj,isa,isb,:,:,ia,iz) &
                                        * (spec(isa)%mass/spec(isb)%mass)**(-3) * (spec(isa)%mass/spec(isb)%mass)**-0.5, iz, psijnorm(ll,jj,isa,isb,iz) )
                                else if (.not.((exact_conservation).or.(exact_conservation_tp)).and.(ll==0).and.(jj==1)) then
                                    ! non-exact conservation of energy
                                    if (spec(isa)%mass/spec(isb)%mass < 1.) then
                                        ! use self-adjointness to avoid resolution problems
                                        call integrate_vmu( (-velvpamu(:,:,iz)**2)*deltaj(ll,jj,isb,isa,:,:,ia,iz), iz, psijnorm(ll,jj,isa,isb,iz) )
                                    else
                                        call integrate_vmu( (-velvpamu(:,:,iz)**2)*deltaj(ll,jj,isa,isb,:,:,ia,iz) * (spec(isa)%mass/spec(isb)%mass)**(-3) &
                                            * (spec(isa)%mass/spec(isb)%mass)**-0.5, iz, psijnorm(ll,jj,isa,isb,iz) )
                                    end if
                                else
                                    ! non-exact conservation of momentum, and higher-order terms
                                    if (spec(isa)%mass/spec(isb)%mass < 1.) then
                                        ! use self-adjointness to avoid resolution problems
                                        call integrate_vmu( vlaguerre_vmu(ll,jj,:,:,ia,iz)*deltaj(ll,jj,isb,isa,:,:,ia,iz), iz, psijnorm(ll,jj,isa,isb,iz) )
                                    else
                                        call integrate_vmu( vlaguerre_vmu(ll,jj,:,:,ia,iz)*deltaj(ll,jj,isa,isb,:,:,ia,iz) * (spec(isa)%mass/spec(isb)%mass)**(-3) &
                                            * (spec(isa)%mass/spec(isb)%mass)**-0.5, iz, psijnorm(ll,jj,isa,isb,iz) )
                                    end if
                                end if

                            end if
                        end do
                    end do
                end do
            end do
        end do

        psijnorm = psijnorm / (4*pi) ! account for theta and phi integrals included in integrate_vmu

        ! to check self-adjointness
        !!print('Integral A =', np.trapz(vlinspace**lll * assoc_laguerre(vlinspace**2, jjj, lll+0.5) * fp0(vlinspace, ma, mb, jjj, lll) * vlinspace**2, x=vlinspace) )
        !!print('Integral B =', np.trapz(ma**0.5/mb**0.5 * ma**3/mb**3 * vlinspace**lll * assoc_laguerre(vlinspace**2, jjj, lll+0.5) * fp0(vlinspace, mb, ma, jjj, lll) * vlinspace**2, x=vlinspace) )
        !if (nspec > 1) then
        !    ia  = 1
        !    ll  = 0
        !    jj  = 1
        !    isa = 1
        !    isb = 2
        !    print*,''
        !    call integrate_vmu( vlaguerre_vmu(ll,jj,:,:,ia,0)*deltaj(ll,jj,isa,isb,:,:,ia,0), 0, psijnorm(ll,jj,isa,isb,0) )
        !    print*,'Integral A j1l0 =', psijnorm(ll,jj,isa,isb,0)
        !    call integrate_vmu( vlaguerre_vmu(ll,jj,:,:,ia,0)*deltaj(ll,jj,isb,isa,:,:,ia,0) * (spec(isa)%mass/spec(isb)%mass)**3.0 * (spec(isa)%mass/spec(isb)%mass)**0.5, 0, psijnorm(ll,jj,isb,isa,0) )
        !    print*,'Integral B j1l0 =', psijnorm(ll,jj,isb,isa,0)
        !    call integrate_vmu( (-velvpamu(:,:,0)**2)*deltaj(ll,jj,isb,isa,:,:,ia,0) * (spec(isa)%mass/spec(isb)%mass)**3.0 * (spec(isa)%mass/spec(isb)%mass)**0.5, 0, psijnorm(ll,jj,isb,isa,0) )
        !    print*,'Integral B j1l0b=', psijnorm(ll,jj,isb,isa,0)
        !    print*,''
        !    ia  = 1
        !    ll  = 1
        !    jj  = 0
        !    isa = 1
        !    isb = 2
        !    print*,''
        !    call integrate_vmu( vlaguerre_vmu(ll,jj,:,:,ia,0)*deltaj(ll,jj,isa,isb,:,:,ia,0), 0, psijnorm(ll,jj,isa,isb,0) )
        !    print*,'Integral A l1j0 =', psijnorm(ll,jj,isa,isb,0)
        !    call integrate_vmu( vlaguerre_vmu(ll,jj,:,:,ia,0)*deltaj(ll,jj,isb,isa,:,:,ia,0) * (spec(isa)%mass/spec(isb)%mass)**3.0 * (spec(isa)%mass/spec(isb)%mass)**0.5, 0, psijnorm(ll,jj,isb,isa,0) )
        !    print*,'Integral B l1j0 =', psijnorm(ll,jj,isb,isa,0)
        !    print*,''
        !end if

        ! to write interspec delt0 to file
        !if (nspec > 1) then
        !    do iv = 1, nvel_local
        !        call calc_delta0 (vel(iv), 1, 0, 1, 2, delt0test1(iv)) ! jj=1, ll=0, ie
        !        call calc_delta0 (vel(iv), 1, 0, 2, 1, delt0test2(iv)) ! jj=1, ll=0, ei
        !        call calc_delta0 (vel(iv), 0, 1, 1, 2, delt0test3(iv)) ! jj=1, ll=0, ie
        !        call calc_delta0 (vel(iv), 0, 2, 2, 2, delt0test4(iv)) ! jj=1, ll=0, ei
        !    end do
        !    call open_output_file (tmpunit,'.delt0test')
        !    do iv = 1, nvel_local
        !      write (tmpunit,'(9es15.4e3)') vel(iv), delt0test1(iv), delt0test2(iv), delt0test3(iv), delt0test4(iv)
        !    end do
        !    write (tmpunit,*)
        !    call close_output_file (tmpunit)
        !    call open_output_file (tmpunit,'.delt0l2test')
        !    do iv = 1, nvel_local
        !      write (tmpunit,'(9es15.4e3)') vel(iv), delt0test4(iv)
        !    end do
        !    write (tmpunit,*)
        !    call close_output_file (tmpunit)
        !end if

        conservative_wgts = .false.
        call set_vpa_weights (conservative_wgts)

        ! save deltajl for inspection
        !call open_output_file (tmpunit,'.deltaj1l1_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 1, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

        !call open_output_file (tmpunit,'.deltaj2l1_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 2, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

        !call open_output_file (tmpunit,'.deltaj3l1_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 3, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)


        !call open_output_file (tmpunit,'.deltaj1l1_ei')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 1, 2, 1, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

        !call open_output_file (tmpunit,'.deltaj2l1_ei')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 2, 2, 1, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

        !call open_output_file (tmpunit,'.deltaj3l1_ei')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 3, 2, 1, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)


        !call open_output_file (tmpunit,'.deltaj0l2_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(2, 0, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

        !print*,'maxval deltaj0l2'
        !print*,maxval(abs(deltaj(2, 0, 2, 2, :, :, 1, 0)))
        !print*,''

        !call open_output_file (tmpunit,'.deltaj1l2_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(2, 1, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

        !call open_output_file (tmpunit,'.deltaj2l2_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(2, 2, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)


        !call open_output_file (tmpunit,'.deltaj0l1_ee')
        !do iv = 1, nvpa
        !    write (tmpunit,*) ( deltaj(1, 0, 2, 2, iv, imu, 1, 0), imu=1, nmu ) ! at z = 0
        !end do
        !write (tmpunit,*)
        !call close_output_file (tmpunit)

    end subroutine init_deltaj_vmu

    subroutine get_testpart_density (isa, isb, g, fld)

        ! if isb==0:
        ! get the field tp_den_isa(g_isa), store it in fld(:,:,:,:,isa)

        ! if isa=0, fix the species indices of the operator tp_den, ie
        ! get the fields tp_den_isb(g_a), tp_den_isb(g_b), tp_den_isb(g_c) ...

        use mp, only: sum_allreduce
        use zgrid, only: nzgrid
        use vpamu_grids, only: integrate_vmu, set_vpa_weights, nvpa, nmu, vpa, mu
        use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, it_idx, is_idx
        use gyro_averages, only: aj0v, aj1v
        use constants, only: pi
        use dist_fn_arrays, only: kperp2
        use species, only: spec, nspec
        use stella_geometry, only: bmag
        use file_utils, only: open_output_file, close_output_file
        use stella_time, only: code_dt

        implicit none

        integer, intent (in) :: isa, isb
        complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
        complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

        integer :: ikxkyz, iky, ikx, iz, it, is, ia, iv, imu, tmpunit, ikxkyz_isb, is_b, iky_b, ikx_b, iz_b, it_b
        complex, dimension (:,:), allocatable :: g0
        complex, dimension (:), allocatable :: ghrs
        real :: clm, j1argnomu, intg
        logical :: conservative_wgts

        allocate (g0(nvpa,nmu))
        allocate (ghrs(nmu*nvpa))

        !conservative_wgts = .false.
        !call set_vpa_weights (conservative_wgts)

        ia = 1

        fld = 0.

        if (isb==0) then
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo,ikxkyz)
               ikx = ikx_idx(kxkyz_lo,ikxkyz)
               iz = iz_idx(kxkyz_lo,ikxkyz)
               is = is_idx(kxkyz_lo,ikxkyz)
               it = it_idx(kxkyz_lo,ikxkyz)
               if (is/=isa) cycle
               do iv = 1, nvpa
                   ghrs(nmu*(iv-1)+1 : nmu*iv) = g(iv,:,ikxkyz)
               end do
               ! blockmatrix_sum contains sum of interspecies test particle operators
               ghrs = matmul(-blockmatrix_sum(:,:,ikxkyz)/code_dt, ghrs)
               do iv = 1, nvpa
                   g0(iv,:) = ghrs(nmu*(iv-1)+1 : nmu*iv)
               end do
               call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
            end do
        else if (isa==isb) then
            ! apply the operator tp_den_isb[] to every species index of g; only used in the calculation of the response matrix
            ! where the species index, is, of g contains the response \delta h_a / \delta \psi^{a,is}; always given on an x_a grid
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo,ikxkyz)
               ikx = ikx_idx(kxkyz_lo,ikxkyz)
               iz = iz_idx(kxkyz_lo,ikxkyz)
               is = is_idx(kxkyz_lo,ikxkyz)
               it = it_idx(kxkyz_lo,ikxkyz)
               ! AVB: to do - cumbersome below, fix
               do iv = 1, nvpa
                   ghrs(nmu*(iv-1)+1 : nmu*iv) = g(iv,:,ikxkyz)
               end do
               do ikxkyz_isb = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                   is_b = is_idx(kxkyz_lo,ikxkyz_isb)
                   iky_b = iky_idx(kxkyz_lo,ikxkyz_isb)
                   ikx_b = ikx_idx(kxkyz_lo,ikxkyz_isb)
                   iz_b = iz_idx(kxkyz_lo,ikxkyz_isb)
                   it_b = it_idx(kxkyz_lo,ikxkyz_isb)
                   if ((is_b/=isb).or.(iky_b/=iky).or.(ikx_b/=ikx).or.(iz_b/=iz).or.(it_b/=it)) cycle
                   ghrs = matmul(-blockmatrix_sum(:,:,ikxkyz_isb)/code_dt,ghrs)
               end do
               do iv = 1, nvpa
                   g0(iv,:) = ghrs(nmu*(iv-1)+1 : nmu*iv)
               end do
               call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
            end do
        else
            fld = 0.
        end if

        deallocate (g0)
        deallocate (ghrs)

        call sum_allreduce (fld)

    end subroutine get_testpart_density

    subroutine init_fp_conserve

    use linear_solve, only: lu_decomposition
    use stella_time, only: code_dt
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu, nmu, nvpa, vpa, vperp2, set_vpa_weights, wgts_vpa, wgts_mu
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx, it_idx
    use dist_fn_arrays, only: gvmu
    use gyro_averages, only: aj0v, aj1v
    use fields, only: get_fields, get_fields_by_spec_idx
    use stella_geometry, only: bmag
    use job_manage, only: time_message, timer_local
    use file_utils, only: open_output_file, close_output_file
    use constants, only: pi

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is, it, iv, imu, idx, ix, ia, idx1, idx2, il, im, ij, mm, ll, jj, tmpunit, cnt
    integer :: il1, im1, ij1, mm1, ll1, jj1, is1, is2, is1a, is1b, is2a, is2b, isb2, isa2
    integer :: il2, im2, ij2, mm2, ll2, jj2, isa, isb
    logical :: conservative_wgts
    real :: dum2, dum3, t1

    complex, dimension (:,:,:,:), allocatable :: dum1
    complex, dimension (:,:,:,:,:), allocatable :: field
    complex, dimension (:,:), allocatable :: sumdelta
    complex, dimension (:,:), allocatable :: gvmutr
    complex, dimension (:), allocatable :: ghrs
    complex, dimension (:,:,:,:,:), allocatable :: response_vpamu

      if (.not.allocated(fp_response)) then
          if (fieldpart) then
              nresponse = 1 + (jmax+1)*(lmax+1)**2 * nspec**2
          else
              nresponse = 1
          end if
          allocate (fp_response(nresponse,nresponse,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); fp_response = 0.
          allocate (diff_idx(nresponse,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
          allocate (response_vpamu(nvpa,nmu,nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      end if

    allocate (dum1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))
    allocate (sumdelta(nvpa,nmu)); sumdelta = 0.
    allocate (gvmutr(nvpa,nmu))
    allocate (ghrs(nmu*nvpa)); ghrs = 0.
    allocate (deltajint(jmax+1,nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

    ia = 1

    ! set wgts to uniform to ensure exact conservation properties
    if (density_conservation) then
        conservative_wgts = .true.
        call set_vpa_weights (conservative_wgts)
    end if
    if (exact_conservation_tp) then
        conservative_wgts = .false.
        call set_vpa_weights (conservative_wgts)
    end if

    ! phi response
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo,ikxkyz)
         ikx = ikx_idx(kxkyz_lo,ikxkyz)
         iz = iz_idx(kxkyz_lo,ikxkyz)
         is = is_idx(kxkyz_lo,ikxkyz)

         ghrs = reshape(transpose(spread(ztmax(:,is),2,nmu)*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(jm0(:,iky,ikx,iz,is),1,nvpa)), shape=(/ nmu*nvpa /))
         call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,iky,ikx,iz,is), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,is), ghrs, nvpa*nmu, info)
         gvmu(:,:,ikxkyz) = transpose(reshape(ghrs, shape=(/nmu,nvpa/)))
    end do

  ! gvmu contains dhs/dphi
  ! for phi equation, need 1-P[dhs/dphi]
  call get_fields (gvmu, field(:,:,:,:,1), dum1, dist='h') ! note that get_fields sums over species, as required in response matrix

  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
      iky = iky_idx(kxkyz_lo,ikxkyz)
      ikx = ikx_idx(kxkyz_lo,ikxkyz)
      iz = iz_idx(kxkyz_lo,ikxkyz)
      it = it_idx(kxkyz_lo,ikxkyz)
      fp_response(1,1,ikxkyz) = 1.0-field(iky,ikx,iz,it,1)
  end do

  ! field particle operator

  ! collect all tp terms for species a in blockmatrix_sum(isa)
  ! required for tp density conservation to machine precision
  ! to do - avoid operations with blockmatrix or blockmatrix_sum, use band storage
  if (density_conservation_tp) then
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          is = is_idx(kxkyz_lo,ikxkyz)
          blockmatrix_sum(:,:,ikxkyz) = blockmatrix(:,:,ikxkyz,is)
          do isb = 1, nspec
              if (isb==is) cycle
              blockmatrix_sum(:,:,ikxkyz) = blockmatrix_sum(:,:,ikxkyz) + blockmatrix(:,:,ikxkyz,isb)
          end do
      end do
  end if

  ! field particle terms
  if (fieldpart) then
      idx1 = 1 ! first row (phi equation)
      do idx2 = 1, (jmax+1)*(lmax+1)**2 * nspec ! column indices

          isa = 1 + int((idx2 - 1)/float((jmax+1)*(lmax+1)**2))
          ij = 1 + mod(1+mod(idx2-1,(jmax+1)*(lmax+1)**2)-1, jmax+1)
          il = 1 + int(sqrt(1.0*(1+mod(idx2-1, (jmax+1)*(lmax+1)**2) - ij)/(jmax+1)))
          im = 1 + (1+mod(idx2-1, (jmax+1)*(lmax+1)**2) - ij)/(jmax+1) - (il-1)**2
          ll = il-1
          mm = -ll + im-1
          jj = ij-1

          ! get psi response, dh/dpsi_{ikn}, need responses:
          ! dh/dpsi_aa, dh/dpsi_ab, ...
          ! dh/dpsi_bb, dh/dpsi_ba, ...

          ! jj=0, ll=0 term is zero because Delta_j^l = 0 for jj=0, ll=0
          ! optionally replace this term with density conserving term
          ! that ensures conservation of density of the test particle operator to machine precision when using non-uniform mu-grid

          if ((jj==0).and.(ll==0).and.(density_conservation_tp)) then

              gvmu = 0.
              ! get testpart_den response for kperp = 0, only need dh_isa / d testpart_den_isa
              ! since this includes all interspecies test-particle operators
              do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                   iky = iky_idx(kxkyz_lo,ikxkyz)
                   ikx = ikx_idx(kxkyz_lo,ikxkyz)
                   iz = iz_idx(kxkyz_lo,ikxkyz)
                   is = is_idx(kxkyz_lo,ikxkyz)

                   if (is/=isa) cycle
                   do iv = 1, nvpa
                       do imu = 1, nmu
                           ghrs(nmu*(iv-1)+imu) = -code_dt*modmw(iv,imu,iz,is)/modmwnorm(iz)
                       end do
                   end do
                   call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,iky,ikx,iz,is), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,is), ghrs, nvpa*nmu, info)
                   do ix = 1, nvpa
                       gvmu(ix,:,ikxkyz) = ghrs(nmu*(ix-1)+1 : nmu*(ix-1)+nmu)
                   end do

              end do

              ! get phi_isa(dh_is2a/dh_is1), phi_isa(dh_is2a/dh_is2) ...
              call get_fields_by_spec_idx (isa, gvmu, field) ! AVB: check - using by_spec_idx instead of by_spec_mod

              do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo,ikxkyz)
                  ikx = ikx_idx(kxkyz_lo,ikxkyz)
                  iz = iz_idx(kxkyz_lo,ikxkyz)
                  it = it_idx(kxkyz_lo,ikxkyz)
                  fp_response(idx1,2+(idx2-1)*nspec:1+idx2*nspec,ikxkyz) = -field(iky,ikx,iz,it,:)
              end do
          else
              ! get the responses dh_is2a/dh_is1, dh_is2a/dh_is2, dh_is2a/dh_is3 ...
              call get_psi_response (ll, mm, jj, isa, gvmu)

              ! get phi_isa(dh_is2a/dh_is1), phi_isa(dh_is2a/dh_is2) ...
              call get_fields_by_spec_idx (isa, gvmu, field) ! AVB: check - using by_spec_idx instead of by_spec_mod

              do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo,ikxkyz)
                  ikx = ikx_idx(kxkyz_lo,ikxkyz)
                  iz = iz_idx(kxkyz_lo,ikxkyz)
                  it = it_idx(kxkyz_lo,ikxkyz)
                  fp_response(idx1,2+(idx2-1)*nspec:1+idx2*nspec,ikxkyz) = -field(iky,ikx,iz,it,:)
              end do

          end if

      end do

      idx2 = 1 ! first column
      do idx1 = 1, (jmax+1)*(lmax+1)**2 * nspec ! row indices

          isa = 1 + int((idx1 - 1)/float((jmax+1)*(lmax+1)**2))
          ij = 1 + mod(1+mod(idx1-1,(jmax+1)*(lmax+1)**2)-1, jmax+1)
          il = 1 + int(sqrt(1.0*(1+mod(idx1-1, (jmax+1)*(lmax+1)**2) - ij)/(jmax+1)))
          im = 1 + (1+mod(idx1-1, (jmax+1)*(lmax+1)**2) - ij)/(jmax+1) - (il-1)**2
          ll = il-1
          mm = -ll + im-1
          jj = ij-1

          ! get phi responses, dh_{isa}/dphi, dh_{isb}/dphi, dh_{isc}/dphi ... store in gvmu
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz  = iz_idx(kxkyz_lo,ikxkyz)
             is  = is_idx(kxkyz_lo,ikxkyz)
             it  = it_idx(kxkyz_lo,ikxkyz)
             do iv = 1, nvpa
                 do imu = 1, nmu
                     ghrs(nmu*(iv-1)+imu) = ztmax(iv,is)*maxwell_mu(1,iz,imu,is)*jm0(imu,iky,ikx,iz,is)
                 end do
             end do
             call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,iky,ikx,iz,is), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,is), ghrs, nvpa*nmu, info)
             do ix = 1, nvpa
                 gvmu(ix,:,ikxkyz) = ghrs(nmu*(ix-1)+1 : nmu*(ix-1)+nmu)
             end do
          end do

          if ((jj==0).and.(ll==0).and.(density_conservation_tp)) then

              ! get C_testpart[isa,isa+isb+isc...][dh_a/dphi], store in field(:,:,:,:,isa)
              field = 0.
              call get_testpart_density (isa, 0, gvmu, field)

              do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo,ikxkyz)
                  ikx = ikx_idx(kxkyz_lo,ikxkyz)
                  iz = iz_idx(kxkyz_lo,ikxkyz)
                  it = it_idx(kxkyz_lo,ikxkyz)
                  fp_response(2+(idx1-1)*nspec:1+idx1*nspec,idx2,ikxkyz) = -field(iky,ikx,iz,it,:)
              end do
          else
              ! get psi_{isa,isa}[dh_{isa}/dphi], psi_{isa,isb}[dh_{isb}/dphi], psi_{isa,isc}[dh_{isc}/dphi] ...
              call get_psi(gvmu,field,isa,0,ll,mm,jj)
              do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo,ikxkyz)
                  ikx = ikx_idx(kxkyz_lo,ikxkyz)
                  iz = iz_idx(kxkyz_lo,ikxkyz)
                  it = it_idx(kxkyz_lo,ikxkyz)
                  fp_response(2+(idx1-1)*nspec:1+idx1*nspec,idx2,ikxkyz) = -field(iky,ikx,iz,it,:)
              end do
          end if

      end do

      ! interior entries
      do idx1 = 1, (jmax+1)*(lmax+1)**2 * nspec ! row indices
          do idx2 = 1, (jmax+1)*(lmax+1)**2 * nspec ! column indices

              isa = 1 + int((idx1 - 1)/float((jmax+1)*(lmax+1)**2))
              ij1 = 1 + mod(1+mod(idx1-1,(jmax+1)*(lmax+1)**2)-1, jmax+1)
              il1 = 1 + int(sqrt(1.0*(1+mod(idx1-1, (jmax+1)*(lmax+1)**2) - ij1)/(jmax+1)))
              im1 = 1 + (1+mod(idx1-1, (jmax+1)*(lmax+1)**2) - ij1)/(jmax+1) - (il1-1)**2
              ll1 = il1-1
              mm1 = -ll1 + im1-1
              jj1 = ij1-1

              isb = 1 + int((idx2 - 1)/float((jmax+1)*(lmax+1)**2))
              ij2 = 1 + mod(1+mod(idx2-1,(jmax+1)*(lmax+1)**2)-1, jmax+1)
              il2 = 1 + int(sqrt(1.0*(1+mod(idx2-1, (jmax+1)*(lmax+1)**2) - ij2)/(jmax+1)))
              im2 = 1 + (1+mod(idx2-1, (jmax+1)*(lmax+1)**2) - ij2)/(jmax+1) - (il2-1)**2
              ll2 = il2-1
              mm2 = -ll2 + im2-1
              jj2 = ij2-1

              if ((jj2==0).and.(ll2==0).and.(density_conservation_tp)) then
                  gvmu = 0.
                  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                       iky = iky_idx(kxkyz_lo,ikxkyz)
                       ikx = ikx_idx(kxkyz_lo,ikxkyz)
                       iz = iz_idx(kxkyz_lo,ikxkyz)
                       is = is_idx(kxkyz_lo,ikxkyz)
                       if (is/=isb) cycle
                       do iv = 1, nvpa
                           do imu = 1, nmu
                               ghrs(nmu*(iv-1)+imu) = -code_dt*modmw(iv,imu,iz,is)/modmwnorm(iz)
                           end do
                       end do
                       call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,iky,ikx,iz,is), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,is), ghrs, nvpa*nmu, info)
                       do ix = 1, nvpa
                           gvmu(ix,:,ikxkyz) = ghrs(nmu*(ix-1)+1 : nmu*(ix-1)+nmu)
                       end do
                  end do
              else
                  gvmu = 0.
                  call get_psi_response (ll2, mm2, jj2, isb, gvmu)
              end if

              if ((jj1==0).and.(ll1==0).and.(density_conservation_tp)) then
                  ! get the fields C_testpart[isa,isa+isb+isc...][dh_a/dpsi^ab], C_testpart[isa,isa+isb+isc...][dh_a/dpsi^ac], C_testpart[isa,isa+isb+isc...][dh_a/dpsi^ad] ...
                  call get_testpart_density (isa, isb, gvmu, field)

                  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                      iky = iky_idx(kxkyz_lo,ikxkyz)
                      ikx = ikx_idx(kxkyz_lo,ikxkyz)
                      iz = iz_idx(kxkyz_lo,ikxkyz)
                      it = it_idx(kxkyz_lo,ikxkyz)
                      ! AVB: check - is index
                      fp_response(2+(idx1-1)*nspec+(isb-1),2+(idx2-1)*nspec:1+idx2*nspec,ikxkyz) = -field(iky,ikx,iz,it,:)
                  end do
              else
                  ! get fields Q_(isa,isb)[dh_{is}/dpsi^{is,isa}], Q_(isa,isb)[dh_{is}/dpsi^{is,isa}], ...
                  call get_psi (gvmu, field, isa, isb, ll1, mm1, jj1)
                  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                      iky = iky_idx(kxkyz_lo,ikxkyz)
                      ikx = ikx_idx(kxkyz_lo,ikxkyz)
                      iz  = iz_idx(kxkyz_lo,ikxkyz)
                      it  = it_idx(kxkyz_lo,ikxkyz)
                      ! AVB: check - is index
                      fp_response(2+(idx1-1)*nspec+(isb-1),2+(idx2-1)*nspec:1+idx2*nspec,ikxkyz) = -field(iky,ikx,iz,it,:)
                  end do
              end if
          end do
      end do

      ! add 1 to diagonal
      do idx1 = 2, 1+(jmax+1)*(lmax+1)**2 * nspec**2
          fp_response(idx1,idx1,:) = fp_response(idx1,idx1,:) + 1
      end do

  end if

  ! save response matrix, to examine contents
  !do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
    ! iky = iky_idx(kxkyz_lo,ikxkyz)
    ! ikx = ikx_idx(kxkyz_lo,ikxkyz)
    ! iz  = iz_idx(kxkyz_lo,ikxkyz)
    ! is  = is_idx(kxkyz_lo,ikxkyz)
    !it  = it_idx(kxkyz_lo,ikxkyz)
    !if ((iky==naky).and.(is==1).and.(iz==0)) then
    !     call open_output_file (tmpunit,'.fp_response')
    !      do idx1 = 1, nresponse
    !          write(tmpunit,*) (real(fp_response(idx1,idx2,ikxkyz)), idx2 = 1,nresponse) ! (12es15.4e3)
    !      end do
    !      write (tmpunit,*)
    !      call close_output_file (tmpunit)
    !  end if
  !end do

  ! LU decomposition for response
  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
      call lu_decomposition (fp_response(:,:,ikxkyz),diff_idx(:,ikxkyz),dum2)
  end do

  conservative_wgts = .false.
  call set_vpa_weights (conservative_wgts)

  deallocate (dum1, field)

end subroutine init_fp_conserve

subroutine get_psi_response (ll, mm, jj, isa, response)

      ! solve for responses dh_a / dpsi^{aa}, dh_a / dpsi^{ab}, dh_a / dpsi^{ac} ... for all species b, c, ...

      use finite_differences, only: tridag
      use linear_solve, only: lu_decomposition
      use stella_time, only: code_dt
      use species, only: nspec, spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu
      use vpamu_grids, only: nmu, nvpa,vpa, vperp2, mu
      use vpamu_grids, only: set_vpa_weights, wgts_vpa, wgts_mu
      use kt_grids, only: naky, nakx
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
      use dist_fn_arrays, only: gvmu
      use gyro_averages, only: aj0v, aj1v
      use fields, only: get_fields, get_fields_by_spec
      use stella_geometry, only: bmag
      use job_manage, only: time_message, timer_local
      use constants, only: pi
      use dist_fn_arrays, only: kperp2
      use file_utils, only: open_output_file, close_output_file

      implicit none

      complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (out) :: response
      integer, intent (in) :: ll, mm, jj, isa
      complex, dimension (:), allocatable :: ghrs
      integer :: ikxkyz, iky, ikx, iz, it, is, ia, iv, imu, ix, tmpunit
      real :: clm, j1arg

      allocate (ghrs(nmu*nvpa))

      ia = 1

      clm = sqrt(((2*ll+1)*gamma(ll-mm+1.))/(4*pi*gamma(ll+mm+1.)))

      ! calculate response dh/dpsi_jlm, for unit impulse to psi_jlm
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo,ikxkyz)
         ikx = ikx_idx(kxkyz_lo,ikxkyz)
         iz = iz_idx(kxkyz_lo,ikxkyz)
         is = is_idx(kxkyz_lo,ikxkyz) ! isb

         ! supply unit impulse to psi_j^(lm)^{isa,isb}
         do iv = 1, nvpa
             do imu = 1, nmu
                 if (mm == 0) then
                     ghrs(nmu*(iv-1)+imu) = code_dt*spec(isa)%vnew(is)*clm*legendre_vpamu(ll,mm,iv,imu,iz)&
                        *jm(imu,mm,iky,ikx,iz,isa)*(spec(isa)%mass/spec(is)%mass)**-1.5*deltaj(ll,jj,isa,is,iv,imu,ia,iz)
                 else if (mm > 0) then
                     ghrs(nmu*(iv-1)+imu) = code_dt*spec(isa)%vnew(is)*clm*legendre_vpamu(ll,mm,iv,imu,iz)&
                        *jm(imu,mm,iky,ikx,iz,isa)*(spec(isa)%mass/spec(is)%mass)**-1.5*deltaj(ll,jj,isa,is,iv,imu,ia,iz)
                 else if (mm < 0) then
                     ghrs(nmu*(iv-1)+imu) = (-1)**mm*code_dt*spec(isa)%vnew(is)*clm*legendre_vpamu(ll,mm,iv,imu,iz)&
                        *jm(imu,abs(mm),iky,ikx,iz,isa)*(spec(isa)%mass/spec(is)%mass)**-1.5*deltaj(ll,jj,isa,is,iv,imu,ia,iz)
                 end if
             end do
         end do

        ! solve for response
        ! need to solve [1 - Deltat C_{test}] dh_a/dpsi^ab = delta_{ab} for dh_a/dpsi^ab. Here, C_{test} includes self collisions and a-b collisions,
        call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, &
            cdiffmat_band(:,:,iky,ikx,iz,isa), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,isa), ghrs, nvpa*nmu, info)

        do iv = 1, nvpa
            response(iv,:,ikxkyz) = ghrs(nmu*(iv-1)+1 : nmu*(iv-1)+nmu)
        end do

        ! to zero l=1,j=1 term:
        if (no_j1l1) then
            if ((ll==1).and.(jj==1)) then
                response(:,:,ikxkyz) = 0.
            end if
        end if

        if (no_j1l2) then
            if ((ll==2).and.(jj==1)) then
                response(:,:,ikxkyz) = 0.
            end if
        end if

        if (no_j0l2) then
            if ((ll==2).and.(jj==0)) then
                response(:,:,ikxkyz) = 0.
            end if
        end if

        if (spitzer_problem) then
            if (.not.((isa==2).and.(is==2))) then
                response(:,:,ikxkyz) = 0.
            end if
        end if

      end do

      deallocate (ghrs)

  end subroutine get_psi_response

  subroutine get_psi (g, fld, isa, isb, ll, mm, jj)

      ! if isb==0:
      ! get the fields psi_aa^lmj(g_a), psi_ab^lmj(g_b), psi_ac^lmj(g_c) ..., for species b, c, ...
      ! if isb/=0, fix the species indices of the operator psi_ab, ie
      ! get the fields psi_ab^lmj(g_a), psi_ab^lmj(g_b), psi_ab^lmj(g_c) ...

      use mp, only: sum_allreduce
      use zgrid, only: nzgrid
      use vpamu_grids, only: integrate_vmu, set_vpa_weights, nvpa, nmu, vpa, mu
      use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use gyro_averages, only: aj0v, aj1v
      use constants, only: pi
      use dist_fn_arrays, only: kperp2
      use species, only: spec, nspec
      use stella_geometry, only: bmag
      use file_utils, only: open_output_file, close_output_file
      use stella_time, only: code_dt

      implicit none

      complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
      complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld
      integer, intent (in) :: isa, isb, ll, mm, jj

      integer :: ikxkyz, iky, ikx, iz, it, is, ia, iv, imu, tmpunit, ikxkyz_isb, is_b, iky_b, ikx_b, iz_b, it_b
      complex, dimension (:,:), allocatable :: g0
      complex, dimension (:), allocatable :: ghrs
      real :: clm, j1argnomu, intg
      logical :: conservative_wgts

      allocate (g0(nvpa,nmu))
      allocate (ghrs(nmu*nvpa))

      if (density_conservation) then
          conservative_wgts = .true.
          call set_vpa_weights (conservative_wgts)
      else
          conservative_wgts = .false.
          call set_vpa_weights (conservative_wgts)
      end if

      clm = (-1)**mm*sqrt(((2*ll+1)*gamma(ll+mm+1.))/(4*pi*gamma(ll-mm+1.)))

      ia = 1

      if (isb==0) then
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             is = is_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             if (mm == 0) then

                 if ( (exact_conservation).and.( ((ll==0).and.(jj==1)).or.((ll==1).and.(jj==0)) ) ) then
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,-mm,iky,ikx,iz,is),1,nvpa)*legendre_vpamu(ll,-mm,:,:,iz)*deltaj_tp(ll,jj,is,isa,:,:,ia,iz)

                 else if ( (exact_conservation_tp).and.((ll==1).and.(jj==0)) ) then
                    do iv = 1, nvpa
                        ghrs(nmu*(iv-1)+1 : nmu*iv) = clm*(spec(is)%mass/spec(isa)%mass)**2*jm(:,-mm,iky,ikx,iz,is)*g(iv,:,ikxkyz)
                    end do
                    ghrs = matmul(-blockmatrix(:,:,ikxkyz,isa)/code_dt/spec(is)%vnew(isa), ghrs)
                    do iv = 1, nvpa
                        g0(iv,:) = -vpa(iv)*ghrs(nmu*(iv-1)+1 : nmu*iv)
                    end do
                    if (spec(is)%vnew(isa)==0) then
                        g0 = 0.
                    end if

                 else if ( (exact_conservation_tp).and.((ll==0).and.(jj==1)) ) then
                     do iv = 1, nvpa
                         ghrs(nmu*(iv-1)+1 : nmu*iv) = clm*(spec(is)%mass/spec(isa)%mass)**1.5*jm(:,-mm,iky,ikx,iz,is)*g(iv,:,ikxkyz)
                     end do
                     ghrs = matmul(-blockmatrix(:,:,ikxkyz,isa)/code_dt/spec(is)%vnew(isa), ghrs)
                     do iv = 1, nvpa
                         g0(iv,:) = velvpamu(iv,:,iz)**2*ghrs(nmu*(iv-1)+1 : nmu*iv)
                     end do
                     if (spec(is)%vnew(isa)==0) then
                         g0 = 0.
                     end if

                 else
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,-mm,iky,ikx,iz,is),1,nvpa)*legendre_vpamu(ll,-mm,:,:,iz)*deltaj(ll,jj,is,isa,:,:,ia,iz)

                 end if
             else if (mm < 0) then
                 if ((exact_conservation).and.( ((ll==0).and.(jj==1)).or.((ll==1).and.(jj==0)) ) ) then
                     g0 = (-1)**mm*clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,abs(mm),iky,ikx,iz,is),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj_tp(ll,jj,is,isa,:,:,ia,iz)
                 else
                     g0 = (-1)**mm*clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,abs(mm),iky,ikx,iz,is),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj(ll,jj,is,isa,:,:,ia,iz)
                 end if

             else if (mm > 0) then
                 if ((exact_conservation).and.( ((ll==0).and.(jj==1)).or.((ll==1).and.(jj==0)) ) ) then
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,mm,iky,ikx,iz,is),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj_tp(ll,jj,is,isa,:,:,ia,iz)
                 else
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,mm,iky,ikx,iz,is),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj(ll,jj,is,isa,:,:,ia,iz)
                 end if
             end if

             call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))

             ! to zero l=1,j=1 term:
             if (no_j1l1) then
                 if ((ll==1).and.(jj==1)) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if
             if (no_j1l2) then
                 if ((ll==2).and.(jj==1)) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if
             if (no_j0l2) then
                 if ((ll==2).and.(jj==0)) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if

             if (spitzer_problem) then
                 if (.not.((isa==2).and.(is==2))) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if

          end do

      else ! isb /= 0:
          ! apply the operator \psi_ab^jlm[] to every species index of g; only used in the calculation of the response matrix
          ! where the species index, is, of g contains the response \delta h_a / \delta \psi^{a,is}
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             is = is_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             if (mm == 0) then
                 if ((exact_conservation).and.( ((ll==0).and.(jj==1)).or.((ll==1).and.(jj==0)) ) ) then
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,mm,iky,ikx,iz,isb),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj_tp(ll,jj,isb,isa,:,:,ia,iz)

                 else if ( (exact_conservation_tp).and.((ll==1).and.(jj==0)) ) then
                    do iv = 1, nvpa
                        ghrs(nmu*(iv-1)+1 : nmu*iv) = clm*(spec(isb)%mass/spec(isa)%mass)**2&
                            *jm(:,mm,iky,ikx,iz,isb)*g(iv,:,ikxkyz)
                    end do
                    do ikxkyz_isb = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                        is_b = is_idx(kxkyz_lo,ikxkyz_isb)
                        iky_b = iky_idx(kxkyz_lo,ikxkyz_isb)
                        ikx_b = ikx_idx(kxkyz_lo,ikxkyz_isb)
                        iz_b = iz_idx(kxkyz_lo,ikxkyz_isb)
                        it_b = it_idx(kxkyz_lo,ikxkyz_isb)
                        if ((is_b/=isb).or.(iky_b/=iky).or.(ikx_b/=ikx).or.(iz_b/=iz).or.(it_b/=it)) cycle
                        ghrs = matmul(-blockmatrix(:,:,ikxkyz_isb,isa)/code_dt/spec(isb)%vnew(isa),ghrs)
                    end do
                    do iv = 1, nvpa
                        g0(iv,:) = -vpa(iv)*ghrs(nmu*(iv-1)+1 : nmu*iv)
                    end do
                    if (spec(isb)%vnew(isa)==0) then
                        g0 = 0.
                    end if

                 else if ( (exact_conservation_tp).and.((ll==0).and.(jj==1)) ) then
                     do iv = 1, nvpa
                         ghrs(nmu*(iv-1)+1 : nmu*iv) = clm*(spec(isb)%mass/spec(isa)%mass)**1.5&
                            *jm(:,mm,iky,ikx,iz,isb)*g(iv,:,ikxkyz)
                     end do
                     do ikxkyz_isb = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                         is_b = is_idx(kxkyz_lo,ikxkyz_isb)
                         iky_b = iky_idx(kxkyz_lo,ikxkyz_isb)
                         ikx_b = ikx_idx(kxkyz_lo,ikxkyz_isb)
                         iz_b = iz_idx(kxkyz_lo,ikxkyz_isb)
                         it_b = it_idx(kxkyz_lo,ikxkyz_isb)
                         if ((is_b/=isb).or.(iky_b/=iky).or.(ikx_b/=ikx).or.(iz_b/=iz).or.(it_b/=it)) cycle
                         ghrs = matmul(-blockmatrix(:,:,ikxkyz_isb,isa)/code_dt/spec(isb)%vnew(isa),ghrs)
                     end do
                     do iv = 1, nvpa
                         g0(iv,:) = velvpamu(iv,:,iz)**2*ghrs(nmu*(iv-1)+1 : nmu*iv)
                     end do
                     if (spec(isb)%vnew(isa)==0) then
                         g0 = 0.
                     end if
                 else
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,mm,iky,ikx,iz,isb),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj(ll,jj,isb,isa,:,:,ia,iz)
                 end if
             else if (mm < 0) then
                 if ((exact_conservation).and.( ((ll==0).and.(jj==1)).or.((ll==1).and.(jj==0)) ) ) then
                     g0 = (-1)**mm*clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,abs(mm),iky,ikx,iz,isb),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj_tp(ll,jj,isb,isa,:,:,ia,iz)
                 else
                     g0 = (-1)**mm*clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,abs(mm),iky,ikx,iz,isb),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj(ll,jj,isb,isa,:,:,ia,iz)
                 end if

             else if (mm > 0) then
                 if ((exact_conservation).and.( ((ll==0).and.(jj==1)).or.((ll==1).and.(jj==0)) ) ) then
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,mm,iky,ikx,iz,isb),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj_tp(ll,jj,isb,isa,:,:,ia,iz)
                 else
                     g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,mm,iky,ikx,iz,isb),1,nvpa)&
                        *legendre_vpamu(ll,-mm,:,:,iz)*deltaj(ll,jj,isb,isa,:,:,ia,iz)
                 end if
             end if

             call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))

             ! to zero l=1,j=1 term:
             if (no_j1l1) then
                 if ((ll==1).and.(jj==1)) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if
             if (no_j1l2) then
                 if ((ll==2).and.(jj==1)) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if
             if (no_j0l2) then
                 if ((ll==2).and.(jj==0)) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if

             if (spitzer_problem) then
                 if (.not.((isa==2).and.(isb==2))) then
                     fld(:,:,:,:,is) = 0.
                 end if
             end if

          end do
      end if

      ! normalise psijs
      if (isb==0) then
          do is = 1, nspec
              do iz = -nzgrid, nzgrid
                 fld(:,:,iz,:,is) = fld(:,:,iz,:,is) / psijnorm(ll,jj,isa,is,iz)
              end do
          end do
      else
          do iz = -nzgrid, nzgrid
              fld(:,:,iz,:,:) = fld(:,:,iz,:,:) / psijnorm(ll,jj,isa,isb,iz)
          end do
      end if

      deallocate (g0)
      deallocate (ghrs)

      conservative_wgts = .false.
      call set_vpa_weights (conservative_wgts)

      call sum_allreduce (fld)

  end subroutine get_psi

  subroutine init_vpadiff_matrix

    use stella_time, only: code_dt
    use species, only: nspec, spec
    use vpamu_grids, only: dvpa, vpa, nvpa
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use stella_geometry, only: bmag
    use dist_fn_arrays, only: kperp2

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is
    integer :: ia

    if (.not.allocated(aa_vpa)) allocate (aa_vpa(nvpa,nspec))
!    if (.not.allocated(bb_vpa)) allocate (bb_vpa(nvpa,nspec))
    if (.not.allocated(bb_vpa)) allocate (bb_vpa(nvpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if (.not.allocated(cc_vpa)) allocate (cc_vpa(nvpa,nspec))

    ! deal with boundary points (BC is f(vpa)=0 beyond +/- vpa_max)
    aa_vpa(1,:) = 0.0 ; cc_vpa(nvpa,:) = 0.0
    ! 2nd order centered differences for d/dvpa (1/2 dh/dvpa + vpa h)
    do is = 1, nspec
       aa_vpa(2:,is) = -code_dt*spec(is)%vnew(is)*0.5*(1.0/dvpa-vpa(:nvpa-1))/dvpa
!       bb_vpa(:,is) = 1.0+code_dt*spec(is)%vnew(is)/dvpa**2
       cc_vpa(:nvpa-1,is) = -code_dt*spec(is)%vnew(is)*0.5*(1.0/dvpa+vpa(2:))/dvpa
    end do

    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       bb_vpa(:,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          *  (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 + 1./dvpa**2)
    end do

  end subroutine init_vpadiff_matrix

  subroutine init_mudiff_matrix

    use stella_time, only: code_dt
    use species, only: nspec, spec
    use zgrid, only: nzgrid
    use stella_geometry, only: bmag
    use vpamu_grids, only: dmu, mu, nmu
    use vpamu_grids, only: dmu_cell, mu_cell, wgts_mu_bare
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use dist_fn_arrays, only: kperp2

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is
    integer :: ia
    ! TMP FOR TESTING -- MAB
!    integer :: imu

    if (.not.allocated(aa_mu)) allocate (aa_mu(-nzgrid:nzgrid,nmu,nspec))
    if (.not.allocated(bb_mu)) allocate (bb_mu(nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if (.not.allocated(cc_mu)) allocate (cc_mu(-nzgrid:nzgrid,nmu,nspec))

    ia = 1

    ! deal with boundary points (BC is f(mu)=0 beyond mu_max and collision operator vanishes for mu -> 0)
    aa_mu(:,1,:) = 0.0 ; cc_mu(:,nmu,:) = 0.0
    ! 2nd order centered differences for dt * nu * d/dmu (mu/B*dh/dmu + 2*mu*h)
    do is = 1, nspec
       do iz = -nzgrid, nzgrid
          aa_mu(iz,2:,is)     = -code_dt*spec(is)%vnew(is)*mu_cell(:nmu-1)*(1.0/(bmag(ia,iz)*dmu)-1.0)/wgts_mu_bare(2:)
          cc_mu(iz,:nmu-1,is) = -code_dt*spec(is)%vnew(is)*mu_cell(:nmu-1)*(1.0/(bmag(ia,iz)*dmu)+1.0)/wgts_mu_bare(:nmu-1)
       end do
    end do

    !1st derivative here is  2 d/dmu( mu h) -> (mu(i+1/2)h(i+1/2) - mu(i-1/2)h(i-1/2))/(mu(i+1/2)-mu(i-1/2))
    !where h(i+1/2) = 0.5*[h(i+1)+h(i)] and mu(i+1/2) = 0.5*(mu(i+1)+mu(i))
    !2nd derivative here is d/dmu ( mu dh/dmu ) = (mu(i+1/2)h'(i+1/2) - mu(i-1/2)h'(i-1/2))/(mu(i+1/2)-mu(i-1/2))
    !where h'(i+1/2) = (h(i+1)-h(i))/(mu(i+1)-mu(i))
    !left  endpoint i=1   -> mu(i-1/2) = 0
    !right endpoint i=nmu -> h(i+1/2) = h'(i+1/2) = 0 (these are the two boundary conditions)
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       bb_mu(1,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          * (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 &
            + dmu_cell(1)*(1.0/(dmu(1)*bmag(ia,iz)) - 1.0)/wgts_mu_bare(1))
       bb_mu(2:nmu-1,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          * (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 &
            + (mu_cell(2:nmu-1)/dmu(2:)+mu_cell(:nmu-2)/dmu(:nmu-2)) &
            /(wgts_mu_bare(2:nmu-1)*bmag(ia,iz)) - dmu_cell(2:nmu-1)/wgts_mu_bare(2:nmu-1))
       bb_mu(nmu,ikxkyz) = 1.0 + code_dt*spec(is)%vnew(is) &
          * (0.25*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2 &
            + mu_cell(nmu-1)*(1.0/(dmu(nmu-1)*bmag(ia,iz)) + 1.0)/wgts_mu_bare(nmu))
    end do

  end subroutine init_mudiff_matrix

  subroutine init_vpadiff_conserve

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_decomposition, lu_inverse
    use stella_time, only: code_dt
    use species, only: nspec, spec, has_electron_species
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu
    use vpamu_grids, only: nmu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx, zonal_mode
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use stella_geometry, only: dl_over_b
    use dist_fn_arrays, only: gvmu
    use gyro_averages, only: aj0v
    use fields, only: get_fields, get_fields_by_spec, efac, gamtot_h
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    implicit none

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    integer :: imu
    integer :: idx
    logical :: conservative_wgts
    real :: dum2
    complex, dimension (:,:,:,:), allocatable :: dum1
    complex, dimension (:,:,:,:,:), allocatable :: field
    complex, dimension (:,:), allocatable :: temp_mat

    ia = 1

    nresponse_vpa = 1
    if (momentum_conservation) nresponse_vpa = nresponse_vpa + nspec
    if (energy_conservation) nresponse_vpa = nresponse_vpa + nspec

    if (.not.allocated(vpadiff_response)) then
       allocate (vpadiff_response(nresponse_vpa,nresponse_vpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       vpadiff_response = 0.
       allocate (vpadiff_idx(nresponse_vpa,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    end if

    if (.not.has_electron_species(spec).and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
      if(.not.allocated(vpadiff_zf_response)) then
        allocate (vpadiff_zf_response(nresponse_vpa,nresponse_vpa,nakx))
        vpadiff_zf_response = 0.
        allocate (vpadiff_zf_idx(nresponse_vpa,nakx))
      endif
    endif

    allocate (dum1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))

    ! set wgts to be equally spaced to ensure exact conservation properties
    conservative_wgts = .true.
    call set_vpa_weights (conservative_wgts)

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          gvmu(:,imu,ikxkyz) = ztmax(:,is)*maxwell_mu(1,iz,imu,is)*aj0v(imu,ikxkyz)
          call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), gvmu(:,imu,ikxkyz))
       end do
    end do

    ! gvmu contains dhs/dphi
    ! for phi equation, need 1-P[dhs/dphi]
    ! for upar equations, need -Us[dhs/dphi]
    ! for energy conservation, need -Qs[dhs/dphi]
    call get_fields (gvmu, field(:,:,:,:,1), dum1, dist='h', skip_fsa=.true.)

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       vpadiff_response(1,1,ikxkyz) = 1.0-field(iky,ikx,iz,it,1)
    end do
    idx = 2
    if (momentum_conservation) then
       call get_upar (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
       idx = idx + nspec
    end if
    if (energy_conservation) then
       call get_temp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
    end if
    idx = 2

    if (momentum_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do imu = 1, nmu
             gvmu(:,imu,ikxkyz) = 2.*code_dt*spec(is)%vnew(is)*vpa*aj0v(imu,ikxkyz)*maxwell_vpa(:,is)*maxwell_mu(1,iz,imu,is)
             call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), gvmu(:,imu,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dupars
       ! need to get -Ps[dhs/dupars] for phi equation
       ! need to get 1-Us[dhs/dupars] for momentum conservation
       ! need to get -Qs[dhs/dupars] for energy conservation
       call get_fields_by_spec (gvmu, field, skip_fsa=.true.)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       call get_upar (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             vpadiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do

       if (energy_conservation) then
          call get_temp (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                vpadiff_response(idx+is+nspec-1,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if
       idx = idx + nspec
    end if

    if (energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do imu = 1, nmu
             gvmu(:,imu,ikxkyz) = 2.*code_dt*spec(is)%vnew(is)*(vpa**2+vperp2(1,iz,imu)-1.5) &
                  *aj0v(imu,ikxkyz)*maxwell_vpa(:,is)*maxwell_mu(1,iz,imu,is)
             call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), gvmu(:,imu,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dQs
       ! need to get -Ps[dhs/dQs] for phi equation
       ! need to get 1-Us[dhs/dQs] for momentum conservation
       ! need to get -Qs[dhs/dQs] for energy conservation
       call get_fields_by_spec (gvmu, field, skip_fsa=.true.)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          vpadiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       if (momentum_conservation) then
          call get_upar (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                vpadiff_response(idx+is-1-nspec,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if

       call get_temp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             vpadiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do
    end if

    ! now get LU decomposition for vpadiff_response
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       call lu_decomposition (vpadiff_response(:,:,ikxkyz),vpadiff_idx(:,ikxkyz),dum2)
    end do

    ! if electrons are adiabatic, compute the matrices for the flux-surface-average
    if (.not.has_electron_species(spec).and.zonal_mode(1) &
        .and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
      allocate (temp_mat(nresponse_vpa,nresponse_vpa))
      vpadiff_zf_response = 0.0
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        it = it_idx(kxkyz_lo,ikxkyz)
        if (iky.ne.1.or.it.ne.1) cycle

        !calculate inverse of vpadiff_response
        call lu_inverse (vpadiff_response(:,:,ikxkyz),vpadiff_idx(:,ikxkyz),temp_mat)

        !calculate -inv(vpadiff_response).Q, where Q has a single entry
        do idx = 1, nresponse_vpa
          vpadiff_zf_response(idx,1,ikx) = vpadiff_zf_response(idx,1,ikx) &
                                           - temp_mat(idx,1)*(efac/gamtot_h)*dl_over_b(ia,iz)
        enddo
      end do

      !finish the flux surface average
      call sum_allreduce(vpadiff_zf_response)

      !calculate 1 - fsa(inv(vpadiff_response).Q)
      do idx = 1, nresponse_vpa
        vpadiff_zf_response(idx,idx,:) = vpadiff_zf_response(idx,idx,:) + 1.0
      enddo

      do ikx = 1,nakx
         call lu_decomposition (vpadiff_zf_response(:,:,ikx),vpadiff_zf_idx(:,ikx),dum2)
      end do

      deallocate (temp_mat)
    endif

    ! reset wgts to default setting
    conservative_wgts = .false.
    call set_vpa_weights (conservative_wgts)

    deallocate (dum1, field)

  end subroutine init_vpadiff_conserve

  subroutine init_mudiff_conserve

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_decomposition, lu_inverse
    use stella_time, only: code_dt
    use species, only: nspec, spec, has_electron_species
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: ztmax, maxwell_vpa, maxwell_mu
    use vpamu_grids, only: nvpa, vpa, vperp2
    use kt_grids, only: naky, nakx, zonal_mode
    use stella_geometry, only: dl_over_b, bmag
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use dist_fn_arrays, only: gvmu, kperp2
    use gyro_averages, only: aj0v, aj1v
    use fields, only: get_fields, get_fields_by_spec, efac, gamtot_h
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    implicit none

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    integer :: iv
    integer :: idx
    real :: dum2
    complex, dimension (:,:), allocatable :: temp_mat
    complex, dimension (:,:,:,:), allocatable :: dum1
    complex, dimension (:,:,:,:,:), allocatable :: field

    nresponse_mu = 1
    if (momentum_conservation) nresponse_mu = nresponse_mu + nspec
    if (energy_conservation) nresponse_mu = nresponse_mu + nspec

    if (.not.allocated(mudiff_response)) then
       allocate (mudiff_response(nresponse_mu,nresponse_mu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       mudiff_response = 0.
       allocate (mudiff_idx(nresponse_mu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    end if

    if (.not.has_electron_species(spec).and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
      if(.not.allocated(mudiff_zf_response)) then
        allocate (mudiff_zf_response(nresponse_mu,nresponse_mu,nakx))
        mudiff_zf_response = 0.
        allocate (mudiff_zf_idx(nresponse_mu,nakx))
      endif
    endif

    allocate (dum1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do iv = 1, nvpa
          gvmu(iv,:,ikxkyz) = ztmax(iv,is)*maxwell_mu(1,iz,:,is)*aj0v(:,ikxkyz)
          call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), gvmu(iv,:,ikxkyz))
       end do
    end do

    ! gvmu contains dhs/dphi
    ! for phi equation, need 1-P[dhs/dphi]
    ! for uperp equations, need -Us[dhs/dphi]
    ! for energy conservation, need -Qs[dhs/dphi]
    call get_fields (gvmu, field(:,:,:,:,1), dum1, dist='h', skip_fsa =.true.)

    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       mudiff_response(1,1,ikxkyz) = 1.0-field(iky,ikx,iz,it,1)
    end do
    idx = 2
    if (momentum_conservation) then
       call get_uperp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
       idx = idx + nspec
    end if
    if (energy_conservation) then
       call get_temp_mu (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(idx:idx+nspec-1,1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do
    end if
    idx = 2

    if (momentum_conservation) then
       ia = 1
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do iv = 1, nvpa
             gvmu(iv,:,ikxkyz) = 2.*code_dt*spec(is)%vnew(is)*kperp2(iky,ikx,ia,iz)*vperp2(ia,iz,:) &
                  *(spec(is)%smz/bmag(ia,iz))**2*aj1v(:,ikxkyz)*maxwell_vpa(iv,is)*maxwell_mu(ia,iz,:,is)
             call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), gvmu(iv,:,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dupars
       ! need to get -Ps[dhs/dupars] for phi equation
       ! need to get 1-Us[dhs/dupars] for momentum conservation
       ! need to get -Qs[dhs/dupars] for energy conservation
       call get_fields_by_spec (gvmu, field, skip_fsa =.true.)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       call get_uperp (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             mudiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do

       if (energy_conservation) then
          call get_temp_mu (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                mudiff_response(idx+is+nspec-1,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if
       idx = idx + nspec
    end if

    if (energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iz = iz_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          do iv = 1, nvpa
             gvmu(iv,:,ikxkyz) = 2.0*code_dt*spec(is)%vnew(is)*(vpa(iv)**2+vperp2(1,iz,:)-1.5) &
                  *aj0v(:,ikxkyz)*maxwell_vpa(iv,is)*maxwell_mu(1,iz,:,is)
             call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), gvmu(iv,:,ikxkyz))
          end do
       end do
       ! gvmu now contains dhs/dQs
       ! need to get -Ps[dhs/dQs] for phi equation
       ! need to get 1-Us[dhs/dQs] for momentum conservation
       ! need to get -Qs[dhs/dQs] for energy conservation
       call get_fields_by_spec (gvmu, field, skip_fsa=.true.)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          mudiff_response(1,idx:idx+nspec-1,ikxkyz) = -field(iky,ikx,iz,it,:)
       end do

       if (momentum_conservation) then
          call get_uperp (gvmu, field)
          do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
             iky = iky_idx(kxkyz_lo,ikxkyz)
             ikx = ikx_idx(kxkyz_lo,ikxkyz)
             iz = iz_idx(kxkyz_lo,ikxkyz)
             it = it_idx(kxkyz_lo,ikxkyz)
             do is = 1, nspec
                mudiff_response(idx+is-1-nspec,idx+is-1,ikxkyz) = -field(iky,ikx,iz,it,is)
             end do
          end do
       end if

       call get_temp_mu (gvmu, field)
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          do is = 1, nspec
             mudiff_response(idx+is-1,idx+is-1,ikxkyz) = 1.0-field(iky,ikx,iz,it,is)
          end do
       end do
    end if

    ! now get LU decomposition for mudiff_response
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       call lu_decomposition (mudiff_response(:,:,ikxkyz),mudiff_idx(:,ikxkyz),dum2)
    end do

    ! if electrons are adiabatic, compute the matrices for the flux-surface-average
    if (.not.has_electron_species(spec).and.zonal_mode(1) &
        .and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
      allocate (temp_mat(nresponse_mu,nresponse_mu))
      mudiff_zf_response = 0.0
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        it = it_idx(kxkyz_lo,ikxkyz)
        if (iky.ne.1.or.it.ne.1) cycle

        !calculate inverse of mudiff_response
        call lu_inverse (mudiff_response(:,:,ikxkyz),mudiff_idx(:,ikxkyz),temp_mat)

        !calculate -inv(mudiff_response).Q, where Q has a single entry
        do idx = 1, nresponse_mu
          mudiff_zf_response(idx,1,ikx) = mudiff_zf_response(idx,1,ikx) &
                                           - temp_mat(idx,1)*(efac/gamtot_h)*dl_over_b(ia,iz)
        enddo
      end do

      !finish the flux surface average
      call sum_allreduce(mudiff_zf_response)

      !calculate 1 - fsa(inv(mudiff_response).Q)
      do idx = 1, nresponse_mu
        mudiff_zf_response(idx,idx,:) = mudiff_zf_response(idx,idx,:) + 1.0
      enddo

      do ikx = 1,nakx
         call lu_decomposition (mudiff_zf_response(:,:,ikx),mudiff_zf_idx(:,ikx),dum2)
      end do

      deallocate (temp_mat)
    endif

    deallocate (dum1, field)

  end subroutine init_mudiff_conserve

  subroutine get_upar (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj0v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       g0 = g(:,:,ikxkyz)*spread(vpa,2,nmu)*spread(aj0v(:,ikxkyz),1,nvpa)
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_upar

  subroutine get_uperp (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vperp2
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj1v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
!       g0 = 2.0*g(:,:,ikxkyz)*spread((vperp2(1,iz,:)-0.5)*aj1v(:,ikxkyz),1,nvpa)
       g0 = g(:,:,ikxkyz)*spread((vperp2(1,iz,:)-0.5)*aj1v(:,ikxkyz),1,nvpa)
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_uperp

  subroutine get_temp (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj0v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       g0 = g(:,:,ikxkyz)*(spread(vpa**2,2,nmu)-0.5) &
          * spread(aj0v(:,ikxkyz),1,nvpa)/1.5
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_temp

  subroutine get_temp_mu (g, fld)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: nvpa, nmu, vperp2
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use gyro_averages, only: aj0v

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,:), intent (out) :: fld

    integer :: ikxkyz, iky, ikx, iz, it, is
    complex, dimension (:,:), allocatable :: g0

    allocate (g0(nvpa,nmu))

    fld = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       g0 = g(:,:,ikxkyz)*(spread(vperp2(1,iz,:),1,nvpa)-1.0) &
          * spread(aj0v(:,ikxkyz),1,nvpa)/1.5
       call integrate_vmu (g0,iz,fld(iky,ikx,iz,it,is))
    end do
    deallocate (g0)

    call sum_allreduce (fld)

  end subroutine get_temp_mu

  subroutine finish_dissipation

    implicit none

    call finish_collisions

  end subroutine finish_dissipation

  subroutine finish_collisions

    implicit none

    if (collisions_implicit) then
       call finish_vpadiff_matrix
       call finish_mudiff_matrix
       call finish_vpadiff_response
       call finish_mudiff_response
    end if

    if (collision_model == "fokker-planck") then
        call finish_nusDpa
        call finish_fp_diffmatrix
        call finish_fp_response
        call finish_deltaj
    end if

    collisions_initialized = .false.

  end subroutine finish_collisions

  subroutine finish_deltaj

        implicit none
        if (allocated(deltaj)) deallocate (deltaj)
        if (allocated(psijnorm)) deallocate (psijnorm)
        if (allocated(mwnorm)) deallocate (mwnorm)

  end subroutine finish_deltaj

  subroutine finish_fp_diffmatrix

    implicit none

    if (allocated(aa_vpa)) deallocate (aa_vpa)
    if (allocated(bb_vpa)) deallocate (bb_vpa)
    if (allocated(cc_vpa)) deallocate (cc_vpa)
    if (allocated(aa_blcs)) deallocate (aa_blcs)
    if (allocated(bb_blcs)) deallocate (bb_blcs)
    if (allocated(cc_blcs)) deallocate (cc_blcs)
    if (allocated(cdiffmat_band)) deallocate (cdiffmat_band)
    if (allocated(blockmatrix)) deallocate (blockmatrix)
    if (allocated(blockmatrix_sum)) deallocate (blockmatrix_sum)

  end subroutine finish_fp_diffmatrix

 subroutine finish_fp_response

    implicit none

    if (allocated(fp_response)) deallocate (fp_response)
    if (allocated(diff_idx)) deallocate (diff_idx)

  end subroutine finish_fp_response

  subroutine finish_vpadiff_matrix

    implicit none

    if (allocated(aa_vpa)) deallocate (aa_vpa)
    if (allocated(bb_vpa)) deallocate (bb_vpa)
    if (allocated(cc_vpa)) deallocate (cc_vpa)

  end subroutine finish_vpadiff_matrix

  subroutine finish_mudiff_matrix

    implicit none

    if (allocated(aa_mu)) deallocate (aa_mu)
    if (allocated(bb_mu)) deallocate (bb_mu)
    if (allocated(cc_mu)) deallocate (cc_mu)

  end subroutine finish_mudiff_matrix

  subroutine finish_vpadiff_response

    implicit none

    if (allocated(vpadiff_response)) deallocate (vpadiff_response)
    if (allocated(vpadiff_idx)) deallocate (vpadiff_idx)

  end subroutine finish_vpadiff_response

  subroutine finish_mudiff_response

    implicit none

    if (allocated(mudiff_response)) deallocate (mudiff_response)
    if (allocated(mudiff_idx)) deallocate (mudiff_idx)

  end subroutine finish_mudiff_response

  subroutine advance_collisions_explicit (g, phi, gke_rhs)

    use mp, only: proc0
    use job_manage, only: time_message
    use redistribute, only: scatter, gather
    use stella_time, only: code_dt
    use zgrid, only: nzgrid, ntubes
    use species, only: spec
    use run_parameters, only: fphi
    use physics_flags, only: radial_variation
    use kt_grids, only: naky, nakx, multiply_by_rho, rho_d_clamped
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: set_vpa_weights
    use stella_geometry, only: bmag, dBdrho
    use stella_layouts, only: vmu_lo, kxkyz_lo
    use stella_layouts, only: is_idx, iky_idx, ikx_idx, iz_idx
    use dist_redistribute, only: kxkyz2vmu
    use dist_fn_arrays, only: gvmu, kperp2, dkperp2dr
    use fields_arrays, only: phi_corr_QN
    use g_tofrom_h, only: g_to_h
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    integer :: is, ikxkyz, imu, iv, ivmu, ikx, iky, iz, ia, it
    logical :: conservative_wgts
    real :: tfac, kfac

    complex, dimension (:), allocatable :: mucoll
    complex, dimension (:,:,:), allocatable :: coll
    complex, dimension (:,:,:,:,:), allocatable :: tmp_vmulo

    complex, dimension (:,:,:), allocatable :: mucoll_fp
    complex, dimension (:,:,:), allocatable :: coll_fp

    complex, dimension (:,:), allocatable :: g0k, g0x

    ia = 1

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

    kfac = 0.0
    if (mu_operator)  kfac = kfac + 0.5
    if (vpa_operator) kfac = kfac + 0.5

    allocate (tmp_vmulo(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    ! want exact conservation properties for collision operator
    conservative_wgts = .true.
    call set_vpa_weights (conservative_wgts)

    if (radial_variation) then
      allocate (g0k(naky,nakx))
      allocate (g0x(naky,nakx))
      !TODO (DSO) - could perhaps operator split the profile variation pieces off the main pieces, and so
      !             this portion of the code could just treat the terms that vary in x

      if (collision_model=="dougherty") then
        ! switch from g = <f> to h = f + Z*e*phi/T * F0
        tmp_vmulo = g
        call g_to_h (tmp_vmulo, phi, fphi, phi_corr_QN)

        !handle gyroviscous term
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          do it = 1, ntubes
            do iz = -nzgrid, nzgrid
              g0k = 0.5*kfac*tmp_vmulo(:,:,iz,it,ivmu)*kperp2(:,:,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2
              gke_rhs(:,:,iz,it,ivmu) = gke_rhs(:,:,iz,it,ivmu) - code_dt*spec(is)%vnew(is)*g0k

              g0k = g0k * (dkperp2dr(:,:,ia,iz) - 2.0*dBdrho(iz)/bmag(ia,iz) - spec(is)%tprim)
              call multiply_by_rho (g0k)
              gke_rhs(:,:,iz,it,ivmu) = gke_rhs(:,:,iz,it,ivmu) - code_dt*spec(is)%vnew(is)*g0k
            enddo
          enddo
        enddo

        !handle the conservation terms
        if (momentum_conservation) call conserve_momentum_vmulo (tmp_vmulo, gke_rhs)
        if (energy_conservation) call conserve_energy_vmulo (tmp_vmulo, gke_rhs)

        !since Bessel functions do not appear under the velocity derivatives, these terms are one-point in x space
        ! and we can simply inverse Fourier transform
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          do it = 1, ntubes
            do iz = -nzgrid, nzgrid
              call transform_kx2x_unpadded (tmp_vmulo(:,:,iz,it,ivmu),g0x)
              tmp_vmulo(:,:,iz,it,ivmu) = g0x
            enddo
          enddo
        enddo

        ! remap so that (vpa,mu) local
        if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
        call scatter (kxkyz2vmu, tmp_vmulo, gvmu)
        if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')

        ! take vpa derivatives
        allocate (coll(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        allocate (mucoll(nmu))
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
           iky = iky_idx(kxkyz_lo,ikxkyz)
           ikx = ikx_idx(kxkyz_lo,ikxkyz)
           iz = iz_idx(kxkyz_lo,ikxkyz)
           is = is_idx(kxkyz_lo,ikxkyz)

           if (vpa_operator) then
             !fix the temperature term
             tfac = (spec(is)%temp/spec(is)%temp_psi0)*(1.0 - rho_d_clamped(ikx)*spec(is)%tprim)
             do imu = 1, nmu
                call vpa_differential_operator (tfac, gvmu(:,imu,ikxkyz), coll(:,imu,ikxkyz))
             end do
           else
             coll(:,:,ikxkyz) = 0.0
           end if

           if (mu_operator) then
             !fix the temperature/bmag term
             tfac = (spec(is)%temp/spec(is)%temp_psi0) &
                    * (1.0 - rho_d_clamped(ikx)*(spec(is)%tprim + dBdrho(iz)/bmag(ia,iz)))
             do iv = 1, nvpa
                call mu_differential_operator (tfac, iz, ia, gvmu(iv,:,ikxkyz), mucoll)
                coll(iv,:,ikxkyz) = coll(iv,:,ikxkyz) + mucoll
             end do
           end if
           gvmu(:,:,ikxkyz) = coll(:,:,ikxkyz)
        end do
        deallocate (coll, mucoll)

        ! remap so that (ky,kx,z,tube) local
        call gather (kxkyz2vmu, gvmu, tmp_vmulo)

        !don't forget to Fourier transform
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          do it = 1, ntubes
            do iz = -nzgrid, nzgrid
              call transform_x2kx_unpadded (tmp_vmulo(:,:,iz,it,ivmu),g0k)
              tmp_vmulo(:,:,iz,it,ivmu) = g0k
            enddo
          enddo
        enddo

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          gke_rhs(:,:,:,:,ivmu) =  gke_rhs(:,:,:,:,ivmu) + code_dt*spec(is)%vnew(is)*tmp_vmulo(:,:,:,:,ivmu)
        end do
      end if

      deallocate (g0k,g0x)
    else
      ! switch from g = <f> to h = f + Z*e*phi/T * F0
      tmp_vmulo = g
      call g_to_h (tmp_vmulo, phi, fphi)

      ! remap so that (vpa,mu) local
      if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
      call scatter (kxkyz2vmu, tmp_vmulo, gvmu)
      if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')

      ia = 1

      ! take vpa derivatives
      if (collision_model=="dougherty") then
        allocate (coll(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
        allocate (mucoll(nmu))
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
           iky = iky_idx(kxkyz_lo,ikxkyz)
           ikx = ikx_idx(kxkyz_lo,ikxkyz)
           iz = iz_idx(kxkyz_lo,ikxkyz)
           is = is_idx(kxkyz_lo,ikxkyz)
           if (vpa_operator) then
              do imu = 1, nmu
                 call vpa_differential_operator (1.0, gvmu(:,imu,ikxkyz), coll(:,imu,ikxkyz))
              end do
           else
              coll(:,:,ikxkyz) = 0.0
           end if
           if (mu_operator) then
              do iv = 1, nvpa
                 call mu_differential_operator (1.0, iz, ia, gvmu(iv,:,ikxkyz), mucoll)
                 coll(iv,:,ikxkyz) = coll(iv,:,ikxkyz) + mucoll
              end do
           end if
           if (momentum_conservation) call conserve_momentum (iky, ikx, iz, is, ikxkyz, gvmu(:,:,ikxkyz), coll(:,:,ikxkyz))
           if (energy_conservation) call conserve_energy (iz, is, ikxkyz, gvmu(:,:,ikxkyz), coll(:,:,ikxkyz))
           ! save memory by using gvmu and deallocating coll below
           ! before re-allocating tmp_vmulo
           gvmu(:,:,ikxkyz) = coll(:,:,ikxkyz) - 0.5*kfac*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*gvmu(:,:,ikxkyz)
        end do
        deallocate (coll, mucoll)

        ! remap so that (ky,kx,z,tube) local
        call gather (kxkyz2vmu, gvmu, tmp_vmulo)

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          gke_rhs(:,:,:,:,ivmu) =  gke_rhs(:,:,:,:,ivmu) + code_dt*spec(is)%vnew(is)*tmp_vmulo(:,:,:,:,ivmu)
        end do
      end if

      if (collision_model=="fokker-planck") then
        allocate (coll_fp(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); coll_fp = 0.0
        allocate (mucoll_fp(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); mucoll_fp = 0.0

        if (density_conservation) then
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                iky = iky_idx(kxkyz_lo,ikxkyz)
                ikx = ikx_idx(kxkyz_lo,ikxkyz)
                iz = iz_idx(kxkyz_lo,ikxkyz)
                is = is_idx(kxkyz_lo,ikxkyz)
                if (vpa_operator) then
                  do imu = 1, nmu
                      call vpa_differential_operator_fp_conservative (gvmu(:,:,ikxkyz), coll_fp(:,:,ikxkyz), imu, iz, is, ia)
                  end do
                end if
                if (mu_operator) then
                  do iv = 1, nvpa
                     call mu_differential_operator_fp_conservative (gvmu(:,:,ikxkyz), mucoll_fp(:,:,ikxkyz), iv, iz, is, ia, iky, ikx, cfac)
                  end do
                end if
                gvmu(:,:,ikxkyz) = coll_fp(:,:,ikxkyz) + mucoll_fp(:,:,ikxkyz)
            end do
        else
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                iky = iky_idx(kxkyz_lo,ikxkyz)
                ikx = ikx_idx(kxkyz_lo,ikxkyz)
                iz = iz_idx(kxkyz_lo,ikxkyz)
                is = is_idx(kxkyz_lo,ikxkyz)
                if (vpa_operator) then
                  do imu = 1, nmu
                      call vpa_differential_operator_fp (gvmu(:,:,ikxkyz), coll_fp(:,:,ikxkyz), imu, iz, is, ia)
                  end do
                end if
                if (mu_operator) then
                  do iv = 1, nvpa
                     call mu_differential_operator_fp (gvmu(:,:,ikxkyz), mucoll_fp(:,:,ikxkyz), iv, iz, is, ia, iky, ikx, cfac)
                  end do
                end if
                gvmu(:,:,ikxkyz) = coll_fp(:,:,ikxkyz) + mucoll_fp(:,:,ikxkyz)
            end do
        end if

        deallocate (coll_fp, mucoll_fp)

        call gather (kxkyz2vmu, gvmu, tmp_vmulo)

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          gke_rhs(:,:,:,:,ivmu) = gke_rhs(:,:,:,:,ivmu) + code_dt*tmp_vmulo(:,:,:,:,ivmu)
        end do
      end if
    endif

    deallocate (tmp_vmulo)

    ! reset to default integration wgts
    conservative_wgts = .false.
    call set_vpa_weights (conservative_wgts)

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

  end subroutine advance_collisions_explicit

  subroutine vpa_differential_operator_fp (h, Dh, imu, iz, is, ia)

      use vpamu_grids, only: nvpa, vpa, dvpa, mu, dmu, nmu, equally_spaced_mu_grid, maxwell_mu, maxwell_vpa
      use stella_geometry, only: bmag
      use constants, only: pi
      use species, only: spec

      implicit none

      complex, dimension (:,:), intent (out) :: Dh
      complex, dimension (:,:), intent (in) :: h
      integer, intent (in) :: imu, iz, ia, is
      integer :: iv
      complex :: Dhmu, Dhmu_u, Dhmu_l, dmuhp, dmuhm, dvpah, dvpahp, dvpahm
      real :: xp, xm, vpap, vpam, nupap, nupam, nuDp, nuDm, mwp, mwm, a, b, c

      iv = 1
      vpap = 0.5*(vpa(iv)+vpa(iv+1))
      xp   = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
      nuDp = spec(is)%vnew(is)*(erf(xp) - (erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
      nupap= spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
      mwp  = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
      dvpahp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa

      if (imu == 1) then
          ! second-order accurate
          !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
          !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
          !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
          !dmuhp = a*h(iv+1,imu)/mw(iv+1,imu,iz,is)+ b*h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) + c*h(iv+1,imu+2)/mw(iv+1,imu+2,iz)
          ! first-order accurate, as in implicit routine:
          dmuhp = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv+1,imu)/mw(iv+1,imu,iz,is))/dmu(imu)
      else if (imu == nmu) then
          ! second-order accurate
          !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
          !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
          !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
          !dmuhp = a*h(iv+1,imu-2)/mw(iv+1,imu-2,iz) + b*h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) + c*h(iv+1,imu)/mw(iv+1,imu,iz,is)
          ! first-order accurate, as in implicit routine:
          dmuhp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))/dmu(imu-1)
      else
          dmuhp = ((h(iv+1,imu)/mw(iv+1,imu,iz,is)-h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
                  + (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is)-h(iv+1,imu)/mw(iv+1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
      end if
      Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp)/(2*dvpa)

      iv = nvpa
      vpam = 0.5*(vpa(iv)+vpa(iv-1))
      xm   = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
      nuDm = spec(is)%vnew(is)*(erf(xm) - (erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
      nupam= spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
      mwm  = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
      dvpahm = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa

      if (imu == 1) then
          ! second-order accurate
          !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
          !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
          !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
          !dmuhm = a*h(iv-1,imu)/mw(iv-1,imu,iz,is) + b*h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) + c*h(iv-1,imu+2)/mw(iv-1,imu+2,iz)
          ! first-order accurate, as in implicit routine:
          dmuhm = (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dmu(imu)
      else if (imu == nmu) then
          ! second-order accurate
          !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
          !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
          !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
          !dmuhm = a*h(iv-1,imu-2)/mw(iv-1,imu-2,iz) + b*h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is) + c*h(iv-1,imu)/mw(iv-1,imu,iz,is)
          ! first-order accurate, as in implicit routine:
          dmuhm = (h(iv-1,imu)/mw(iv-1,imu,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dmu(imu-1)
      else
          dmuhm = ((h(iv-1,imu)/mw(iv-1,imu,iz,is)-h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
            + (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is)-h(iv-1,imu)/mw(iv-1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
      end if
      Dh(iv,imu) = (-2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)

      do iv = 2, nvpa-1
          ! AVB: interior nodes:
          ! quantities at half-grid-points:
          vpap = 0.5*(vpa(iv)+vpa(iv+1))
          vpam = 0.5*(vpa(iv)+vpa(iv-1))
          xp   = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
          xm   = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
          nuDp = spec(is)%vnew(is)*(erf(xp) - (erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
          nuDm = spec(is)%vnew(is)*(erf(xm) - (erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
          nupap= spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
          nupam= spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
          mwp  = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
          mwm  = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
          dvpahp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv  ,imu)/mw(iv  ,imu,iz,is))/dvpa
          dvpahm = (h(iv  ,imu)/mw(iv  ,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa

          if (imu == 1) then
              ! second-order accurate:
              !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
              !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
              !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
              !dmuhp = a*h(iv+1,imu)/mw(iv+1,imu,iz,is) + b*h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) + c*h(iv+1,imu+2)/mw(iv+1,imu+2,iz)
              !dmuhm = a*h(iv-1,imu)/mw(iv-1,imu,iz,is) + b*h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) + c*h(iv-1,imu+2)/mw(iv-1,imu+2,iz)
              ! or first-order accurate, as in implicit routine:
              dmuhp = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv+1,imu)/mw(iv+1,imu,iz,is))/dmu(imu)
              dmuhm = (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dmu(imu)
          else if (imu == nmu) then
              ! second-order accurate:
              !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
              !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
              !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
              !dmuhp = a*h(iv+1,imu-2)/mw(iv+1,imu-2,iz) + b*h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) + c*h(iv+1,imu)/mw(iv+1,imu,iz,is)
              !dmuhm = a*h(iv-1,imu-2)/mw(iv-1,imu-2,iz) + b*h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is) + c*h(iv-1,imu)/mw(iv-1,imu,iz,is)
              ! or first-order accurate, as in implicit routine:
              dmuhp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))/dmu(imu-1)
              dmuhm = (h(iv-1,imu)/mw(iv-1,imu,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dmu(imu-1)
          else
              dmuhp = ((h(iv+1,imu)/mw(iv+1,imu,iz,is)-h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
                + (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is)-h(iv+1,imu)/mw(iv+1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
              dmuhm = ((h(iv-1,imu)/mw(iv-1,imu,iz,is)-h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
                + (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is)-h(iv-1,imu)/mw(iv-1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
          end if
          Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp &
            + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp &
            - 2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm &
            - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)
      end do

  end subroutine vpa_differential_operator_fp

  subroutine mu_differential_operator_fp (h, Dh, iv, iz, is, ia, iky, ikx, cfac)

      use vpamu_grids, only: nmu, mu, dmu, vpa, dvpa, nvpa, maxwell_mu, maxwell_vpa, equally_spaced_mu_grid
      use stella_geometry, only: bmag
      use species, only: spec
      use dist_fn_arrays, only: kperp2
      use constants, only: pi
      use job_manage, only: timer_local, time_message
      use mp, only: proc0

      implicit none

      complex, dimension (:,:), intent (in) :: h
      complex, dimension (:,:), intent (out) :: Dh
      integer, intent (in) :: iv, iz, is, ia, iky, ikx
      real, intent (in) :: cfac
      complex ::  Dvpah, Dvpah_p, Dvpah_m, Dmuh, Dmuh_m, Dmuh_p, Dmuh1, Dmuh2
      real :: nuDp, nuDm, nupap, nupam, mup, mum, nusp, nuxp, xp, xm, mwm, mwp, a, b, c
      integer :: imu

      imu = 1
      ! vpa-differential terms:
      if (iv == 1) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa
          Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv,imu+1)/mw(iv,imu+1,iz,is))/dvpa
      else if (iv == nvpa) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa
          Dvpah_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/dvpa
      else
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/(2*dvpa)
          Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/(2*dvpa)
      end if

      ! first mu-derivative term at mu_{i}:
      ! use ghost cell at mu_{0} = 0, where mu*vpa*nux(vpa,mu)*F0 vanishes, dmu(0) = mu(1).
      Dmuh1 = ((vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,is)*mw(iv,imu,iz,is)*Dvpah)*dmu(imu)/mu(imu) &
              +(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,is)*mw(iv,imu+1,iz,is)*Dvpah_p - vpa(iv)*mu(imu)*nux(iv,imu,iz,is,is)*mw(iv,imu,iz,is)*Dvpah)*mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu))
      ! first derivative of h, at mu_{i+1/2}, and at mu_i:
      Dmuh  = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dmu(imu) ! first-order accurate, as used in implicit routine
      ! for second-order accuracy:
      !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
      !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
      !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
      !Dmuh   = a*h(iv,imu)/mw(iv,imu,iz,is) + b*h(iv,imu+1)/mw(iv,imu+1,iz,is) + c*h(iv,imu+2)/mw(iv,imu+2,iz) ! second order accurate
      Dmuh_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dmu(imu) ! second-order accurate
      ! quantities at mu_{i+1/2}:
      mup = 0.5*(mu(imu)+mu(imu+1))
      mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
      xp  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
      nuDp  = spec(is)%vnew(is)*( erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
      nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
      ! second mu-derivative term at mu_{i}:
      ! use d/dmu[...]_{1} = ([...]_{1+1/2} - [...]_{0})/(dmu_{1}/2+mu(1)), where [...]_{0} is a ghost cell at mu_{0} = 0, with [...]_{0} = 0.
      Dmuh2 = ( (2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu)/2./mu(imu) &
               +(2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp*Dmuh_p &
               - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh)*mu(imu)/(dmu(imu)/2.) )/(mu(imu)+dmu(imu)/2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2&
                *(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)

      imu = nmu
      ! vpa-differential terms:
      if (iv == 1) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa
          Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dvpa
      else if (iv == nvpa) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa
          Dvpah_m = (h(iv,imu-1)/mw(iv,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dvpa
      else
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/(2*dvpa)
          Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/(2*dvpa)
      end if

      ! first mu-derivative term at mu_{nmu}:
      Dmuh1 = (( vpa(iv)*mu(imu )*nux(iv,imu  ,iz,is,is)*mw(iv,imu,iz,is)*Dvpah &
        - vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,is)*mw(iv,imu-1,iz,is)*Dvpah_m)*dmu(imu-1)/dmu(imu-1) &
        +( -vpa(iv)*mu(imu )*nux(iv,imu  ,iz,is,is)*mw(iv,imu,iz,is)*Dvpah )*dmu(imu-1)/dmu(imu-1)) / (dmu(imu-1)+dmu(imu-1))
      ! first derivative of h, at mu_{nmu} and mu_{nmu-1/2}:
      Dmuh_m = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dmu(imu-1)
      Dmuh   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dmu(imu-1) ! first-order accurate
      ! for second-order accuracy:
      !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
      !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
      !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
      !Dmuh = a*h(iv,imu-2)/mw(iv,imu-2,iz) + b*h(iv,imu-1)/mw(iv,imu-1,iz,is) + c*h(iv,imu)/mw(iv,imu,iz,is)
      ! quantities at mu_{nmu-1/2}:
      mum = 0.5*(mu(imu)+mu(imu-1))
      mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
      xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
      nuDm  = spec(is)%vnew(is)*( erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
      nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
      ! second mu-derivative term at mu_{nmu}:
      ! use d/dmu[...]_{nmu} = ([...]_{nmu+1} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)), where [...]_{nmu+1} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1), with [...]_{nmu+1} = 0.
      Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))&
        *mu(imu))*mw(iv,imu,iz,is)*Dmuh - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m)*dmu(imu-1)/(dmu(imu-1)/2) &
               +(-2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))&
               *mu(imu))*mw(iv,imu,iz,is)*Dmuh)*dmu(imu-1)/2./dmu(imu-1) )/(dmu(imu-1)/2.+dmu(imu-1))
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2&
            *(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)

      do imu = 2, nmu-1
          ! vpa-differential terms:
          if (iv == 1) then
              ! AVB: first order accurate:
              Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa
              Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv,imu+1)/mw(iv,imu+1,iz,is))/dvpa
              Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dvpa
          else if (iv == nvpa) then
              ! AVB: first order accurate:
              Dvpah   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa
              Dvpah_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/dvpa
              Dvpah_m = (h(iv,imu-1)/mw(iv,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dvpa
          else
              Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/(2*dvpa)
              Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/(2*dvpa)
              Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/(2*dvpa)
          end if
          ! first mu-derivative of vpa-derivative term, at mu_{i}:
          Dmuh1 = ( (vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,is)*mw(iv,imu  ,iz,is)*Dvpah   &
            - vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,is)*mw(iv,imu-1,iz,is)*Dvpah_m)*dmu(imu)/dmu(imu-1) &
            +(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,is)*mw(iv,imu+1,iz,is)*Dvpah_p&
            - vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,is)*mw(iv,imu  ,iz,is)*Dvpah  )*dmu(imu-1)/dmu(imu) ) / (dmu(imu-1)+dmu(imu))
          ! first mu-derivatives of h, at mu_i, mu_{i+1/2} and mu_{i-1/2}:
          Dmuh   = ( (h(iv,imu)/mw(iv,imu,iz,is)-h(iv,imu-1)/mw(iv,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
            + (h(iv,imu+1)/mw(iv,imu+1,iz,is)-h(iv,imu)/mw(iv,imu,iz,is))*dmu(imu-1)/dmu(imu) ) / (dmu(imu-1)+dmu(imu))
          Dmuh_m = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dmu(imu-1)
          Dmuh_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dmu(imu)
          ! quantities at mu_{i+1/2} and mu_{i-1/2}:
          mup = 0.5*(mu(imu)+mu(imu+1))
          mum = 0.5*(mu(imu)+mu(imu-1))
          mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
          mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
          xp  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
          xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
          nuDp  = spec(is)%vnew(is)*( erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2) ) / (2*xp**2) )/xp**3
          nuDm  = spec(is)%vnew(is)*( erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2) ) / (2*xm**2) )/xm**3
          nupap = spec(is)%vnew(is)*2*( erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2) ) / (2*xp**2) / xp**3
          nupam = spec(is)%vnew(is)*2*( erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2) ) / (2*xm**2) / xm**3
          ! second mu-derivative term at mu_{i}:
          Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh &
            - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m )*dmu(imu)/dmu(imu-1) &
            +( 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp*Dmuh_p - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2&
            +nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu-1)/dmu(imu) )*2/(dmu(imu-1)+dmu(imu))
          ! add differential terms and gyro-diffusive term:
          Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
            + nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)
      end do

  end subroutine mu_differential_operator_fp

  subroutine vpa_differential_operator_fp_conservative (h, Dh, imu, iz, is, ia)

      use vpamu_grids, only: nvpa, vpa, dvpa, mu, dmu, nmu, equally_spaced_mu_grid, maxwell_mu, maxwell_vpa
      use stella_geometry, only: bmag
      use constants, only: pi
      use species, only: spec

      implicit none

      complex, dimension (:,:), intent (out) :: Dh
      complex, dimension (:,:), intent (in) :: h
      integer, intent (in) :: imu, iz, ia, is
      integer :: iv
      complex :: Dhmu, Dhmu_u, Dhmu_l, dmuhp, dmuhm, dvpah, dvpahp, dvpahm
      real :: xp, xm, vpap, vpam, nupap, nupam, nuDp, nuDm, mwp, mwm, a, b, c

      iv = 1
      vpap = 0.5*(vpa(iv)+vpa(iv+1))
      xp   = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
      nuDp = spec(is)%vnew(is)*(erf(xp) - (erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
      nupap= spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
      mwp  = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
      dvpahp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa

      if (imu == 1) then
          ! second-order accurate
          !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
          !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
          !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
          !dmuhp = a*h(iv+1,imu)/mw(iv+1,imu,iz,is)+ b*h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) + c*h(iv+1,imu+2)/mw(iv+1,imu+2,iz)
          ! first-order accurate, as in implicit routine:
          dmuhp = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv+1,imu)/mw(iv+1,imu,iz,is))/dmu(imu)
      else if (imu == nmu) then
          ! second-order accurate
          !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
          !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
          !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
          !dmuhp = a*h(iv+1,imu-2)/mw(iv+1,imu-2,iz) + b*h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) + c*h(iv+1,imu)/mw(iv+1,imu,iz,is)
          ! first-order accurate, as in implicit routine:
          dmuhp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))/dmu(imu-1)
      else
          dmuhp = ((h(iv+1,imu)/mw(iv+1,imu,iz,is)-h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
                  + (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is)-h(iv+1,imu)/mw(iv+1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
      end if
      Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp)/(2*dvpa)

      iv = nvpa
      vpam = 0.5*(vpa(iv)+vpa(iv-1))
      xm   = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
      nuDm = spec(is)%vnew(is)*(erf(xm) - (erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
      nupam= spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
      mwm  = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
      dvpahm = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa

      if (imu == 1) then
          ! second-order accurate
          !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
          !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
          !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
          !dmuhm = a*h(iv-1,imu)/mw(iv-1,imu,iz,is) + b*h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) + c*h(iv-1,imu+2)/mw(iv-1,imu+2,iz)
          ! first-order accurate, as in implicit routine:
          dmuhm = (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dmu(imu)
      else if (imu == nmu) then
          ! second-order accurate
          !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
          !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
          !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
          !dmuhm = a*h(iv-1,imu-2)/mw(iv-1,imu-2,iz) + b*h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is) + c*h(iv-1,imu)/mw(iv-1,imu,iz,is)
          ! first-order accurate, as in implicit routine:
          dmuhm = (h(iv-1,imu)/mw(iv-1,imu,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dmu(imu-1)
      else
          dmuhm = ((h(iv-1,imu)/mw(iv-1,imu,iz,is)-h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
            + (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is)-h(iv-1,imu)/mw(iv-1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
      end if
      Dh(iv,imu) = (-2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)

      do iv = 2, nvpa-1
          ! AVB: interior nodes:
          ! quantities at half-grid-points:
          vpap = 0.5*(vpa(iv)+vpa(iv+1))
          vpam = 0.5*(vpa(iv)+vpa(iv-1))
          xp  = sqrt(vpap**2 + 2*bmag(ia,iz)*mu(imu))
          xm  = sqrt(vpam**2 + 2*bmag(ia,iz)*mu(imu))
          nuDp = spec(is)%vnew(is)*(erf(xp) - (erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
          nuDm = spec(is)%vnew(is)*(erf(xm) - (erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
          nupap= spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
          nupam= spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
          mwp = exp(-vpap**2)*maxwell_mu(1,iz,imu,is)
          mwm = exp(-vpam**2)*maxwell_mu(1,iz,imu,is)
          dvpahp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv  ,imu)/mw(iv  ,imu,iz,is))/dvpa
          dvpahm = (h(iv  ,imu)/mw(iv  ,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa

          if (imu == 1) then
              ! second-order accurate:
              !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
              !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
              !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
              !dmuhp = a*h(iv+1,imu)/mw(iv+1,imu,iz,is) + b*h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) + c*h(iv+1,imu+2)/mw(iv+1,imu+2,iz)
              !dmuhm = a*h(iv-1,imu)/mw(iv-1,imu,iz,is) + b*h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) + c*h(iv-1,imu+2)/mw(iv-1,imu+2,iz)
              ! or first-order accurate, as in implicit routine:
              dmuhp = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv+1,imu)/mw(iv+1,imu,iz,is))/dmu(imu)
              dmuhm = (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dmu(imu)
          else if (imu == nmu) then
              ! second-order accurate:
              !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
              !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
              !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
              !dmuhp = a*h(iv+1,imu-2)/mw(iv+1,imu-2,iz) + b*h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) + c*h(iv+1,imu)/mw(iv+1,imu,iz,is)
              !dmuhm = a*h(iv-1,imu-2)/mw(iv-1,imu-2,iz) + b*h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is) + c*h(iv-1,imu)/mw(iv-1,imu,iz,is)
              ! or first-order accurate, as in implicit routine:
              dmuhp = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))/dmu(imu-1)
              dmuhm = (h(iv-1,imu)/mw(iv-1,imu,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dmu(imu-1)
          else
              dmuhp = ((h(iv+1,imu)/mw(iv+1,imu,iz,is)-h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
                + (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is)-h(iv+1,imu)/mw(iv+1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
              dmuhm = ((h(iv-1,imu)/mw(iv-1,imu,iz,is)-h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))*dmu(imu)/dmu(imu-1) &
                + (h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is)-h(iv-1,imu)/mw(iv-1,imu,iz,is))*dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu))
          end if

          if (iv==2) then
              ! assume vpa*mu*nux*F0*dh/dmu vanishes at iv=1, to ensure density conservation
              Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp &
                + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp &
                - 2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm &
                - 0*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)
          else if (iv==nvpa-1) then
              ! assume vpa*mu*nux*F0*dh/dmu vanishes at iv=nvpa, to ensure density conservation
              Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp &
                + 0*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp &
                - 2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm &
                - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)
          else
              Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp &
                + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp &
                - 2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm &
                - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)
          end if

      end do

  end subroutine vpa_differential_operator_fp_conservative

  subroutine mu_differential_operator_fp_conservative (h, Dh, iv, iz, is, ia, iky, ikx, cfac)

      use vpamu_grids, only: nmu, mu, dmu, vpa, dvpa, nvpa, maxwell_mu, maxwell_vpa, equally_spaced_mu_grid
      use stella_geometry, only: bmag
      use species, only: spec
      use dist_fn_arrays, only: kperp2
      use constants, only: pi
      use job_manage, only: timer_local, time_message
      use mp, only: proc0

      implicit none

      complex, dimension (:,:), intent (in) :: h
      complex, dimension (:,:), intent (out) :: Dh
      integer, intent (in) :: iv, iz, is, ia, iky, ikx
      real, intent (in) :: cfac
      complex ::  Dvpah, Dvpah_p, Dvpah_m, Dmuh, Dmuh_m, Dmuh_p, Dmuh1, Dmuh2
      real :: nuDp, nuDm, nupap, nupam, mup, mum, nusp, nuxp, xp, xm, mwm, mwp, a, b, c
      integer :: imu

      imu = 1
      ! vpa-differential terms:
      if (iv == 1) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa
          Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv,imu+1)/mw(iv,imu+1,iz,is))/dvpa
      else if (iv == nvpa) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa
          Dvpah_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/dvpa
      else
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/(2*dvpa)
          Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/(2*dvpa)
      end if
      ! first mu-derivative term at mu_{i}:
      ! use ghost cell at mu_{0} = 0, where mu*vpa*nux(vpa,mu)*F0 vanishes, dmu(0) = mu(1)
      ! to ensure conservation of density we approximate as follows:
      Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz,is,is)*mw(iv,imu,iz,is)*Dvpah + vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,is)*mw(iv,imu+1,iz,is)*Dvpah_p) / (2.*dmu(imu))
      !!Dmuh1 = ((vpa(iv)*mu(imu  )*nux(iv,imu  ,iz)*mw(iv,imu  ,iz,is)*Dvpah)*dmu(imu)/mu(imu) &
      !!      +(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz,is)*Dvpah_p - vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah)*mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu))

      ! first derivative of h, at mu_{i+1/2}, and at mu_i:
      Dmuh  = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dmu(imu) ! first-order accurate, as used in implicit routine
      ! for second-order accuracy:
      !a = -(2.*dmu(1) + dmu(2))/(dmu(1)*(dmu(1)+dmu(2)))
      !b = (dmu(1)+dmu(2))/(dmu(1)*dmu(2))
      !c = -dmu(1)/(dmu(2)*(dmu(1)+dmu(2)))
      !Dmuh   = a*h(iv,imu)/mw(iv,imu,iz,is) + b*h(iv,imu+1)/mw(iv,imu+1,iz,is) + c*h(iv,imu+2)/mw(iv,imu+2,iz) ! second order accurate
      Dmuh_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dmu(imu) ! second-order accurate
      ! quantities at mu_{i+1/2}:
      mup = 0.5*(mu(imu)+mu(imu+1))
      mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
      xp  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
      nuDp  = spec(is)%vnew(is)*( erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) )/xp**3
      nupap = spec(is)%vnew(is)*2*(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2)) / (2*xp**2) / xp**3
      ! second mu-derivative term at mu_{i}:
      ! use d/dmu[...]_{1} = ([...]_{1+1/2} - [...]_{0})/(dmu_{1}/2+mu(1)), where [...]_{0} is a ghost cell at mu_{0} = 0, with [...]_{0} = 0.
      Dmuh2 = ( (2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu)/2./mu(imu) &
        +(2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp*Dmuh_p &
        - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh)*mu(imu)/(dmu(imu)/2.) )/(mu(imu)+dmu(imu)/2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2&
        *(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) + nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)

      do imu = 2, nmu-1
          ! vpa-differential terms:
          if (iv == 1) then
              ! AVB: first order accurate:
              Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa
              Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv,imu+1)/mw(iv,imu+1,iz,is))/dvpa
              Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dvpa
          else if (iv == nvpa) then
              ! AVB: first order accurate:
              Dvpah   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa
              Dvpah_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/dvpa
              Dvpah_m = (h(iv,imu-1)/mw(iv,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dvpa
          else
              Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/(2*dvpa)
              Dvpah_p = (h(iv+1,imu+1)/mw(iv+1,imu+1,iz,is) - h(iv-1,imu+1)/mw(iv-1,imu+1,iz,is))/(2*dvpa)
              Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/(2*dvpa)
          end if

          ! first mu-derivative of vpa-derivative term, at mu_{i}:
          if (imu==nmu-1) then
              ! to ensure conservation of density we assume mu*nux*F0*d(h/F0)/dvpa vanishes at nmu
              Dmuh1 = ( -vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,is)*mw(iv,imu-1,iz,is)*Dvpah_m)  / (dmu(imu-1)+dmu(imu-1))
          else
              Dmuh1 = ( -vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,is)*mw(iv,imu-1,iz,is)*Dvpah_m &
                        +vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,is)*mw(iv,imu+1,iz,is)*Dvpah_p ) / (dmu(imu-1)+dmu(imu-1))
          end if

          ! first mu-derivatives of h, at mu_i, mu_{i+1/2} and mu_{i-1/2}:
          Dmuh   = ( (h(iv,imu)/mw(iv,imu,iz,is)-h(iv,imu-1)/mw(iv,imu-1,iz,is))*dmu(imu)/dmu(imu-1) + (h(iv,imu+1)/mw(iv,imu+1,iz,is)-h(iv,imu)/mw(iv,imu,iz,is))&
            *dmu(imu-1)/dmu(imu) ) / (dmu(imu-1)+dmu(imu))
          Dmuh_m = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dmu(imu-1)
          Dmuh_p = (h(iv,imu+1)/mw(iv,imu+1,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dmu(imu)
          ! quantities at mu_{i+1/2} and mu_{i-1/2}:
          mup = 0.5*(mu(imu)+mu(imu+1))
          mum = 0.5*(mu(imu)+mu(imu-1))
          mwp = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mup)
          mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
          xp = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mup)
          xm = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
          nuDp = spec(is)%vnew(is)*( erf(xp)-(erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2) ) / (2*xp**2) )/xp**3
          nuDm = spec(is)%vnew(is)*( erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2) ) / (2*xm**2) )/xm**3
          nupap = spec(is)%vnew(is)*2*( erf(xp)-xp*(2/sqrt(pi))*exp(-xp**2) ) / (2*xp**2) / xp**3
          nupam = spec(is)%vnew(is)*2*( erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2) ) / (2*xm**2) / xm**3
          ! second mu-derivative term at mu_{i}:
          Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh &
            - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m )*dmu(imu)/dmu(imu-1) &
            +( 2*(nupap*mup**2+nuDp*vpa(iv)**2/(2*bmag(ia,iz))*mup)*mwp*Dmuh_p &
            - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu-1)/dmu(imu) )*2/(dmu(imu-1)+dmu(imu))
          ! add differential terms and gyro-diffusive term:
          Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
            + nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)
      end do


      imu = nmu
      ! vpa-differential terms:
      if (iv == 1) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv,imu)/mw(iv,imu,iz,is))/dvpa
          Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dvpa
      else if (iv == nvpa) then
          ! AVB: first order accurate:
          Dvpah   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/dvpa
          Dvpah_m = (h(iv,imu-1)/mw(iv,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/dvpa
      else
          Dvpah   = (h(iv+1,imu)/mw(iv+1,imu,iz,is) - h(iv-1,imu)/mw(iv-1,imu,iz,is))/(2*dvpa)
          Dvpah_m = (h(iv+1,imu-1)/mw(iv+1,imu-1,iz,is) - h(iv-1,imu-1)/mw(iv-1,imu-1,iz,is))/(2*dvpa)
      end if
      ! first mu-derivative term at mu_{nmu}:
      ! to ensure conservation of density, assume that term is zero beyond nmu
      Dmuh1 = ( -vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,is)*mw(iv,imu-1,iz,is)*Dvpah_m ) / (dmu(imu-1)+dmu(imu-1))

      ! first derivative of h, at mu_{nmu} and mu_{nmu-1/2}:
      Dmuh_m = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dmu(imu-1)
      Dmuh   = (h(iv,imu)/mw(iv,imu,iz,is) - h(iv,imu-1)/mw(iv,imu-1,iz,is))/dmu(imu-1) ! first-order accurate
      ! for second-order accuracy:
      !a = dmu(nmu-1)/(dmu(nmu-2)*(dmu(nmu-2)+dmu(nmu-1)))
      !b = -(dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-2)*dmu(nmu-1))
      !c = (2.*dmu(nmu-1)+dmu(nmu-2))/(dmu(nmu-1)*(dmu(nmu-1)+dmu(nmu-2)))
      !Dmuh = a*h(iv,imu-2)/mw(iv,imu-2,iz) + b*h(iv,imu-1)/mw(iv,imu-1,iz,is) + c*h(iv,imu)/mw(iv,imu,iz,is)
      ! quantities at mu_{nmu-1/2}:
      mum = 0.5*(mu(imu)+mu(imu-1))
      mwm = maxwell_vpa(iv,is)*exp(-2*bmag(ia,iz)*mum)
      xm  = sqrt(vpa(iv)**2 + 2*bmag(ia,iz)*mum)
      nuDm  = spec(is)%vnew(is)*( erf(xm)-(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) )/xm**3
      nupam = spec(is)%vnew(is)*2*(erf(xm)-xm*(2/sqrt(pi))*exp(-xm**2)) / (2*xm**2) / xm**3
      ! second mu-derivative term at mu_{nmu}:
      ! to ensure density conservation
      ! use d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2), where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
      Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh &
            - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m)*dmu(imu-1)/dmu(imu-1) &
            +(-2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))&
            *mw(iv,imu,iz,is)*Dmuh)*dmu(imu-1)/dmu(imu-1) )/(dmu(imu-1)/2.+dmu(imu-1)/2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - cfac*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
            + nuD(iv,imu,iz,is,is)*(vpa(iv)**2 + bmag(ia,iz)*mu(imu)))*h(iv,imu)

  end subroutine mu_differential_operator_fp_conservative

  subroutine vpa_differential_operator (tfac, h, Dh)

    use vpamu_grids, only: nvpa, vpa, dvpa

    implicit none

    real, intent (in) :: tfac
    complex, dimension (:), intent (in) :: h
    complex, dimension (:), intent (out) :: Dh

    integer :: iv

    ! use h = 0 at ghost cells beyond +/- vpa_max
    iv = 1
    Dh(iv) = (0.5*h(iv+1)*(tfac/dvpa+vpa(iv+1))-tfac*h(iv)/dvpa)/dvpa
    iv = nvpa
    Dh(iv) = (-tfac*h(iv)/dvpa+0.5*h(iv-1)*(tfac/dvpa-vpa(iv-1)))/dvpa
    do iv = 2, nvpa-1
       Dh(iv) = (0.5*h(iv+1)*(tfac/dvpa+vpa(iv+1))-tfac*h(iv)/dvpa+0.5*h(iv-1)*(tfac/dvpa-vpa(iv-1)))/dvpa
    end do

  end subroutine vpa_differential_operator

  subroutine mu_differential_operator (tfac, iz, ia, h, Dh)

    use vpamu_grids, only: nmu, mu, dmu
    use vpamu_grids, only: mu_cell, dmu_cell, wgts_mu_bare
    use vpamu_grids, only: equally_spaced_mu_grid
    use stella_geometry, only: bmag

    implicit none

    real, intent (in) :: tfac
    integer, intent (in) :: iz, ia
    complex, dimension (:), intent (in) :: h
    complex, dimension (:), intent (out) :: Dh

    integer :: imu
    real :: mm, m0, mp

    ! the following finite difference method is explained in init_mudiff_matrix

    imu = 1
    m0 = dmu_cell(imu)*(1.0 - tfac/(dmu(imu)*bmag(ia,iz)))/wgts_mu_bare(imu)
    mp = mu_cell(imu)*(tfac/(bmag(ia,iz)*dmu(imu)) + 1.0)/wgts_mu_bare(imu)
    Dh(imu) = m0*h(imu) + mp*h(imu+1)

    imu = nmu
    mm =  mu_cell(imu-1)*(tfac/(bmag(ia,iz)*dmu(imu-1)) - 1.0)/wgts_mu_bare(imu)
    m0 = -mu_cell(imu-1)*(tfac/(dmu(imu-1)*bmag(ia,iz)) + 1.0)/wgts_mu_bare(imu)
    Dh(imu) = mm*h(imu-1) + m0*h(imu)

    do imu = 2, nmu-1
      mm =   mu_cell(imu-1)*(tfac/(bmag(ia,iz)*dmu(imu-1)) - 1.0)/wgts_mu_bare(imu)
      m0 = -(mu_cell(imu)/dmu(imu) + mu_cell(imu-1)/dmu(imu-1))*tfac/(wgts_mu_bare(imu)*bmag(ia,iz)) + dmu_cell(imu)/wgts_mu_bare(imu)
      mp =   mu_cell(imu  )*(tfac/(bmag(ia,iz)*dmu(imu  )) + 1.0)/wgts_mu_bare(imu)
      Dh(imu) = mm*h(imu-1) + m0*h(imu) + mp*h(imu+1) 
    enddo

  end subroutine mu_differential_operator

  subroutine conserve_momentum (iky, ikx, iz, is, ikxkyz, h, Ch)

    use species, only: spec
    use stella_geometry, only: bmag
    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
!   use vpamu_grids, only: int_vpa2
    use dist_fn_arrays, only: kperp2
    use gyro_averages, only: aj0v, aj1v

    implicit none

    integer, intent (in) :: iky, ikx, iz, is, ikxkyz
    complex, dimension (:,:), intent (in) :: h
    complex, dimension (:,:), intent (in out) :: Ch

    complex, dimension (:,:), allocatable :: u_fac
    complex :: integral
    integer :: ia

    real :: norm

    allocate (u_fac(nvpa,nmu))

    ia = 1

!   norm = 1.0/int_vpa2(ia,iz,is)
    norm = 2.0

    if (vpa_operator) then 
      u_fac = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)
      call integrate_vmu (u_fac*h,iz,integral)

      Ch = Ch + norm*u_fac*integral*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)
    endif

    if (mu_operator) then 
      u_fac = spread(vperp2(ia,iz,:)*aj1v(:,ikxkyz),1,nvpa)*sqrt(kperp2(iky,ikx,ia,iz))*spec(is)%smz/bmag(ia,iz)
      call integrate_vmu (u_fac*h,iz,integral)

      Ch = Ch + norm*u_fac*integral*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)
    endif

    deallocate (u_fac)

  end subroutine conserve_momentum

  subroutine conserve_energy (iz, is, ikxkyz, h, Ch)

    use vpamu_grids, only: integrate_vmu
    use vpamu_grids, only: vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
!   use vpamu_grids, only: int_unit, int_vpa2, int_vperp2, int_vfrth
    use gyro_averages, only: aj0v

    implicit none

    integer, intent (in) :: iz, is, ikxkyz
    complex, dimension (:,:), intent (in) :: h
    complex, dimension (:,:), intent (in out) :: Ch

    complex, dimension (:,:), allocatable :: T_fac
    complex :: integral
    real :: norm
    
    integer :: ia

    allocate (T_fac(nvpa,nmu))

    ia = 1
    T_fac = 0.0
    if (vpa_operator) T_fac = spread(aj0v(:,ikxkyz),1,nvpa)*(spread(vpa**2,2,nmu))-0.5
    if (mu_operator)  T_fac = T_fac + spread(vperp2(ia,iz,:),1,nvpa)-1.0
    call integrate_vmu (T_fac*h,iz,integral)

!   norm = 1.0/(int_vfrth(ia,iz,is) - (int_vperp2(ia,iz,is) + int_vpa2(ia,iz,is))**2/int_unit(ia,iz,is))
    norm = 2.0/3.0

    Ch = Ch + 2.0*norm*T_fac*integral*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)

    deallocate (T_fac)

  end subroutine conserve_energy

  subroutine conserve_momentum_vmulo (h, gke_rhs)

    use mp, only: sum_allreduce
    use stella_time, only: code_dt
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, iv_idx, is_idx
    use species, only: spec
    use physics_flags, only: radial_variation
    use stella_geometry, only: bmag, dBdrho
    use kt_grids, only: nakx, naky, multiply_by_rho
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: integrate_species, mu, vpa, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: h
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    complex, dimension (:,:), allocatable :: g0k, g1k
    complex, dimension (:,:,:), allocatable :: gyro_g
    complex, dimension (:,:,:,:), allocatable :: field1, field2
    real :: prefac, energy
    integer :: it, iz, ivmu, imu, iv, ia, is

    ia = 1

    allocate (g0k(naky,nakx))
    allocate (g1k(naky,nakx))
    allocate (gyro_g(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    allocate (field1(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (field2(naky,nakx,-nzgrid:nzgrid,ntubes))

    !component from vpa
    do it = 1, ntubes
      do iz = -nzgrid, nzgrid
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          iv = iv_idx(vmu_lo,ivmu)

          call gyro_average (h(:,:,iz,it,ivmu), iz, ivmu, gyro_g(:,:,ivmu))
          gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu)*vpa(iv)
          g0k = 0.0
          if(radial_variation) then
            g0k(:,:) = gyro_g(:,:,ivmu) &
                * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 &
                * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                + dBdrho(iz)/bmag(ia,iz))

            call multiply_by_rho(g0k)
          endif
          gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu) + g0k

        end do
        call integrate_species (gyro_g, iz, spec%dens_psi0*spec%temp_psi0, field1(:,:,iz,it), reduce_in=.false.)
      end do
    end do
    call sum_allreduce(field1)


    !component from vperp
    do it = 1, ntubes
      do iz = -nzgrid, nzgrid
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          iv = iv_idx(vmu_lo,ivmu)

          !component from vpa
          call gyro_average_j1 (h(:,:,iz,it,ivmu), iz, ivmu, gyro_g(:,:,ivmu))
          gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu)*vperp2(ia,iz,imu)*sqrt(kperp2(:,:,ia,iz))*spec(is)%smz_psi0/bmag(ia,iz)
          g0k = 0.0
          if(radial_variation) then
            g0k = gyro_g(:,:,ivmu)*(dBdrho(iz)/bmag(ia,iz) + 0.5*dkperp2dr(:,:,ia,iz)) &
                  + h(:,:,iz,it,ivmu)*(0.5*aj0x(:,:,iz,ivmu) - aj1x(:,:,iz,ivmu)) &
                     *(dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz))*vperp2(ia,iz,imu) &
                     *sqrt(kperp2(:,:,ia,iz))*spec(is)%smz_psi0/bmag(ia,iz)

            call multiply_by_rho(g0k)
          endif
          gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu) + g0k

        end do
        call integrate_species (gyro_g, iz, spec%dens_psi0*spec%temp_psi0, field2(:,:,iz,it), reduce_in=.false.)

      end do
    end do
    call sum_allreduce(field2)
    deallocate (gyro_g)


    do it = 1, ntubes
      do iz = -nzgrid, nzgrid
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          iv = iv_idx(vmu_lo,ivmu)

          prefac = 2.0/(spec(is)%dens*spec(is)%temp)*code_dt*spec(is)%vnew(is) &
                    * maxwell_mu(ia,iz,imu,is)*maxwell_vpa(iv,is)*maxwell_fac(is)

          g0k = aj0x(:,:,iz,ivmu)*vpa(iv)*field1(:,:,iz,it) &
              + aj1x(:,:,iz,ivmu)*vperp2(ia,iz,imu)*field2(:,:,iz,it) &
                *spec(is)%smz_psi0*sqrt(kperp2(:,:,ia,iz))/bmag(ia,iz)


          gke_rhs(:,:,iz,it,ivmu) = gke_rhs(:,:,iz,it,ivmu) + prefac*g0k

          if(radial_variation) then
            energy = (vpa(iv)**2 + vperp2(ia,iz,imu))*(spec(is)%temp_psi0/spec(is)%temp)

            g1k = field1(:,:,iz,it)*vpa(iv)*(-0.5*aj1x(:,:,iz,ivmu)*(spec(is)%smz)**2 &
                   * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                   * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz))) &
                + field2(:,:,iz,it)*spec(is)%smz_psi0*vperp2(ia,iz,imu)/bmag(ia,iz)*sqrt(kperp2(:,:,ia,iz)) &
                   * (0.5*aj1x(:,:,iz,ivmu)*dkperp2dr(:,:,ia,iz) + (0.5*aj0x(:,:,iz,ivmu) - aj1x(:,:,iz,ivmu))  &
                       *(dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)))
            g1k = g1k + g0k*(spec(is)%tprim*(energy-2.5) + 2*mu(imu)*dBdrho(iz))

            call multiply_by_rho(g1k)

            gke_rhs(:,:,iz,it,ivmu) = gke_rhs(:,:,iz,it,ivmu) + prefac*g1k
          endif
        enddo
      end do
    end do

    deallocate (g0k,g1k)
    deallocate (field1,field2)

  end subroutine conserve_momentum_vmulo

  subroutine conserve_energy_vmulo (h, gke_rhs)

    use mp, only: sum_allreduce
    use stella_time, only: code_dt
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, iv_idx, is_idx
    use species, only: spec
    use physics_flags, only: radial_variation
    use stella_geometry, only: bmag, dBdrho
    use kt_grids, only: nakx, naky, multiply_by_rho
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: integrate_species
    use vpamu_grids, only: mu, vpa, nvpa, nmu, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: h
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    complex, dimension (:,:), allocatable :: g0k, g1k
    complex, dimension (:,:,:), allocatable :: gyro_g
    complex, dimension (:,:,:,:), allocatable :: field
    complex :: integral
    real :: prefac, energy
    integer :: it, iz, ivmu, imu, iv, ia, is

    ia = 1

    allocate (g0k(naky,nakx))
    allocate (g1k(naky,nakx))
    allocate (gyro_g(naky,nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes))

    !component from vpa
    do it = 1, ntubes
      do iz = -nzgrid, nzgrid
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          iv = iv_idx(vmu_lo,ivmu)

          call gyro_average (h(:,:,iz,it,ivmu), iz, ivmu, gyro_g(:,:,ivmu))
          gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu)*(vpa(iv)**2 + vperp2(ia,iz,imu)-1.5*spec(is)%temp/spec(is)%temp_psi0)
          g0k = 0.0
          if(radial_variation) then
            g0k(:,:) = gyro_g(:,:,ivmu) &
                        * (-0.5*aj1x(:,:,iz,ivmu)/aj0x(:,:,iz,ivmu)*(spec(is)%smz)**2 &
                        * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                        * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                           + dBdrho(iz)/bmag(ia,iz)) &
                     + h(:,:,iz,it,ivmu)*aj0x(:,:,iz,ivmu)*(vperp2(ia,iz,imu)*dBdrho(iz)/bmag(ia,iz) &
                                                            + 1.5*spec(is)%tprim)
            call multiply_by_rho(g0k)
          endif
          gyro_g(:,:,ivmu) = gyro_g(:,:,ivmu) + g0k

        end do
        call integrate_species (gyro_g, iz, spec%dens_psi0*spec%temp_psi0**2, field(:,:,iz,it), reduce_in=.false.)
      end do
    end do
    call sum_allreduce(field)

    deallocate (gyro_g)

    do it = 1, ntubes
      do iz = -nzgrid, nzgrid
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo,ivmu)
          imu = imu_idx(vmu_lo,ivmu)
          iv = iv_idx(vmu_lo,ivmu)

          prefac = 4.0/(3.0*spec(is)%dens*spec(is)%temp**2)*code_dt*spec(is)%vnew(is) &
                    * maxwell_mu(ia,iz,imu,is)*maxwell_vpa(iv,is)*maxwell_fac(is)

          g0k = aj0x(:,:,iz,ivmu)*field(:,:,iz,it) &
                  *(vpa(iv)**2 + vperp2(ia,iz,imu) - 1.5*spec(is)%temp/spec(is)%temp_psi0)

          gke_rhs(:,:,iz,it,ivmu) = gke_rhs(:,:,iz,it,ivmu) + prefac*g0k

          if(radial_variation) then
            energy = (vpa(iv)**2 + vperp2(ia,iz,imu))*(spec(is)%temp_psi0/spec(is)%temp)

            g1k = field(:,:,iz,it)*(energy - 1.5)*(-0.5*aj1x(:,:,iz,ivmu)*(spec(is)%smz)**2 &
                   * (kperp2(:,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                   * (dkperp2dr(:,:,ia,iz) - dBdrho(iz)/bmag(ia,iz))) &
                + field(:,:,iz,it)*aj0x(:,:,iz,ivmu)*(vperp2(ia,iz,imu)*dBdrho(iz)/bmag(ia,iz) + 1.5*spec(is)%tprim)


            g1k = g1k + g0k*(spec(is)%tprim*(energy-3.5) + 2*mu(imu)*dBdrho(iz))

            call multiply_by_rho(g1k)

            gke_rhs(:,:,iz,it,ivmu) = gke_rhs(:,:,iz,it,ivmu) + prefac*g1k
          endif
        enddo
      end do
    end do

    deallocate (g0k,g1k)
    deallocate (field)

  end subroutine conserve_energy_vmulo

  subroutine advance_collisions_implicit (mirror_implicit, phi, apar, g)

    use mp, only: proc0
    use redistribute, only: gather, scatter
    use dist_redistribute, only: kxkyz2vmu
    use job_manage, only: time_message
    use zgrid, only: nzgrid
    use vpamu_grids, only: set_vpa_weights
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: gvmu

    implicit none

    logical, intent (in) :: mirror_implicit
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    logical :: conservative_wgts

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

    conservative_wgts = .true.
    call set_vpa_weights (conservative_wgts)

    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
    call scatter (kxkyz2vmu, g, gvmu)
    if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')

    if (collision_model == "dougherty") then
        if (vpa_operator) call advance_vpadiff_implicit (phi, apar, gvmu)
        if (mu_operator) call advance_mudiff_implicit (phi, apar, gvmu)
    end if
    if (collision_model == "fokker-planck") then

        if (density_conservation) then
            conservative_wgts = .true.
            call set_vpa_weights (conservative_wgts)
        end if
        if (exact_conservation_tp) then
            conservative_wgts = .false.
            call set_vpa_weights (conservative_wgts)
        end if

        call advance_implicit_fp (phi, apar, gvmu)
    end if

    if (.not.mirror_implicit) then
       ! then take the results and remap again so ky,kx,z local.
       if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
       call gather (kxkyz2vmu, gvmu, g)
       if (proc0) call time_message(.false.,time_collisions(:,2),' coll_redist')
    end if

    conservative_wgts = .false.
    call set_vpa_weights (conservative_wgts)

    if (proc0) call time_message(.false.,time_collisions(:,1),' collisions')

  end subroutine advance_collisions_implicit

  subroutine advance_implicit_fp (phi, apar, g)

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_back_substitution
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use species, only: nspec, spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: mu, nmu, nvpa, integrate_vmu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
    use g_tofrom_h, only: g_to_h
    use gyro_averages, only: aj0v, aj1v
    use fields, only: get_fields, fields_updated
    use constants, only: pi
    use dist_fn_arrays, only: kperp2
    use stella_geometry, only: bmag
    use stella_time, only: code_dt

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g

    real, dimension (:,:), allocatable :: tmp
    complex, dimension (:,:,:,:,:), allocatable :: flds
    complex, dimension (:,:,:), allocatable :: g_in
    complex, dimension (:,:), allocatable :: gvmutr
    complex, dimension (:), allocatable :: ghrs
    complex, dimension (:,:,:,:,:), allocatable :: field

    complex, dimension (:,:), allocatable :: g0spitzer


    integer :: ikxkyz, iky, ikx, iz, is, imu, iv, it, ia
    integer :: idx, idx1, ij, il, im, jj, ll, mm, il1, im1, ij1, ll1 , mm1, jj1, isa, isb
    real :: clm, j1argnomu

    real :: spitzer_i1, spitzer_i2, applied_Epar, gradpar_lnp0, gradpar_lnT0
    complex :: spitzer_int

    ! store input g for use later, as g will be overwritten below
    allocate (g_in(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    allocate (g0spitzer(nvpa,nmu))

    ia = 1

    if (spitzer_problem) then
        ! to solve the Spitzer problem, we add a source term associated with a constant
        ! applied electric field, and pressure and temperature gradients here
        ! all other non-collisional terms are disabled, and Delta t --> \infty

        applied_Epar = 0.01
        gradpar_lnp0 = 0
        gradpar_lnT0 = 0.01

        spitzer_i1 = (applied_Epar - gradpar_lnp0)*i1fac
        spitzer_i2 = gradpar_lnT0*i2fac

        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
           iky = iky_idx(kxkyz_lo,ikxkyz)
           ikx = ikx_idx(kxkyz_lo,ikxkyz)
           is = is_idx(kxkyz_lo,ikxkyz)
           iz = iz_idx(kxkyz_lo,ikxkyz)
           if (is==1) cycle
           g(:,:,ikxkyz) = g(:,:,ikxkyz) - 1./sqrt(spec(is)%mass)*code_dt*( spread(vpa,2,nmu)*spitzer_i1 &
            + (spread(vpa,2,nmu)*velvpamu(:,:,iz)**2 - 5./2.*spread(vpa,2,nmu))*spitzer_i2 )*mw(:,:,iz,is)
        end do
    end if

    g_in = g

    allocate (gvmutr(nvpa,nmu))
    allocate (ghrs(nmu*nvpa))
    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes,nspec))
    allocate (flds(naky,nakx,-nzgrid:nzgrid,ntubes,nresponse))

    ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
    ! invert above equation to get h_inh^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)

       do iv = 1, nvpa
           ghrs(nmu*(iv-1)+1 : nmu*iv) = g(iv,:,ikxkyz)
       end do

       call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,iky,ikx,iz,is), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,is), ghrs, nvpa*nmu, info)

       do iv = 1, nvpa
           g(iv,:,ikxkyz) = ghrs(nmu*(iv-1)+1 : nmu*iv)
       end do
   end do

    ! obtain phi^{n+1} and conservation terms using response matrix approach

    ! first get phi_inh^{n+1}
    if (advfield_coll) then
        call get_fields (g, phi, apar, dist='h')
        flds(:,:,:,:,1) = phi
    end if

    ! next get the psi^{s1s2,jlm}_inh^{n+1}
    ! note g contains h_s1_inh, h_s2_inh, ... ie all species
    if (fieldpart) then
        ! layout of field(,,,,:) is phi; jlm0 psi_aa, jlm0 psi_ab, jlm0 psi_ba,  jlm0 psi_bb; jlm1 psi_aa, etc, because we want species to be contiguous
        do idx1 = 1, (jmax+1)*(lmax+1)**2*nspec

            isa = 1 + int((idx1-1)/((jmax+1)*(lmax+1)**2)) !1 + mod(idx1-1,(jmax+1)*(lmax+1)**2)
            ij = 1 + mod(1+mod(idx1-1,(jmax+1)*(lmax+1)**2)-1,jmax+1)
            il = 1 + int(sqrt(1.0*(1+mod(idx1-1,(jmax+1)*(lmax+1)**2) - ij)/(jmax+1)))
            im = (1+mod(idx1-1,(jmax+1)*(lmax+1)**2) - ij)/(jmax+1) - (il-1)**2 + 1
            ll = il-1
            mm = -ll + im-1
            jj = ij-1

            if (density_conservation_tp.and.(jj==0).and.(ll==0)) then
                ! get the density produced by the combined test particle operator C_{test,isa} = C_{test,isa,isa} + C_{test,isa,isb} + ...
                ! this is stored in field(:,:,:,:,isa), all other entries are zero
                field = 0.
                call get_testpart_density (isa, 0, g, field)
            else
                ! get psi^{isa[is1....isN],jlm}_inh^{n+1}
                call get_psi(g, field, isa, 0, ll, mm, jj)
            end if

            ! add to rhs vector
            flds(:,:,:,:,2 + (idx1-1)*nspec : 1 + idx1*nspec) = field(:,:,:,:,:)

        end do
    end if

    ! AVB: obtain phi^{n+1} and psijlm^{n+1} from response matrix
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ! all is indices inside ikxkyz super-index have same info
       ! no need to compute multiple times
       is = is_idx(kxkyz_lo,ikxkyz); if (is /= 1) cycle
       call lu_back_substitution (fp_response(:,:,ikxkyz), diff_idx(:,ikxkyz), flds(iky,ikx,iz,it,:))
    end do

    if (advfield_coll) then
        phi(:,:,:,:) = flds(:,:,:,:,1)
        call sum_allreduce (phi)
    end if

    g = g_in

    ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + sum_jlm psi_jlm^{n+1}*delta_jl
    ! first two terms added via g_to_h subroutine
    if (advfield_coll) then
        call g_to_h (g, phi, fphi)
    end if

    ! add field particle contribution to RHS:
    if (fieldpart) then
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
           iky = iky_idx(kxkyz_lo,ikxkyz)
           ikx = ikx_idx(kxkyz_lo,ikxkyz)
           iz = iz_idx(kxkyz_lo,ikxkyz)
           it = it_idx(kxkyz_lo,ikxkyz)
           is = is_idx(kxkyz_lo,ikxkyz)

           do idx1 = 1, (jmax+1)*(lmax+1)**2
               ij  = 1 + mod(idx1-1,jmax+1)
               il  = 1 + int(sqrt(1.0*(idx1 - ij)/(jmax+1)))
               im  = (idx1 - ij)/(jmax+1) - (il-1)**2 + 1
               ll1 = il-1
               mm1 = -ll1 + im-1
               jj1 = ij-1

               if (density_conservation_tp.and.(jj1==0).and.(ll1==0)) then
                   isb = is
                   g(:,:,ikxkyz) = g(:,:,ikxkyz) - code_dt*flds(iky,ikx,iz,it, 2 + (is-1)*((jmax+1)*(lmax+1)**2*nspec) + (idx1-1)*nspec + (isb-1))&
                                                          *modmw(:,:,iz,is)/modmwnorm(iz)
               else
                   clm = sqrt(((2*ll1+1)*gamma(ll1-mm1+1.))/(4*pi*gamma(ll1+mm1+1.)))
                   do isb = 1, nspec
                       if (mm1 == 0) then
                           g(:,:,ikxkyz) = g(:,:,ikxkyz) + code_dt*spec(is)%vnew(isb)*clm*legendre_vpamu(ll1,mm1,:,:,iz)*spread(jm(:,mm1,iky,ikx,iz,is),1,nvpa) &
                                *flds(iky,ikx,iz,it, 2 + (is-1)*((jmax+1)*(lmax+1)**2*nspec) + (idx1-1)*nspec + (isb-1))&
                                *(spec(is)%mass/spec(isb)%mass)**(-1.5)*deltaj(ll1,jj1,is,isb,:,:,ia,iz)
                       else if (mm1 > 0) then
                           g(:,:,ikxkyz) = g(:,:,ikxkyz) + code_dt*spec(is)%vnew(isb)*clm*legendre_vpamu(ll1,mm1,:,:,iz)*spread(jm(:,mm1,iky,ikx,iz,is),1,nvpa) &
                                *flds(iky,ikx,iz,it, 2 + (is-1)*((jmax+1)*(lmax+1)**2*nspec) + (idx1-1)*nspec + (isb-1))&
                                *(spec(is)%mass/spec(isb)%mass)**(-1.5)*deltaj(ll1,jj1,is,isb,:,:,ia,iz)
                       else if (mm1 < 0) then
                           g(:,:,ikxkyz) = g(:,:,ikxkyz) + (-1)**mm1*code_dt*spec(is)%vnew(isb)*clm*legendre_vpamu(ll1,mm1,:,:,iz)*spread(jm(:,abs(mm1),iky,ikx,iz,is),1,nvpa) &
                                *flds(iky,ikx,iz,it, 2 + (is-1)*((jmax+1)*(lmax+1)**2*nspec) + (idx1-1)*nspec + (isb-1))&
                                *(spec(is)%mass/spec(isb)%mass)**(-1.5)*deltaj(ll1,jj1,is,isb,:,:,ia,iz)
                       end if
                   end do

               end if

           end do
        end do

    end if

    deallocate (flds)

    ! invert system to get h^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       do iv = 1, nvpa
           ghrs(nmu*(iv-1)+1 : nmu*iv) = g(iv,:,ikxkyz)
       end do
       call zgbtrs('No transpose', nvpa*nmu, nmu+1, nmu+1, 1, cdiffmat_band(:,:,iky,ikx,iz,is), 3*(nmu+1)+1, ipiv(:,iky,ikx,iz,is), ghrs, nvpa*nmu, info)
       do iv = 1, nvpa
           g(iv,:,ikxkyz) = ghrs(nmu*(iv-1)+1 : nmu*iv)
       end do
    end do

    ! get g^{n+1} from h^{n+1} and phi^{n+1}
    if (advfield_coll) then
        call g_to_h (g, phi, -fphi)
    end if

    !fields_updated = .false.

    deallocate (g_in)
    deallocate (field)
    deallocate (gvmutr)
    deallocate (ghrs)

  end subroutine advance_implicit_fp

  subroutine advance_vpadiff_implicit (phi, apar, g)

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_back_substitution
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use species, only: nspec, spec, has_electron_species
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx, zonal_mode
    use stella_geometry, only: dl_over_b
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use g_tofrom_h, only: g_to_h
    use gyro_averages, only: aj0v
    use fields, only: get_fields, efac, gamtot_h
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    integer :: imu
    integer :: idx
    real, dimension (:,:), allocatable :: tmp
    complex, dimension (:,:,:,:,:), allocatable :: flds
    complex, dimension (:), allocatable :: tmp2
    complex, dimension (:,:,:), allocatable :: flds_zf, g_in

    ia = 1

    ! store input g for use later, as g will be overwritten below
    allocate (g_in(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    g_in = g

    ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
    ! g = g^{***}.  tridag below inverts above equation to get h_inh^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), g(:,imu,ikxkyz))
       end do
    end do

    allocate (flds(naky,nakx,-nzgrid:nzgrid,ntubes,nresponse_vpa))
    allocate (flds_zf(nakx,ntubes,nresponse_vpa)); flds_zf = 0.

    ! need to obtain phi^{n+1} and conservation terms using response matrix approach
    ! first get phi_inh^{n+1}
    call get_fields (g, phi, apar, dist='h', skip_fsa=.true.)
    flds(:,:,:,:,1) = phi

    idx = 2
    ! get upar_inh^{n+1}
    if (momentum_conservation) then
       call get_upar (g, flds(:,:,:,:,idx:idx+nspec-1))
       idx = idx + nspec
    end if

    ! get temp_inh^{n+1}
    if (energy_conservation) call get_temp (g, flds(:,:,:,:,idx:idx+nspec-1))

    phi = 0.0
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ! all is indices inside ikxkyz super-index have same info
       ! no need to compute multiple times
       is = is_idx(kxkyz_lo,ikxkyz) ; if (is /= 1) cycle
       call lu_back_substitution (vpadiff_response(:,:,ikxkyz), vpadiff_idx(:,ikxkyz), &
                                  flds(iky,ikx,iz,it,:))
       if (.not.has_electron_species(spec).and.zonal_mode(iky) &
        .and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
          flds_zf(ikx,it,:) = flds_zf(ikx,it,:) + dl_over_b(ia,iz)*flds(iky,ikx,iz,it,:)
       endif
       phi(iky,ikx,iz,it) = flds(iky,ikx,iz,it,1)
    end do

    if (.not.has_electron_species(spec).and.zonal_mode(1) &
        .and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
      !complete flux-surface average
      call sum_allreduce(flds_zf)
      do it=1, ntubes
        do ikx=1, nakx
          call lu_back_substitution (vpadiff_zf_response(:,:,ikx), vpadiff_zf_idx(:,ikx), &
                                     flds_zf(ikx,it,:))
          !multiply by Q, which has a single non-zero component
          flds_zf(ikx,it,1) = (efac/gamtot_h)*flds_zf(ikx,it,1)
          flds_zf(ikx,it,2:) = 0.
        enddo
      enddo

      allocate (tmp2(nresponse_vpa))
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        it = it_idx(kxkyz_lo,ikxkyz)
        is = is_idx(kxkyz_lo,ikxkyz)

        if(iky.ne.1.or.is.ne.1) cycle

        tmp2 = flds_zf(ikx,it,:)
        call lu_back_substitution (vpadiff_response(:,:,ikxkyz), vpadiff_idx(:,ikxkyz), &
                                  tmp2)

        phi(1,ikx,iz,it) = phi(1,ikx,iz,it) + tmp2(1)
      enddo
      deallocate (tmp2)

    endif

    call sum_allreduce (phi)

    g = g_in

    ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + 2*dt*nu*J0*F0*(vpa*upar+(v^2-3/2)*temp)
    ! first two terms added via g_to_h subroutine
    call g_to_h (g, phi, fphi)

    allocate (tmp(nvpa,nmu))

    if (momentum_conservation .or. energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          tmp = 2.0*code_dt*spec(is)%vnew(is) &
               *spread(maxwell_vpa(:,is),2,nmu)*spread(aj0v(:,ikxkyz)*maxwell_mu(1,iz,:,is),1,nvpa)
          if (momentum_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) + tmp*spread(vpa*flds(iky,ikx,iz,it,is+1),2,nmu)
          if (energy_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) &
               + tmp*(spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa)-1.5)*flds(iky,ikx,iz,it,idx+is-1)
       end do
    end if

    deallocate (tmp, flds, flds_zf)

    ! now invert system to get h^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          call tridag (1, aa_vpa(:,is), bb_vpa(:,ikxkyz), cc_vpa(:,is), g(:,imu,ikxkyz))
       end do
    end do

    ! now get g^{n+1} from h^{n+1} and phi^{n+1}
    call g_to_h (g, phi, -fphi)

    deallocate (g_in)

  end subroutine advance_vpadiff_implicit

  subroutine advance_mudiff_implicit (phi, apar, g)

    use mp, only: sum_allreduce
    use finite_differences, only: tridag
    use linear_solve, only: lu_back_substitution
    use stella_time, only: code_dt
    use run_parameters, only: fphi
    use species, only: nspec, spec, has_electron_species
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nmu, nvpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
    use vpamu_grids, only: set_vpa_weights
    use kt_grids, only: naky, nakx, zonal_mode
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use dist_fn_arrays, only: kperp2
    use gyro_averages, only: aj0v, aj1v
    use g_tofrom_h, only: g_to_h
    use fields, only: get_fields, efac, gamtot_h
    use stella_geometry, only: bmag, dl_over_b
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg

    ! TMP FOR TESTING
!    use vpamu_grids, only: mu

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in out) :: g

    integer :: ikxkyz, iky, ikx, iz, it, is, ia
    integer :: iv
    integer :: idx

    ! TMP FOR TESTING
!    integer :: imu

    real, dimension (:,:), allocatable :: tmp
    complex, dimension (:), allocatable :: tmp2
    complex, dimension (:,:,:,:,:), allocatable :: flds
    complex, dimension (:,:,:), allocatable :: flds_zf, g_in

    ia = 1

    ! store input g for use later, as g will be overwritten below
    allocate (g_in(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    g_in = g

    ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
    ! g = g^{***}.  tridag below inverts above equation to get h_inh^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       ! TMP FOR TESTING
!       do imu = 1, nmu
!          g(:,imu,ikxkyz) = maxwell_vpa*maxwell_mu(1,iz,imu)
!       end do
       do iv = 1, nvpa
          call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), g(iv,:,ikxkyz))
       end do
       ! TMP FOR TESTING
!       iv = nvpa/2
!       do imu = 1, nmu
!          write (*,*) 'ggg', mu(imu), real(g(iv,imu,ikxkyz)), aimag(g(iv,imu,ikxkyz)), maxwell_vpa(iv,is)*maxwell_mu(1,iz,imu)
!       end do
    end do

    allocate (flds(naky,nakx,-nzgrid:nzgrid,ntubes,nresponse_mu))
    allocate (flds_zf(nakx,ntubes,nresponse_mu)); flds_zf = 0.

    ! need to obtain phi^{n+1} and conservation terms using response matrix approach
    ! first get phi_inh^{n+1}
    call get_fields (g, phi, apar, dist='h', skip_fsa=.true.)
    flds(:,:,:,:,1) = phi

    idx = 2
    ! get upar_inh^{n+1}
    if (momentum_conservation) then
       call get_uperp (g, flds(:,:,:,:,idx:idx+nspec-1))
       idx = idx + nspec
    end if

    ! get temp_inh^{n+1}
    if (energy_conservation) call get_temp_mu (g, flds(:,:,:,:,idx:idx+nspec-1))

    phi = 0.0
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ! all is indices inside ikxkyz super-index have same info
       ! no need to compute multiple times
       is = is_idx(kxkyz_lo,ikxkyz) ; if (is /= 1) cycle
       call lu_back_substitution (mudiff_response(:,:,ikxkyz), mudiff_idx(:,ikxkyz), &
            flds(iky,ikx,iz,it,:))
       if (.not.has_electron_species(spec).and.zonal_mode(iky) &
        .and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
          flds_zf(ikx,it,:) = flds_zf(ikx,it,:) + dl_over_b(ia,iz)*flds(iky,ikx,iz,it,:)
       endif
       phi(iky,ikx,iz,it) = flds(iky,ikx,iz,it,1)
    end do

    if (.not.has_electron_species(spec).and.zonal_mode(1) &
        .and.adiabatic_option_switch==adiabatic_option_fieldlineavg) then
      !complete flux-surface average
      call sum_allreduce(flds_zf)
      do it=1, ntubes
        do ikx=1, nakx
          call lu_back_substitution (mudiff_zf_response(:,:,ikx), mudiff_zf_idx(:,ikx), &
                                     flds_zf(ikx,it,:))
          !multiply by Q, which has a single non-zero component
          flds_zf(ikx,it,1) = (efac/gamtot_h)*flds_zf(ikx,it,1)
          flds_zf(ikx,it,2:) = 0.
        enddo
      enddo

      allocate (tmp2(nresponse_mu))
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iky = iky_idx(kxkyz_lo,ikxkyz)
        ikx = ikx_idx(kxkyz_lo,ikxkyz)
        iz = iz_idx(kxkyz_lo,ikxkyz)
        it = it_idx(kxkyz_lo,ikxkyz)
        is = is_idx(kxkyz_lo,ikxkyz)

        if(iky.ne.1.or.is.ne.1) cycle

        tmp2 = flds_zf(ikx,it,:)
        call lu_back_substitution (mudiff_response(:,:,ikxkyz), mudiff_idx(:,ikxkyz), &
                                  tmp2)

        phi(1,ikx,iz,it) = phi(1,ikx,iz,it) + tmp2(1)
      enddo
      deallocate (tmp2)

    endif

    call sum_allreduce (phi)

    g = g_in

    ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + ...
    ! first two terms added via g_to_h subroutine
    call g_to_h (g, phi, fphi)

    allocate (tmp(nvpa,nmu))

    if (momentum_conservation .or. energy_conservation) then
       do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
          iky = iky_idx(kxkyz_lo,ikxkyz)
          ikx = ikx_idx(kxkyz_lo,ikxkyz)
          iz = iz_idx(kxkyz_lo,ikxkyz)
          it = it_idx(kxkyz_lo,ikxkyz)
          is = is_idx(kxkyz_lo,ikxkyz)
          tmp = 2.0*code_dt*spec(is)%vnew(is) &
               *spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(1,iz,:,is),1,nvpa)
          if (momentum_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) + tmp*kperp2(iky,ikx,ia,iz) &
             * spread(vperp2(ia,iz,:)*aj1v(:,ikxkyz),1,nvpa)*(spec(is)%smz/bmag(ia,iz))**2 &
             * flds(iky,ikx,iz,it,is+1)
          if (energy_conservation) &
               g(:,:,ikxkyz) = g(:,:,ikxkyz) &
               + flds(iky,ikx,iz,it,idx+is-1)*tmp*spread(aj0v(:,ikxkyz),1,nvpa) &
             * (spread(vpa**2,2,nmu)+spread(vperp2(1,iz,:),1,nvpa)-1.5)
       end do
    end if

    deallocate (tmp, flds, flds_zf)

    ! now invert system to get h^{n+1}
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do iv = 1, nvpa
          call tridag (1, aa_mu(iz,:,is), bb_mu(:,ikxkyz), cc_mu(iz,:,is), g(iv,:,ikxkyz))
       end do
    end do

    ! now get g^{n+1} from h^{n+1} and phi^{n+1}
    call g_to_h (g, phi, -fphi)

    deallocate (g_in)

  end subroutine advance_mudiff_implicit

  subroutine advance_hyper_dissipation (g)

    use stella_time, only: code_dt
    use zgrid, only: nzgrid, ntubes, zed
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: kperp2
    use kt_grids, only: ikx_max, naky, nakx
    use kt_grids, only: aky, akx, theta0, zonal_mode
    use stella_geometry, only: geo_surf, q_as_x

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: ia, ivmu, iz, it, iky
    real :: tfac
    real :: k2max

    if (.not.use_physical_ksqr) then
       ! avoid spatially dependent kperp

       !get k2max at outboard midplane
       k2max = akx(ikx_max)**2 + aky(naky)**2
       tfac= geo_surf%shat**2
       if(q_as_x) tfac = 1.0

       ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1,ntubes
           do iz = -nzgrid, nzgrid
             do iky = 1, naky
               if(zonal_mode(iky)) then
                  g(iky,:,iz,it,ivmu) = g(iky,:,iz,it,ivmu)/(1.+code_dt*(akx(:)**2/k2max)**2*D_hyper)
               else
                  g(iky,:,iz,it,ivmu) = g(iky,:,iz,it,ivmu)/(1.+code_dt*(aky(iky)**2 &
                                       * (1.0+ tfac*(zed(iz) - theta0(iky,:))**2)/k2max)**2*D_hyper)
               endif
             enddo
           enddo
         enddo
!        g(:,:,:,:,ivmu) = g(:,:,:,:,ivmu)/(1.+code_dt &
!           * (spread(spread(spread(akx**2,1,naky)+spread(aky**2,2,nakx),3,nztot),4,ntubes)/k2max)**2*D_hyper)
       end do
    else
       k2max = maxval(kperp2)
       ia = 1
       ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          g(:,:,:,:,ivmu) = g(:,:,:,:,ivmu)/(1.+code_dt*(spread(kperp2(:,:,ia,:),4,ntubes)/k2max)**2*D_hyper)
       end do
    end if

  end subroutine advance_hyper_dissipation

end module dissipation
