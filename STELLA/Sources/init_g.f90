!> This module contains the subroutines which set the initial value of the
!! fields and the distribution function.

module init_g
  implicit none

  public :: ginit
  public :: init_init_g, finish_init_g
  public :: width0
  public :: scale_to_phiinit, phiinit
  public :: tstart
  public :: reset_init

  private

  ! knobs
  integer :: ginitopt_switch
  integer, parameter :: ginitopt_default = 1,  &
       ginitopt_noise = 2, ginitopt_restart_many = 3, &
       ginitopt_kpar = 4, ginitopt_nltest = 5, &
       ginitopt_kxtest = 6, ginitopt_rh = 7, &
       ginitopt_remap = 8

  real :: width0, phiinit, imfac, refac, zf_init
  real :: den0, upar0, tpar0, tperp0
  real :: den1, upar1, tpar1, tperp1
  real :: den2, upar2, tpar2, tperp2
  real :: tstart, scale, kxmax, kxmin
  logical :: chop_side, left, even, scale_to_phiinit
  character(300), public :: restart_file
  character (len=150) :: restart_dir

  logical :: initialized = .false.
  logical :: exist

contains

  subroutine init_init_g

    use stella_save, only: init_save, read_many
    use stella_layouts, only: init_stella_layouts
    use system_fortran, only: systemf
    use mp, only: proc0, broadcast

    implicit none

    integer :: ind_slash

    if (initialized) return
    initialized = .true.

    call init_stella_layouts

    if (proc0) call read_parameters

    ! prepend restart_dir to restart_file
    ! append trailing slash if not exists
    if(restart_dir(len_trim(restart_dir):) /= "/") &
         restart_dir=trim(restart_dir)//"/"

    if (proc0) call systemf('mkdir -p '//trim(restart_dir))

    !Determine if restart file contains "/" if so split on this point to give DIR//FILE
    !so restart files are created in DIR//restart_dir//FILE
    ind_slash=index(restart_file,"/",.True.)
    if (ind_slash.EQ.0) then !No slash present
       restart_file=trim(restart_dir)//trim(restart_file)
    else !Slash present
       restart_file=trim(restart_file(1:ind_slash))//trim(restart_dir)//trim(restart_file(ind_slash+1:))
    endif

    call broadcast (ginitopt_switch)
    call broadcast (width0)
    call broadcast (refac)
    call broadcast (imfac)
    call broadcast (den0)
    call broadcast (upar0)
    call broadcast (tpar0)
    call broadcast (tperp0)
    call broadcast (den1)
    call broadcast (upar1)
    call broadcast (tpar1)
    call broadcast (tperp1)
    call broadcast (den2)
    call broadcast (upar2)
    call broadcast (tpar2)
    call broadcast (tperp2)
    call broadcast (phiinit)
    call broadcast (zf_init)
    call broadcast (kxmax)
    call broadcast (kxmin)
    call broadcast (tstart)
    call broadcast (chop_side)
    call broadcast (even)
    call broadcast (left)
    call broadcast (restart_file)
    call broadcast (read_many)
    call broadcast (scale_to_phiinit)
    call broadcast (scale)

    call init_save (restart_file)

  end subroutine init_init_g

  subroutine ginit (restarted,istep0)

    use stella_save, only: init_tstart
    logical, intent (out) :: restarted
    integer, intent (out) :: istep0
    integer :: istatus

    restarted = .false.
    istep0 = 0
    select case (ginitopt_switch)
    case (ginitopt_default)
       call ginit_default
    case (ginitopt_noise)
       call ginit_noise
    case (ginitopt_kpar)
       call ginit_kpar
     case (ginitopt_rh)
        call ginit_rh
     case (ginitopt_remap)
        call ginit_remap
    case (ginitopt_restart_many)
       call ginit_restart_many 
       call init_tstart (tstart, istep0, istatus)
       restarted = .true.
       scale = 1.
!    case (ginitopt_nltest)
!       call ginit_nltest
!    case (ginitopt_kxtest)
!       call ginit_kxtest
    end select
    
  end subroutine ginit

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, run_name, input_unit_exist
    use text_options, only: text_option, get_option_value
    use stella_save, only: read_many

    implicit none

    type (text_option), dimension (8), parameter :: ginitopts = &
         (/ text_option('default', ginitopt_default), &
            text_option('noise', ginitopt_noise), &
            text_option('many', ginitopt_restart_many), &
            text_option('nltest', ginitopt_nltest), &
            text_option('kxtest', ginitopt_kxtest), &
            text_option('kpar', ginitopt_kpar), &
            text_option('rh', ginitopt_rh), &
            text_option('remap', ginitopt_remap) &
            /)
    character(20) :: ginit_option
    namelist /init_g_knobs/ ginit_option, width0, phiinit, chop_side, &
         restart_file, restart_dir, read_many, left, scale, tstart, zf_init, &
         den0, upar0, tpar0, tperp0, imfac, refac, even, &
         den1, upar1, tpar1, tperp1, &
         den2, upar2, tpar2, tperp2, &
         kxmax, kxmin, scale_to_phiinit

    integer :: ierr, in_file

    tstart = 0.
    scale = 1.0
    ginit_option = "default"
    width0 = -3.5
    refac = 1.
    imfac = 0.
    den0 = 1.
    upar0 = 0.
    tpar0 = 0.
    tperp0 = 0.
    den1 = 0.
    upar1 = 0.
    tpar1 = 0.
    tperp1 = 0.
    den2 = 0.
    upar2 = 0.
    tpar2 = 0.
    tperp2 = 0.
    phiinit = 1.0
    zf_init = 1.0
    kxmax = 1.e100
    kxmin = 0.
    chop_side = .false.
    scale_to_phiinit = .false.
    left = .true.
    even = .true.

    restart_file = trim(run_name)//".nc"
    restart_dir = "./"
    in_file = input_unit_exist ("init_g_knobs", exist)
!    if (exist) read (unit=input_unit("init_g_knobs"), nml=init_g_knobs)
    if (exist) read (unit=in_file,nml=init_g_knobs)

    ierr = error_unit()
    call get_option_value &
         (ginit_option, ginitopts, ginitopt_switch, &
         ierr, "ginit_option in ginit_knobs",stop_on_error=.true.)
  end subroutine read_parameters

  subroutine ginit_default

    use constants, only: zi
    use species, only: spec
    use zgrid, only: nzgrid, zed
    use kt_grids, only: naky, nakx
    use kt_grids, only: theta0
    use kt_grids, only: reality, zonal_mode
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: gvmu
    use stella_layouts, only: kxkyz_lo, iz_idx, ikx_idx, iky_idx, is_idx
    use ran, only: ranf

    implicit none

    complex, dimension (naky,nakx,-nzgrid:nzgrid) :: phi
    logical :: right
    integer :: ikxkyz
    integer :: iz, iky, ikx, is, ia

    right = .not. left

    do iz = -nzgrid, nzgrid
       phi(:,:,iz) = exp(-((zed(iz)-theta0)/width0)**2)*cmplx(1.0,1.0)
    end do

    ! this is a messy way of doing things
    ! could tidy it up a bit
    if (sum(cabs(phi)) < epsilon(0.)) then
       do iz = -nzgrid, nzgrid
          phi(:,:,iz) = exp(-(zed(iz)/width0)**2)*cmplx(1.0,1.0)
       end do
    end if

    if (chop_side) then
       if (left) phi(:,:,:-1) = 0.0
       if (right) phi(:,:,1:) = 0.0
    end if

    if (reality .and. zonal_mode(1)) phi(1,:,:) = 0.0

    ia = 1

    gvmu = 0.
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       gvmu(:,:,ikxkyz) = phiinit*phi(iky,ikx,iz)/abs(spec(is)%z) &
            * (den0 + 2.0*zi*spread(vpa,2,nmu)*upar0) &
            * spread(maxwell_mu(ia,iz,:,is),1,nvpa)*spread(maxwell_vpa(:,is),2,nmu)*maxwell_fac(is)
    end do

  end subroutine ginit_default

  ! initialize two kys and kx=0
!   subroutine ginit_nltest

!     use mp, only: proc0
!     use species, only: spec
!     use zgrid, only: nzgrid, bmag
!     use kt_grids, only: naky, ntheta0
!     use vpamu_grids, only: nvgrid, vpa, mu
!     use dist_fn_arrays, only: gnew, gold
!     use stella_layouts, only: gxyz_lo, iv_idx, is_idx, imu_idx

!     implicit none

!     complex, dimension (-nzgrid:nzgrid,ntheta0,naky) :: phi
!     logical :: right
!     integer :: iglo
!     integer :: ig, ik, it, is, iv, imu

!     right = .not. left

!     if (naky < 4 .or. ntheta0 < 2) then
!        if (proc0) write (*,*) 'must have at least 2 kxs and 4 kys to use nltest init option. aborting.'
!        stop
!     end if

!     phi = 0.0
!     do ig = -nzgrid, nzgrid
!        phi(ig,2,2) = 1.0!exp(-((theta(ig)-theta0(2,2))/width0)**2)*cmplx(1.0,1.0)
!     end do
    
!     gnew = 0.0
!     do iglo = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
!        iv = iv_idx(gxyz_lo,iglo)
!        is = is_idx(gxyz_lo,iglo)
!        imu = imu_idx(gxyz_lo,iglo)
!        it = 1
!        do ik = 2, 3
!           gnew(:,it,ik,iglo) = exp(-2.0*mu(imu)*bmag)*phi(:,it,ik) &
!                *spec(is)%z*phiinit*exp(-vpa(iv)**2)
!        end do
!     end do
!     gold = gnew

!   end subroutine ginit_nltest

!   subroutine ginit_kxtest

!     use constants, only: zi
!     use species, only: spec
!     use zgrid, only: itor_over_b
!     use kt_grids, only: ntheta0, akx, naky
!     use vpamu_grids, only: nvgrid, energy, vpa
!     use dist_fn_arrays, only: gnew, gold
!     use stella_layouts, only: gxyz_lo, iv_idx, is_idx, imu_idx

!     implicit none

!     integer :: iglo
!     integer :: ik, it, is, imu, iv

!     do iglo = gxyz_lo%llim_proc, gxyz_lo%ulim_proc
!        iv = iv_idx(gxyz_lo,iglo)
!        is = is_idx(gxyz_lo,iglo)
!        imu = imu_idx(gxyz_lo,iglo)
!        do it = 1, ntheta0
!           do ik = 1, naky
!              gnew(:,it,ik,iglo) = exp(-zi*akx(it)*itor_over_B*vpa(iv)/spec(is)%zstm) &
!                   *exp(-energy(:,iv,imu))*spec(is)%z*phiinit
!           end do
!        end do
!     end do
!     gold = gnew

!   end subroutine ginit_kxtest

!   !> Initialise with only the kparallel = 0 mode.
  
!   subroutine single_initial_kx(phi)
!     use zgrid, only: nzgrid 
!     use kt_grids, only: naky, ntheta0
!     use mp, only: mp_abort
!     implicit none
!     complex, dimension (-nzgrid:nzgrid,ntheta0,naky), intent(inout) :: phi
!     real :: a, b
!     integer :: ig, ik, it

!     if (ikx_init  < 2 .or. ikx_init > (ntheta0+1)/2) then
!       call mp_abort("The subroutine single_initial_kx should only be called when 1 < ikx_init < (ntheta0+1)/2")
!     end if

!     do it = 1, ntheta0
!       if (it .ne. ikx_init) then 
!          do ik = 1, naky
!             do ig = -nzgrid, nzgrid
!                a = 0.0
!                b = 0.0 
!                phi(ig,it,ik) = cmplx(a,b)
!              end do
!          end do
!        end if
!     end do
!   end subroutine single_initial_kx

  !> Initialise the distribution function with random noise. This is the default
  !! initialisation option. Each different mode is given a random amplitude
  !! between zero and one.

  subroutine ginit_noise

    use mp, only: proc0, broadcast
    use dist_fn_arrays, only: kperp2
    use species, only: spec
    use zgrid, only: nzgrid, ntubes
    use extended_zgrid, only: ikxmod, nsegments, neigen
    use extended_zgrid, only: it_right
    use extended_zgrid, only: periodic
    use kt_grids, only: naky, nakx, reality, zonal_mode
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: gvmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
    use mp, only: proc0, broadcast, max_allreduce
    use mp, only: scope, crossdomprocs, subprocs
    use file_utils, only: runtype_option_switch, runtype_multibox
    use ran

    implicit none

    complex, dimension (naky,nakx,-nzgrid:nzgrid,ntubes) :: phi
    real :: a, b, kmin
    integer :: ikxkyz, iz, it, iky, ikx, is, ie, iseg, ia
    integer :: itmod

    if (naky == 1 .and. nakx==1) then
       if (proc0) then
          write (*,*) 'noise initialization option is not suited for single mode simulations.'
          write (*,*) 'using default initialization option'
       end if
       call ginit_default
       return
    else
       ! zero out ky=kx=0 mode
       phi(1,1,:,:) = 0.0
    end if

    ia = 1
    if (proc0) then
       phi(1,1,:,:) = 0.0
       kmin = 1.e6
       if (naky > 1) kmin = minval(kperp2(2,1,ia,:))
       if (nakx > 1) kmin = min(kmin,minval(kperp2(1,2,ia,:)))

       if(runtype_option_switch == runtype_multibox) then
         call scope(crossdomprocs)
         call max_allreduce (kmin)
         call scope(subprocs)
       end if

       ! keep old (ikx, iky) loop order to get old results exactly: 
       !Fill phi with random (complex) numbers between -0.5 and 0.5
       do ikx = 1, nakx
          do iky = 1, naky
             do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                   a = ranf()-0.5
                   b = ranf()-0.5
                   ! do not populate high k modes with large amplitudes
                   if ((ikx > 1 .or. iky > 1) .and. (kperp2(iky,ikx,ia,iz) .ge. kmin))  then
                     !the following as an extra factor of kmin to offset the Gamma-1 in quasineutrality
                     phi(iky,ikx,iz,it) =cmplx(a,b)*kmin*kmin/kperp2(iky,ikx,ia,iz)
                   else
                     phi(iky,ikx,iz,it) = 0.0
                   endif
                end do
                if (chop_side) then
                   if (left) then
                      phi(iky,ikx,:-1,it) = 0.0
                   else
                      phi(iky,ikx,1:,it) = 0.0
                   endif
                end if
             end do
          end do
       end do

       ! enforce periodicity where required
       do iky = 1, naky
          if (periodic(iky)) then
             phi(1,:,nzgrid,:) = phi(1,:,-nzgrid,:)
          end if
       end do

       ! zero out the kx=ky=0 mode and apply optional
       ! scaliing factor to all zonal modes
       if (zonal_mode(1)) then
          !Apply scaling factor
          phi(1,:,:,:) = phi(1,:,:,:)*zf_init
          
          !Set ky=kx=0.0 mode to zero in amplitude
          phi(1,1,:,:) = 0.0
       end if
       
       !Apply reality condition (i.e. -kx mode is conjugate of +kx mode)
       if (reality) then
          do ikx = nakx/2+2, nakx
             phi(1,ikx,:,:) = conjg(phi(1,nakx-ikx+2,:,:))
          enddo
       end if
       
    end if

    do iky = 1, naky
       do ie = 1, neigen(iky)
          ! enforce zero BC at ends of domain, unless periodic
          if (.not.periodic(iky)) then
             phi(iky,ikxmod(1,ie,iky),-nzgrid,:) = 0.0
             phi(iky,ikxmod(nsegments(ie,iky),ie,iky),nzgrid,:) = 0.0
          end if
          ! enforce equality of g values at duplicate zed points
          if (nsegments(ie,iky) > 1) then
             do it = 1, ntubes
                itmod = it
                do iseg = 2, nsegments(ie,iky)
                   phi(iky,ikxmod(iseg,ie,iky),-nzgrid,it_right(itmod)) = phi(iky,ikxmod(iseg-1,ie,iky),nzgrid,itmod)
                   itmod = it_right(itmod)
                end do
             end do
          end if
       end do
    end do
    
    call broadcast (phi)

    !Now set g using data in phi
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iz = iz_idx(kxkyz_lo,ikxkyz)
       it = it_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iky = iky_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       gvmu(:,:,ikxkyz) = spec(is)%z*phiinit*phi(iky,ikx,iz,it) &
            * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is)
    end do

  end subroutine ginit_noise

  subroutine ginit_kpar

!    use species, only: spec, has_electron_species
    use zgrid, only: nzgrid, zed
    use kt_grids, only: naky, nakx, theta0
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: vpa, vperp2
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use dist_fn_arrays, only: gvmu
    use stella_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx
    use constants, only: zi

    implicit none

    complex, dimension (naky,nakx,-nzgrid:nzgrid) :: phi, odd
    real, dimension (-nzgrid:nzgrid) :: dfac, ufac, tparfac, tperpfac
    integer :: ikxkyz
    integer :: iz, iky, ikx, imu, iv, ia, is
    
    phi = 0.
    odd = 0.
    if (width0 > 0.) then
       do iz = -nzgrid, nzgrid
          phi(:,:,iz) = exp(-((zed(iz)-theta0)/width0)**2)*cmplx(refac, imfac)
       end do
    else
       do iz = -nzgrid, nzgrid
          phi(:,:,iz) = cmplx(refac, imfac)
       end do
    end if
    if (chop_side) then
       if (left) then
          phi(:,:,:-1) = 0.0
       else
          phi(:,:,1:) = 0.0
       end if
    end if
    
    odd = zi * phi
    
    dfac     = den0   + den1 * cos(zed) + den2 * cos(2.*zed) 
    ufac     = upar0  + upar1* sin(zed) + upar2* sin(2.*zed) 
    tparfac  = tpar0  + tpar1* cos(zed) + tpar2* cos(2.*zed) 
    tperpfac = tperp0 + tperp1*cos(zed) + tperp2*cos(2.*zed) 
    
    ia = 1
    ! charge dependence keeps initial Phi from being too small
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       do imu = 1, nmu
          do iv = 1, nvpa
             gvmu(iv,imu,ikxkyz) = phiinit &
                  * ( dfac(iz)*phi(iky,ikx,iz) &
                  + 2.0*vpa(iv)*ufac(iz)*odd(iky,ikx,iz) &
                  + (vpa(iv)**2-0.5)*tparfac(iz)*phi(iky,ikx,iz) &
                  + tperpfac(iz)*(vperp2(ia,iz,imu)-1.)*phi(iky,ikx,iz) )
          end do
       end do
       gvmu(:,:,ikxkyz) = gvmu(:,:,ikxkyz) &
            * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is)
    end do

! FLAG -- should be uncommented, which means I need to fix flae
!    if (has_electron_species(spec)) then
!       call flae (gold, gnew)
!       gold = gold - gnew
!    end if   
!    gnew = gold

  end subroutine ginit_kpar

  subroutine ginit_rh
    
    use species, only: spec
    use dist_fn_arrays, only: gvmu, kperp2
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use vpamu_grids, only: nvpa, nmu
    use kt_grids, only: akx

    implicit none
    
    integer :: ikxkyz, iky, ikx, iz, is, ia
    
    ! initialize g to be a Maxwellian with a constant density perturbation

    gvmu = 0.
    
    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz) ; if (iky /= 1) cycle
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)
       
       if(abs(akx(ikx)) < kxmax .and. abs(akx(ikx)) > kxmin) then
         gvmu(:,:,ikxkyz) = spec(is)%z*0.5*phiinit*kperp2(iky,ikx,ia,iz) &
            * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is)
       end if
    end do

  end subroutine ginit_rh

  subroutine ginit_remap

    use species, only: spec
    use dist_fn_arrays, only: gvmu
    use stella_layouts, only: kxkyz_lo
    use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use vpamu_grids, only: nvpa, nmu

    implicit none

    integer :: ikxkyz, iky, ikx, iz, is, ia

    ! initialize g to be a Maxwellian with a constant density perturbation

    gvmu = 0.

    ia = 1
    do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
       iky = iky_idx(kxkyz_lo,ikxkyz)
       ikx = ikx_idx(kxkyz_lo,ikxkyz)
       iz = iz_idx(kxkyz_lo,ikxkyz)
       is = is_idx(kxkyz_lo,ikxkyz)

       !if((ikx.eq.15.and.iky.eq.5).or.((ikx-nakx).eq.-12.and.iky.eq.3)) then
       if((ikx.eq.1.and.iky.eq.2)) then
         gvmu(:,:,ikxkyz) = spec(is)%z*phiinit &
            * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is)
       endif
    end do

  end subroutine ginit_remap

  subroutine ginit_restart_many

    use dist_fn_arrays, only: gvmu
    use stella_save, only: stella_restore
    use mp, only: proc0
    use file_utils, only: error_unit
    
    implicit none

    integer :: istatus, ierr

    ! should really check if profile_variation=T here but need
    ! to move profile_variation to module that is accessible here
    call stella_restore (gvmu, scale, istatus)

    if (istatus /= 0) then
       ierr = error_unit()
       if (proc0) write(ierr,*) "Error reading file: ", trim(restart_file)
       gvmu = 0.
    end if

  end subroutine ginit_restart_many

  subroutine reset_init

    ginitopt_switch = ginitopt_restart_many

  end subroutine reset_init

!   subroutine flae (g, gavg)

!     use species, only: spec, electron_species 
!     use zgrid, only: nzgrid, delthet, jacob
!     use kt_grids, only: aky, ntheta0
!     use vpamu_grids, only: nvgrid
!     use stella_layouts, only: gxyz_lo, is_idx
!     complex, dimension (-nzgrid:,:,:,gxyz_lo%llim_proc:), intent (in) :: g
!     complex, dimension (-nzgrid:,:,:,gxyz_lo%llim_proc:), intent (out) :: gavg

!     real :: wgt
!     integer :: iglo, it, ik
    
!     gavg = 0.
!     wgt = 1./sum(delthet*jacob)

!     do iglo = gxyz_lo%llim_proc, gxyz_lo%ulim_proc       
!        if (spec(is_idx(gxyz_lo, iglo))%type /= electron_species) cycle
!        ik = 1
!        if (aky(ik) > epsilon(0.)) cycle
!        do it = 1, ntheta0
!           gavg(:,it,ik,iglo) = sum(g(:,it,ik,iglo)*delthet*jacob)*wgt
!        end do
!     end do
    
!   end subroutine flae

  subroutine finish_init_g

    use stella_save, only: finish_save

    implicit none

    initialized = .false.

    call finish_save

  end subroutine finish_init_g

end module init_g
