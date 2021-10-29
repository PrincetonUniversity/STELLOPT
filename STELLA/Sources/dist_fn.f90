module dist_fn

  implicit none

  public :: init_gxyz
  public :: init_dist_fn, finish_dist_fn
  public :: adiabatic_option_switch
  public :: adiabatic_option_fieldlineavg

  private
  
  logical :: dist_fn_initialized = .false.
  logical :: gxyz_initialized = .false.
  logical :: kp2init = .false.
  logical :: vp2init = .false.
!  logical :: bessinit = .false.
  logical :: readinit = .false.

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3

  logical :: debug = .false.

contains

  subroutine init_gxyz (restarted)

    use dist_fn_arrays, only: gvmu, gold, gnew
    use redistribute, only: gather, scatter
    use dist_redistribute, only: kxkyz2vmu
    use physics_flags, only: radial_variation
    use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
    use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
    use kt_grids, only: nalpha, nakx, naky, multiply_by_rho
    use vpamu_grids, only: mu, vpa, vperp2
    use zgrid, only: nzgrid, ntubes
    use species, only: spec, pfac
    use stella_geometry, only: dBdrho, gfac


    implicit none

    real :: corr
    integer :: ivmu, is, imu, iv, it, iz, ia
    real, dimension (:,:), allocatable :: energy
    complex, dimension (:,:), allocatable :: g0k
    logical, intent(in) :: restarted


    if (gxyz_initialized) return
    gxyz_initialized = .false.


    ! get version of g that has ky,kx,z local
    call gather (kxkyz2vmu, gvmu, gnew)

    ia = 1

    !calculate radial corrections to F0 for use in Krook operator, as well as g1 from initialization
    if(radial_variation) then
      !init_g uses maxwellians, so account for variation in temperature, density, and B

      allocate (energy(nalpha,-nzgrid:nzgrid))
      allocate (g0k(naky,nakx))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        is  = is_idx(vmu_lo, ivmu)
        imu = imu_idx(vmu_lo, ivmu)
        iv  = iv_idx(vmu_lo, ivmu)
        energy = (vpa(iv)**2 + vperp2(:,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid

            corr = -( pfac*(spec(is)%fprim+spec(is)%tprim*(energy(ia,iz)-1.5)) &
                     + 2*gfac*mu(imu)*dBdrho(iz))
         
            if(.not.restarted) then
              g0k = corr*gnew(:,:,iz,it,ivmu)
              call multiply_by_rho(g0k)
              gnew(:,:,iz,it,ivmu) = gnew(:,:,iz,it,ivmu) + g0k
            endif
          enddo
        enddo
      enddo
      deallocate(energy, g0k)

      if(.not.restarted) call scatter(kxkyz2vmu,gnew,gvmu)
    endif

    gold = gnew


  end subroutine init_gxyz

  subroutine init_dist_fn

    use mp, only: proc0
    use stella_layouts, only: init_dist_fn_layouts
    use gyro_averages, only: init_bessel

    implicit none

    if (dist_fn_initialized) return
    dist_fn_initialized = .true.

    debug = debug .and. proc0
    
    if (debug) write (*,*) 'dist_fn::init_dist_fn::read_parameters'
    call read_parameters
    if (debug) write (*,*) 'dist_fn::init_dist_fn::allocate_arrays'
    call allocate_arrays
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_kperp2'
    call init_kperp2
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_vperp2'
    call init_vperp2
    if (debug) write (*,*) 'dist_fn::init_dist_fn::init_bessel'
    call init_bessel

  end subroutine init_dist_fn

  subroutine read_parameters

    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast

    implicit none

    logical :: dfexist

    type (text_option), dimension (6), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg) /)
    character(30) :: adiabatic_option
            
    namelist /dist_fn_knobs/ adiabatic_option

    integer :: ierr, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       adiabatic_option = 'default'

       in_file = input_unit_exist("dist_fn_knobs", dfexist)
       if (dfexist) read (unit=in_file, nml=dist_fn_knobs)

       ierr = error_unit()
       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")
    end if

    call broadcast (adiabatic_option_switch)

  end subroutine read_parameters 

  subroutine init_kperp2

    use dist_fn_arrays, only: kperp2, dkperp2dr
    use stella_geometry, only: gds2, gds21, gds22
    use stella_geometry, only: dgds2dr, dgds21dr
    use stella_geometry, only: dgds22dr
    use stella_geometry, only: geo_surf, q_as_x
    use zgrid, only: nzgrid
    use kt_grids, only: naky, nakx, theta0
    use kt_grids, only: akx, aky
    use kt_grids, only: zonal_mode
    use kt_grids, only: nalpha


    implicit none

    integer :: iky, ikx

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(naky,nakx,nalpha,-nzgrid:nzgrid))
    allocate (dkperp2dr(naky,nakx,nalpha,-nzgrid:nzgrid))
    do iky = 1, naky
       if (zonal_mode(iky)) then
          do ikx = 1, nakx
             if(q_as_x) then
               kperp2(iky,ikx,:,:) = akx(ikx)*akx(ikx)*gds22
               where (kperp2(iky,ikx,:,:) .gt. epsilon(0.0))
                 dkperp2dr(iky,ikx,:,:) = akx(ikx)*akx(ikx)*dgds22dr/kperp2(iky,ikx,:,:)
               elsewhere
                 dkperp2dr(iky,ikx,:,:) = 0.0
               endwhere
             else
               kperp2(iky,ikx,:,:) = akx(ikx)*akx(ikx)*gds22/(geo_surf%shat**2)
               where (kperp2(iky,ikx,:,:) .gt. epsilon(0.0))
                 dkperp2dr(iky,ikx,:,:) = akx(ikx)*akx(ikx)*dgds22dr/(geo_surf%shat**2*kperp2(iky,ikx,:,:))
               elsewhere
                 dkperp2dr(iky,ikx,:,:) = 0.0
               endwhere
             endif
          end do
       else
          do ikx = 1, nakx
             kperp2(iky,ikx,:,:) = aky(iky)*aky(iky) &
                  *(gds2 + 2.0*theta0(iky,ikx)*gds21 &
                  + theta0(iky,ikx)*theta0(iky,ikx)*gds22)
             dkperp2dr(iky,ikx,:,:) = aky(iky)*aky(iky) &
                  *(dgds2dr + 2.0*theta0(iky,ikx)*dgds21dr &
                  + theta0(iky,ikx)*theta0(iky,ikx)*dgds22dr)
             dkperp2dr(iky,ikx,:,:)=dkperp2dr(iky,ikx,:,:)/kperp2(iky,ikx,:,:)
             if(any(kperp2(iky,ikx,:,:) .lt. epsilon(0.))) dkperp2dr(iky,ikx,:,:) = 0.
          end do
       end if
    end do
    
    call enforce_single_valued_kperp2

!   filename=trim(run_name)//".kperp2"
!   open (1232,file=trim(filename),status='unknown')
!   ia=1
!   do iz=-nzgrid,nzgrid
!     do ikx=1, naky
!       do iky = 1, 1
!         write(1232,'(2e15.8)') kperp2(iky,ikx,ia,iz), dkperp2dr(iky,ikx,ia,iz)
!       enddo
!     enddo
!   enddo
!   close (1232)

    
  end subroutine init_kperp2

  subroutine enforce_single_valued_kperp2

    use dist_fn_arrays, only: kperp2
    use kt_grids, only: naky, nalpha
    use zgrid, only: nzgrid
    use extended_zgrid, only: neigen, nsegments, ikxmod

    implicit none

    integer :: iky, ie, iseg
    real, dimension (:), allocatable :: tmp

    allocate (tmp(nalpha)) ; tmp = 0.0

    do iky = 1, naky
       do ie = 1, neigen(iky)
          if (nsegments(ie,iky) > 1) then
             do iseg = 2, nsegments(ie,iky)
                tmp = 0.5*(kperp2(iky,ikxmod(iseg-1,ie,iky),:,nzgrid) + kperp2(iky,ikxmod(iseg,ie,iky),:,-nzgrid))
                kperp2(iky,ikxmod(iseg,ie,iky),:,-nzgrid) = tmp
                kperp2(iky,ikxmod(iseg-1,ie,iky),:,nzgrid) = tmp
             end do
          end if
       end do
    end do

    deallocate (tmp)

  end subroutine enforce_single_valued_kperp2

  subroutine allocate_arrays

    use stella_layouts, only: kxkyz_lo, vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use vpamu_grids, only: nvpa, nmu
    use dist_fn_arrays, only: gnew, gold
    use dist_fn_arrays, only: gvmu

    implicit none

    if (.not.allocated(gnew)) &
         allocate (gnew(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gnew = 0.
    if (.not.allocated(gold)) &
         allocate (gold(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    gold = 0.
    if (.not.allocated(gvmu)) &
         allocate (gvmu(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    gvmu = 0.


  end subroutine allocate_arrays

  subroutine init_vperp2

    use stella_geometry, only: bmag
    use zgrid, only: nzgrid
    use vpamu_grids, only: vperp2
    use vpamu_grids, only: nmu, mu
    use kt_grids, only: nalpha

    implicit none

    integer :: imu
    
    if (vp2init) return
    vp2init = .true.

    if (.not.allocated(vperp2)) allocate (vperp2(nalpha,-nzgrid:nzgrid,nmu)) ; vperp2 = 0.
    
    do imu = 1, nmu
       vperp2(:,:,imu) = 2.0*mu(imu)*bmag
    end do

  end subroutine init_vperp2

  subroutine finish_dist_fn

    use gyro_averages, only: finish_bessel

    implicit none

    call finish_bessel
    call finish_kperp2
    call finish_vperp2
    call deallocate_arrays

    dist_fn_initialized = .false.
    readinit = .false.
    gxyz_initialized = .false.

  end subroutine finish_dist_fn

  subroutine deallocate_arrays

    use dist_fn_arrays, only: gnew, gold, gvmu

    implicit none

    if (allocated(gnew)) deallocate (gnew)
    if (allocated(gold)) deallocate (gold)
    if (allocated(gvmu)) deallocate (gvmu)

  end subroutine deallocate_arrays

  subroutine finish_kperp2

    use dist_fn_arrays, only: kperp2, dkperp2dr

    implicit none

    if (allocated(kperp2)) deallocate (kperp2)
    if (allocated(dkperp2dr)) deallocate (dkperp2dr)

    kp2init = .false.

  end subroutine finish_kperp2

  subroutine finish_vperp2

    use vpamu_grids, only: vperp2

    implicit none

    if (allocated(vperp2)) deallocate (vperp2)

    vp2init = .false.
    
  end subroutine finish_vperp2

end module dist_fn
