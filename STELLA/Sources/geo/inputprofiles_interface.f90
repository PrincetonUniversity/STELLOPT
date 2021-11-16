module inputprofiles_interface

  implicit none

  public :: read_inputprof_geo, read_inputprof_spec

  private

  integer :: n_exp

  real, dimension (:), allocatable :: rhotor, rmin, rmaj_in, qinp, kappa
  real, dimension (:), allocatable :: delta, Te, ne, z_eff, omega0
  real, dimension (:), allocatable :: ni, Ti
  real, dimension (:), allocatable :: dr
  real, dimension (:), allocatable :: rhoc, rmaj, shift
  real, dimension (:), allocatable :: shat, d2qdr2
  real, dimension (:), allocatable :: tri, triprim
  real, dimension (:), allocatable :: kapprim
  real, dimension (:), allocatable :: neprim, Teprim
  real, dimension (:), allocatable :: niprim, Tiprim
  real, dimension (:), allocatable :: nedbprim, Tedbprim
  real, dimension (:), allocatable :: nidbprim, Tidbprim
  real, dimension (:), allocatable :: betaprim, betadbprim
  real, dimension (:), allocatable :: psitor
  real, dimension (:), allocatable :: drhotordrho, dpsitordrho, d2psitordrho2
  real, dimension (:), allocatable :: pres_tot
  real, dimension (:), allocatable :: loglam

  real :: bref, omega_ref, rho_ref
  real :: bunit

contains

  subroutine read_inputprof_geo (surf)

    use constants, only: pi
    use common_types, only: flux_surface_type
    use finite_differences, only: fd3pt, d2_3pt
    use splines, only: geo_spline
    use millerlocal, only: local

    implicit none

    type (flux_surface_type), intent (in out) :: surf

    integer :: in_unit = 101
    character (10) :: dum
    character (500) :: line

    real :: bt_exp
    real :: arho_exp

    real :: aref
    real :: mu0

    integer :: ir

    open (unit=in_unit, file='input.profiles', status='old', action='read')

    ! read in header and ignore
    read (in_unit,*) line
    read (in_unit,*) line
    read (in_unit,*) line
    read (in_unit,*) line
    read (in_unit,*) line
    
    read (in_unit,*) dum, n_exp
    read (in_unit,*) dum, bt_exp
    read (in_unit,*) dum, arho_exp
    
    call allocate_arrays_geo

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) rhotor(ir), rmin(ir), rmaj_in(ir), qinp(ir), kappa(ir)
    end do

    ! aref is stella reference length (device minor radius)
    aref = rmin(n_exp)

    rhoc = rmin/aref
    rmaj = rmaj_in/aref

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) delta(ir), Te(ir), ne(ir), line
    end do
    
    ! stella redefines delat to be ArcSin(delta_miller)
    tri = asin(delta)

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line

    ! read in some stuff we don't need to use at the moment
    do ir = 1, n_exp
       read (in_unit,*) line
    end do
    
    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    ! read in some stuff we don't need to use at the moment
    do ir = 1, n_exp
       read (in_unit,*) line
    end do
    
    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) ni(ir), line
    end do

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) Ti(ir), line
    end do
    
    close (in_unit)
    
    dr = rhoc(2:)-rhoc(:n_exp-1)
    
    ! obtain s_hat
    call fd3pt (qinp, shat, dr)
    shat = rhoc*shat/qinp

    ! obtain d2q/dr2
    call d2_3pt (qinp, d2qdr2, dr)

    ! obtain d (kappa) / drho
    call fd3pt (kappa, kapprim, dr)

    ! obtain d (ArcSin(delta_miller)) / drho
    call fd3pt (tri, triprim, dr)

    ! obtain dR/drho
    call fd3pt (rmaj, shift, dr)

    pres_tot = ne*Te+ni*Ti
    call fd3pt (pres_tot, betaprim, dr)
    call d2_3pt (pres_tot, betadbprim, dr)

    mu0 = 4.*pi*1.e-7
    betaprim = -mu0*betaprim*1.6022e3/bt_exp**2
    betadbprim = -mu0*betadbprim*1.6022e3/bt_exp**2

    ! this is psi_toroidal/(Bref aref**2)/2pi
    ! assumption here is that Bref = BT_EXP
    psitor = 0.5*rhotor*rhotor*(arho_exp/aref)**2

    ! get drhotordr
    call fd3pt (rhotor, drhotordrho, dr)
    call fd3pt (psitor, dpsitordrho, dr)
    call d2_3pt (psitor, d2psitordrho2, dr)

    ! this sets the reference B-field to be bt_exp
!    dpsitordrho = rhotor*drhotordrho*(arho_exp/aref)**2

    ! next need to pick out the correct flux surface
    ! and assign various local% values
    call geo_spline (rhoc, rmaj, local%rhoc, local%rmaj)
    call geo_spline (rhoc, dpsitordrho, local%rhoc, local%dpsitordrho)
    call geo_spline (rhoc, d2psitordrho2, local%rhoc, local%d2psitordrho2)
    call geo_spline (rhoc, shift, local%rhoc, local%shift)
    call geo_spline (rhoc, kappa, local%rhoc, local%kappa)
    call geo_spline (rhoc, kapprim, local%rhoc, local%kapprim)
    call geo_spline (rhoc, qinp, local%rhoc, local%qinp)
    call geo_spline (rhoc, shat, local%rhoc, local%shat)
    call geo_spline (rhoc, d2qdr2, local%rhoc, local%d2qdr2)
    call geo_spline (rhoc, tri, local%rhoc, local%tri)
    call geo_spline (rhoc, triprim, local%rhoc, local%triprim)
    call geo_spline (rhoc, betaprim, local%rhoc, local%betaprim)
    call geo_spline (rhoc, betadbprim, local%rhoc, local%betadbprim)
    call geo_spline (rhoc, rhotor, local%rhoc, local%rhotor)
    call geo_spline (rhoc, drhotordrho, local%rhoc, local%drhotordrho)
    local%psitor_lcfs = psitor(n_exp)
    local%zed0_fac = 1.0

    surf = local

    call deallocate_arrays_geo

  end subroutine read_inputprof_geo

  subroutine read_inputprof_spec (nspec, spec)

    use mp, only: mp_abort
    use finite_differences, only: fd3pt, d2_3pt
    use splines, only: geo_spline
    use common_types, only: spec_type
    use millerlocal, only: local
    use physics_parameters, only: vnew_ref, rhostar

    implicit none

    integer, intent (in) :: nspec
    type (spec_type), dimension (:), intent (in out) :: spec

    integer, parameter :: electron_species = 2

    integer :: in_unit = 102
    character (10) :: dum
    character (500) :: line

    real :: bt_exp
    real :: arho_exp

    real :: aref
    real :: mref, nref, tref
    real :: vtref, local_loglam

    integer :: ir, is

    open (unit=in_unit, file='input.profiles', status='old', action='read')

    ! read in header and ignore
    read (in_unit,*) line
    read (in_unit,*) line
    read (in_unit,*) line
    read (in_unit,*) line
    read (in_unit,*) line
    
    read (in_unit,*) dum, n_exp
    read (in_unit,*) dum, bt_exp
    read (in_unit,*) dum, arho_exp
    
    call allocate_arrays_spec

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) rhotor(ir), rmin(ir), line
    end do

    ! aref is stella reference length (device minor radius)
    aref = rmin(n_exp)

    rhoc = rmin/aref

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) dum, Te(ir), ne(ir), z_eff(ir), omega0(ir)
    end do
    
    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line

    ! read in some stuff we don't need to use at the moment
    do ir = 1, n_exp
       read (in_unit,*) line
    end do
    
    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    ! read in some stuff we don't need to use at the moment
    do ir = 1, n_exp
       read (in_unit,*) line
    end do
    
    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) ni(ir), line
    end do

    ! read in more headers
    read (in_unit,*) line
    read (in_unit,*) line
    
    do ir = 1, n_exp
       read (in_unit,*) Ti(ir), line
    end do
    
    close (in_unit)
    
    dr = rhoc(2:)-rhoc(:n_exp-1)
    
    ! obtain -d ln(ne) / drho
    call fd3pt (ne, neprim, dr)
    call d2_3pt (ne, nedbprim, dr)
    neprim = -neprim/ne
    nedbprim = nedbprim/ne

    ! obtain -d ln(Te) / drho
    call fd3pt (Te, Teprim, dr)
    call d2_3pt (Te, Tedbprim, dr)
    Teprim = -Teprim/Te
    Tedbprim = Tedbprim/Te

    ! obtain -d ln(ni) / drho
    call fd3pt (ni, niprim, dr)
    call d2_3pt (ni, nidbprim, dr)
    niprim = -niprim/ni
    nidbprim = nidbprim/ni

    ! obtain -d ln(Ti) / drho
    call fd3pt (Ti, Tiprim, dr)
    call d2_3pt (Ti, Tidbprim, dr)
    Tiprim = -Tiprim/Ti
    Tidbprim = Tidbprim/Ti

    ! next need to pick out the correct flux surface
    ! and assign various local% values

    ! choose first species as reference species
    is = 1
    spec(is)%dens = 1.0
    spec(is)%temp = 1.0
    ! get reference density and temperature
    if (spec(is)%type == electron_species) then
       call geo_spline (rhoc, Te, local%rhoc, tref)
       call geo_spline (rhoc, ne, local%rhoc, nref)
       ! only needed for collisions -- reference mass in units of proton mass
       ! and is only included for clarity -- factors of mref cancel in the end
       mref = 5.44625e-4
    else
       call geo_spline (rhoc, Ti, local%rhoc, tref)
       call geo_spline (rhoc, ni, local%rhoc, nref)
       ! only needed for collisions -- reference mass in units of proton mass
       ! and is only included for clarity -- factors of mref cancel in the end
       mref = 2.0
    end if
    ! next get the normalized density and temperature for all other species
    if (nspec == 2) then
       do is = 2, nspec
          if (spec(is)%type == electron_species) then
             call geo_spline (rhoc, Te/tref, local%rhoc, spec(is)%temp)
             call geo_spline (rhoc, ne/nref, local%rhoc, spec(is)%dens)
          else
             call geo_spline (rhoc, Ti/tref, local%rhoc, spec(is)%temp)
             call geo_spline (rhoc, ni/tref, local%rhoc, spec(is)%dens)
          end if
       end do
    else if (nspec > 2) then
       call mp_abort ('multiple ion species not currently supported for input.profiles. aborting.')
    end if

    ! now get the density and temperature gradients at the requested flux surface
    do is = 1, nspec
       if (spec(is)%type == electron_species) then
          call geo_spline (rhoc, Teprim, local%rhoc, spec(is)%tprim)
          call geo_spline (rhoc, Tedbprim, local%rhoc, spec(is)%d2Tdr2)
          call geo_spline (rhoc, neprim, local%rhoc, spec(is)%fprim)
          call geo_spline (rhoc, nedbprim, local%rhoc, spec(is)%d2ndr2)
       else
          call geo_spline (rhoc, Tiprim, local%rhoc, spec(is)%tprim)
          call geo_spline (rhoc, Tidbprim, local%rhoc, spec(is)%d2Tdr2)
          call geo_spline (rhoc, niprim, local%rhoc, spec(is)%fprim)
          call geo_spline (rhoc, nidbprim, local%rhoc, spec(is)%d2ndr2)
       end if
    end do
       
    ! get collisionalities for stella
    loglam = 24.0 - log(1e4*sqrt(0.1*ne)/te)
    call geo_spline (rhoc, loglam, local%rhoc, local_loglam)

    ! vtref = sqrt(2*Tref/mref), with Tref and mref in SI units
    ! so vtref has dimensions of meters / second
    ! note that tref below is T in units of keV and mref is in units of proton mass
    vtref = 3.09497e5*sqrt(2.*tref/mref)
    
    ! reference collision frequency for stella
    ! uses the mass, density and temperature of the reference species,
    ! along with the proton charge in the expression
    ! vnew_ref = (aref/vtref)*(4/3)sqrt(2pi)/(4pi*epsilon_0)**2 * nref * e**4 * loglam / sqrt(mref) / Tref**1.5
    ! note that all quantities are given in SI units and epsilon_0 is permittivity of vacuum
    vnew_ref = 28.5134*(aref/vtref)*local_loglam*nref/(sqrt(mref)*tref**1.5)

    bref = bt_exp
    omega_ref = 9.5791e7*bref/mref

    rho_ref = vtref/omega_ref

    rhostar = rho_ref/aref
    ! I think this is incorrect -- MAB
    bunit = (bref*arho_exp**2)/(local%rhotor/aref/local%rhoc)*local%drhotordrho

!    vnewki = 9.21e-5*aref*zi**4/sqrt(2.)*loglam*ni/ti**2
!    vnewke = 3.95e-3*aref*zi**2*sqrt(0.5*mi)*loglam*ne &
!         /(sqrt(ti)*te**1.5)

    call deallocate_arrays_spec

  end subroutine read_inputprof_spec

  subroutine allocate_arrays_geo
    
    implicit none
    
    allocate (rhotor(n_exp))
    allocate (psitor(n_exp))
    allocate (rmin(n_exp))
    allocate (rmaj_in(n_exp))
    allocate (qinp(n_exp))
    allocate (kappa(n_exp))
    allocate (rhoc(n_exp))
    allocate (rmaj(n_exp))
    allocate (delta(n_exp))
    allocate (tri(n_exp))
    allocate (ne(n_exp))
    allocate (Te(n_exp))
    allocate (ni(n_exp))
    allocate (Ti(n_exp))
    allocate (dr(n_exp-1))
    allocate (shat(n_exp))
    allocate (d2qdr2(n_exp))
    allocate (kapprim(n_exp))
    allocate (triprim(n_exp))
    allocate (shift(n_exp))
    allocate (betaprim(n_exp))
    allocate (betadbprim(n_exp))
    allocate (pres_tot(n_exp))
    allocate (drhotordrho(n_exp))
    allocate (dpsitordrho(n_exp))
    allocate (d2psitordrho2(n_exp))

  end subroutine allocate_arrays_geo
  
  subroutine allocate_arrays_spec
    
    implicit none
    
    allocate (rhotor(n_exp))
    allocate (rmin(n_exp))
    allocate (rhoc(n_exp))
    allocate (Te(n_exp))
    allocate (ne(n_exp))
    allocate (z_eff(n_exp))
    allocate (omega0(n_exp))
    allocate (ni(n_exp))
    allocate (Ti(n_exp))
    allocate (dr(n_exp-1))
    allocate (neprim(n_exp))
    allocate (nedbprim(n_exp))
    allocate (Teprim(n_exp))
    allocate (Tedbprim(n_exp))
    allocate (niprim(n_exp))
    allocate (nidbprim(n_exp))
    allocate (Tiprim(n_exp))
    allocate (Tidbprim(n_exp))
    allocate (loglam(n_exp))
!    allocate (vnewki(n_exp))
!    allocate (vnewke(n_exp))
    
  end subroutine allocate_arrays_spec
  
  subroutine deallocate_arrays_geo
    
    implicit none
    
    deallocate (rhotor)
    deallocate (rmin)
    deallocate (rmaj_in)
    deallocate (rmaj)
    deallocate (qinp)
    deallocate (kappa)
    deallocate (rhoc)
    deallocate (delta)
    deallocate (Te)
    deallocate (ne)
    deallocate (tri)
    deallocate (ni)
    deallocate (Ti)
    deallocate (dr)
    deallocate (shat)
    deallocate (d2qdr2)
    deallocate (kapprim)
    deallocate (triprim)
    deallocate (shift)
    deallocate (betaprim)
    deallocate (betadbprim)
    deallocate (pres_tot)
    deallocate (drhotordrho)
    deallocate (dpsitordrho)
    deallocate (d2psitordrho2)
    
  end subroutine deallocate_arrays_geo

  subroutine deallocate_arrays_spec
    
    implicit none
    
    deallocate (rhotor)
    deallocate (rmin)
    deallocate (rhoc)
    deallocate (Te)
    deallocate (ne)
    deallocate (z_eff)
    deallocate (omega0)
    deallocate (ni)
    deallocate (Ti)
    deallocate (dr)
    deallocate (neprim)
    deallocate (nedbprim)
    deallocate (Teprim)
    deallocate (Tedbprim)
    deallocate (niprim)
    deallocate (nidbprim)
    deallocate (Tiprim)
    deallocate (Tidbprim)
    deallocate (loglam)
!    deallocate (vnewki)
!    deallocate (vnewke)
    
  end subroutine deallocate_arrays_spec

end module inputprofiles_interface
