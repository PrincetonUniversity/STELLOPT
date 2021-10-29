module euterpe_interface

  implicit none

contains

  subroutine read_species_euterpe (nspec, spec)

    use mp, only: mp_abort
    use finite_differences, only: fd3pt, d2_3pt
    use common_types, only: spec_type
    use splines, only: geo_spline
    use physics_parameters, only: vnew_ref, rhostar, tite, nine
    use stella_geometry, only: geo_surf, aref, bref

    implicit none

    integer, intent (in) :: nspec
    type (spec_type), dimension (:), intent (in out) :: spec

    integer, parameter :: electron_species = 2
    integer :: euterpe_unit = 1099, out_unit = 1098
    integer :: nradii_euterpe
    integer :: is, irad
    character (1000) :: euterpe_infile

    real :: mref, tref, nref, local_loglam, vtref, omega_ref, rho_ref
    real :: vnew_ref_euterpe, rhostar_euterpe

    real, dimension (:), allocatable :: dr, rhotor, psitor
    real, dimension (:), allocatable :: ni, Ti
    real, dimension (:), allocatable :: ne, Te
    real, dimension (:), allocatable :: dlnneds, dlnTeds
    real, dimension (:), allocatable :: dlnnids, dlnTids

    real, dimension (:), allocatable :: neprim, Teprim
    real, dimension (:), allocatable :: nedbprim, Tedbprim
    real, dimension (:), allocatable :: niprim, Tiprim
    real, dimension (:), allocatable :: nidbprim, Tidbprim

    real, dimension (:), allocatable :: loglam

    call read_euterpe_parameters (nradii_euterpe, euterpe_infile)

    open (unit=euterpe_unit, file=trim(euterpe_infile), status='old', action='read')

    allocate (dr(nradii_euterpe))
    allocate (psitor(nradii_euterpe))
    allocate (rhotor(nradii_euterpe))
    allocate (Ti(nradii_euterpe))
    allocate (ni(nradii_euterpe))
    allocate (Te(nradii_euterpe))
    allocate (ne(nradii_euterpe))
    allocate (dlnTids(nradii_euterpe))
    allocate (dlnTeds(nradii_euterpe))
    allocate (dlnnids(nradii_euterpe))
    allocate (dlnneds(nradii_euterpe))

    allocate (neprim(nradii_euterpe))
    allocate (nedbprim(nradii_euterpe))
    allocate (niprim(nradii_euterpe))
    allocate (nidbprim(nradii_euterpe))
    allocate (Teprim(nradii_euterpe))
    allocate (Tedbprim(nradii_euterpe))
    allocate (Tiprim(nradii_euterpe))
    allocate (Tidbprim(nradii_euterpe))

    allocate (loglam(nradii_euterpe))
    
    ! column 1 is s=psitor/psitor_LCFS, 2 is dlog(Ti)/ds, 3 is Ti in eV, 4 is dlog(Te)/ds, 5 is Te in eV
    ! 6 is dlog(ni)/ds, 7 is ni in 1/m^3, 8 is dlog(ne)/ds, 9 is ne in 1/m^3
    do irad = 1, nradii_euterpe
       read (euterpe_unit,*) psitor(irad), dlnTids(irad), Ti(irad), dlnTeds(irad), &
            Te(irad), dlnnids(irad), ni(irad), dlnneds(irad), ne(irad)
    end do

    close (euterpe_unit)

    rhotor = sqrt(psitor)
    dr = rhotor(2:)-rhotor(:nradii_euterpe-1)

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
    ! get reference density and temperature at local surface
    if (spec(is)%type == electron_species) then
       call geo_spline (rhotor, Te, geo_surf%rhotor, tref)
       call geo_spline (rhotor, ne, geo_surf%rhotor, nref)
    else
       call geo_spline (rhotor, Ti, geo_surf%rhotor, tref)
       call geo_spline (rhotor, ni, geo_surf%rhotor, nref)
    end if
    ! next get the normalized density and temperature for all other species
    if (nspec == 2) then
       do is = 2, nspec
          if (spec(is)%type == electron_species) then
             call geo_spline (rhotor, Te/tref, geo_surf%rhotor, spec(is)%temp)
             call geo_spline (rhotor, ne/nref, geo_surf%rhotor, spec(is)%dens)
          else
             call geo_spline (rhotor, Ti/tref, geo_surf%rhotor, spec(is)%temp)
             call geo_spline (rhotor, ni/tref, geo_surf%rhotor, spec(is)%dens)
          end if
       end do
    else if (nspec > 2) then
       call mp_abort ('multiple ion species not currently supported for euterpe option. aborting.')
    end if

    ! assume mass in stella input file given in units of proton mass
    mref = spec(1)%mass
    ! convert from eV to keV
    tref = tref*0.001
    ! convert from 1/m^3 to 10^19/m^3
    nref = nref*1.e-19

    ! now get the density and temperature gradients at the requested flux surface
    do is = 1, nspec
       if (spec(is)%type == electron_species) then
          if (spec(is)%tprim < -999.0) call geo_spline (rhotor, Teprim, geo_surf%rhotor, spec(is)%tprim)
          call geo_spline (rhotor, Tedbprim, geo_surf%rhotor, spec(is)%d2Tdr2)
          if (spec(is)%fprim < -999.0) call geo_spline (rhotor, neprim, geo_surf%rhotor, spec(is)%fprim)
          call geo_spline (rhotor, nedbprim, geo_surf%rhotor, spec(is)%d2ndr2)
       else
          if (spec(is)%tprim < -999.0) call geo_spline (rhotor, Tiprim, geo_surf%rhotor, spec(is)%tprim)
          call geo_spline (rhotor, Tidbprim, geo_surf%rhotor, spec(is)%d2Tdr2)
          if (spec(is)%fprim < -999.0) call geo_spline (rhotor, niprim, geo_surf%rhotor, spec(is)%fprim)
          call geo_spline (rhotor, nidbprim, geo_surf%rhotor, spec(is)%d2ndr2)
       end if
    end do

    ! get quantities needed for runs with Boltzmann electrons
    call geo_spline (rhotor, Ti/Te, geo_surf%rhotor, tite)
    call geo_spline (rhotor, ni/ne, geo_surf%rhotor, nine)

    ! get collisionalities for stella
    loglam = 24.0 - log(1e4*sqrt(1.e-20*ne)/(te*0.001))
    call geo_spline (rhotor, loglam, geo_surf%rhotor, local_loglam)

    ! vtref = sqrt(2*Tref/mref), with Tref and mref in SI units
    ! so vtref has dimensions of meters / second
    ! note that tref below is T in units of keV and mref is in units of proton mass
    vtref = 3.09497e5*sqrt(2.*tref/mref)

    ! reference collision frequency for stella
    ! uses the mass, density and temperature of the reference species,
    ! along with the proton charge in the expression
    ! vnew_ref = (aref/vtref)*(4/3)sqrt(2pi)/(4pi*epsilon_0)**2 * nref * e**4 * loglam / sqrt(mref) / Tref**1.5
    ! note that all quantities are given in SI units and epsilon_0 is permittivity of vacuum
    vnew_ref_euterpe = 28.5134*(aref/vtref)*local_loglam*nref/(sqrt(mref)*tref**1.5)

    omega_ref = 9.5791e7*bref/mref
    rho_ref = vtref/omega_ref

    rhostar_euterpe = rho_ref/aref

    if (rhostar < 0.0) rhostar = rhostar_euterpe
    if (vnew_ref < 0.0) vnew_ref = vnew_ref_euterpe

    open (unit=out_unit, file='euterpe.input', status='replace', action='write')

    write (out_unit,*) 'aref: ', aref, 'mref: ', mref, 'nref: ', nref, 'tref: ', tref
    write (out_unit,*) 'loglam: ', local_loglam, 'vnew_ref_euterpe: ', vnew_ref_euterpe, 'vnew_ref: ', vnew_ref
    write (out_unit,*) 'omega_ref: ', omega_ref, 'rho_ref: ', rho_ref
    write (out_unit,*) 'rhostar_euterpe: ', rhostar_euterpe, 'rhostar: ', rhostar
    write (out_unit,*) 'nine: ', nine, 'tite: ', tite, 'fprim: ', spec(1)%fprim, 'tprim: ', spec(1)%tprim
    write (out_unit,*) 'd2ndr2: ', spec(1)%d2ndr2, 'd2Tdr2: ', spec(1)%d2Tdr2

    close (out_unit)

    deallocate (dr, rhotor, psitor)
    deallocate (ni, ne, Ti, Te)
    deallocate (dlnTids, dlnTeds, dlnnids, dlnneds)

    deallocate (niprim, neprim, nidbprim, nedbprim)
    deallocate (Tiprim, Teprim, Tidbprim, Tedbprim)

    deallocate (loglam)
    
  end subroutine read_species_euterpe

  subroutine read_euterpe_parameters (nradii_out, data_file_out)

    use file_utils, only: input_unit_exist

    implicit none
    
    integer :: in_file
    logical :: exist

    integer, intent (out) :: nradii_out
    character (*), intent (out) :: data_file_out

    integer :: nradii
    character (1000) :: data_file

    namelist /euterpe_parameters/ nradii, data_file

    nradii = 1000
    data_file = 'euterpe.dat'

    in_file = input_unit_exist("euterpe_parameters", exist)
    if (exist) read (unit=in_file, nml=euterpe_parameters)

    nradii_out = nradii
    data_file_out = data_file

  end subroutine read_euterpe_parameters

end module euterpe_interface
