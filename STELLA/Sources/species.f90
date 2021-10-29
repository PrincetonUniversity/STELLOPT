module species

  use common_types, only: spec_type

  implicit none

  public :: init_species, finish_species
  public :: read_species_knobs
  public :: reinit_species
  public :: communicate_species_multibox
  !public :: init_trin_species
  public :: nspec, spec, pfac
  public :: ion_species, electron_species, slowing_down_species, tracer_species
  public :: has_electron_species, has_slowing_down_species
  public :: ions, electrons, impurity

  private

  integer, parameter :: ion_species = 1
  integer, parameter :: electron_species = 2 ! for collision operator
  integer, parameter :: slowing_down_species = 3 ! slowing-down distn
  integer, parameter :: tracer_species = 4 ! for test particle diffusion studies

  integer :: species_option_switch
  integer, parameter :: species_option_stella = 1
  integer, parameter :: species_option_inputprofs = 2
  integer, parameter :: species_option_euterpe = 3
  integer, parameter :: species_option_multibox = 4

  integer :: nspec
  logical :: read_profile_variation, write_profile_variation
  logical :: ecoll_zeff

  type (spec_type), dimension (:), allocatable :: spec

  integer :: ions, electrons, impurity
  real :: pfac
!  integer :: ntspec_trin
!  real, dimension (:), allocatable :: dens_trin, temp_trin, fprim_trin, tprim_trin, nu_trin

  character (20) :: species_option


  logical :: initialized = .false.

contains

  subroutine init_species

!    use mp, only: trin_flag
    use mp, only: proc0, broadcast
    use physics_parameters, only: vnew_ref, zeff
    use physics_flags, only: include_pressure_variation
    use inputprofiles_interface, only: read_inputprof_spec
    use euterpe_interface, only: read_species_euterpe

    implicit none

    integer :: is, is2

    if (initialized) return
    initialized = .true.

    allocate (spec(nspec))
    if (proc0) then
       select case (species_option_switch)
       case (species_option_stella)
          call read_species_stella
       case (species_option_inputprofs)
          call read_species_stella
          call read_inputprof_spec (nspec, spec)
       case (species_option_euterpe)
          call read_species_stella
          call read_species_euterpe (nspec, spec)
       case (species_option_multibox)
          call read_species_stella
          !this will be called by the central box in stella.f90 after
          !ktgrids is set up as we need to know the radial box size
          call communicate_species_multibox
       end select


       if (ecoll_zeff) then
           ! AVB: only intra-species collisions, account for e-i and e-impurity collisions using zeff
           do is = 1, nspec
              ! initialize nu_ss' = 0 for all s'
              spec(is)%vnew = 0.
             ! FLAG -- only contains self-collisions at the moment
              spec(is)%vnew(is) = vnew_ref*spec(is)%dens*spec(is)%z**4 &
                   / (sqrt(spec(is)%mass)*spec(is)%temp**1.5)
              ! include electron-ion collisions
              if (spec(is)%type == electron_species) then
                 spec(is)%vnew(is) = spec(is)%vnew(is)*(1.+zeff)
              end if
           end do
       else
           print*,'using full inter-species collisions'
           ! AVB: full intra- and inter-species collision frequencies
           do is = 1, nspec
               do is2 = 1, nspec
                   if (spec(is)%type == electron_species) then
                       spec(is)%vnew(is2) = vnew_ref*spec(is2)%dens*spec(is)%z**2*spec(is2)%z**2&
                               / (sqrt(spec(is)%mass)*spec(is)%temp**1.5)
                   else if ((spec(is)%type == ion_species).and.(spec(is2)%type == ion_species)) then
                       spec(is)%vnew(is2) = vnew_ref*spec(is2)%dens*spec(is)%z**2*spec(is2)%z**2&
                               / (sqrt(spec(is)%mass)*spec(is)%temp**1.5)
                   else if ((spec(is)%type == ion_species).and.(spec(is2)%type == electron_species)) then
                       spec(is)%vnew(is2) = vnew_ref*spec(is2)%dens*spec(is)%z**2*spec(is2)%z**2&
                               / (sqrt(spec(is)%mass)*spec(is)%temp**1.5)
                   end if
                end do
           end do
       end if

       call dump_species_input

    end if

    pfac = 1.0
    if(.not.include_pressure_variation) pfac = 0.0

    call broadcast_parameters

!    if (trin_flag) call reinit_species (ntspec_trin, dens_trin, &
!         temp_trin, fprim_trin, tprim_trin, nu_trin)
  end subroutine init_species

  subroutine read_species_knobs

    use mp, only: proc0, job, broadcast
    use file_utils, only: error_unit, input_unit_exist
    use file_utils, only: runtype_option_switch, runtype_multibox
    use physics_flags, only: radial_variation
    use text_options, only: text_option, get_option_value

    implicit none

    integer :: ierr, in_file
    logical :: exist

    namelist /species_knobs/ nspec, species_option, &
                             read_profile_variation, &
                             write_profile_variation, &
                             ecoll_zeff


    type (text_option), dimension (4), parameter :: specopts = (/ &
         text_option('default', species_option_stella), &
         text_option('stella', species_option_stella), &
         text_option('input.profiles', species_option_inputprofs), &
         text_option('euterpe', species_option_euterpe) /)

    if (proc0) then
       nspec = 2
       read_profile_variation  = .false.
       write_profile_variation = .false.
       species_option = 'stella'

       ecoll_zeff = .false.

       in_file = input_unit_exist("species_knobs", exist)
       if (exist) read (unit=in_file, nml=species_knobs)

       ierr = error_unit()
       call get_option_value (species_option, specopts, species_option_switch, &
            ierr, "species_option in species_knobs")

       if (runtype_option_switch.eq.runtype_multibox.and.(job.ne.1).and.radial_variation) then
         !will need to readjust the species parameters in the left/right boxes
         species_option_switch = species_option_multibox
       endif

       if (nspec < 1) then
          ierr = error_unit()
          write (unit=ierr, &
               fmt="('Invalid nspec in species_knobs: ', i5)") nspec
          stop
       end if
    end if
    call broadcast (nspec)
    call broadcast (read_profile_variation)
    call broadcast (write_profile_variation)
    call broadcast (ecoll_zeff)

  end subroutine read_species_knobs

  subroutine read_species_stella

    use file_utils, only: error_unit, get_indexed_namelist_unit
    use text_options, only: text_option, get_option_value
    use stella_geometry, only: geo_surf

    implicit none

    real :: z, mass, dens, temp, tprim, fprim, d2ndr2, d2Tdr2, dr, bess_fac
    integer :: ierr, unit, is
    character(len=128) :: filename

    character(20) :: type
    type (text_option), dimension (9), parameter :: typeopts = (/ &
         text_option('default', ion_species), &
         text_option('ion', ion_species), &
         text_option('electron', electron_species), &
         text_option('e', electron_species), &
         text_option('beam', slowing_down_species), &
         text_option('fast', slowing_down_species), &
         text_option('alpha', slowing_down_species), &
         text_option('slowing-down', slowing_down_species), &
         text_option('trace', tracer_species) /)

    namelist /species_parameters/ z, mass, dens, temp, &
         tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type

    do is = 1, nspec
       call get_indexed_namelist_unit (unit, "species_parameters", is)
       z = 1
       mass = 1.0
       dens = 1.0
       temp = 1.0
       tprim = -999.9
       fprim = -999.9
       d2ndr2 = 0.0
       d2Tdr2 = 0.0
       bess_fac = 1.0
       type = "default"
       read (unit=unit, nml=species_parameters)
       close (unit=unit)

       spec(is)%z = z
       spec(is)%mass = mass
       spec(is)%dens = dens
       spec(is)%temp = temp
       spec(is)%tprim = tprim
       spec(is)%fprim = fprim
       ! this is (1/n_s)*d^2 n_s / drho^2
       spec(is)%d2ndr2 = d2ndr2
       ! this is (1/T_s)*d^2 T_s / drho^2
       spec(is)%d2Tdr2 = d2Tdr2

       spec(is)%dens_psi0 = dens
       spec(is)%temp_psi0 = temp

       spec(is)%bess_fac = bess_fac

       if(write_profile_variation) then
         write (filename,"(A,I1)") "specprof_", is
         open  (1002, file=filename,status='unknown')
         write (1002,'(6e13.5)') dens, temp, fprim, tprim, d2ndr2, d2Tdr2
         close (1002)
       endif
       if(read_profile_variation) then
         write (filename,"(A,I1)") "specprof_", is
         open  (1002, file=filename,status='unknown')
         read  (1002,'(6e13.5)') dens, temp, fprim, tprim, d2ndr2, d2Tdr2
         close (1002)

         dr = geo_surf%rhoc - geo_surf%rhoc_psi0
         spec(is)%dens = dens*(1.0 - dr*fprim)! + 0.5*dr**2*d2ndr2)
         spec(is)%temp = temp*(1.0 - dr*tprim)! + 0.5*dr**2*d2Tdr2)
         spec(is)%fprim = (fprim - dr*d2ndr2)*(dens/spec(is)%dens)
         spec(is)%tprim = (tprim - dr*d2Tdr2)*(temp/spec(is)%temp)
         !spec(is)%dens = 1.0
         !spec(is)%temp = 1.0
       endif

       ierr = error_unit()
       call get_option_value (type, typeopts, spec(is)%type, ierr, "type in species_parameters_x")


    end do

  end subroutine read_species_stella

  subroutine broadcast_parameters

    use mp, only: broadcast

    implicit none

    integer :: is

    do is = 1, nspec
       call broadcast (spec(is)%z)
       call broadcast (spec(is)%mass)
       call broadcast (spec(is)%dens)
       call broadcast (spec(is)%temp)
       call broadcast (spec(is)%tprim)
       call broadcast (spec(is)%fprim)
       call broadcast (spec(is)%vnew)
       call broadcast (spec(is)%d2ndr2)
       call broadcast (spec(is)%d2Tdr2)
       call broadcast (spec(is)%dens_psi0)
       call broadcast (spec(is)%temp_psi0)
       call broadcast (spec(is)%bess_fac)
       call broadcast (spec(is)%type)

       spec(is)%stm  = sqrt(spec(is)%temp/spec(is)%mass)
       spec(is)%zstm = spec(is)%z/sqrt(spec(is)%temp*spec(is)%mass)
       spec(is)%tz   = spec(is)%temp/spec(is)%z
       spec(is)%zt   = spec(is)%z/spec(is)%temp
       spec(is)%smz  = abs(sqrt(spec(is)%temp*spec(is)%mass)/spec(is)%z)

       spec(is)%stm_psi0  = sqrt(spec(is)%temp_psi0/spec(is)%mass)
       spec(is)%zstm_psi0 = spec(is)%z/sqrt(spec(is)%temp_psi0*spec(is)%mass)
       spec(is)%tz_psi0   = spec(is)%temp_psi0/spec(is)%z
       spec(is)%zt_psi0   = spec(is)%z/spec(is)%temp_psi0
       spec(is)%smz_psi0  = abs(sqrt(spec(is)%temp_psi0*spec(is)%mass)/spec(is)%z)
    end do

  end subroutine broadcast_parameters

  pure function has_electron_species (spec)
    use common_types, only: spec_type
    implicit none
    type (spec_type), dimension (:), intent (in) :: spec
    logical :: has_electron_species
    has_electron_species = any(spec%type == electron_species)
  end function has_electron_species

  pure function has_slowing_down_species (spec)
    use common_types, only: spec_type
    implicit none
    type (spec_type), dimension (:), intent (in) :: spec
    logical :: has_slowing_down_species
    has_slowing_down_species = any(spec%type == slowing_down_species)
  end function has_slowing_down_species

  subroutine finish_species

    implicit none

    deallocate (spec)

    initialized = .false.

  end subroutine finish_species

  subroutine reinit_species (ntspec, dens, temp, fprim, tprim, bess_fac)

   use mp, only: broadcast, proc0

     implicit none

     integer, intent (in) :: ntspec
     real, dimension (:), intent (in) :: dens, fprim, temp, tprim, bess_fac


     integer :: is
     logical, save :: first = .true.

     if (first) then
        if (nspec == 1) then
            ions = 1
            electrons = 0
            impurity = 0
         else
            ! if 2 or more species in GS2 calculation, figure out which is main ion
            ! and which is electron via mass (main ion mass assumed to be one)
            do is = 1, nspec
               if (abs(spec(is)%mass-1.0) <= epsilon(0.0)) then
                  ions = is
               else if (spec(is)%mass < 0.3) then
                  ! for electrons, assuming electrons are at least a factor of 3 less massive
                  ! than main ion and other ions are no less than 30% the mass of the main ion
                  electrons = is
               else if (spec(is)%mass > 1.0 + epsilon(0.0)) then
                  impurity = is
               else
                  if (proc0) write (*,*) &
                       "Error: TRINITY requires the main ions to have mass 1", &
                       "and the secondary ions to be impurities (mass > 1)"
                  stop
               end if
            end do
         end if
         first = .false.
      end if

      if (proc0) then

         nspec = ntspec

         ! Species are passed in following order: main ion, electron, impurity (if present)
         if (nspec == 1) then
            spec(1)%dens = dens(1)
            spec(1)%temp = temp(1)
            spec(1)%fprim = fprim(1)
            spec(1)%tprim = tprim(1)
            spec(1)%bess_fac = bess_fac(1)
         else
            spec(ions)%dens = dens(1)
            spec(ions)%temp = temp(1)
            spec(ions)%fprim = fprim(1)
            spec(ions)%tprim = tprim(1)
            spec(ions)%bess_fac = bess_fac(1)

            spec(electrons)%dens = dens(2)
            spec(electrons)%temp = temp(2)
            spec(electrons)%fprim = fprim(2)
            spec(electrons)%tprim = tprim(2)
            spec(electrons)%bess_fac = bess_fac(2)

            if (nspec > 2) then
               spec(impurity)%dens = dens(3)
               spec(impurity)%temp = temp(3)
               spec(impurity)%fprim = fprim(3)
               spec(impurity)%tprim = tprim(3)
               spec(impurity)%bess_fac = bess_fac(3)
            end if
         end if

         do is = 1, nspec
            spec(is)%stm = sqrt(spec(is)%temp/spec(is)%mass)
            spec(is)%zstm = spec(is)%z/sqrt(spec(is)%temp*spec(is)%mass)
            spec(is)%tz = spec(is)%temp/spec(is)%z
            spec(is)%zt = spec(is)%z/spec(is)%temp
            spec(is)%smz = abs(sqrt(spec(is)%temp*spec(is)%mass)/spec(is)%z)

  !          write (*,100) 'reinit_species', rhoc_ms, spec(is)%temp, spec(is)%fprim, &
  !               spec(is)%tprim, spec(is)%vnewk, real(is)
         end do

         call dump_species_input

      end if

  !100 format (a15,9(1x,1pg18.11))

      call broadcast (nspec)

      do is = 1, nspec
         call broadcast (spec(is)%dens)
         call broadcast (spec(is)%temp)
         call broadcast (spec(is)%fprim)
         call broadcast (spec(is)%tprim)
         call broadcast (spec(is)%bess_fac)
         call broadcast (spec(is)%stm)
         call broadcast (spec(is)%zstm)
         call broadcast (spec(is)%tz)
         call broadcast (spec(is)%zt)
         call broadcast (spec(is)%smz)
      end do

    end subroutine reinit_species


    subroutine communicate_species_multibox(dr_m,dr_p)
      use job_manage, only: njobs
      use mp, only: job, scope, mp_abort,  &
                  crossdomprocs, subprocs,  &
                  send, receive

      implicit none

      real, optional, intent (in) :: dr_m, dr_p

      real, dimension (:), allocatable :: dens, ldens, ltemp, lfprim, ltprim
      real, dimension (:), allocatable :: temp, rdens, rtemp, rfprim, rtprim

      integer :: i

      allocate(dens(nspec))
      allocate(temp(nspec))
      allocate(ldens(nspec))
      allocate(ltemp(nspec))
      allocate(lfprim(nspec))
      allocate(ltprim(nspec))
      allocate(rdens(nspec))
      allocate(rtemp(nspec))
      allocate(rfprim(nspec))
      allocate(rtprim(nspec))


      if(job == 1) then

        ! recall that fprim and tprim are the negative gradients
        ldens  =  spec%dens*(1.0 - dr_m*spec%fprim)! + 0.5*dr_m**2*spec%d2ndr2)
        ltemp  =  spec%temp*(1.0 - dr_m*spec%tprim)! + 0.5*dr_m**2*spec%d2Tdr2)
        lfprim = (spec%fprim - dr_m*spec%d2ndr2)*(spec%dens/ldens)
        ltprim = (spec%tprim - dr_m*spec%d2Tdr2)*(spec%temp/ltemp)

        rdens  =  spec%dens*(1.0 - dr_p*spec%fprim)! + 0.5*dr_p**2*spec%d2ndr2)
        rtemp  =  spec%temp*(1.0 - dr_p*spec%tprim)! + 0.5*dr_p**2*spec%d2Tdr2)
        rfprim = (spec%fprim - dr_p*spec%d2ndr2)*(spec%dens/rdens)
        rtprim = (spec%tprim - dr_p*spec%d2Tdr2)*(spec%temp/rtemp)

        do i=1,nspec
          if(ldens(i) < 0 .or. ltemp(i) < 0 .or. &
            rdens(i) < 0 .or. rtemp(i) < 0) then
            call mp_abort('Negative n/T encountered. Try reducing rhostar.')
          endif
        enddo
      endif

      call scope(crossdomprocs)

      if(job==1) then
        call send(ldens    ,0,120)
        call send(ltemp    ,0,121)
        call send(lfprim   ,0,122)
        call send(ltprim   ,0,123)
        call send(spec%dens,0,124)
        call send(spec%temp,0,125)
        call send(rdens    ,njobs-1,130)
        call send(rtemp    ,njobs-1,131)
        call send(rfprim   ,njobs-1,132)
        call send(rtprim   ,njobs-1,133)
        call send(spec%dens,njobs-1,134)
        call send(spec%temp,njobs-1,135)
      elseif(job == 0) then
        call receive(ldens, 1,120)
        call receive(ltemp, 1,121)
        call receive(lfprim,1,122)
        call receive(ltprim,1,123)
        call receive(dens,  1,124)
        call receive(temp,  1,125)
        spec%dens       = ldens
        spec%temp       = ltemp
        spec%fprim      = lfprim
        spec%tprim      = ltprim
        spec%dens_psi0  = dens
        spec%temp_psi0  = temp
      elseif(job== njobs-1) then
        call receive(rdens, 1,130)
        call receive(rtemp, 1,131)
        call receive(rfprim,1,132)
        call receive(rtprim,1,133)
        call receive(dens,  1,134)
        call receive(temp,  1,135)
        spec%dens  = rdens
        spec%temp  = rtemp
        spec%fprim = rfprim
        spec%tprim = rtprim
        spec%dens_psi0  = dens
        spec%temp_psi0  = temp
      endif

      call scope(subprocs)

      deallocate(dens)
      deallocate(temp)
      deallocate(ldens)
      deallocate(ltemp)
      deallocate(lfprim)
      deallocate(ltprim)
      deallocate(rdens)
      deallocate(rtemp)
      deallocate(rfprim)
      deallocate(rtprim)

    end subroutine communicate_species_multibox

    subroutine dump_species_input

      use file_utils, only: run_name

      implicit none

      integer :: is
      character (300) :: filename

      filename = trim(trim(run_name)//'.species.input')
      open (1003,file=filename,status='unknown')
      write (1003,'(9a12,a9)') '#1.z', '2.mass', '3.dens', &
            '4.temp', '5.tprim','6.fprim', '7.vnewss',  &
           '8.dens_psi0', '9.temp_psi0', '11.type'
      do is = 1, nspec
        write (1003,'(9e12.4,i9)') spec(is)%z, spec(is)%mass, &
               spec(is)%dens, spec(is)%temp, spec(is)%tprim, &
               spec(is)%fprim, spec(is)%vnew(is), &
               spec(is)%dens_psi0, spec(is)%temp_psi0, &
               spec(is)%type
      end do
      close (1003)

    end subroutine dump_species_input

!   subroutine init_trin_species (ntspec_in, dens_in, temp_in, fprim_in, tprim_in, nu_in)

!     implicit none

!     integer, intent (in) :: ntspec_in
!     real, dimension (:), intent (in) :: dens_in, fprim_in, temp_in, tprim_in, nu_in

!     if (.not. allocated(temp_trin)) then
!        allocate (dens_trin(size(dens_in)))
!        allocate (fprim_trin(size(fprim_in)))
!        allocate (temp_trin(size(temp_in)))
!        allocate (tprim_trin(size(tprim_in)))
!        allocate (nu_trin(size(nu_in)))
!     end if

!     ntspec_trin = ntspec_in
!     dens_trin = dens_in
!     temp_trin = temp_in
!     fprim_trin = fprim_in
!     tprim_trin = tprim_in
!     nu_trin = nu_in

!   end subroutine init_trin_species

end module species
