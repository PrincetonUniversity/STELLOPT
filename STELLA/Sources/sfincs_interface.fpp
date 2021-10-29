module sfincs_interface

  implicit none
  
  public :: get_neo_from_sfincs

  private

# ifdef USE_SFINCS
  integer :: nproc_sfincs
  integer :: irad_min, irad_max
  real :: Er_window
  logical :: includeXDotTerm
  logical :: includeElectricFieldTermInXiDot
!  logical :: includeRadialExBDrive
  integer :: magneticDriftScheme
  logical :: includePhi1
  logical :: includePhi1InKineticEquation
!  logical :: includePhi1InCollisionOperator
  integer :: geometryScheme
  integer :: VMECRadialOption
  integer :: coordinateSystem
  integer :: inputRadialCoordinate
  integer :: inputRadialCoordinateForGradients
  character (200) :: equilibriumFile
  real :: aHat, psiAHat, Delta
  real :: nu_n, dPhiHatdrN
  integer :: nxi, nx, ntheta, nzeta
  logical :: calculate_radial_electric_field
  logical :: read_sfincs_output_from_file
  logical :: sfincs_finished = .true.
  real, dimension (:), allocatable :: fprim_local, tprim_local
# endif

contains

  subroutine get_neo_from_sfincs (nradii, drho, f_neoclassical, phi_neoclassical, &
       dfneo_dalpha, dphineo_dalpha)

# ifdef USE_SFINCS
    use mp, only: proc0, iproc
    use mp, only: comm_split, comm_free
    use stella_geometry, only: geo_surf
    use species, only: spec, nspec
# endif
    use mp, only: mp_abort
    use zgrid, only: nzgrid
    
    implicit none

    integer, intent (in) :: nradii
    real, intent (in) :: drho
    real, dimension (:,-nzgrid:,:,:,:,-nradii/2:), intent (out) :: f_neoclassical
    real, dimension (:,-nzgrid:,-nradii/2:), intent (out) :: phi_neoclassical
    real, dimension (:,-nzgrid:,:,:,:), intent (out) :: dfneo_dalpha
    real, dimension (:,-nzgrid:), intent (out) :: dphineo_dalpha

# ifdef USE_SFINCS
    integer :: sfincs_comm
    integer :: color, ierr
    integer :: irad
    real :: dPhiHatdrN_best_guess
    logical :: Er_converged
    integer :: nsfincs_calls

    if (.not.allocated(fprim_local)) allocate (fprim_local(nspec))
    if (.not.allocated(tprim_local)) allocate (tprim_local(nspec))

    if (proc0) call read_sfincs_parameters (nradii)
    call broadcast_sfincs_parameters
    if (iproc < nproc_sfincs) then
       color = 0
    else
       color = 1
    end if
    call comm_split (color, sfincs_comm, ierr)
    if (iproc < nproc_sfincs) then
       do irad = irad_min, irad_max
          ! get local values of -dlog(ns)/drho and -dlog(Ts)/drho
          ! using dlog(n)/drho = dlog(n0)/drho + delrho*d/drho(dlog(n)/drho)
          fprim_local = 1.0/geo_surf%drhotordrho*(spec%fprim &
               + irad*drho*(spec%fprim**2-spec%d2ndr2)/geo_surf%drhotordrho)
          tprim_local = 1.0/geo_surf%drhotordrho*(spec%tprim &
               + irad*drho*(spec%tprim**2-spec%d2Tdr2)/geo_surf%drhotordrho)
!          fprim_local = 1.0/geo_surf%drhotordrho*(spec%fprim - irad*drho*spec%d2ndr2)
!          tprim_local = 1.0/geo_surf%drhotordrho*(spec%tprim - irad*drho*spec%d2Tdr2)
          if (calculate_radial_electric_field) then

             ! get best guess at radial electric field
             ! using force balance with radial pressure gradient
             if (dPhiHatdrN > -9999.0) then
                dPhiHatdrN_best_guess = dPhiHatdrN
             else
                dPhiHatdrN_best_guess = fprim_local(1)+tprim_local(1)
             end if
             call iterate_sfincs_until_electric_field_converged (sfincs_comm, &
                  irad, drho, irad_max, dPhiHatdrN_best_guess, &
                  Er_converged, nsfincs_calls)

             if (proc0) then
                write (*,*)
                write (*,*) 'At irad= ', irad, 'Er_converged= ', Er_converged, &
                  'nsfincs_calls_required= ', nsfincs_calls, 'dPhiHatdrN= ', dPhiHatdrN
                write (*,*)
             end if

             ! write_and_finish_sfincs manipulates sfincs output
             ! to get the neoclassical distribution function and potential
             ! on the stella (zed,alpha,vpa,mu) grid; it then
             ! deallocates sfincs arrays, etc. to make it ready
             ! for running again if need be
             call write_and_finish_sfincs (f_neoclassical(:,:,:,:,:,irad), &
                  phi_neoclassical(:,:,irad), dfneo_dalpha, dphineo_dalpha, irad)
          else
             ! init_and_run_sfincs initializes sfincs,
             ! including passing geometry info if necessary;
             ! and runs sfincs (if requested)
             call init_and_run_sfincs (sfincs_comm, irad, drho, irad_max)
             ! write_and_finish_sfincs manipulates sfincs output
             ! to get the neoclassical distribution function and potential
             ! on the stella (zed,alpha,vpa,mu) grid; it then
             ! deallocates sfincs arrays, etc. to make it ready
             ! for running again if need be
             call write_and_finish_sfincs (f_neoclassical(:,:,:,:,:,irad), &
                  phi_neoclassical(:,:,irad), dfneo_dalpha, dphineo_dalpha, irad)
          end if
       end do
    end if

    call comm_free (sfincs_comm, ierr)

    ! NB: NEED TO CHECK THIS BROADCAST OF SFINCS RESULTS 
    do irad = irad_min, irad_max
       call broadcast_sfincs_output &
            (f_neoclassical(:,:,:,:,:,irad), phi_neoclassical(:,:,irad))
    end do
    call broadcast_sfincs_output (dfneo_dalpha, dphineo_dalpha)

    deallocate (fprim_local, tprim_local)

    if (proc0) then
       write (*,*)
       write (*,*) 'maxval(fneo): ', maxval(f_neoclassical(:,:,:,:,:,irad_min:irad_max)), &
            'maxval(phineo): ', maxval(phi_neoclassical(:,:,irad_min:irad_max))
       write (*,*)
    end if

    if ((irad_min /= -nradii/2) .or. (irad_max /= nradii/2)) &
         call mp_abort ('WARNING: irad_min must equal -nradii/2 and irad_max must equal &
         & nradii/2 to proceed to stella calculation.  aborting.')

# else
    real :: dum
    ! this pointless dum assignment only here to avoid 
    ! annoying warning messages during compilation
    ! about unused variable drho 
    dum = drho
    f_neoclassical = 0. ; phi_neoclassical = 0.
    dfneo_dalpha = 0. ; dphineo_dalpha = 0.
    call mp_abort ("to run with include_neoclassical_terms=.true., &
         & USE_SFINCS must be defined at compilation time.  Aborting.")
# endif

  end subroutine get_neo_from_sfincs

# ifdef USE_SFINCS
  subroutine iterate_sfincs_until_electric_field_converged (sfincs_comm, irad, drho, &
       nrad_max, dPhiHatdrN_best_guess, dphiHatdrN_is_converged, number_of_sfincs_calls_for_convergence)

    use mp, only: proc0, iproc

    implicit none

    integer, intent (in) :: sfincs_comm, irad, nrad_max
    real, intent (in) :: drho, dPhiHatdrN_best_guess
    logical, intent (out) :: dPhiHatdrN_is_converged
    integer, intent (out) :: number_of_sfincs_calls_for_convergence

    integer :: itmax_bracket = 10
    integer :: itmax_root = 10
    real :: window = 0.3
    real :: tol = 0.1

    integer :: it
    real :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,eps
    real :: converged_dPhiHatdrN

    dPhiHatdrN_is_converged = .false.

    a = dPhiHatdrN_best_guess*(1.0-Er_window)
    b = dPhiHatdrN_best_guess*(1.0+Er_window)
    ! initialize sfincs, run it, and return the total charge flux as fa
    call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, a, fa)
    call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, b, fb)
    number_of_sfincs_calls_for_convergence = 2
    do it = 1, itmax_bracket
       eps = epsilon(a)
       if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
          if (proc0) then
             write (*,*)
             write (*,*) 'dPhiHatdrN values ', a, ' and ', b, ' do not bracket root.'
             write (*,*) 'flux at ', a, ' is ', fa, '.'
             write (*,*) 'flux at ', b, ' is ', fb, '.'
          end if
          a = a*(1.0-Er_window)
          b = b*(1.0+Er_window)
          if (proc0) then
             write (*,*) 'Trying again with values ', a, ' and ', b, ' .'
             write (*,*)
          end if
          call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, a, fa)
          call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, b, fb)
          number_of_sfincs_calls_for_convergence = number_of_sfincs_calls_for_convergence + 2

!           ! eliminate the endpoint corresonding to the flux that is furthest from zero in magnitude
!           if (abs(fa) > abs(fb)) then
!              ! keep b as an endpoint and eliminate a
!              a = b ; fa = fb
!              b = a*(1.0+window)
!              if (proc0) then
!                 write (*,*) 'Trying again with values ', a, ' and ', b, ' .'
!                 write (*,*)
!              end if
!              call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, b, fb)
!           else
!              ! keep a as an endpoint and eliminate b
!              b = a ; fb = fa
!              a = b*(1.0-window)
!              if (proc0) then
!                 write (*,*) 'Trying again with values ', a, ' and ', b, ' .'
!                 write (*,*)
!              end if
!              call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, a, fa)
!           end if
       else
          exit
       end if
    end do

    c=b
    fc=fb
    do it = 1, itmax_root
       if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
          c=a
          fc=fa
          d=b-a
          e=d
       end if
       if (abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       end if
       tol1=2.0*eps*abs(b)+0.5*tol
       xm=0.5*(c-b)
       if (abs(xm) <= tol1 .or. fb == 0.0) then
          converged_dPhiHatdrN = b
          dPhiHatdrN_is_converged = .true.
!          number_of_sfincs_calls_for_convergence = it+1
          exit
       end if
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa
          if (a==c) then
             p=2.0*xm*s
             q=1.0-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
             q=(q-1.0)*(r-1.0)*(s-1.0)
          end if
          if (p > 0.0) q=-q
          p=abs(p)
          if (2.0*p < min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          end if
       else
          d=xm
          e=d
       end if
       a=b
       fa=fb
       b=b+merge(d,sign(tol1,xm), abs(d) > tol1)
       call get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, b, fb)
       number_of_sfincs_calls_for_convergence = number_of_sfincs_calls_for_convergence + 1
    end do

  end subroutine iterate_sfincs_until_electric_field_converged

  subroutine get_total_charge_flux (sfincs_comm, irad, drho, nrad_max, &
       dPhiHatdrN_in, total_charge_flux)

    use mp, only: iproc, broadcast_with_comm
    use sfincs_main, only: finish_sfincs
    use globalVariables, only: Zs, particleFlux_vd_psiHat
    use species, only: nspec

    implicit none

    integer, intent (in) :: sfincs_comm, irad, nrad_max
    real, intent (in) :: drho
    real, intent (in) :: dPhiHatdrN_in
    real, intent (out) :: total_charge_flux

    dPhiHatdrN = dPhiHatdrN_in
    call init_and_run_sfincs (sfincs_comm, irad, drho, nrad_max)

    call broadcast_with_comm (Zs, sfincs_comm)
    call broadcast_with_comm (particleFlux_vd_psiHat, sfincs_comm)

    total_charge_flux = sum(Zs(:nspec)*particleFlux_vd_psiHat(:nspec))

  end subroutine get_total_charge_flux

  subroutine init_and_run_sfincs (sfincs_comm, irad, drho, nrad_max)

    use mp, only: proc0, iproc
    use sfincs_main, only: init_sfincs, prepare_sfincs, run_sfincs, finish_sfincs
    use globalVariables, only: Zs, particleFlux_vd_psiHat
    use species, only: nspec

    implicit none

    integer, intent (in) :: sfincs_comm, irad, nrad_max
    real, intent (in) :: drho

    if (.not.sfincs_finished) then
       call finish_sfincs
       sfincs_finished = .true.
    end if
    call init_sfincs (sfincs_comm)
    call pass_inputoptions_to_sfincs (irad*drho)
    call pass_outputoptions_to_sfincs
    call prepare_sfincs
    ! if geometryScheme = 5, then sfincs will read in equilibrium
    ! parameters from vmec file separately
    ! otherwise, assume system is axisymmetric and pass geometry
    ! from stella (miller local equilibrium or similar)
    if (geometryScheme /= 5) call pass_geometry_to_sfincs (irad*drho)
    if (read_sfincs_output_from_file) then
       if (proc0) call read_sfincs_output (irad, nrad_max)
    else
       if (proc0) then
          write (*,*)
          write (*,*) 'Running sfincs at irad= ', irad, ', with dPhiHatdrN= ', dPhiHatdrN
       end if
       call run_sfincs
       if (proc0) then
          write (*,*) 'sfincs finished running.  total charge flux= ', sum(Zs(:nspec)*particleFlux_vd_psiHat(:nspec))
          write (*,*)
       end if
       ! write Phi1Hat and delta_f to file
       ! so we have the option of using it
       ! again without re-running sfincs
       if (proc0) call write_sfincs (irad, nrad_max)
    end if
    sfincs_finished = .false.

  end subroutine init_and_run_sfincs
  
  subroutine write_and_finish_sfincs (fneo, phineo, dfneo, dphineo, irad)

    use mp, only: proc0
    use sfincs_main, only: finish_sfincs
    use zgrid, only: nzgrid

    implicit none

    real, dimension (:,-nzgrid:,:,:,:), intent (out) :: fneo
    real, dimension (:,-nzgrid:), intent (out) :: phineo
    real, dimension (:,-nzgrid:,:,:,:), intent (in out) :: dfneo
    real, dimension (:,-nzgrid:), intent (in out) :: dphineo
    integer, intent (in) :: irad

    if (proc0) then
       ! only need to compute dfneo_dalpha and dphineo_dalpha
       ! for central radius and for stellarator calculation
       if (irad == 0 .and. geometryScheme == 5) then
          call get_sfincs_output (fneo, phineo, dfneo, dphineo)
       else
          call get_sfincs_output (fneo, phineo)
       end if
    end if

    call finish_sfincs

  end subroutine write_and_finish_sfincs

  subroutine read_sfincs_parameters (nradii)

    use constants, only: pi
    use mp, only: nproc
    use file_utils, only: input_unit_exist
    use species, only: nspec
    use physics_parameters, only: rhostar, vnew_ref
    use stella_geometry, only: geo_surf, aref, bref

    implicit none
    
    integer, intent (in) :: nradii

    namelist /sfincs_input/ nproc_sfincs, &
         calculate_radial_electric_field, &
         includeXDotTerm, &
         includeElectricFieldTermInXiDot, &
         irad_min, irad_max, &
         magneticDriftScheme, &
         includePhi1, &
         includePhi1InKineticEquation, &
!         includePhi1InCollisionOperator, &
         geometryScheme, &
         VMECRadialOption, &
         equilibriumFile, &
         coordinateSystem, &
         inputRadialCoordinate, &
         inputRadialCoordinateForGradients, &
         aHat, psiAHat, nu_N, nxi, nx, Delta, &
         dPhiHatdrN, &
         ntheta, nzeta, &
         read_sfincs_output_from_file, Er_window

    logical :: exist
    integer :: in_file

    ! if read_sfincs_output_from_file=.true.,
    ! will try to read in Phi1Hat and delta_f
    ! from pre-saved file named sfincs.output
    ! otherwise, run sfincs to compute these
    ! quantities on the fly
    read_sfincs_output_from_file = .false.
    ! number of processors to use for sfincs calculation
    nproc_sfincs = 1
    ! minimum and maximum radial index (irad=0 corresponds to central radius)
    irad_min = -nradii/2 ; irad_max = nradii/2
    ! if calculate_radial_electric_field, then
    ! will scan in radial electric field to find value
    ! for which ambipolarity is satisfied, and then
    ! use this value to obtain neoclassical fluxes,
    ! distribution function, and potential
    calculate_radial_electric_field = .true.
    ! do not include radial electric field term if set to .false.
    includeXDotTerm = .true.
    includeElectricFieldTermInXiDot = .true.
    ! include v_E . grad r term
!    includeRadialExBDrive = .true.
    ! no poloidal or toroidal magnetic drifts
    magneticDriftScheme = 0
    ! combo of next two variables means
    ! phi1 will be calculated via quasineutrality
    includePhi1 = .true.
    includePhi1InKineticEquation = .false.
!    includePhi1InCollisionOperator = .false.
    ! will be overridden by direct input of geometric quantities
    ! unless geometryScheme = 5 (vmec equilibrium)
    geometryScheme = 1
    ! only relevant if geometryScheme = 5
    ! radial option to use for vmec equilibrium
    ! 0 corresponds to using radial interpolation to get desired surface
    ! 1 corresponds to using nearest surface on VMEC HALF grid
    ! 2 corresponds to using nearest surface on VMEC FULL grid
    ! should not change this unless self-consistently change in the
    ! vmec input namelist
    VMECRadialOption = 0
    ! path of vmec equilibrium file
    equilibriumFile = 'wout_161s1.nc'
    ! seems to be a nonsensical option
    coordinateSystem = 3
    ! option 3 corresponds to using sqrt of toroidal flux
    ! normalized by toroidal flux enclosed by the LCFS
    inputRadialCoordinate = 3
    ! option 3 corresponds to same choice
    ! when calculating gradients of density, temperature, and potential
    inputRadialCoordinateForGradients = 3
    ! corresponds to r_LCFS as reference length in sfincs
    ! only used in sfincs when geometryScheme=1
    aHat = 1.0
    ! psitor_LCFS / (B_ref * a_ref^2)
    psiAHat = geo_surf%psitor_lcfs
    ! Delta is rho* = mref*vt_ref/(e*Bref*aref), with reference
    ! quantities given in SI units
    ! unless geometryScheme = 5, in which case Bref=1T
    ! and aref = 1m (these are hardwired in sfincs)
    ! set negative to allow check later to see if any value given in input file
    Delta = -1.0
    ! nu_n = nu_ref * aref/vt_ref
    ! nu_ref = 4*sqrt(2*pi)*nref*e**4*loglam/(3*sqrt(mref)*Tref**3/2)
    ! (with nref, Tref, and mref in Gaussian units)
    ! set negative to allow check later to see if any value given in input file
    nu_n = -1.0
    ! radial derivative of normalized phi
    dPhiHatdrN = -9999.9
    Er_window = 0.3
    ! number of spectral coefficients in pitch angle
    nxi = 48
    ! number of speeds
    nx = 12
    ! number of poloidal angles
    Ntheta = 65
    ! number of toroidal angles, 1 is appropriate for tokamak
    Nzeta = 1

    in_file = input_unit_exist("sfincs_input", exist)
    if (exist) read (unit=in_file, nml=sfincs_input)

    if (nproc_sfincs > nproc) then
       write (*,*) 'requested number of processors for sfincs is greater &
            & than total processor count.'
       write (*,*) 'allocating ', nproc, ' processors for sfincs.'
    end if

    if (Delta < 0.0) then
       Delta = rhostar
       ! if geometryScheme=5, Bref=1T and aref=1m are hard-wired in sfincs
       ! but these are not the values used in stella to define rhostar
       if (geometryScheme == 5) Delta = rhostar*bref*aref
    end if

    if (nu_n < 0.0) then
       nu_n = vnew_ref*(4./(3.*sqrt(pi)))
       ! if geometryScheme=5, aref=1m is hard-wired in sfincs
       ! but this is not the value used in stella
       if (geometryScheme == 5) nu_n = nu_n/aref
    end if

! FLAG -- NOT YET SURE IF THIS SHOULD BE HERE
!    if (nspec == 1 .and. includePhi1) then
!       write (*,*) 'includePhi1 = .true. is incompatible with a single-species run.'
!       write (*,*) 'forcing includePhi1 = .false.'
!       includePhi1 = .false.
!    end if
    
    ! ensure that ntheta and nzeta are odd for SFINCS
    ntheta = 2*(ntheta/2)+1
    nzeta = 2*(nzeta/2)+1

  end subroutine read_sfincs_parameters

  subroutine broadcast_sfincs_parameters

    use mp, only: broadcast
    use physics_parameters, only: rhostar

    implicit none

    call broadcast (read_sfincs_output_from_file)
    call broadcast (nproc_sfincs)
    call broadcast (irad_min)
    call broadcast (irad_max)
    call broadcast (calculate_radial_electric_field)
    call broadcast (includeXDotTerm)
    call broadcast (includeElectricFieldTermInXiDot)
!    call broadcast (includeRadialExBDrive)
    call broadcast (magneticDriftScheme)
    call broadcast (includePhi1)
    call broadcast (includePhi1InKineticEquation)
!    call broadcast (includePhi1InCollisionOperator)
    call broadcast (geometryScheme)
    call broadcast (VMECRadialOption)
    call broadcast (equilibriumFile)
    call broadcast (coordinateSystem)
    call broadcast (inputRadialCoordinate)
    call broadcast (inputRadialCoordinateForGradients)
    call broadcast (aHat)
    call broadcast (psiAHat)
    call broadcast (Delta)
    call broadcast (nu_N)
    call broadcast (dPhiHatdrN)
    call broadcast (Er_window)
    call broadcast (nxi)
    call broadcast (nx)
    call broadcast (ntheta)
    call broadcast (nzeta)

    write (*,*) 'Delta stella', Delta, rhostar

  end subroutine broadcast_sfincs_parameters

  subroutine pass_inputoptions_to_sfincs (delrho)

    use mp, only: mp_abort
    use stella_geometry, only: geo_surf
    use species, only: spec, nspec
    use zgrid, only: nzed
    use physics_parameters, only: nine, tite
    use globalVariables, only: includeXDotTerm_sfincs => includeXDotTerm
    use globalVariables, only: includeElectricFieldTermInXiDot_sfincs => includeElectricFieldTermInXiDot
!    use globalVariables, only: includeRadialExBDrive_sfincs => includeRadialExBDrive
    use globalVariables, only: magneticDriftScheme_sfincs => magneticDriftScheme
    use globalVariables, only: includePhi1_sfincs => includePhi1
    use globalVariables, only: includePhi1InKineticEquation_sfincs => includePhi1InKineticEquation
!    use globalVariables, only: includePhi1InCollisionOperator_sfincs => includePhi1InCollisionOperator
    use globalVariables, only: geometryScheme_sfincs => geometryScheme
    use globalVariables, only: equilibriumFile_sfincs => equilibriumFile
    use globalVariables, only: VMECRadialOption_sfincs => VMECRadialOption
    use globalVariables, only: coordinateSystem_sfincs => coordinateSystem
    use globalVariables, only: RadialCoordinate => inputRadialCoordinate
    use globalVariables, only: RadialCoordinateForGradients => inputRadialCoordinateForGradients
    use globalVariables, only: rN_wish
    use globalVariables, only: Nspecies, nHats, THats, MHats, Zs
    use globalVariables, only: adiabaticNHat, adiabaticTHat, adiabaticZ
    use globalVariables, only: nxi_sfincs => Nxi
    use globalVariables, only: nx_sfincs => Nx
    use globalVariables, only: ntheta_sfincs => Ntheta
    use globalVariables, only: nzeta_sfincs => Nzeta
    use globalVariables, only: dnHatdrNs, dTHatdrNs
    use globalVariables, only: aHat_sfincs => aHat
    use globalVariables, only: psiAHat_sfincs => psiAHat
    use globalVariables, only: Delta_sfincs => Delta
    use globalVariables, only: nu_n_sfincs => nu_n
    use globalVariables, only: dPhiHatdrN_sfincs => dPhiHatdrN
    use globalVariables, only: withAdiabatic

    implicit none

    real, intent (in) :: delrho

    includeXDotTerm_sfincs = includeXDotTerm
    includeElectricFieldTermInXiDot_sfincs = includeElectricFieldTermInXiDot
!    includeRadialExBDrive_sfincs = includeRadialExBDrive
    magneticDriftScheme_sfincs = magneticDriftScheme
    includePhi1_sfincs = includePhi1
    includePhi1InKineticEquation_sfincs = includePhi1InKineticEquation
!    includePhi1InCollisionOperator_sfincs = includePhi1InCollisionOperator
    geometryScheme_sfincs = geometryScheme
    VMECRadialOption_sfincs = VMECRadialOption
    equilibriumFile_sfincs = trim(equilibriumFile)
    coordinateSystem_sfincs = coordinateSystem
    RadialCoordinate = inputRadialCoordinate
    RadialCoordinateForGradients = inputRadialCoordinateForGradients
    Nspecies = nspec
    nHats(:nspec) = spec%dens*(1.0-delrho*spec%fprim)
    THats(:nspec) = spec%temp*(1.0-delrho*spec%tprim)
    mHats(:nspec) = spec%mass
    Zs(:nspec) = spec%z
    nzeta_sfincs = nzeta
    ntheta_sfincs = ntheta
    nx_sfincs = nx
    nxi_sfincs = nxi
    aHat_sfincs = aHat
    psiAHat_sfincs = psiAHat
    Delta_sfincs = Delta
    nu_n_sfincs = nu_n
    dPhiHatdrN_sfincs = dPhiHatdrN
    if (nspec == 1) then
       withAdiabatic = .true.
       adiabaticNHat = nHats(1)/nine
       adiabaticTHat = THats(1)/tite
       adiabaticZ = -1
    end if

    if (inputRadialCoordinate == 3) then
       rN_wish = geo_surf%rhotor + delrho*geo_surf%drhotordrho
    else
       call mp_abort ('only inputRadialCoordinate=3 currently supported. aborting.')
    end if
    if (inputRadialCoordinateForGradients == 3) then
       ! radial density gradient with respect to rhotor = sqrt(psitor/psitor_LCFS)
       ! normalized by reference density (not species density)
       dnHatdrNs(:nspec) = -spec%dens*fprim_local
       ! radial temperature gradient with respect to rhotor = sqrt(psitor/psitor_LCFS)
       ! normalized by reference tmperatures (not species temperature)
       dTHatdrNs(:nspec) = -spec%temp*tprim_local
    else
       call mp_abort ('only inputRadialCoordinateForGradients=3 currently supported. aborting.')
    end if

  end subroutine pass_inputoptions_to_sfincs

  subroutine pass_outputoptions_to_sfincs
    use export_f, only: export_f_theta_option
    use export_f, only: export_f_zeta_option
    use export_f, only: export_f_xi_option
    use export_f, only: export_f_x_option
    use export_f, only: export_delta_f, export_full_f
    use export_f, only: export_f_stella
    implicit none
    export_f_theta_option = 0
    export_f_zeta_option = 0
    export_f_xi_option = 0
    export_f_x_option = 0
    export_delta_f = .true.
    export_full_f = .false.
    export_f_stella = .true.
  end subroutine pass_outputoptions_to_sfincs

  ! if this subroutine is being called, then
  ! using sfincs in tokamak geometry
  ! so zed in stella is theta
  subroutine pass_geometry_to_sfincs (delrho)

    use constants, only: pi
    use splines, only: linear_interp_periodic
    use zgrid, only: nz2pi, zed
    use stella_geometry, only: bmag, dbdzed, gradpar
    use stella_geometry, only: dBdrho, d2Bdrdth, dgradpardrho, dIdrho
    use stella_geometry, only: geo_surf
    use globalVariables, only: BHat
    use globalVariables, only: dBHatdtheta
    use globalVariables, only: iota
    use globalVariables, only: DHat
    use globalVariables, only: BHat_sup_theta
    use globalVariables, only: BHat_sub_zeta
    use export_f, only: export_f_theta

    implicit none

    real, intent (in) :: delrho

    integer :: nzeta = 1
    integer :: nzpi
    real :: q_local
    real, dimension (:), allocatable :: B_local, dBdz_local, gradpar_local
    real, dimension (:), allocatable :: zed_stella
    real, dimension (:), allocatable :: theta_sfincs

    nzpi = nz2pi/2
    allocate (B_local(-nzpi:nzpi))
    allocate (dBdz_local(-nzpi:nzpi))
    allocate (gradpar_local(-nzpi:nzpi))
    allocate (theta_sfincs(ntheta))
    allocate (zed_stella(-nzpi:nzpi))

    call init_zero_arrays

    ! first get some geometric quantities at this radius
    ! for theta from -pi to pi
    q_local = geo_surf%qinp*(1.0+delrho*geo_surf%shat/geo_surf%rhoc)
    B_local = bmag(1,-nzpi:nzpi) + delrho*dBdrho(-nzpi:nzpi)
    dBdz_local = dbdzed(1,-nzpi:nzpi) + delrho*d2Bdrdth(-nzpi:nzpi)
    gradpar_local = gradpar(-nzpi:nzpi) + delrho*dgradpardrho(-nzpi:nzpi)

    zed_stella = zed(-nzpi:nzpi)+pi
    theta_sfincs = export_f_theta(:ntheta)

    iota = 1./q_local

    ! interpolate from stella zed-grid to sfincs theta grid
    ! point at -pi (stella) is same as point at 0 (sfincs)
    BHat(1,1) = B_local(-nzpi)
    call linear_interp_periodic (zed_stella, B_local, theta_sfincs(2:), BHat(2:,1))
    ! FLAG -- needs to be changed for stellarator runs
    BHat = spread(BHat(:,1),2,nzeta)
    
    dBHatdtheta(1,1) = dBdz_local(-nzpi)
    call linear_interp_periodic (zed_stella, dBdz_local, theta_sfincs(2:), dBHatdtheta(2:,1))
    dBHatdtheta = spread(dBHatdtheta(:,1),2,nzeta)

    ! this is bhat . grad theta
    BHat_sup_theta(1,1) = B_local(-nzpi)*gradpar_local(-nzpi)
    call linear_interp_periodic (zed_stella, B_local*gradpar_local, theta_sfincs(2:), BHat_sup_theta(2:,1))
    BHat_sup_theta = spread(BHat_sup_theta(:,1),2,nzeta)
    ! this is I(psi) / (aref*Bref)
    BHat_sub_zeta = geo_surf%rgeo + delrho*dIdrho
    ! this is grad psitor . (grad theta x grad zeta)
    ! note that + sign below relies on B = I grad zeta + grad zeta x grad psi
    DHat = q_local*BHat_sup_theta

    deallocate (B_local, dBdz_local, gradpar_local)
    deallocate (theta_sfincs, zed_stella)

  end subroutine pass_geometry_to_sfincs

  subroutine init_zero_arrays
    use globalVariables, only: dBHatdzeta
    use globalVariables, only: dBHatdpsiHat
    use globalVariables, only: BHat_sup_zeta
    use globalVariables, only: BHat_sub_psi
    use globalVariables, only: BHat_sub_theta
    use globalVariables, only: dBHat_sub_psi_dtheta
    use globalVariables, only: dBHat_sub_psi_dzeta
    use globalVariables, only: dBHat_sub_theta_dpsiHat
    use globalVariables, only: dBHat_sub_theta_dzeta
    use globalVariables, only: dBHat_sub_zeta_dpsiHat
    use globalVariables, only: dBHat_sub_zeta_dtheta
    use globalVariables, only: dBHat_sup_theta_dpsiHat
    use globalVariables, only: dBHat_sup_theta_dzeta
    use globalVariables, only: dBHat_sup_zeta_dpsiHat
    use globalVariables, only: dBHat_sup_zeta_dtheta
    implicit none
    dBHatdzeta = 0.
    dBHatdpsiHat = 0.
    BHat_sup_zeta = 0.
    BHat_sub_psi = 0.
    BHat_sub_theta = 0.
    dBHat_sub_psi_dtheta = 0.
    dBHat_sub_psi_dzeta = 0.
    dBHat_sub_theta_dpsiHat = 0.
    dBHat_sub_theta_dzeta = 0.
    dBHat_sub_zeta_dpsiHat = 0.
    dBHat_sub_zeta_dtheta = 0.
    dBHat_sup_theta_dpsiHat = 0.
    dBHat_sup_theta_dzeta = 0.
    dBHat_sup_zeta_dpsiHat = 0.
    dBHat_sup_zeta_dtheta = 0.
  end subroutine init_zero_arrays

  subroutine get_sfincs_output (f_neoclassical, phi_neoclassical, &
       dfneo_dalpha, dphineo_dalpha)
    
    use constants, only: pi
    use sort, only: sort_array_ascending, unsort_array_ascending
    use species, only: nspec
    use zgrid, only: nzgrid, nz2pi
    use export_f, only: h_sfincs => delta_f
    use globalVariables, only: Phi1Hat
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,-nzgrid:,:,:,:), intent (out) :: f_neoclassical
    real, dimension (:,-nzgrid:), intent (out) :: phi_neoclassical
    real, dimension (:,-nzgrid:,:,:,:), intent (out), optional :: dfneo_dalpha
    real, dimension (:,-nzgrid:), intent (out), optional :: dphineo_dalpha

    integer :: i, j
    integer :: ialpha, iz, is

    real, dimension (:), allocatable :: zed_stella
    integer :: nfp, nfp_stella
    integer, dimension (:), allocatable :: nzed_per_field_period
    real, dimension (:,:), allocatable :: zed_stella_by_field_period
    real, dimension (:,:), allocatable :: alpha_like_stella
    integer, dimension (:,:), allocatable :: alpha_sort_map

    integer :: nzed_sfincs, nalpha_sfincs
    real, dimension (:), allocatable :: zed_sfincs, alpha_like_sfincs
    real, dimension (:,:), allocatable :: phi_sfincs

    real, dimension (:,:), allocatable :: phi_stella_zgrid, phi_stella
    real, dimension (:,:), allocatable :: dphi_dalpha_stella_zgrid
    real, dimension (:,:), allocatable :: tmp_sfincs, tmp_stella_zgrid, tmp_stella
    real, dimension (:,:), allocatable :: dh_stella_zgrid
    real, dimension (:,:,:,:), allocatable :: h_stella, dh_stella

    ! zed coordinate in stella is zeta when simulating stellarators (using vmec)
    ! and theta otherwise.  this leads to some complications, treated below

    allocate (zed_stella(nz2pi))
    allocate (alpha_like_stella(nalpha,nz2pi))
    allocate (alpha_sort_map(nalpha,nz2pi))
    ! obtain theta and zeta grids used in stella, and assign them
    ! to the zed and alpha-like coordinates.
    ! for stellarator, zed=zeta and alpha_like=theta.
    ! for tokamak, zed=theta and alpha_like = zeta.
    call get_stella_theta_zeta_grids (alpha_like_stella, zed_stella)
!    do iz = 1, nz2pi
!       do ialpha = 1, size(alpha_like_stella,1)
!          write (*,*) 'unsorted_alpha_like_stella', alpha_like_stella(ialpha,iz)
!       end do
!       write (*,*)
!    end do
    ! rearrange alpha_like_stella to be in ascending order and store map
    ! so that sorting can be undone later
    do iz = 1, nz2pi
       call sort_array_ascending (alpha_like_stella(:,iz), alpha_sort_map(:,iz))
    end do
!    do iz = 1, nz2pi
!       do ialpha = 1, size(alpha_like_stella,1)
!          write (*,*) 'sorted_alpha_like_stella', alpha_like_stella(ialpha,iz)
!       end do
!       write (*,*)
!    end do
    ! obtain the number of alpha-like and zed grid points to use in sfincs theta-zeta grid
    ! also obtain the number of field periods per 2*pi segment in zed
    call get_sfincs_theta_zeta_grid_sizes (nalpha_sfincs, nzed_sfincs, nfp)
    ! obtain the alpha-like and zed coordinate grids
    ! note that additional points are added at periodic points
    ! that are not sampled in sfincs
    allocate (zed_sfincs(nzed_sfincs))
    allocate (alpha_like_sfincs(nalpha_sfincs))
    call get_sfincs_theta_zeta_grids (alpha_like_sfincs, zed_sfincs)

    allocate (phi_sfincs(nalpha_sfincs,nzed_sfincs))
    call get_sfincs_field_theta_zeta (Phi1Hat, phi_sfincs)    

    ! this is the number of field periods included in stella
    ! simulation domain
    nfp_stella = int((zed_stella(nz2pi)*nfp-100.*epsilon(0.))/(2.*pi)) + 1

    ! obtain the number of zed grid points in each field period within the zed domain
    allocate (nzed_per_field_period(nfp_stella))
    call get_nzed_per_field_period (zed_stella, nfp, nzed_per_field_period)
    ! obtain the zed grid within each field period
    allocate (zed_stella_by_field_period(maxval(nzed_per_field_period),nfp_stella))
    call sort_zed_by_field_period (zed_stella, nzed_per_field_period, nfp, zed_stella_by_field_period)

    ! interpolate phi from sfincs zed grid to stella zed grid
    allocate (phi_stella_zgrid(nalpha_sfincs,nz2pi))
    call get_field_on_stella_zed_grid (phi_sfincs, nfp_stella, nfp, nzed_per_field_period, &
         nalpha_sfincs, zed_sfincs, zed_stella_by_field_period, phi_stella_zgrid)

    allocate (phi_stella(nalpha,nz2pi))
    ! interpolate onto stella (sorted) alpha grid
    call get_field_stella (phi_stella_zgrid, alpha_like_sfincs, alpha_like_stella, phi_stella)
    ! need to remap from ascending (sorted) alpha grid to original ordering
    do iz = 1, nz2pi
       call unsort_array_ascending (phi_stella(:,iz), alpha_sort_map(:,iz))
    end do
    call get_field_on_extended_zed (phi_stella, phi_neoclassical)

    if (present(dphineo_dalpha)) then
       allocate (dphi_dalpha_stella_zgrid(nalpha_sfincs,nz2pi))
       call get_dfield_dalpha (phi_stella_zgrid, alpha_like_sfincs, dphi_dalpha_stella_zgrid)
       call get_field_stella (dphi_dalpha_stella_zgrid, alpha_like_sfincs, alpha_like_stella, phi_stella)
       do iz = 1, nz2pi
          call unsort_array_ascending (phi_stella(:,iz), alpha_sort_map(:,iz))
       end do
       call get_field_on_extended_zed (phi_stella, dphineo_dalpha)
       deallocate (dphi_dalpha_stella_zgrid)
    end if

    deallocate (phi_stella, phi_stella_zgrid, phi_sfincs)

    allocate (tmp_sfincs(nalpha_sfincs,nzed_sfincs))
    allocate (tmp_stella_zgrid(nalpha_sfincs,nz2pi))
    allocate (tmp_stella(nalpha,nz2pi))
    allocate (h_stella(nalpha,nz2pi,size(h_sfincs,4),size(h_sfincs,5)))
    if (present(dfneo_dalpha)) then
       allocate (dh_stella_zgrid(nalpha_sfincs,nz2pi))
       allocate (dh_stella(nalpha,nz2pi,size(h_sfincs,4),size(h_sfincs,5)))
    end if

    do is = 1, nspec
       do i = 1, size(h_sfincs,5)
          do j = 1, size(h_sfincs,4)
             ! re-order theta and zeta indices for sfincs h to ensure alpha-like coordinate
             ! appears before zed coordinate
             call get_sfincs_field_theta_zeta (h_sfincs(is,:,:,j,i), tmp_sfincs)
             ! interpolate onto stella zed grid
             call get_field_on_stella_zed_grid (tmp_sfincs, nfp_stella, nfp, nzed_per_field_period, &
                  nalpha_sfincs, zed_sfincs, zed_stella_by_field_period, tmp_stella_zgrid)
             ! interpolate onto (sorted) stella alpha-like coordinate
             call get_field_stella (tmp_stella_zgrid, alpha_like_sfincs, alpha_like_stella, tmp_stella)
             do iz = 1, nz2pi
                call unsort_array_ascending (tmp_stella(:,iz), alpha_sort_map(:,iz))
             end do
             ! use periodicity to copy onto extended zed grid if nperiod > 1
             call get_field_on_extended_zed (tmp_stella, h_stella(:,:,j,i))

             if (present(dfneo_dalpha)) then
                call get_dfield_dalpha (tmp_stella_zgrid, alpha_like_sfincs, dh_stella_zgrid)
                call get_field_stella (dh_stella_zgrid, alpha_like_sfincs, &
                     alpha_like_stella, tmp_stella)
                do iz = 1, nz2pi
                   call unsort_array_ascending (tmp_stella(:,iz), alpha_sort_map(:,iz))
                end do
                call get_field_on_extended_zed (tmp_stella, dh_stella(:,:,j,i))
             end if
          end do
       end do

       do ialpha = 1, nalpha
          do iz = -nzgrid, nzgrid
             call sfincs_vspace_to_stella_vspace (ialpha, iz, is, h_stella(ialpha,iz+nzgrid+1,:,:), &
                  phi_neoclassical(ialpha,iz), f_neoclassical(ialpha,iz,:,:,is))
             if (present(dfneo_dalpha)) call sfincs_vspace_to_stella_vspace (ialpha, iz, is, &
                  dh_stella(ialpha,iz+nzgrid+1,:,:), dphineo_dalpha(ialpha,iz), &
                  dfneo_dalpha(ialpha,iz,:,:,is))
          end do
       end do

    end do

    deallocate (tmp_sfincs, tmp_stella_zgrid, tmp_stella)
    deallocate (zed_stella, alpha_like_stella)
    deallocate (alpha_sort_map)
    deallocate (zed_sfincs, alpha_like_sfincs)
    deallocate (nzed_per_field_period, zed_stella_by_field_period)
    deallocate (h_stella)
    if (present(dfneo_dalpha)) deallocate (dh_stella, dh_stella_zgrid)

  end subroutine get_sfincs_output

  subroutine get_stella_theta_zeta_grids (alpha_like_stella, zed_stella)
    
    use constants, only: pi
    use zgrid, only: nz2pi, zed
    use stella_geometry, only: zed_scalefac
    use stella_geometry, only: alpha
    use kt_grids, only: nalpha
    use globalVariables, only: iota
    
    implicit none
    
    real, dimension (:,:), intent (out) :: alpha_like_stella
    real, dimension (:), intent (out) :: zed_stella
    
    integer :: nzpi
    
    nzpi = nz2pi/2
    
    ! convert from scaled zed grid on [-pi,pi]
    ! to un-scaled grid with lower bound of zero
    ! note that zed_scalefac=1 unless geometryScheme=5 (VMEC)
    ! for geometryScheme=5, this will get extended zeta domain
    ! with lower bound of 0
    zed_stella = (zed(-nzpi:nzpi) + pi)/zed_scalefac
    
    ! if geometryScheme is 5, then using vmec geo
    ! and thus zed in stella is scaled zeta
    ! otherwise zed in stella is theta
    if (geometryScheme == 5) then
       ! alpha_like coordinate is theta = alpha + iota*zeta
       alpha_like_stella = spread(alpha,2,nz2pi) + spread(iota*zed_stella,1,nalpha)
    else
       ! alpha_like coordinate is zeta = (theta-alpha)/iota
       alpha_like_stella = (spread(zed_stella,1,nalpha)-spread(alpha,2,nz2pi))/iota
    end if

    ! restrict alpha to [0,2*pi]
    alpha_like_stella = modulo(alpha_like_stella,2.*pi)
    
  end subroutine get_stella_theta_zeta_grids

  subroutine get_nzed_per_field_period (zed_stella, nfp, nzed_per_field_period)

    use constants, only: pi
    use zgrid, only: nz2pi

    implicit none

    real, dimension (:), intent (in) :: zed_stella
    integer, intent (in) :: nfp
    integer, dimension (:), intent (out)  :: nzed_per_field_period

    integer :: ifp, iz

    nzed_per_field_period = 0

    ifp = 1 ; iz = 1
    do while (iz <= nz2pi)
       if (zed_stella(iz) <= ifp*2.*pi/nfp) then
          nzed_per_field_period(ifp) = nzed_per_field_period(ifp) + 1
          iz = iz + 1
       else
          ifp = ifp + 1
       end if
    end do

  end subroutine get_nzed_per_field_period

  subroutine sort_zed_by_field_period (zed_ext, nzed_per_fp, nfp, zed_by_fp)

    use constants, only: pi

    implicit none

    real, dimension (:), intent (in) :: zed_ext
    integer, dimension (:), intent (in) :: nzed_per_fp
    integer, intent (in) :: nfp
    real, dimension (:,:), intent (out) :: zed_by_fp
    
    integer :: ifp, llim, ulim

    ulim = 0
    do ifp = 1, size(zed_by_fp,2)
       llim = ulim + 1
       ulim = llim + nzed_per_fp(ifp) - 1
       zed_by_fp(:nzed_per_fp(ifp),ifp) = modulo(zed_ext(llim:ulim)-100.*epsilon(0.),2.*pi/nfp)
    end do

    ! avoid special case of setting zed = 0 to zed=2*pi/nfp
    zed_by_fp(1,1) = 0.0

  end subroutine sort_zed_by_field_period

  subroutine get_sfincs_theta_zeta_grid_sizes (nalpha_sfincs, nzed_sfincs, nfp)
    
    use constants, only: pi
    use export_f, only: export_f_zeta

    implicit none
    
    integer, intent (out) :: nalpha_sfincs, nzed_sfincs, nfp
    
    ! if geometryScheme is 5, then using vmec geo
    ! and thus zed in stella is scaled zeta
    ! otherwise zed in stella is theta
    if (geometryScheme == 5) then
       ! note that zeta grid in sfincs only covers one field period of the stellarator
       ! get the number of field periods from the sfincs zeta grid
       nfp = nint(2.*pi/(export_f_zeta(nzeta)+export_f_zeta(2)))
       ! zed coordinate is zeta
       nzed_sfincs = nzeta+1
       ! alpha_like coordinate is theta = alpha + iota*zeta
       nalpha_sfincs = ntheta+1
    else
       nfp = 1
       ! zed coordinate is theta
       nzed_sfincs = ntheta+1
       ! alpha_like coordinate is zeta = (theta-alpha)*q
       nalpha_sfincs = nzeta
    end if

  end subroutine get_sfincs_theta_zeta_grid_sizes

  subroutine get_sfincs_theta_zeta_grids (alpha_like_sfincs, zed_sfincs)

    use export_f, only: export_f_theta, export_f_zeta
    
    implicit none

    real, dimension (:), intent (out) :: zed_sfincs, alpha_like_sfincs

    if (geometryScheme == 5) then
       ! zed is zeta.  it goes from 0 to 2*pi/nfp - dzeta in sfincs,
       ! where nfp is the number of field periods
       zed_sfincs(:nzeta) = export_f_zeta(:nzeta)
       ! add in point at 2*pi/nfp using periodicity
       zed_sfincs(nzeta+1) = zed_sfincs(nzeta) + zed_sfincs(2)

       ! alpha-like coordinate is theta.  it goes from 0 to 2*pi-dtheta in sfincs
       alpha_like_sfincs(:ntheta) = export_f_theta(:ntheta)
       ! add in point at 2*pi using periodicity
       alpha_like_sfincs(ntheta+1) = export_f_theta(ntheta)+export_f_theta(2)
    else
       ! zed is theta.  it goes from 0 to 2*pi-dtheta in sfincs
       zed_sfincs(:ntheta) = export_f_theta(:ntheta)
       ! add in point at 2*pi using periodicity
       zed_sfincs(ntheta+1) = export_f_theta(ntheta) + export_f_theta(2)

       ! alpha-like coordinate is zeta.  there should only be 1 zeta for tokamak calculation
       alpha_like_sfincs(:nzeta) = export_f_zeta(:nzeta)
    end if

  end subroutine get_sfincs_theta_zeta_grids

  subroutine get_sfincs_field_theta_zeta (field_in, field_out)

    implicit none

    real, dimension (:,:), intent (in) :: field_in
    real, dimension (:,:), intent (out) :: field_out

    integer :: itheta, izeta

    ! want phi from sfincs such that alpha-like coordinate
    ! appears in first index and zed coordinate in second index
    if (geometryScheme == 5) then
       ! get phi on sfincs (theta,zeta) grid
       field_out(:ntheta,:nzeta) = field_in
       ! use periodicity in theta to add point at theta=2*pi
       field_out(ntheta+1,:nzeta) = field_out(1,:nzeta)
       ! use periodicity in zeta to add point at zeta = 2*pi/nfp
       field_out(:,nzeta+1) = field_out(:,1)
    else
       ! get phi on sfincs (zeta,theta) grid
       do izeta = 1, nzeta
          do itheta = 1, ntheta
             field_out(izeta,itheta) = field_in(itheta,izeta)
          end do
       end do
       ! use periodicity in theta to add point at 2*pi
       field_out(:,ntheta+1) = field_out(:,1)
    end if

  end subroutine get_sfincs_field_theta_zeta

  subroutine get_field_on_stella_zed_grid (field_sfincs, nfp_stella, nfp, nzed_per_field_period, &
       nalpha_sfincs, zed_sfincs, zed_stella_by_field_period, field_stella_zgrid)

    use constants, only: pi
    use splines, only: linear_interp_periodic

    implicit none

    real, dimension (:,:), intent (in) :: field_sfincs
    integer, intent (in) :: nfp_stella, nfp, nalpha_sfincs
    integer, dimension (:), intent (in) :: nzed_per_field_period
    real, dimension (:), intent (in) ::  zed_sfincs
    real, dimension (:,:), intent (in) :: zed_stella_by_field_period
    real, dimension (:,:), intent (out) :: field_stella_zgrid

    integer :: ialpha, ifp
    integer :: llim, ulim

    do ialpha = 1, nalpha_sfincs
       ulim = 0
       do ifp = 1, nfp_stella
          llim = ulim + 1
          ulim = llim + nzed_per_field_period(ifp) - 1
          call linear_interp_periodic (zed_sfincs, field_sfincs(ialpha,:), &
               zed_stella_by_field_period(:nzed_per_field_period(ifp),ifp), &
               field_stella_zgrid(ialpha,llim:ulim), 2.*pi/nfp)
       end do
    end do

  end subroutine get_field_on_stella_zed_grid

  subroutine get_field_stella (field_stella_zgrid, alpha_like_sfincs, alpha_like_stella, field_stella)

    use splines, only: linear_interp_periodic
    use zgrid, only: nz2pi

    implicit none

    real, dimension (:,:), intent (in) :: field_stella_zgrid
    real, dimension (:), intent (in) :: alpha_like_sfincs
    real, dimension (:,:), intent (in) :: alpha_like_stella
    real, dimension (:,:), intent (out) :: field_stella

    integer :: iz

    do iz = 1, nz2pi
       call linear_interp_periodic (alpha_like_sfincs, field_stella_zgrid(:,iz), &
            alpha_like_stella(:,iz), field_stella(:,iz))
    end do

  end subroutine get_field_stella

  subroutine get_field_on_extended_zed (field_stella, field_neoclassical)

    use zgrid, only: nzgrid, nz2pi, nperiod
    use kt_grids, only: nalpha

    implicit none

    real, dimension (:,:), intent (in) :: field_stella
    real, dimension (:,-nzgrid:), intent (out) :: field_neoclassical

    integer :: ialpha
    integer :: iz_low, iz_up
    integer :: ip

    ! need to account for cases with nperiod > 1
    do ialpha = 1, nalpha
       iz_low = -nzgrid
       iz_up = -nzgrid+nz2pi-1
       field_neoclassical(ialpha,iz_low:iz_up) = field_stella(ialpha,:)
       ! if nperiod > 1 need to make copies of
       ! neoclassical potential for other 2pi segments
       if (nperiod > 1) then
          do ip = 2, 2*nperiod-1
             iz_low = iz_up + 1
             iz_up = iz_low + nz2pi -2
             field_neoclassical(ialpha,iz_low:iz_up) = field_stella(ialpha,2:)
          end do
       end if
    end do
    
  end subroutine get_field_on_extended_zed

  subroutine sfincs_vspace_to_stella_vspace (ialpha, iz, is, h_stella, phi_neoclassical, f_neoclassical)

    use constants, only: pi
    use species, only: spec
    use vpamu_grids, only: nvpa, nvgrid, nmu
    use vpamu_grids, only: vpa, vperp2
    use vpamu_grids, only: maxwell_mu, maxwell_vpa
    use globalVariables, only: nxi_sfincs => nxi
    use globalVariables, only: nx_sfincs => nx
    use globalVariables, only: x_sfincs => x
    use xGrid, only: xGrid_k

    implicit none

    integer, intent (in) :: ialpha, iz, is
    real, dimension (:,:), intent (in) :: h_stella
    real, intent (in) :: phi_neoclassical
    real, dimension (:,:), intent (out) :: f_neoclassical

    integer :: iv, imu, ixi
    real, dimension (1) :: x_stella
    integer, dimension (2) :: sgnvpa
    integer :: nxi_stella
    real, dimension (:), allocatable :: xi_stella, hstella
    real, dimension (:,:), allocatable :: xsfincs_to_xstella, legpoly
    real, dimension (:), allocatable :: htmp

    ! each (vpa,mu) pair in stella specifies a speed
    ! on each speed arc, there are two (vpa,mu) pairs:
    ! one each corresponding to a given +/- vpa
    nxi_stella = 2
    allocate (xi_stella(nxi_stella))
    allocate (hstella(nxi_stella))
    allocate (legpoly(nxi_stella,0:nxi_sfincs-1))

    allocate (htmp(nxi_sfincs))
    allocate (xsfincs_to_xstella(1,nx_sfincs))

    sgnvpa(1) = 1 ; sgnvpa(2) = -1

    ! h_stella is on the sfincs energy grid
    ! but is spectral in pitch-angle
    do imu = 1, nmu
       ! loop over positive vpa values
       ! negative vpa values will be treated inside loop using symmetry
       do iv = nvgrid+1, nvpa
          ! x_stella is the speed 
          ! corresponding to this (vpa,mu) grid point
          x_stella = sqrt(vpa(iv)**2+vperp2(ialpha,iz,imu))
          ! xi_stella contains the two pitch angles (+/-)vpa/v
          ! corresponding to this (vpa,mu) grid point
          xi_stella = sgnvpa*vpa(iv)/x_stella(1)
          
          ! set up matrix that interpolates from sfincs speed grid
          ! to the speed corresponding to this (vpa,mu) grid point
          call polynomialInterpolationMatrix (nx_sfincs, 1, &
               x_sfincs, x_stella, exp(-x_sfincs*x_sfincs)*(x_sfincs**xGrid_k), &
               exp(-x_stella*x_stella)*(x_stella**xGrid_k), xsfincs_to_xstella)
                   
          ! do the interpolation
          do ixi = 1, nxi_sfincs
             htmp(ixi) = sum(h_stella(ixi,:)*xsfincs_to_xstella(1,:))
          end do
                   
          ! next need to Legendre transform in pitch-angle
          ! first evaluate Legendre polynomials at requested pitch angles
          call legendre (xi_stella, legpoly)
          
          ! then do the transforms
          call legendre_transform (legpoly, htmp, hstella)
          
          f_neoclassical(nvpa-iv+1,imu) = hstella(2)
          f_neoclassical(iv,imu) = hstella(1)
          
       end do
       ! h_sfincs is H_nc / (nref/vt_ref^3), with H_nc the non-Boltzmann part of F_nc
       ! NB: n_ref, etc. is fixed in stella to be the reference density
       ! at the central sfincs simulation; i.e., it does not vary with radius
       ! to be consistent with stella distribution functions,
       ! want H_nc / (n_s / vt_s^3 * pi^(3/2)) / exp(-v^2/vts^2)
       f_neoclassical(:,imu) = f_neoclassical(:,imu) &
            * pi**1.5 * spec(is)%stm**3/spec(is)%dens
       
       ! phi_sfincs is e phi / Tref as long as alpha=1 (default)
       ! need to multiply by Z_s * Tref/T_s
       f_neoclassical(:,imu) = f_neoclassical(:,imu) &
            - phi_neoclassical*spec(is)%zt*maxwell_vpa*maxwell_mu(ialpha,iz,imu)
    end do

    deallocate (xi_stella, hstella, legpoly)
    deallocate (xsfincs_to_xstella)
    deallocate (htmp)

  end subroutine sfincs_vspace_to_stella_vspace

  subroutine get_dfield_dalpha (field, alpha_like_sfincs, dfield_dalpha)

    use zgrid, only: nz2pi

    implicit none

    real, dimension (:,:), intent (in) :: field
    real, dimension (:), intent (in) :: alpha_like_sfincs
    real, dimension (:,:), intent (out) :: dfield_dalpha

    integer :: iz

    do iz = 1, nz2pi
       call get_periodic_derivative (field(:,iz),alpha_like_sfincs,dfield_dalpha(:,iz))
    end do

  end subroutine get_dfield_dalpha

  ! 4th order accurate, centered differences, assumes 
  ! first and last elements of f are equal (periodic)
  subroutine get_periodic_derivative (f, x, dfdx)

    implicit none

    real, dimension (:), intent (in) :: f, x
    real, dimension (:), intent (out) :: dfdx

    integer :: i, n
    real, dimension (:), allocatable :: fp

    n = size(x)

    allocate (fp(-1:n+2))
    ! extend f using periodicity, as these additional points 
    ! required near boundaries for finite differences
    fp(1:n) = f
    fp(-1:0) = f(n-1:n)
    fp(n+1:n+2) = f(1:2)

    do i = 1, n
       dfdx(i) = (0.25*(fp(i-2)-fp(i+2)) + 2.0*(fp(i+1)-fp(i-1)))/3.0
    end do

    deallocate (fp)

  end subroutine get_periodic_derivative

!   subroutine bilinear_interpolation (alpha_in, zed_in, phi_in, alpha_out, zed_out, phi_out)
    
!     use zgrid, only: nz2pi
    
!     implicit none
    
!     real, dimension (:,:), intent (in) :: alpha_sfincs, phi_sfincs
!     real, dimension (:), intent (in) :: zed_sfincs
!     real, dimension (:), intent (in) :: alpha_stella, zed_stella
!     real, dimension (:,:), intent (out) :: phi_stella
    
!     integer :: ialpha, iz
    
!     do ialpha = 1, nalpha
!        do iz = 1, nz2pi
!           call find_sfincs_cell (alpha_stella(ialpha), zed_stella(iz), alpha_sfincs, zed_sfincs, &
!                alpha_grids, zed_grids)
!        end do
!     end do
    
!   contains
    
!     subroutine find_sfincs_cell (alpha_target, zed_target, alpha, zed, ialpha_out, ized_out)
      
!       use zgrid, only: nz2pi

!       implicit none

!       real, intent (in) :: alpha_target, zed_target
!       real, dimension (:,:), intent (in) :: alpha
!       real, dimension (:), intent (in) :: zed
!       integer, dimension (2), intent (out) :: ialpha_out, ized_out

!       integer :: iz, ia

!       do iz = 1, nz2pi
!          do ia = 1, size(alpha,1)
!             if (
!          end do
!       end do

!     end subroutine find_sfincs_cell
    
!   end subroutine bilinear_interpolation

  ! returns the Legendre polynomials (legp)
  ! on requested grid (x)
  subroutine legendre (x, legp)
    
    implicit none
    
    real, dimension (:), intent (in) :: x
    real, dimension (:,0:), intent (out) :: legp
    
    integer :: n, idx
    
    n = size(legp,2)-1
    
    legp(:,0) = 1.0
    legp(:,1) = x
    
    do idx = 2, n
       legp(:,idx) = ((2.*idx-1.)*x*legp(:,idx-1) + (1.-idx)*legp(:,idx-2))/idx
    end do
    
  end subroutine legendre
  
  subroutine legendre_transform (legp, coefs, func)

    implicit none
    
    real, dimension (:,0:), intent (in) :: legp
    real, dimension (:), intent (in) :: coefs
    real, dimension (:), intent (out) :: func
    
    integer :: i

    func = 0.
    do i = 1, size(coefs)
       func = func + legp(:,i-1)*coefs(i)
    end do

  end subroutine legendre_transform

  subroutine broadcast_sfincs_output (fneo, phineo)
    use mp, only: broadcast
    use zgrid, only: nzgrid
    implicit none
    real, dimension (-nzgrid:,:,:,:,:), intent (in out) :: fneo
    real, dimension (-nzgrid:,:), intent (in out) :: phineo
    call broadcast (fneo)
    call broadcast (phineo)
  end subroutine broadcast_sfincs_output

  subroutine write_sfincs (irad, nrad_max)

    use species, only: nspec
    use globalVariables, only: Phi1Hat
    use export_f, only: export_f_zeta, export_f_theta
    use export_f, only: delta_f

    implicit none

    integer, intent (in) :: irad, nrad_max
    
    integer :: unit = 999
    integer :: izeta, itheta, is, i, j
    character (1) :: irad_str

    write (irad_str,'(i0)') irad+min(nrad_max,1)
    open (unit=unit,file='sfincs.output.rad'//irad_str,status='replace',action='write')

    do izeta = 1, nzeta
       do itheta = 1, ntheta
          write (unit,'(a8,3e13.5,i3)') 'Phi1Hat', export_f_theta(itheta), export_f_zeta(izeta), Phi1Hat(itheta,izeta), irad
       end do
       write (unit,*)
    end do
    write (unit,*)

    do is = 1, nspec
       do i = 1, size(delta_f,5)
          do j = 1, size(delta_f,4)
             do izeta = 1, nzeta
                do itheta = 1, ntheta
                   write (unit,'(a8,3e13.5,4i5)') 'delta_f', export_f_theta(itheta), export_f_zeta(izeta), &
                        delta_f(is,itheta,izeta,j,i), i, j, is, irad
                end do
                write (unit,*)
             end do
          end do
       end do
    end do
    write (unit,*)

    close (unit)

  end subroutine write_sfincs

  subroutine read_sfincs_output (irad, nrad_max)

    use species, only: nspec
    use globalVariables, only: Phi1Hat
    use export_f, only: export_f_zeta, export_f_theta
    use export_f, only: delta_f

    implicit none

    integer, intent (in) :: irad, nrad_max
    
    integer :: unit = 999
    integer :: izeta, itheta, is, i, j
    character (8) :: dum
    character (1) :: irad_str

    write (irad_str,'(i0)') irad+nrad_max
    open (unit=unit,file='sfincs.output.rad'//irad_str,status='old',action='read')

    do izeta = 1, nzeta
       do itheta = 1, ntheta
          read (unit,*) dum, export_f_theta(itheta), export_f_zeta(izeta), Phi1Hat(itheta,izeta), dum
       end do
       read (unit,*)
    end do
    read (unit,*)
    
    do is = 1, nspec
       do i = 1, size(delta_f,5)
          do j = 1, size(delta_f,4)
             do izeta = 1, nzeta
                do itheta = 1, ntheta
                   read (unit,*) dum, export_f_theta(itheta), export_f_zeta(izeta), &
                        delta_f(is,itheta,izeta,j,i), dum, dum, dum, dum
                end do
                read (unit,*)
             end do
          end do
       end do
    end do
    read (unit,*)
    
    close (unit)

  end subroutine read_sfincs_output

# endif

end module sfincs_interface
