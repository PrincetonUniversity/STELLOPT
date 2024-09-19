!-----------------------------------------------------------------------------
!
!  -THIS IS A BETA VERSION OF THE PENTA CODE-
!
! -- PENTA3 v1.3
!
! AUTHORS: J.LORE AND D.SPONG
!
!-----------------------------------------------------------------------------
!
!  PENTA calculates the neoclassical parallel flows, radial particle and 
!   energy fluxes, and the radial electric field for a surface given the
!   plasma profiles (n, T), surface geometry information (from VMEC) and
!   the monoenergetic transport coefficients (from DKES)
!
! To run,type: PENTA3 [command line arguments]
!
! where the command line arguments are:
!
! 1  unique identifier text on DKES database files
!     i.e., for D11_star_***, this would be
!     the text that replaces *** for the specific surface (ex: hsx_s10)
! 2  Er,min
! 3  Er,max  - args 2 and 3 specify the interval in Er over which
!     the search is done for the ambipolar electric field root.
!    Er,min and Er,max are dimensionless normalized electric
!     field variables, i.e.: e<a>Er,min/max/kT_e
!   - Note if the logical variable input_is_Er below is true then 
!     the search interval is set by args 2 and 3 as Er in (V/cm)
! 4  flux surface number (index from VMEC, typically)
! 5  parameter that determines whether output data should be
!     written into new files (=0) or appended to
!     existing files (=1). Typically,when a script is run for a
!     sequence of flux surfaces, this is set to 0 for the first
!     surface and to 1 for subsequent surfaces
! 6  extension of the profile_data_*** file (i.e., the *** text)
! 7  unique identifier on plasma profiles file, i.e., the ***
!    text in the filename plasma_profiles_***.dat. This allows
!    multiple profiles to be run without having to rename files
! 8  <B*E||> in T*V/m  [real] (Used if a parallel electric field is applied
!     or determined from equilibrium calculations)
! 9  Smax -- The upper limit on the summation of Laguerre polynomial terms.
!     Since the summation is from 0 to Smax the number of terms used is
!     Smax+1.  This affects the number of parallel flow moments calculated and
!     output.
!
!  Input files:
!
!   [D11,D13,D33]_star_lijs_*** - where *** is arg1 above.
!       These files contain the values of efield and cmul and the
!       corresponding normalized monoenergetic coefficients from DKES.
!
!   ion_params - A file containing the namelist "ion_params" which
!       defines the variables num_ion_species, Z_ion_init and 
!       miomp_init which are the number of ion species (integer),
!       the corresponding ion charge numbers (real, array) and 
!       the ion to proton mass ratio (real, array) for each 
!       non-electron species.
!  
!  run_params - A file containing the namelist "run_params" which contains
!    run settings.  These settings are described below.
!
!   profile_data_*** - where *** is arg1 above.  Contains geometry
!       variables, most likely from a VMEC run.
!
!   plasma_profiles_XXX.dat - A file containing the plasma profile
!       data, where XXX is arg7 above.  This file has the format:
!       row 1: number of radial points for which data is specified.
!       all other rows: r/a, ne, Te, Ti1, ni1, Ti2, ni2 ...
!       where the densities are in units of 10**12/cc and the 
!       temperatures in eV and i1, i2... are the ion species.
!
!   beam_force.dat - Same format as plasma_profiles. Units of Newtons
!       (Torque/major radius). Columns : r/a, F.  This is force on the
!       plasma mass, and is divided by species density inside the code.
!
!   Utilde2_profile - Contains the quantity <U**2>, where U is the
!     PS flow function as defined by Sugama and Nishimura.  The
!     first row is the number of points, then r/a points and the
!     corresponding <U**2> value.  Note that if read_U2_file is /false.
!     then <U**2> is calculated from D11 at high collisionality
!      
!   
!  Output files:
!
!   "fluxes_vs_roa" - Contains the radial particle and energy fluxes for 
!     each run surface for each ambipolar root.  The particle fluxes are
!     listed for each species [units of particles/m**2/s] and the energy
!     fluxes as Q/T [units of m**-2s**-1].  The format is given as
!     r/a  Er e<a>Er/kTe Gamma_e Q_e/T_e Gamma_i1 Q_i1/T_i1 ...
!     .
!     .
!     .
!     Where the horizontal ... indicates that additional ion species
!     will be listed horizontally, and the vertical ... indicates that
!     additional ambipolar root results will be listed vertically. 
!     Each additional surface for a given run will then follow the same
!     format and can be appended to this file. 
!
!  "fluxes_vs_Er" - Contains the radial paricle fluxes for each surface
!    for each value used in the ambipolar Er search loop. The format of 
!    this file is:
!    r/a Er e<a>Er/kTe Gamma_e Gamma_i1 Gamma_i2 ...
!    .
!    .
!    .
!    Where the horizontal ... indicates that additional ion species
!    will be listed horizontally, and the vertical ... indicates that
!    each additional Er in the search range will be listed.
!
!  "flows_vs_roa" - Contains the parallel flow moments for each species
!    evaluated at the ambipolar radial electric field. The format of the
!    file is:
!    r/a Er e<a>Er/kTe <B*u_||ke>/<B**2>  <B*u_||ki>/<B**2> ...
!    .
!    .
!    .
!    where < > indicates a flux surface average, u_||ks is the Sonine
!    polynomial of order "k" weighted parallel flow moment for the 
!    species "s". So, after the Er and e<a>Er/kTe first each Sonine
!    weighted parallel flow moment (from k=0 to k=Smax) is listed for 
!    the electrons, then each moment for the first ion species, each
!    moment for the second ion species, etc. 
!    The Sonine polynomial weighted flow moments are closely related 
!    to physical flows for low orders.  u_||0s = u_||s, where u_||s
!    is the parallel particle flow for species s. 
!    q_||s = -(5/2)*p_s*u_||1s, where p_s is the pressure for species
!    s, and q_||s is the parallel heat flow.
!    As in the "fluxes_vs_roa" file, the results are given for each 
!    surface, for each ambipolar root vertically.
!
!  "flows_vs_Er" - Contains the Sonine weighted parallel flow moments
!    for each Er in the ambipolar search range. Similar to the 
!    "fluxes_vs_Er" file, and the flow moments are described under
!    "flows_vs_roa".  
!  
!  "plasma_profiles_check" - Contains the plasma profile values actually
!    used at each surface (T,n,gradients).  Allows for the fits to the
!    provided data points to be evaluated.
!
!  "Jprl_vs_roa" - Contains the parallel current densities in a similar
!    manner as fluxes_vs_roa above. The format of the file is
!    r/a Er e<a>Er/kTe J_prl_e J_prl_i1 J_prl_i2 ... Jprl
!    .
!    .
!    .
!    where the units of the current densities are [A/m**2] and 
!    J_prl_s = ns*qs*<b*u_||s> where ns is the density, qs
!    the charge and b is the normalized bfield b=<B>/<B**2>**1/2.
!    The total parallel current Jprl is then the sum over species.
!    
!  "ucontra_vs_roa" - Contains the contravarient Boozer coordinate poloidal
!    and toroidal components of the flow. These flows are given as 
!    <upol> = <u_vector dot grad theta> (toroidal case is similar). The units
!    are thus 1/m.  The format of the file is
!    r/a Er e<a>Er/kTe <upol_e> <utor_e> <upol_i1> <utor_i1> ...
!    .
!    .
!    . 
!    Similarly to the above files with _vs_roa suffix.
!
!  Versioning:
!   See file DOC/version.txt
!
!   
!-----------------------------------------------------------------------------
!   NOTES
!    - based on 
!         Sugama, Nishimura PoP 9 4637 (2002), 
!         Sugama, Nishimura PoP 15, 042502 (2008),
!         Maassberg, et al, PoP 15, 072504 (2009),
!         Taguchi, Phys. Fluids B, 4 3638 (1992)
!         Spong, PoP 12, 056114 (2005)
!         Note that the code is written from PhD thesis of J. Lore,
!          from which notation typically follows.
!    - All units are SI, except T is [eV] (unless otherwise noted)
!    - Some code used from orignal PENTA code, written by Don Spong and 
!        modified by J. Lore
!
!  7/2009-9/22/2011 JL
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
!+ The main PENTA program
!
Program penta3
 
! Description: 
!  PENTA calculates the neoclassical parallel flows, radial particle and 
!   energy fluxes, and the radial electric field for a surface given the
!   plasma profiles (n, T), surface geometry information (from VMEC) and
!   the monoenergetic transport coefficients (from DKES).  See above for
!   input/output and references.
!
! History:
! Version   Date        Comment
! -------   ----        -------
! ...     .........    Older versions did not follow this scheme
! 1.0     05/24/2010   Original code for "PENTA3". JL 
! 1.1     01/06/2011   Added flow vs Er output, fixed PS flux bug, check
!                      for monotonically increasing cmul, efield in D* 
!                      files, parameter esmall in fit_coeffs.
! 1.2     02/01/2011   Fixed bug where flows(num_Er_test) was used instead
!                      of Flows_ambi(:,iroot) when evaluating ambipolar
!                      fluxes. Updated documentation. Added parallel current
!                      and contravarient flow output files.
! 1.3     09/27/2011   Major update, see version.txt. JL
! 
! Author(s): J. Lore 7/2009 - current
!            D. Spong - pre 7/2009
!
! Code Description:
!   Language:   Fortran 90 with use of Fortran 2003 COMMAND_ARGUMENT_COUNT
!   Standard:   Based on "European Standards for Writing and Documenting
!               Exchangeable Fortran 90 Code"  v1.1
! 
! External subroutines used:
!   define_friction_coeffs: Contained in file penta_subroutines.f90
!   find_Er_roots: Contained in file penta_subroutines.f90
!   form_Xvec: Contained in file penta_subroutines.f90
!   fit_coeffs: Contained in file penta_subroutines.f90

! Modules used:
Use penta_kind_mod                     ! Import rknd, iknd specifications
Use io_unit_spec, Only :             &
  iu_nl,                             & ! Ion parameter namelist i/o unit #
  iu_flux_out,                       & ! flux vs r/a i/o unit #
  iu_pprof_out,                      & ! plasma profile check i/o unit #
  iu_fvEr_out,                       & ! Flux vs Er i/o unit #
  iu_QoTvEr_out,                     & ! Energy flux vs Er i/o unit #
  iu_flows_out,                      & ! flows vs r/a i/o unit #
  iu_flowvEr_out,                    & ! flows vs Er i/o unit #
  iu_Jprl_out,                       & ! Parallel current den vs r/a i/o unit #
  iu_contraflows_out                   ! Contravariant flows vs roa
Use read_input_file_mod, Only :      &
  ! Imported Subroutines
  read_vmec_file,                    & ! Reads VMEC data file
  read_pprof_file,                   & ! Reads plasma profile data file
  read_dkes_star_files,              & ! Reads DKES data files
  read_Utilde2_file,                 & ! Reads Utilde2 file (<U**2> data)
  read_beam_file
Use vmec_var_pass , Only :           & ! VMEC variables
  ! Imported scalar variables
    arad, roa_surf, Bsq, B0,         &
    chip, psip, vol_p, btheta, bzeta
Use pprof_pass, Only :               & ! Plasma profile variables
  ! Imported scalar variables
    ne, Te, dnedr, dTedr,            &
    beam_force,                      & 
  ! Imported array variables (1D)
    ni, Ti, dnidr, dTidr
Use coeff_var_pass, Only :           & ! DKES coeff. variables
  ! Imported scalar variables
    num_c, num_e, &
  ! Imported array variables (1D)
    cmul,efield, &
  ! Imported array variables (2D)
    D11_mat,D13_mat,D31_mat,         &
    D33_mat
Use phys_const, Only :               & ! Physical constants
  ! Imported parameters (scalar)
    p_mass,                          & ! proton mass
    e_mass,                          & ! electron mass
    elem_charge,                     & ! elementary charge
    pi
Use penta_functions_mod, Only :      & ! Functions
  ! Imported functions
     calc_flows_T,                   & ! Taguchi method of parallel flows
     calc_flows_SN,                  & ! S-N method of parallel flow
     calc_fluxes_SN,                 & ! S-N method of radial part. fluxes
     calc_QoTs_SN,                   & ! S-N method of radial energy fluxes
     calc_fluxes_MBT,                & ! M-B-T method of radial part. fluxes
     calc_QoTs_MBT,                  & ! M-B-T method of radial energy fluxes
     calc_flows_DKES,                & ! DKES method of parallel flows
     calc_fluxes_DKES,               & ! DKES method of radial part. fluxes
     calc_QoTs_DKES                    ! DKES method of radial part. fluxes

Use penta_math_routines_mod, Only :  & 
  ! Imported functions
     rlinspace                         ! Creates linearly spaced arrays

Use PENTA_subroutines, Only : find_Er_roots, fit_coeffs, &
     define_friction_coeffs, form_xvec

Implicit None

! Default values for variables set in run_params namelist file:
Logical ::    &
  input_is_Er    = .true.,     & ! If true, Er range is (V/cm) else e<a>Er/kT_e
  log_interp     = .true.,     & ! If true, log. interp. of DKES coeffs is used
  use_quanc8     = .false.,    & ! If false, rect. approx. to convolution used
  read_U2_file   = .true.,     & ! If false <U**2> is calculated from D11*
  flux_cap       = .true.,     & ! If true, min(L_radial) = 0
  output_QoT_vs_Er = .false.,  & ! If true Q/T vs Er output file is written
  Add_Spitzer_to_D33 = .true., & ! If true collisional portion of D33* is added
                                 ! else it is assumed to be included in-file
  use_beam       = .false.       ! If true, beam force in Newtons (Torque/Major radius) is read
                                 ! as a function of 

Integer(iknd) ::    &
  num_Er_test = 50_iknd,        & ! Number of Er points in search range
  numKsteps   = 10000_iknd,     & ! Number of K points (linear) for convolution 
                                 !  (used if use_quanc8=.false.)
  kord_pprof  = 3_iknd,        & ! Spline order for plasma profile fitting 
                                 !  (and <U**2> file)
  keord       = 2_iknd,        & ! Spline order for DKES coeff fitting (efield)
  kcord       = 2_iknd         ! Spline order for DKES coeff fitting (cmul)
Real(rknd) ::    &
  Kmin   = 1.e-5_rknd,         & ! Minimum K in energy convolution
  Kmax   = 20._rknd,           & ! Maximum K in energy convolution
  epsabs = 1.e-8_rknd,         & ! Absolute tolerance for quanc8  
                                 !  (used if use_quanc8=.true.)
  epsrel = 1.e-6_rknd            ! Relative tolerance for quanc8  
                                 !  (used if use_quanc8=.true.)
Character(Len=10) ::     &
  Method  = 'DKES'               ! Which algorithm to use.  Options are
                                 !  'T'    = Taguchi
                                 !  'SN'   = Sugama-Nishimura
                                 !  'MBT'  = Maassberg-Beidler-Turkin
                                 !  'DKES' = Direct energy convolution

! Local parameters used to set max array sizes
Integer(iknd), Parameter :: num_roots_max = 10_iknd
Integer(iknd), Parameter :: num_ion_max = 20_iknd

! Local variables (scalar)
Integer(iknd) :: numargs          ! Number of command line args.
Integer(iknd) :: js               ! Booz_xform (VMEC) index of the surface
Integer(iknd) :: i_append         ! Set to 1 to append to output files
Integer(iknd) :: num_species      ! Total number of plasma species
Integer(iknd) :: num_ion_species  ! Number of ion species
Integer(iknd) :: Smax             ! Max index of Sonine poly. expansion
Integer(iknd) :: iocheck          ! Used to check for file open errors
Integer(iknd) :: ie               ! Er loop index
Integer(iknd) :: ind_X, ind_A     ! Thermodynamic force vector indices
Integer(iknd) :: ispec1           ! Primary species loop index
Integer(iknd) :: min_ind          ! Index for finding Er = 0
Integer(iknd) :: iroot            ! Ambipolar root index
Integer(iknd) :: num_roots        ! Number of ambipolar roots
Real(rknd)    :: Er_min, Er_max   ! Er search range min and max
Real(rknd)    :: B_Eprl           ! Input parallel electric field
Real(rknd)    :: U2               ! <U**2> Related to Pfirsch-Schlueter factor
!Real(rknd)    :: G2               ! <G**2> Related to Pfirsch-Schlueter factor
Real(rknd)    :: vth_e            ! Electron thermal velocity
Real(rknd)    :: loglambda        ! Coulomb logarithm
Real(rknd)    :: Er_test, abs_Er  ! Current value of Er and |Er| used in loop
Real(rknd)    :: min_Er           ! Min value of Er_test_vals
Real(rknd)    :: eaEr_o_kTe       ! Normalized Er (used as output)
Real(rknd) ::  cmin,            & ! Coefficient axes limits
  cmax, emin, emax

! Local variables (array)
Character(Len=100) ::           & ! Command line args
  arg1, arg2, arg3, arg4,       & 
  arg5, arg6, arg7, arg8,       &
  arg9  
Character(Len=100) :: coeff_ext   ! Identifier for DKES coeff. data files
Character(Len=100) :: run_ident   ! Identifier for VMEC data file
Character(Len=20)  :: pprof_char  ! Identifier for plasma profile file
Character(Len=100) :: fpos        ! Status for writing to files (position)
Character(Len=100) :: fstatus     ! Status for writing to files
character(Len=100) :: str_num     ! Used for converting numbers to strings
Real(rknd)         ::           &
  Z_ion_init(num_ion_max),      & ! Ion charge numbers
  miomp_init(num_ion_max)         ! Ion mass ratios

! Local allocatable arrays (1D)
Real(rknd), Allocatable ::      &
  ion_mass(:),                  & ! Ion masses
  Z_ion(:),                     & ! Ion charges
  vth_i(:),                     & ! Ion thermal velocities
  charges(:),                   & ! All species charges
  dens(:),                      & ! All species densities
  masses(:),                    & ! All species masses
  temps(:),                     & ! All species temperatures (eV)
  vths(:),                      & ! All species thermal velocities  
  dTdrs(:),                     & ! All species temperature gradients (eV/m)
  dndrs(:),                     & ! All species density gradients (1/m^2)
  Xvec(:), Avec(:),             & ! Thermodynamic force vectors
  Flows(:),                     & ! Parallel flow moment array
  Gammas(:),                    & ! Radial flux array
  QoTs(:),                      & ! Radial flux array
  xt_c(:),xt_e(:),              & ! Spline knots for DKES coeffs
  Gamma_e_vs_Er(:),             & ! Electron particle flux vs Er
  QoT_e_vs_Er(:),               & ! Electron Energy flux vs Er
  Er_test_vals(:),              & ! Er values used in ambipolar search
  Jprl_ambi(:)                    ! Parallel current density (roots)    

! Local allocatable arrays (2D)
Real(rknd), Allocatable ::      &
  lmat(:,:),                    & ! Classical friction coefficients
  Dspl_D11(:,:),Dspl_D13(:,:),  & ! DKES coeff spline coeffs
  Dspl_D31(:,:),Dspl_D33(:,:),  & ! DKES coeff spline coeffs
  Dspl_Dex(:,:),                & ! Extra radial flux coefficient for SN & T
  Dspl_DUa(:,:),                & ! Extra viscosity coefficient for SN method
  Dspl_Drat(:,:),               & ! Extra coefficient for SN method
  Dspl_Drat2(:,:),              & ! Extra coefficient for SN method
  Dspl_logD11(:,:),             &
  Dspl_logD33(:,:),             &
  cmesh(:,:),                   & ! Repeated 2D array of cmul
  Gamma_i_vs_Er(:,:),           & ! Ion particle flux vs Er
  QoT_i_vs_Er(:,:)                ! Ion Energy flux vs Er
Real(rknd), Allocatable ::      & 
  Flows_ambi(:,:),              & ! Parallel flow moments (j*species,roots)
  Gammas_ambi(:,:),             & ! Radial particle fluxes (species,roots)
  QoTs_ambi(:,:),               & ! Radial energy fluxes (species,roots)
  Jprl_parts(:,:),              & ! Parallel current density (species,roots)
  upol(:,:),                    & ! fsa contravariant poloidal flow (species,roots)
  utor(:,:)                       ! fsa contravariant toroidal flow (species,roots)

Real(rknd) :: Er_roots(num_roots_max)   ! The maximum number of roots allowed

! Namelist files
Namelist / ion_params / num_ion_species, Z_ion_init, miomp_init
Namelist / run_params / input_is_Er, log_interp, use_quanc8, read_U2_file, &
  Add_Spitzer_to_D33, num_Er_test, numKsteps, kord_pprof, keord, kcord,    &
  Kmin, Kmax, epsabs, epsrel, Method, flux_cap, output_QoT_vs_Er, use_beam

!- End of header -------------------------------------------------------------

!-----------------------------------------------------------------------------
![1.0] Initialize: read command line args, allocate space, read input files
!                  display intro text, calculate thermal velocities and
!                  log(lambda)
!-----------------------------------------------------------------------------

! Read namelist file for ion parameters
Open(iu_nl,file="ion_params",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening ion_params namelist file'
  Stop 'Exiting: I/O Error in penta.f90'
Endif
Read(iu_nl,nml=ion_params)
Close(iu_nl)

! Read namelist file for run parameters
Open(iu_nl,file="run_params",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening run_params namelist file'
  Stop 'Exiting: I/O Error in penta.f90'
Endif
Read(iu_nl,nml=run_params)
Close(iu_nl)

! Get command line arguments
numargs=command_argument_count()
If ( numargs /= 9 ) Then
  Write(*,*) 'Incorrect number of input arguments, see penta.f90 for details'
  Stop 'Exiting: Input arguments error in penta.f90'
Endif
Call Getarg(1, arg1)
Call Getarg(2, arg2)
Call Getarg(3, arg3)
Call Getarg(4, arg4)
Call Getarg(5, arg5)
Call Getarg(6, arg6)
Call Getarg(7, arg7)
Call Getarg(8, arg8)
Call Getarg(9, arg9)

! Store command line args
coeff_ext = Trim(Adjustl(arg1))
Read(arg2,*) Er_min
Read(arg3,*) Er_max
Read(arg4,*) js
Read(arg5,*) i_append
run_ident = Trim(Adjustl(arg6))
pprof_char = Trim(Adjustl(arg7))
Read(arg8,*) B_Eprl
Read(arg9,*) Smax
  
! Allocate variables according to number of ion species defined
num_species = num_ion_species + 1_iknd
Allocate(ni(num_ion_species),Ti(num_ion_species))          ! Ion profile info
Allocate(dnidr(num_ion_species),dTidr(num_ion_species))
Allocate(Z_ion(num_ion_species),ion_mass(num_ion_species)) ! Ion parameters
Allocate(vth_i(num_ion_species))
Allocate(charges(num_species))                             ! Parameters for
Allocate(dens(num_species))                                !  all species
Allocate(vths(num_species)) 
Allocate(masses(num_species))
Allocate(Temps(num_species))
Allocate(dTdrs(num_species))
Allocate(dndrs(num_species))
Allocate(lmat((Smax+1)*num_species,(Smax+1)*num_species))  ! Clas. fric.coeffs
Allocate(Xvec(num_species*2+1),Avec(num_species*3))        ! Thermo. Force vecs.
Allocate(Flows((Smax+1)*num_species))                      ! Prl flow moments
Allocate(Gammas(num_species))                              ! Rad fluxes
Allocate(QoTs(num_species))                                ! Rad energy fluxes
Allocate(Gamma_i_vs_Er(num_Er_test,num_ion_species))       ! Ion flux vs Er
Allocate(Gamma_e_vs_Er(num_Er_test))                       ! Electron flux vs Er
Allocate(Er_test_vals(num_Er_test))                        ! Er to loop over
If ( output_QoT_vs_Er .EQV. .true. ) Then
  Allocate(QoT_i_vs_Er(num_Er_test,num_ion_species))       ! Ion flux vs Er
  Allocate(QoT_e_vs_Er(num_Er_test))                       ! Electron flux vs Er
Endif

! Read input files
Call read_vmec_file(js,run_ident)
Call read_pprof_file(pprof_char,num_ion_species,roa_surf,arad,kord_pprof)
Call read_dkes_star_files(coeff_ext,Add_Spitzer_to_D33,Bsq)
If (use_beam) Then
  Call read_beam_file(roa_surf,kord_pprof)
Else
  beam_force = 0._rknd
Endif

! Allocate DKES coefficients arrays
Allocate(xt_c(num_c + kcord))
Allocate(xt_e(num_e + keord))
Allocate(Dspl_D11(num_c,num_e))
Allocate(Dspl_D13(num_c,num_e))
Allocate(Dspl_D31(num_c,num_e))
Allocate(Dspl_D33(num_c,num_e))
Allocate(Dspl_Dex(num_c,num_e))
Allocate(Dspl_Drat(num_c,num_e)) 
Allocate(Dspl_Drat2(num_c,num_e)) 
Allocate(Dspl_DUa(num_c,num_e))  
Allocate(Dspl_logD11(num_c,num_e))
Allocate(Dspl_logD33(num_c,num_e))
Allocate(cmesh(num_c,num_e))

! Optionally read file containing <U**2> info.  Else this is 
! calculated from the D11* coefficient at high nu/v and Er=0.
If ( read_U2_file ) Then
  Call read_Utilde2_file(roa_surf,U2,kord_pprof)
Else
  U2=1.5d0*D11_mat(num_c,1)/cmul(num_c);
Endif

! Change Er test range to V/cm if necessary
If ( input_is_Er .EQV. .false.)  Then
  Er_min = Er_min * Te / arad
  Er_max = Er_max * Te / arad
Endif

! Assign ion parameters
Z_ion    = Z_ion_init(1:num_ion_species)
ion_mass = miomp_init(1:num_ion_species) * p_mass

! Calculate fitting parameters to the D##* coefficients
Call fit_coeffs(cmul,efield,num_c,num_e,D11_mat,log_interp, &
  kcord,keord,xt_c,xt_e,Dspl_D11,cmin,cmax,emin,emax)
Call fit_coeffs(cmul,efield,num_c,num_e,D13_mat,log_interp, &
  kcord,keord,xt_c,xt_e,Dspl_D13,cmin,cmax,emin,emax)
Call fit_coeffs(cmul,efield,num_c,num_e,D31_mat,log_interp, &
  kcord,keord,xt_c,xt_e,Dspl_D31,cmin,cmax,emin,emax)
Call fit_coeffs(cmul,efield,num_c,num_e,D33_mat,log_interp, &
  kcord,keord,xt_c,xt_e,Dspl_D33,cmin,cmax,emin,emax)

! Fit log(D*) for D11 and D33
Call fit_coeffs(cmul,efield,num_c,num_e,Log(D11_mat),log_interp, &
  kcord,keord,xt_c,xt_e,Dspl_logD11,cmin,cmax,emin,emax)
Call fit_coeffs(cmul,efield,num_c,num_e,Log(D33_mat),log_interp, &
  kcord,keord,xt_c,xt_e,Dspl_logD33,cmin,cmax,emin,emax)

! Display run data to screen (if i_append==0)
If ( i_append == 0 ) Then
  Write(*,*) 
  Write(*,*) "Welcome to PENTA3, please note the following settings:"
  Write(*,*)
  Write(*,'(a,i3)') ' Number of ion species: ',num_ion_species
  If ( input_is_Er .EQV. .true. ) Then
    Write(*,*) 'Interpreting input range as Er (V/cm)'
  Else
    Write(*,*) 'Interpreting input range as e<a>Er/kTe'
  Endif
  If ( log_interp .EQV. .true. ) Then
    Write(*,*) 'Performing logarithmic interpolation in Er, cmul'
  Else
    Write(*,*) 'Performing linear interpolation in Er,cmul'
  Endif
  If ( use_quanc8 .EQV. .true. ) Then
    Write(*,'(a,2(a,e10.4))')                           &
      ' Using quanc8 integrator with tolerances: ',     &
      'abs: ',epsabs,' rel: ', epsrel
  Else
    Write(*,'(a,i6,a)') ' Using ',numKsteps,            &
      ' point integral approximation'
  Endif
  Write(*,'(a,2(" ",e15.4))') ' K range on convolution integral: ', &
    Kmin, Kmax
  If ( Add_Spitzer_to_D33 .EQV. .true. ) Then
    Write(*,*) 'Adding collisional (Spitzer) portion to D33 coefficient'
  Endif
  If ( flux_cap .EQV. .true. ) Then
    Write(*,*) 'Enforcing minimum radial diffusion coefficient = 0'
  Endif
  Write(*,'(a,i2)') ' Number of terms in Sonine expansion: ', Smax+1
  Select Case (Method)
    Case ('T')
      Write(*,*) 'Using Taguchi Method'
    Case ('SN')
      Write(*,*) 'Using Sugama-Nishimura Method'
    Case ('MBT')
      Write(*,*) 'Using Maassberg-Beidler-Turkin Method'
    Case ('DKES')
      Write(*,*) 'Using DKES Method'
    Case Default
      Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
        ''' is not a valid Method'
      Stop 'Error: Exiting, method select error in penta.f90 (1)'
  EndSelect
      
  Write(*,*)
  Write(*,*) " <r>/<a>","   Er root(s) (V/cm)"
Endif ! intro text

! Set specifiers for opening output files
If (i_append == 0) Then
  fstatus = "unknown"
  fpos = "asis"
ElseIf (i_append == 1) Then
  fstatus = "old"
  fpos = "append"
Else
  Write(*,*) 'Bad value for i_append (0 or 1 expected)'
  Stop 'Error: Exiting, i_append error in penta.f90'
EndIf

! Open output files
Open(unit=iu_flux_out, file="fluxes_vs_roa", &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_pprof_out, file="plasma_profiles_check",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_fvEr_out, file="fluxes_vs_Er",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_flows_out, file="flows_vs_roa",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_flowvEr_out, file="flows_vs_Er",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_Jprl_out,file="Jprl_vs_roa",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_contraflows_out,file="ucontra_vs_roa",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))

If ( output_QoT_vs_Er .EQV. .true. ) Then
  Open(unit=iu_QoTvEr_out, file="QoTs_vs_Er",  &
    position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
  Write(iu_QoTvEr_out,'("*",/,"r/a   Er[V/cm]   Q_e/T_e [m**-2s**-1] ",&
    & "   Q_i/T_i [m**-2s**-1]")')
Endif

! Write legends (if i_append = 0) for most files
If ( i_append == 0 ) Then
  ! Fluxes vs r/a
  Write(iu_flux_out,'("*",/,"r/a    Er[V/cm]    e<a>Er/kTe    ",  &
    & "Gamma_e [m**-2s**-1]   Q_e/T_e [m**-2s**-1]     ",         &
    & "Gamma_i [m**-2s**-1]   Q_i/T_i [m**-2s**-1]")')
  ! Flows vs r/a
  Write(iu_flows_out,'("*",/,"r/a   Er[V/cm]    e<a>Er/kTe    ",  &
    & " <B*u_||ke>/<B**2> [m/sT]   <B*u_||ki>/<B**2> [m/sT]")')
  ! Plasma profile check
  Write(iu_pprof_out,'("*",/,"r/a    Te [eV]   ne [m**-3]     ",  & 
    & "dnedr [m**-4]   dTedr [eV/m]  Ti [eV]   ni [m**-3]     ",  &
    & "dnidr [m**-4]   dTidr [eV/m]")')
  Write(iu_Jprl_out,'("*",/,"r/a    Er [V/cm]    e<a>Er/kTe    ",  &
    & "Jprl_e [A/m**2]    Jprli [A/m**2]    Jprl [A/m**2]")')
  Write(iu_contraflows_out,'("*",/,"r/a    Er [V/cm]    e<a>Er/kTe    ",  &
    & "<ue^pol_contra> [1/s]     <ue^tor_contra> [1/s]       ",  &
    & " <ui^pol_contra> [1/s]     <ui^tor_contra> [1/s]")')
EndIf

! Legend for fluxes vs Er is written for each surface
Write(iu_fvEr_out,'("*",/,"r/a   Er[V/cm]   Gamma_e [m**-2s**-1] ",&
  & "   Gamma_i [m**-2s**-1]")')
! Legend for flows vs Er is written for each surface
Write(iu_flowvEr_out,'("*",/,"r/a   Er[V/cm]  ", &
  & "    <B*u_||ke>/<B**2> [m/sT]  <B*u_||ki>/<B**2> [m/sT]")')

! Calculate thermal velocities 
vth_i = Dsqrt(2._rknd*Ti*elem_charge/ion_mass)
vth_e = Dsqrt(2._rknd*Te*elem_charge/e_mass)

! Calculate Coulomb logarithm
If ( Te > 50._rknd ) Then
  loglambda = 25.3_rknd - 1.15_rknd*Dlog10(ne/1.e6_rknd) + 2.3_rknd*Dlog10(Te)
Else
  loglambda = 23.4_rknd - 1.15_rknd*Dlog10(ne/1.e6_rknd) + 3.45_rknd*Dlog10(Te)
Endif

! Assign arrays of parameters for all species (charge, mass, n, T, v_th, dTdr)
charges=elem_charge*(/-1._rknd,Z_ion/)
dens=(/ne, ni/)
masses=(/e_mass, ion_mass/)
Temps=(/Te,Ti/)
vths=(/vth_e,vth_i/)
dTdrs=(/dTedr,dTidr/)
dndrs=(/dnedr,dnidr/)

! Define matrix of friction coefficients (lmat)
Call define_friction_coeffs(masses,charges,vths,Temps,dens,loglambda, &
                            num_species,Smax,lmat)

! Fit radial transport coefficients specific to different methods
SelectCase (Method)
  Case ('T', 'MBT')
    cmesh = Spread(cmul,2,num_e)
    ! Calculate the D11 coefficient minus the P-S contribution  (Dex)
    ! Also, do not allow for negative coefficients
    Call fit_coeffs(cmul,efield,num_c,num_e, &
      Max(D11_mat-(2._rknd/3._rknd)*cmesh*U2,0._rknd), &
      log_interp,kcord,keord,xt_c,xt_e,Dspl_Dex,     &
      cmin,cmax,emin,emax)

  Case ('SN')
    ! Calculate fits to D31*/D33*  (Drat)
    Call fit_coeffs(cmul,efield,num_c,num_e, &
      D31_mat/D33_mat, &
      log_interp,kcord,keord,xt_c,xt_e,Dspl_Drat,     &
      cmin,cmax,emin,emax)

    ! Calculate fits to (D31*)**2/D33*   (Drat2)
    Call fit_coeffs(cmul,efield,num_c,num_e, &
      D31_mat*D31_mat/D33_mat, &
      log_interp,kcord,keord,xt_c,xt_e,Dspl_Drat2,     &
      cmin,cmax,emin,emax)

    cmesh = Spread(cmul,2,num_e)
    ! Calculate coefficient for Ua term  (DUa)
    Call fit_coeffs(cmul,efield,num_c,num_e, &
      (2._rknd*Bsq/(3._rknd*D33_mat) - cmesh), &
      log_interp,kcord,keord,xt_c,xt_e,Dspl_DUa,     &
      cmin,cmax,emin,emax)

    ! Calculate coefficient for radial flux  (Capped term)  (Dex)
    ! Also, do not allow for negative coefficients
    Call fit_coeffs(cmul,efield,num_c,num_e, &
      Max(D11_mat-(2._rknd/3._rknd)*cmesh*U2+D31_mat*D31_mat/D33_mat, &
      0._rknd),log_interp,kcord,keord,xt_c,xt_e,Dspl_Dex,     &
      cmin,cmax,emin,emax)

  Case ('DKES')
  Case Default
    Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
      ''' is not a valid Method'
    Stop 'Error: Exiting, method select error in penta.f90 (2)'
EndSelect

!-----------------------------------------------------------------------------
![2.0] Efield loop: Loop over test values of Er and calculate flows and fluxes
!-----------------------------------------------------------------------------

! Define array of Er values to test [V/m]
Er_test_vals = rlinspace(Er_min,Er_max,num_Er_test)*100._rknd

! Check for Er=0, doesn't work for log interpolation
min_Er = Minval(Dabs(Er_test_vals),DIM=1) 
If ((log_interp .EQV. .true. ) .AND. ( Dabs(min_Er) <= elem_charge ))  Then
  min_ind = Minloc(Dabs(Er_test_vals),DIM=1)
  If ( min_ind == Num_Er_test ) Then 
    Er_test_vals(min_ind) = Er_test_vals(min_ind - 1)/2._rknd
  Else
    Er_test_vals(min_ind) = Er_test_vals(min_ind + 1)/2._rknd
  EndIf
  Write(*,'(a,i4,a,f10.3)') 'Cannot use Er=0 with log_interp, using Er(',  &
     min_ind, ') = ', Er_test_vals(min_ind)
EndIf

! Loop over Er to get fluxes as a function of Er
Do ie = 1,num_Er_test

  Er_test = Er_test_vals(ie)
  abs_Er = Abs(Er_test)

  ! Form thermodynamic force vector (Xvec)
  Call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

  ! Form alternate thermodynamic force vector (Avec)
  Do ispec1 = 1,num_species
    ind_X = (ispec1-1)*2 + 1
    ind_A = (ispec1-1)*3 + 1

    Avec(ind_A)   = -Xvec(ind_X) / (Temps(ispec1)*elem_charge) &
         - 2.5_rknd*dTdrs(ispec1)/Temps(ispec1)
    Avec(ind_A+1)   = -Xvec(ind_X+1) / (Temps(ispec1)*elem_charge)
    Avec(ind_A+2)   = Xvec(num_species*2+1)*charges(ispec1) &
         * B0/(Temps(ispec1)*elem_charge*Sqrt(Bsq)) + &
         beam_force/(Temps(ispec1)*elem_charge*dens(ispec1))
  Enddo

  ! Select the appropriate algorithm and calculate the flows and fluxes
  SelectCase (Method)

    Case ('T', 'MBT')
      ! Calculate array of parallel flow moments
      Flows = calc_flows_T(num_species,Smax,abs_Er,Temps,dens,vths,charges,   &
        masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,   &
        cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c,num_e,kcord,      &
        keord,Avec,Bsq,lmat)
      Gammas = calc_fluxes_MBT(num_species,Smax,abs_Er,Temps,dens,vths,       &
        charges,masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,  &
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,        &
        Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,Flows,U2,B0,flux_cap)   
      If ( output_QoT_vs_Er .EQV. .true. ) Then
        QoTs = calc_QoTs_MBT(num_species,Smax,abs_Er,Temps,dens,vths,charges, &
          masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
          log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,      &
          Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,Flows,U2,B0,flux_cap)   
      Endif    
    Case ('SN')
                           
      Flows = calc_flows_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,       &
         cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_DUa,num_c,num_e,kcord,  &
         keord,Avec,lmat)                                                
      Gammas = calc_fluxes_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,&
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,cmax, &
        emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,Dspl_Dex,Dspl_logD11,        &
        Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,lmat,Flows,U2,dTdrs,        &
        dndrs,flux_cap)  
      If ( output_QoT_vs_Er .EQV. .true. ) Then
        QoTs = calc_QoTs_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
          masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,    &
          cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,Dspl_Dex,Dspl_logD11, &
          Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,lmat,Flows,U2,dTdrs,      &
          dndrs,flux_cap)  
      Endif    
    Case ('DKES')
      Flows = calc_flows_DKES(num_species,Smax,abs_Er,Temps,dens,vths,charges,&
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,  &
         cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c,num_e,kcord,     &
         keord,Avec)
      Gammas = calc_fluxes_DKES(num_species,abs_Er,Temps,dens,vths,charges,   &
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,cmax, &
        emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,keord,     &
        Avec,B0)   
      If ( output_QoT_vs_Er .EQV. .true. ) Then
        QoTs = calc_QoTs_DKES(num_species,abs_Er,Temps,dens,vths,charges,     &
          masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,    &
          cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,    &
          keord,Avec,B0)  
      Endif
    Case Default
      Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
        ''' is not a valid Method'
      Stop 'Error: Exiting, method select error in penta.f90 (3)'
  EndSelect

  Gamma_e_vs_Er(ie)   = Gammas(1)
  Gamma_i_vs_Er(ie,:) = Gammas(2:num_species)

  ! Write fluxes vs Er
  Write(str_num,*) num_ion_species + 2  ! Convert num to string
  Write(iu_fvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,Gamma_e_vs_Er(ie),Gamma_i_vs_Er(ie,:)

  If ( output_QoT_vs_Er .EQV. .true. ) Then
    QoT_e_vs_Er(ie)   = QoTs(1)
    QoT_i_vs_Er(ie,:) = QoTs(2:num_species)
    Write(iu_QoTvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
      roa_surf,Er_test/100._rknd,QoT_e_vs_Er(ie),QoT_i_vs_Er(ie,:)
  Endif

  ! Write flows vs Er
  Write(str_num,*) (Smax+1)*num_species + 2  ! Convert num to string
  Write(iu_flowvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,Flows

Enddo !efield loop


!-----------------------------------------------------------------------------
![3.0] Find ambipolar roots and evaluate quantities at these roots.
!-----------------------------------------------------------------------------

! Check for only one Er test value -- this is then used to evaluate the ambipolar fluxes QQ
!If ( num_Er_test  == 1 ) Then
!  Er_roots = Er_test_vals

! Find the ambipolar root(s) from gamma_e = sum(Z*gamma_i)
Call find_Er_roots(gamma_e_vs_Er,gamma_i_vs_Er,Er_test_vals,Z_ion, &
  num_Er_test,num_ion_species,Er_roots,num_roots)

! Allocate arrays according to number of ambipolar roots
Allocate(Flows_ambi((Smax+1)*num_species,num_roots)) ! Parallel flow moments
Allocate(Gammas_ambi(num_species,num_roots))         ! Rad. particle fluxes
Allocate(QoTs_ambi(num_species,num_roots))           ! Rad. energy fluxes
Allocate(Jprl_ambi(num_roots))                       ! Parallel current density
Allocate(Jprl_parts(num_species,num_roots))          ! Par. curr. dens. per spec.
Allocate(upol(num_species,num_roots))                ! fsa contra pol flow
Allocate(utor(num_species,num_roots))                ! fsa contra tor flow

! Evaluate fluxes and flows at the ambipolar Er
Do iroot = 1_iknd, num_roots

  Er_test = Er_roots(iroot)
  abs_Er = Dabs(Er_test)

  ! Form thermodynamic force vector (Xvec)
  Call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

  ! Form alternate thermodynamic force vector (Avec)
  Do ispec1 = 1,num_species
    ind_X = (ispec1-1)*2 + 1
    ind_A = (ispec1-1)*3 + 1

    Avec(ind_A)   = -Xvec(ind_X) / (Temps(ispec1)*elem_charge) &
         - 2.5_rknd*dTdrs(ispec1)/Temps(ispec1)
    Avec(ind_A+1)   = -Xvec(ind_X+1) / (Temps(ispec1)*elem_charge)
    Avec(ind_A+2)   = Xvec(num_species*2+1)*charges(ispec1) &
         * B0/(Temps(ispec1)*elem_charge*Sqrt(Bsq))
  Enddo

  ! Select the appropriate algorithm and calculate the flows and fluxes
  SelectCase (Method)


    Case ('T', 'MBT')
      ! Calculate array of parallel flow moments 
        ! Note: Flow methods are the same for T and MBT
      Flows_ambi(:,iroot) = calc_flows_T(num_species,Smax,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,     &
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c, &
        num_e,kcord,keord,Avec,Bsq,lmat)
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_MBT(num_species,Smax,abs_Er,Temps,  &
        dens,vths,charges,masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax, &
        numKsteps,log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,      &
        Dspl_D31,Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,                 &
        Flows_ambi(:,iroot),U2,B0,flux_cap)   
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_MBT(num_species,Smax,abs_Er,Temps,dens, &
        vths,charges,masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,      &
        numKsteps,log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,      &
        Dspl_D31,Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,                 &
        Flows_ambi(:,iroot),U2,B0,flux_cap)   

    Case ('SN')
      ! Calculate array of parallel flow moments
                                              
      Flows_ambi(:,iroot) = calc_flows_SN(num_species,Smax,abs_Er,Temps,dens,&
         vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,    &
         log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_DUa,num_c,  &
         num_e,kcord,keord,Avec,lmat)                                                
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_SN(num_species,Smax,abs_Er,Temps,   &
        dens,vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,   &
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,       &
        Dspl_Dex,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,      &
        lmat,Flows_ambi(:,iroot),U2,dTdrs,dndrs,flux_cap)  
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_SN(num_species,Smax,abs_Er,Temps,dens,  &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,       &
        Dspl_Dex,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,      &
        lmat,Flows_ambi(:,iroot),U2,dTdrs,dndrs,flux_cap)  

    Case ('DKES')

      ! Calculate array of parallel flow moments 
      Flows_ambi(:,iroot) = calc_flows_DKES(num_species,Smax,abs_Er,Temps,   &
        dens,vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,&
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c, &
        num_e,kcord,keord,Avec)
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_DKES(num_species,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c, &
        num_e,kcord,keord,Avec,B0)  
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_DKES(num_species,abs_Er,Temps,dens,     &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
        log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c, &
        num_e,kcord,keord,Avec,B0)  
        

    Case Default
      Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
        ''' is not a valid Method'
      Stop 'Error: Exiting, method select error in penta.f90 (4)'
  EndSelect


  ! Calculate parallel current density
  Jprl_parts(:,iroot) = dens*charges*Sqrt(Bsq) *         &
    Flows_ambi(1:(num_species-1)*(Smax+1)+1:Smax+1,iroot)
  Jprl_ambi(iroot) = Sum(Jprl_parts(:,iroot))

  ! Calculate flow components
  upol(:,iroot) = chip*(4._rknd*pi*pi/vol_p)*(                 &
    Flows_ambi(1:(num_species-1)*(Smax+1)+1:Smax+1,iroot) - &
    bzeta*Xvec(1:(num_species-1)*2+1:2)/(charges*chip*Bsq))
  utor(:,iroot) = psip*(4._rknd*pi*pi/vol_p)*(                 &
    Flows_ambi(1:(num_species-1)*(Smax+1)+1:Smax+1,iroot) + &
    btheta*Xvec(1:(num_species-1)*2+1:2)/(charges*psip*Bsq))

EndDo ! Ambipolar root loop



!-----------------------------------------------------------------------------
![4.0] Write output files
!-----------------------------------------------------------------------------

! Loop over ambipolar Er for writing output files
Do iroot = 1_iknd, num_roots

  Er_test = Er_roots(iroot)
  eaEr_o_kTe = arad*Er_test/Te

  ! Write fluxes to file "fluxes_vs_roa"
  Write(str_num,*) 2*num_species + 2
  Write(iu_flux_out,'(f7.3,' // Trim(Adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,eaEr_o_kTe,Gammas_ambi(1,iroot),  &
    QoTs_ambi(1,iroot),Gammas_ambi(2:num_species,iroot),  &
    QoTs_ambi(2:num_species,iroot)

  ! Write flows to file "flows_vs_roa"
  Write(str_num,*) (Smax+1)*num_species + 2
  Write(iu_flows_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,eaEr_o_kTe,Flows_ambi(:,iroot)

  ! Write current densities to file "Jprl_vs_roa"
  Write(str_num,*) num_species + 3
  Write(iu_Jprl_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))')  & 
    roa_surf,Er_test/100._rknd,eaEr_o_kTe,Jprl_parts(:,iroot),Jprl_ambi(iroot)

  ! Write contravariant flows to file "ucontra_vs_roa"
  Write(str_num,*) 2*num_species + 2
  Write(iu_contraflows_out,'(f7.3,' // trim(adjustl(str_num))//'(" ",e15.7))') & 
    roa_surf,Er_test/100._rknd,eaEr_o_kTe,upol(1,iroot),utor(1,iroot),         &
    upol(2:num_species,iroot),utor(2:num_species,iroot)

EndDo ! Ambipolar root loop

! Write plasma profile information to "plasma_profiles_check"
Write(str_num,*) 4*num_species
Write(iu_pprof_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') & 
  roa_surf,Te,ne,dnedr,dTedr,Ti,ni,dnidr,dTidr

! QQ write file with number of roots per surface!

! Write screen output
write(str_num,*) num_roots
write(*,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.4))') & 
  roa_surf,er_roots(1:num_roots)/100._rknd

!-----------------------------------------------------------------------------
![4.0] Cleanup -- Deallocate variables and close files
!-----------------------------------------------------------------------------
! Deallocate variables
Deallocate(ni,Ti,dnidr,dTidr)                     ! Ion profile info
Deallocate(Z_ion,ion_mass)                        ! Ion parameters
Deallocate(cmul,efield)                   ! DKES coeff info
Deallocate(D11_mat,D13_mat,D31_mat,D33_mat)
Deallocate(Temps,masses,vths,charges,dens,dTdrs,dndrs)  ! All species parameters
Deallocate(lmat)                                  ! Class. fric. coeffs.
Deallocate(Xvec,Avec)                             ! Thermo. Force Vecs
Deallocate(Flows)                                 ! Prl flows
Deallocate(Gammas,QoTs)                                ! Rad fluxes
Deallocate(xt_c,xt_e)   ! DKES spline fitting info
Deallocate(Dspl_D11)
Deallocate(Dspl_D13)
Deallocate(Dspl_D31)
Deallocate(Dspl_D33)
Deallocate(Dspl_Dex)
Deallocate(Dspl_Drat)
Deallocate(Dspl_DUa)
Deallocate(cmesh)     ! Spread array of cmul data
Deallocate(Gamma_i_vs_Er) ! Ion particle flux vs Er
Deallocate(Gamma_e_vs_Er) ! Electron particle flux vs Er
If ( output_QoT_vs_Er .EQV. .true. ) Then
  Deallocate(QoT_i_vs_Er)
  Deallocate(QoT_e_vs_Er)
Endif
Deallocate(Gammas_ambi)   ! Ambipolar quantities
Deallocate(QoTs_ambi)
Deallocate(Flows_ambi)
Deallocate(Er_test_vals)  ! Er to loop over
Deallocate(Jprl_ambi,Jprl_parts) ! Parallel current densities
Deallocate(utor,upol) ! Contravariant fsa flows

! Close output files
Close(iu_flux_out)
Close(iu_pprof_out)
Close(iu_fvEr_out)
Close(iu_QoTvEr_out)
Close(iu_flows_out)
Close(iu_flowvEr_out)
Close(iu_Jprl_out)
Close(iu_contraflows_out)

End program penta3

