!     SPH0101408: REPLACED INTEGER(iprec) WITH INTEGER
!*******************************************************************************
!  File eq_T.f
!  Contains module eq_T
!  Defines derived-types: eq_param_fix, eq_param_var, eq_state

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  It deals with the interface with the equilibrium code.
!  Right now (6/2004), the VMEC code is the only equilibrium code.

!  Issues: (As of 9-21-20004) JDH
!  1) presence of pointers to other eq_ structures within the structures.
!     (Does an eq_param_var need to have a pointer to an eq_param_fix?)
!  2) 'Use as allocatable-array' or 'Use as Pointer' status of various arrays.
!  3) Special usage for the eq_state, since it can point to vmec storage space.
!  4) Assignment

!  Other Modules: Module eq_interface deals with the interface between V3FIT 
!  and VMEC.

!*******************************************************************************
!  MODULE eq_T
!    (Equilibrium types)
! SECTION I.      VARIABLE DECLARATIONS
! SECTION II.     DERIVED TYPE DECLARATIONS
! SECTION III.    INTERFACE BLOCKS
! SECTION IV.     CONSTRUCTION SUBROUTINES
! SECTION V.      DESTRUCTION SUBROUTINES
! SECTION VI.     ASSIGNMENT SUBROUTINES
! SECTION VII.    AUX DEFINITION SUBROUTINES
! SECTION VIII.   OTHER SUBROUTINES
!*******************************************************************************
      MODULE eq_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------
      USE stel_kinds

!-------------------------------------------------------------------------------
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_constants

!-------------------------------------------------------------------------------
!  Use Statements for  V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities
      
      IMPLICIT NONE

!*******************************************************************************
! SECTION II. DERIVED TYPE (STRUCTURE) DECLARATIONS
!
!   Derived Type to carry information about the fixed equilibrium parameters
!   These are discrete variables like logicals, array lengths, etc.
!   They are necessary for the initialization, and CAN NOT be changed
!   (without reinitializing the equilibrium solver)
!         eq_param_fix
!
!   Derived Type to carry information about the variable equilibrium parameters
!         eq_param_var
!
!  Derived Type to carry information about the first set of auxiliary
!    equilibrium variables. These are variables that 
!    1) are written in the wout file, or are computed in module read_wout_mod
!    2) are needed for the signal computation for magnetic diagnostics
!         eq_aux_1
!
!  Derived Type to carry information about the second set of auxiliary
!    equilibrium variables. These are physical variables, on an suv 
!    (VMEC spatial variable) grid.
!    They are computed from eq_aux_1 variables, and they are needed for
!    computation of magnetic diagnostic signals.
!         eq_aux_2
!
!   Derived Type to carry information about the equilibrium state
!         eq_state
!
!*******************************************************************************
!
!  
!-------------------------------------------------------------------------------
!  Declare type eq_param_fix
!    Variables for code='VMEC2000'
!                      VMEC:'module declared in':'file name'                                                     
!       lasym      logical, .true., stellarator asymmetric 
!                           .false. stellarator symmetric
!			           VMEC:vmec_input:vmec_input.f
!       lfreeb     logical, .true., use mgrid_file to do free-bdy calc
!                           .false. fixed-bdy
!			           VMEC:vmec_input:vmec_input.f
!                        Overridden (to .false.) if mgrid_file does not exist
!       lforbal    logical, .true., improved avg force bal 
!                           .false. variational force bal
!			           VMEC:vmec_input:vmec_input.f
!       lRFP      logical, .true., Reversed Field Pinch       ! Added 2010-12-04
!                           .false., tokamaks, stellarators
!			           VMEC:vmec_input:vmec_input.f
!       mgrid_file character, name for mgrid file.
! 			           VMEC:vmec_input:vmec_input.f
!       nfp        integer, number of toroidal field periods ( =1 for Tokamak)
! 			           VMEC:vmec_input:vmec_input.f
!       ncurr      integer, flux conserving (=0) or 
! 			                prescribed toroidal current (=1)
! 			           VMEC:vmec_input:vmec_input.f
!       nvacskip   integer, iterations skipped between full update of vacuum 
!                  solution
! 			           VMEC:vmec_input:vmec_input.f
!       mpol       integer, number of poloidal modes
! 			           VMEC:vmec_input:vmec_input.f
!       ntor       integer, number of toroidal modes
! 			           VMEC:vmec_input:vmec_input.f
!       ntheta     integer, number of poloidal points, 0 - 2pi
! 			           VMEC:vmec_input:vmec_input.f or computed in read_indata
!       nzeta      integer, number of toroidal points per field period
! 			           VMEC:vmec_input:vmec_input.f
!                        VMEC:vmec_input:vmec_input.f or computed in read_indata
!       mns        integer, length of Xmn arrays 
!			           VMEC:computed in vmec_dim:vmec_dim.f
!       ns_array   integer array of radial (multi-) grids
!			           VMEC:vmec_dim:vmec_dim.f
!       nextcur    integer - number of external currents to use (extcur array)
!                      VMEC:mgrid_mod.f
!       pcurr_type Character, specify type of current profile specification
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       piota_type Character, specify type of iota or q profile specification
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       pmass_type Character, specify type of pressure profile specification
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!
!    Run control parameters
!
!       niter      integer, maximum number of iterations for convergence (VMEC 
!                  will retry once if ftol convergence level not achieved yet)
!			           VMEC:vmec_dim:vmec_dim.f
!       nstep      integer, iterations between printout of equilibrium data to 
!                      threed1 file/screen
!			           VMEC:vmec_dim:vmec_dim.f
!       delt       real(rprec), initial guess (~1) for time step
!			           VMEC:vmec_dim:vmec_dim.f
!       ftol_array real(rprec) array of force convergence tolerances on radial 
!                      grids defined by ns_array
!			           VMEC:vmec_dim:vmec_dim.f
!       gamma      real, value of compressibility index (gamma=0 =>                                                      
! 			           pressure prescribed)
! 			           VMEC:vmec_input:vmec_input.f
!       tcon0      real, optional parameter weight for spectral constraint 
!			           force (~1)
!			           VMEC:vmec_dim:vmec_dim.f
!
!-------------------------------------------------------------------------------
      TYPE eq_param_fix
!  Variables for code='VMEC2000'
         LOGICAL ::  lasym
         LOGICAL ::  lfreeb
         LOGICAL ::  lforbal
         LOGICAL ::  lRFP                                    ! Added 2010-12-04
         CHARACTER (len=100) :: mgrid_file
         INTEGER ::  nfp
         INTEGER ::  ncurr
         INTEGER ::  nvacskip
         INTEGER ::  mpol
         INTEGER ::  ntor
         INTEGER ::  ntheta
         INTEGER ::  nzeta
         INTEGER ::  mns
         INTEGER, DIMENSION(100) ::  ns_array
         INTEGER ::  nextcur
         CHARACTER(len=20) :: pcurr_type                      ! Added 2010-12-04
         CHARACTER(len=20) :: piota_type                      ! Added 2010-12-04
         CHARACTER(len=20) :: pmass_type                      ! Added 2010-12-04
         INTEGER ::  niter
         INTEGER ::  nstep
         REAL(rprec) ::  delt
         REAL(rprec), DIMENSION(100) ::  ftol_array
         REAL(rprec) ::  gamma
         REAL(rprec) ::  tcon0

      END TYPE eq_param_fix
!-------------------------------------------------------------------------------
!  Declare type eq_param_var
!    Variables for code='VMEC2000'
!     Integer Variables
!       ns_index   integer - index to the ns_array and ftol_array arrays.
!                      VMEC:runvmec
!                      Called multi_ns_grid(vmec_main.f) in readin.f and 
!                                                           eqsolve.f
!     Real Variables
!       ac         real array, (0:20)
!                      surface-averaged toroidal current density (ncurr=1)
!                      expansion coefficients (power series in s)
! 			           VMEC:vmec_input:vmec_input.f
!       ac_aux_s   real array (ndatafmax) (101)
!                      Auxiliary array _s_ for ac specification (splines)
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       ac_aux_f   real array (ndatafmax) (101)
!                      Auxiliary array _f_ for ac specification (splines)
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       ai         real array, (0:20)
!                      dimensionless iota (ncurr=0) expansion coefficients 
! 			           (power series in s)
!                      VMEC:vmec_input:vmec_input.f
!       ai_aux_s   real array (ndatafmax) (101)
!                      Auxiliary array _s_ for ai specification (splines)
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       ai_aux_f   real array (ndatafmax) (101)
!                      Auxiliary array _f_ for ai specification (splines)
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       am         real array, (0:20)
!                      mass or pressure (gamma=0) expansion coefficients 
!                      (power series in s) in MKS units (NWT/M**2) 
! 			           VMEC:vmec_input:vmec_input.f
!       am_aux_s   real array (ndatafmax) (101)
!                      Auxiliary array _s_ for am specification (splines)
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       am_aux_f   real array (ndatafmax) (101)
!                      Auxiliary array _f_ for am specification (splines)
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       bloat      real, used to expand the computation region.
! 			           VMEC:vmec_input:vmec_input.f           ! Added 2010-12-04
!       rbc        real array, (-ntord:ntord,0:mpol1d)
!                      boundary coefficients of COS(m*theta-n*zeta) for R
! 			           VMEC:vmec_input:vmec_input.f
!     ^ USE AS ALLOCATABLE ARRAY
!       zbs        real array, (-ntord:ntord,0:mpol1d)
!                      boundary coefficients of SIN(m*theta-n*zeta) for Z
! 			           VMEC:vmec_input:vmec_input.f
!     ^ USE AS ALLOCATABLE ARRAY
!       rbs        real array, (-ntord:ntord,0:mpol1d)
!                      boundary coefficients of SIN(m*theta-n*zeta) for R
! 			           VMEC:vmec_input:vmec_input.f
!     ^ USE AS ALLOCATABLE ARRAY
!       zbc        real array, (-ntord:ntord,0:mpol1d)
!                      boundary coefficients of COS(m*theta-n*zeta) for Z
! 			           VMEC:vmec_input:vmec_input.f
!     ^ USE AS ALLOCATABLE ARRAY
!       extcur     real array, (nigroup) [nigroup set to 100 in vsvd0.f]
!                      array of currents in each external current group. Used to
!                      multiply Green''s function for fields and loops read in 
!                      from MGRID file. Should use real current units (A).
! 			           VMEC:vmec_input:vmec_input.f
!       curtor     real, net toroidal current (used if ncurr=1) in MKS [A]
!			           VMEC:vmec_dim:vmec_dim.f
!       phiedge    real, value of toroidal flux at plasma edge (defines plasma 
!                      volume in free-bdy case)
! 			           VMEC:vmec_input:vmec_input.f
!       pres_scale real
!                      factor used to scale pressure profile (default value = 1)
!                      useful so user can fix profile and change beta without 
!                      having to change all AM coefficients separately
! 			           VMEC:vmec_input:vmec_input.f
!
! Modified 12/20/05 JMS :  Changed am, ai, ac size from (0:10), (0:20) to conform with
!                          updated version of VMEC2000
!
!-------------------------------------------------------------------------------
      TYPE eq_param_var
!  Variables for code='VMEC2000'
         INTEGER ::  ns_index
         REAL(rprec), DIMENSION(0:20) :: ac
         REAL(rprec), DIMENSION(101)  :: ac_aux_s             ! Added 2010-12-04
         REAL(rprec), DIMENSION(101)  :: ac_aux_f             ! Added 2010-12-04
         REAL(rprec), DIMENSION(0:20) :: ai
         REAL(rprec), DIMENSION(101)  :: ai_aux_s             ! Added 2010-12-04
         REAL(rprec), DIMENSION(101)  :: ai_aux_f             ! Added 2010-12-04
         REAL(rprec), DIMENSION(0:20) :: am
         REAL(rprec), DIMENSION(101)  :: am_aux_s             ! Added 2010-12-04
         REAL(rprec), DIMENSION(101)  :: am_aux_f             ! Added 2010-12-04
         REAL(rprec)                  :: bloat                ! Added 2010-12-04
         REAL(rprec), DIMENSION(:,:), POINTER :: rbc => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: zbs => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: zbc => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: rbs => null()
         REAL(rprec), DIMENSION(100) :: extcur
         REAL(rprec) :: curtor
         REAL(rprec) :: phiedge
         REAL(rprec) :: pres_scale
      END TYPE eq_param_var
!-------------------------------------------------------------------------------
!  Declare type eq_aux_1
!    Variables for code='VMEC'
! These are variables that were written to the wout file, or computed
!       in read_wout_file. Many are mn coefficients of physical variables.
!       Some are arrays on the full mesh.
!       They are needed for computation of signals.
!       
!  mgrid_mode     Character, VMECV:mgrid_mod.f
!  rmnc           cosmn component of cylindrical R, full mesh
!  zmns           sinmn component of cylindrical Z, full mesh
!  lmns           sinmn component of lambda, half mesh
!  gmnc           cosmn component of jacobian, half mesh
!  currumnc           
!  currvmnc           
!  bsubumnc       cosmn covariant u-component of B, half mesh
!  bsubvmnc       cosmn covariant v-component of B, half mesh
!  bsupumnc       cosmn contravariant u-component of B, half mesh
!  bsupvmnc       cosmn contravariant v-component of B, half mesh
!  rmns           sinmn component of cylindrical R, full mesh
!  zmnc           cosmn component of cylindrical Z, full mesh
!  lmnc           cosmn component of lambda, half mesh
!  gmns           sinmn component of jacobian, half mesh
!  currumns           
!  currvmns           
!  bsubumns       sinmn covariant u-component of B, half mesh
!  bsubvmns       sinmn covariant v-component of B, half mesh 
!  bsupumns       sinmn covariant u-component of B, half mesh
!  bsupvmns       sinmn covariant v-component of B, half mesh 
!  xm             Poloidal mode numbers
!  xn             Toroidal mode numbers
!  xm_nyq         Poloidal mode numbers (Nyquist)
!  xn_nyq         Toroidal mode numbers (Nyquist)
!  phi            Toroidal flux on full mesh
!  iotaf          iota on the full mesh
!  mnmax          max size of mn indices for arrays
!  mnmax_nyq      max size of Nyquist mn indices for arrays
!  ns
!  
!-------------------------------------------------------------------------------
      TYPE eq_aux_1
!  Variables for code='VMEC'
         CHARACTER(LEN=1)                     :: mgrid_mode
!    Symmetric
         REAL(rprec), DIMENSION(:,:), POINTER ::  rmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  zmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  lmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  gmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  currumnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  currvmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsubumnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsubvmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupumnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupvmnc => null()
!   Asymmetric (lasym = .TRUE.)
         REAL(rprec), DIMENSION(:,:), POINTER ::  rmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  zmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  lmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  gmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  currumns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  currvmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsubumns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsubvmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupumns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupvmns => null()
!    m and n value arrays
         REAL(rprec), DIMENSION(:), POINTER   ::  xm => null() 
         REAL(rprec), DIMENSION(:), POINTER   ::  xn => null() 
         REAL(rprec), DIMENSION(:), POINTER   ::  xm_nyq => null() 
         REAL(rprec), DIMENSION(:), POINTER   ::  xn_nyq => null()
!    Full Mesh Variables
         REAL(rprec), DIMENSION(:), POINTER   ::  phi => null()
         REAL(rprec), DIMENSION(:), POINTER   ::  iotaf => null()
!     max sizes of mn indices for arrays (added 8/3/07 by JS)
         INTEGER                          :: mnmax
         INTEGER                          :: mnmax_nyq
         INTEGER                          :: ns

!   Contravariant B field arrays (added 8/3/07)

      END TYPE eq_aux_1

!-------------------------------------------------------------------------------
!  Declare type eq_aux_2
!    Variables for code='VMEC'
!       FILL IN HERE
!       These are physical variables, on an suv (VMEC spatial variable) grid.
!       They are computed from eq_aux_1 variables, and they are needed for
!       computation of signals.
!  phitot      total toroidal flux. Defined from the phi array in aux_1.
!     ? Ask Steve Hirshman, is this ever different from phiedge in varp?
!-------------------------------------------------------------------------------
      TYPE eq_aux_2
!  Variables for code='VMEC'
!    Variables at plasma edge, for surface computation
!    First dimension should be of length 2 for these
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  bsubu => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  bsubv => null()
!   Variables throughout the volume
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  rsuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  zsuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  gsuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  currusuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  currvsuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  rusuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  zusuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  rvsuv => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER ::  zvsuv => null()
! Total toroidal flux
         REAL(rprec)                            :: phitot 

! Nyquist indices used by read_wout_mod.  Added by JMS 8/3/07
         INTEGER                            :: mnyq
         INTEGER                            :: nnyq

      END TYPE eq_aux_2
!-------------------------------------------------------------------------------
!  Declare type eq_state
!     code            character, name of the equilibrium code
!     version         character, version number of the equilibrium code
!     s_id            character, short identification of the equilibrium state
!     l_id            character, long identification of the equilibrium state
!     fixp            derived type eq_param_fix. Fixed Parameters of eq_state
!     varp            derived type eq_param_var. Variable Parameters of eq_state
!     l_def_aux1      logical, .true. if aux1 component is defined
!     aux1            derived type eq_aux_1. Auxiliary Variables
!     l_def_aux2      logical, .true. if aux2 component is defined
!     aux2            derived type eq_aux_2. Auxiliary Variables
!     xc              Packed array of geometry coefficients
!     ^ USE AS ALLOCATABLE ARRAY
!     xcdot           Time derivative of xc
!     ^ USE AS ALLOCATABLE ARRAY
!     wout_filename   name of wout file, used for computing eq_aux structures
!     fsq_max          current convergence parameter, compared with ftol array
!         10-4-04, JDH. Currently VMEC checks for convergence by requiring
!            EACH of fsqr, fsqz, and fsql to be less than the ftol_array value.
!            (See subroutine evolve.f). These values are printed out in the
!            threed file as FSQR, FSQZ, and FSQL. Here I store the maximum
!            of those three values
!-------------------------------------------------------------------------------
      TYPE eq_state
         CHARACTER (len=8)   :: code
         CHARACTER (len=20)  :: version
         CHARACTER (len=30)  :: s_id                                 
         CHARACTER (len=80)  :: l_id
!  Component Structures
         TYPE (eq_param_fix) :: fixp
         TYPE (eq_param_var) :: varp
         LOGICAL             :: l_def_aux1
         TYPE (eq_aux_1)     :: aux1
         LOGICAL             :: l_def_aux2
         TYPE (eq_aux_2)     :: aux2
!  Variables for code='VMEC2000'
         REAL(rprec), DIMENSION(:), POINTER :: xc => null()
         REAL(rprec), DIMENSION(:), POINTER :: xcdot => null()
         CHARACTER(LEN=80)   :: wout_filename
         REAL(rprec)         :: fsq_max

      END TYPE eq_state

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for eq_param_var
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE eq_param_var_assign, eq_aux_1_assign,                &
     &               eq_aux_2_assign, eq_state_assign
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE eq_construct
         MODULE PROCEDURE eq_param_fix_construct,                              &
     &               eq_param_var_construct, eq_aux_1_construct,               &
     &               eq_aux_2_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE eq_destroy
         MODULE PROCEDURE eq_param_fix_destroy,                                &
     &               eq_param_var_destroy, eq_aux_1_destroy,                   &
     &               eq_aux_2_destroy
      END INTERFACE

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct an eq_param_fix
!
! Needs to be cleaned up when the contents of an eq_param_fix becomes more
! certain.
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_param_fix_construct(this,                                  &
     &   lasym, lfreeb, lforbal, lRFP, mgrid_file, nfp, ncurr,                 &
     &   nvacskip, mpol, ntor, ntheta, nzeta, mns, ns_array, nextcur,          &
     &   pcurr_type, piota_type, pmass_type,                                   &
     &   niter, nstep, delt, ftol_array, gamma, tcon0)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_param_fix), INTENT (inout)   :: this
      LOGICAL, INTENT(in)                   :: lasym
      LOGICAL, INTENT(in)                   :: lfreeb
      LOGICAL, INTENT(in)                   :: lforbal
      LOGICAL, INTENT(in)                   :: lRFP
      CHARACTER (len=*), INTENT(in)         :: mgrid_file
      INTEGER, INTENT(in)                   :: nfp
      INTEGER, INTENT(in)                   :: ncurr
      INTEGER, INTENT(in)                   :: nvacskip
      INTEGER, INTENT(in)                   :: mpol
      INTEGER, INTENT(in)                   :: ntor
      INTEGER, INTENT(in)                   :: ntheta
      INTEGER, INTENT(in)                   :: nzeta
      INTEGER, INTENT(in)                   :: mns
      INTEGER, INTENT(in), DIMENSION(:)     :: ns_array
      INTEGER, INTENT(in)                   :: nextcur
      CHARACTER(len=20), INTENT(in)         :: pcurr_type
      CHARACTER(len=20), INTENT(in)         :: piota_type
      CHARACTER(len=20), INTENT(in)         :: pmass_type
      INTEGER, INTENT(in)                   :: niter
      INTEGER, INTENT(in)                   :: nstep
      REAL(rprec), INTENT(in)               :: delt
      REAL(rprec), INTENT(in), DIMENSION(:) :: ftol_array   
      REAL(rprec), INTENT(in)               :: gamma
      REAL(rprec), INTENT(in)               :: tcon0

!  Declare local variables
      INTEGER :: nsd, nfd
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_param_fix_construct: '

!  Start of executable code

!      WRITE(*,*) ' Executing eq_construct_fixp'

      this % lasym   = lasym
      this % lfreeb  = lfreeb
      this % lforbal = lforbal
      this % lRFP = lRFP
      this % mgrid_file = TRIM(ADJUSTL(mgrid_file))
      this % nfp     = nfp
      this % ncurr   = ncurr
      this % nvacskip = nvacskip
      this % mpol    = mpol
      this % ntor    = ntor
      this % ntheta  = ntheta 
      this % nzeta   = nzeta
      this % mns   = mns
      CALL assert_eq(SIZE(ns_array),SIZE(this % ns_array),
     &   sub_name // 'ns array sizes unequal','Warn')
      nsd = MIN(SIZE(ns_array),SIZE(this % ns_array))
      this % ns_array(1:nsd) = ns_array(1:nsd)
      this % nextcur   = nextcur
      this % pcurr_type = pcurr_type 
      this % piota_type = piota_type 
      this % pmass_type = pmass_type 
      this % niter   = niter
      this % nstep   = nstep
      this % delt    = delt
      CALL assert_eq(SIZE(ftol_array),SIZE(this % ftol_array),
     &   sub_name // 'ftol array sizes unequal','Warn')
      nfd = MIN(SIZE(ftol_array),SIZE(this % ftol_array))
      this % ftol_array(1:nfd) = ftol_array(1:nfd)
      this % tcon0   = tcon0
      this % gamma   = gamma
      
      END SUBROUTINE eq_param_fix_construct

!-------------------------------------------------------------------------------
!  Construct an eq_param_var
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_param_var_construct(this, ns_index,                        &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_param_var), INTENT (inout)     :: this
      INTEGER, INTENT(in)                     :: ns_index
      REAL(rprec), DIMENSION(:), INTENT(in)   :: ac
      REAL(rprec), DIMENSION(:), INTENT(in)   :: ac_aux_s
      REAL(rprec), DIMENSION(:), INTENT(in)   :: ac_aux_f
      REAL(rprec), DIMENSION(:), INTENT(in)   :: ai
      REAL(rprec), DIMENSION(:), INTENT(in)   :: ai_aux_s
      REAL(rprec), DIMENSION(:), INTENT(in)   :: ai_aux_f
      REAL(rprec), DIMENSION(:), INTENT(in)   :: am
      REAL(rprec), DIMENSION(:), INTENT(in)   :: am_aux_s
      REAL(rprec), DIMENSION(:), INTENT(in)   :: am_aux_f
      REAL(rprec), INTENT(in)                 :: bloat
      REAL(rprec), DIMENSION(:,:), INTENT(in) :: rbc
      REAL(rprec), DIMENSION(:,:), INTENT(in) :: zbs
      REAL(rprec), DIMENSION(:,:), INTENT(in) :: zbc
      REAL(rprec), DIMENSION(:,:), INTENT(in) :: rbs
      REAL(rprec), DIMENSION(:), INTENT(in)   :: extcur
      REAL(rprec), INTENT(in)                 :: curtor
      REAL(rprec), INTENT(in)                 :: phiedge
      REAL(rprec), INTENT(in)                 :: pres_scale

!  Declare local variables
      INTEGER :: nd1, nd2, ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_param_var_construct: '

!  Start of executable code

!      WRITE(*,*) ' Executing eq_construct_varp'

!  The eq_param_var type contains pointers to arrays. Therefore,
!  to deallocate space, need to destroy this
      CALL eq_param_var_destroy(this)

      this % ns_index = ns_index

!  For ac, ai, and am do vector assignment. Should get correct lower limit.
      CALL assert_eq(SIZE(ac),SIZE(this % ac),                                 &
     &   sub_name // 'ac array sizes unequal','Warn')
      this % ac = ac

      CALL assert_eq(SIZE(ac_aux_s),SIZE(this % ac_aux_s),                     &
     &   sub_name // 'ac_aux_s array sizes unequal','Warn')
      nd1 = MIN(SIZE(ac_aux_s),SIZE(this % ac_aux_s))
      this % ac_aux_s(1:nd1) = ac_aux_s(1:nd1)

      CALL assert_eq(SIZE(ac_aux_f),SIZE(this % ac_aux_f),                     &
     &   sub_name // 'ac_aux_f array sizes unequal','Warn')
      nd1 = MIN(SIZE(ac_aux_f),SIZE(this % ac_aux_f))
      this % ac_aux_f(1:nd1) = ac_aux_f(1:nd1)

      CALL assert_eq(SIZE(ai),SIZE(this % ai),                                 &
     &   sub_name // 'ai array sizes unequal','Warn')
      this % ai = ai

      CALL assert_eq(SIZE(ai_aux_s),SIZE(this % ai_aux_s),                     &
     &   sub_name // 'ai_aux_s array sizes unequal','Warn')
      nd1 = MIN(SIZE(ai_aux_s),SIZE(this % ai_aux_s))
      this % ai_aux_s(1:nd1) = ai_aux_s(1:nd1)

      CALL assert_eq(SIZE(ai_aux_f),SIZE(this % ai_aux_f),                     &
     &   sub_name // 'ai_aux_f array sizes unequal','Warn')
      nd1 = MIN(SIZE(ai_aux_f),SIZE(this % ai_aux_f))
      this % ai_aux_f(1:nd1) = ai_aux_f(1:nd1)

      CALL assert_eq(SIZE(am),SIZE(this % am),                                 &
     &   sub_name // 'am array sizes unequal','Warn')
      this % am = am

      CALL assert_eq(SIZE(am_aux_s),SIZE(this % am_aux_s),                     &
     &   sub_name // 'am_aux_s array sizes unequal','Warn')
      nd1 = MIN(SIZE(am_aux_s),SIZE(this % am_aux_s))
      this % am_aux_s(1:nd1) = am_aux_s(1:nd1)

      CALL assert_eq(SIZE(am_aux_f),SIZE(this % am_aux_f),                     &
     &   sub_name // 'am_aux_f array sizes unequal','Warn')
      nd1 = MIN(SIZE(am_aux_f),SIZE(this % am_aux_f))
      this % am_aux_f(1:nd1) = am_aux_f(1:nd1)
      
      this % bloat = bloat

      ALLOCATE(this % rbc(SIZE(rbc,1),SIZE(rbc,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rbc')
      this % rbc = rbc

      ALLOCATE(this % zbs(SIZE(zbs,1),SIZE(zbs,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc zbs')
      this % zbs = zbs

      ALLOCATE(this % zbc(SIZE(zbc,1),SIZE(zbc,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc zbc')
      this % zbc = zbc

      ALLOCATE(this % rbs(SIZE(rbs,1),SIZE(rbs,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rbs')
      this % rbs = rbs

      CALL assert_eq(SIZE(extcur),SIZE(this % extcur),                         &
     &   sub_name // 'extcur array sizes unequal','Warn')
      nd1 = MAX(SIZE(extcur),SIZE(this % extcur))
      this % extcur(1:nd1) = extcur(1:nd1)
         
      this % curtor = curtor
      this % phiedge = phiedge
      this % pres_scale = pres_scale
      
      END SUBROUTINE eq_param_var_construct

!-------------------------------------------------------------------------------
!  Construct an eq_aux_1
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_1_construct(this,state)
!
!  See also subroutine eq_aux_1_define, which calls this subroutine.
!
!  08-03-07  Added contravariant arrays bsupumnc, bsupvmnc, etc.  JMS
!
!  10-15-2004
!  At this time, I'll get the information by reading the wout file. As long as
!  the file writing and reading overhead is small, this will save having to 
!  duplicate the coding in VMEC that converts from the state variables
      
      USE read_wout_mod, only : read_wout_file, rmnc, zmns, lmns,              &
     &   gmnc, currumnc, currvmnc, bsubumnc, bsubvmnc, lasym, rmns,            &
     &   zmnc, lmnc, gmns, currumns, currvmns, bsubumns, bsubvmns,             &
     &   xm, xn, xm_nyq, xn_nyq, ns, mnmax, mnmax_nyq,                         &
     &   phi, iotaf, bsupumnc, bsupvmnc, bsupumns, bsupvmns
     
      USE mgrid_mod, only : mgrid_mode
     
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_aux_1), INTENT (inout)       :: this
      TYPE (eq_state), INTENT (in)          :: state

!  Declare local variables
      INTEGER nwout, ierr
      INTEGER, DIMENSION(3)   :: dims
! 11-29-04      INTEGER :: ns, ntor, mpol1, ntmax
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_aux_1_construct: '

!  Start of executable code

!  Read the wout file
      ierr = 0
      CALL read_wout_file(state % wout_filename, ierr)
      CALL assert_eq(0,ierr,sub_name // 'read_wout error')

!  Clear out this
      CALL eq_aux_1_destroy(this)

!  Non-array variables
      this % mgrid_mode = mgrid_mode
      this % ns         = ns
      this % mnmax      = mnmax
      this % mnmax_nyq  = mnmax_nyq
      
!  Copy arrays in to this
      dims(1:2) = shape(rmnc)
      ALLOCATE(this % rmnc(dims(1),dims(2)))
      this % rmnc = rmnc
      dims(1:2) = shape(zmns)
      ALLOCATE(this % zmns(dims(1),dims(2)))
      this % zmns = zmns
      dims(1:2) = shape(lmns)
      ALLOCATE(this % lmns(dims(1),dims(2)))
      this % lmns = lmns
      dims(1:2) = shape(gmnc)
      ALLOCATE(this % gmnc(dims(1),dims(2)))
      this % gmnc = gmnc
      dims(1:2) = shape(currumnc)
      ALLOCATE(this % currumnc(dims(1),dims(2)))
      this % currumnc = currumnc
      dims(1:2) = shape(currvmnc)
      ALLOCATE(this % currvmnc(dims(1),dims(2)))
      this % currvmnc = currvmnc
      dims(1:2) = shape(bsubumnc)
      ALLOCATE(this % bsubumnc(dims(1),dims(2)))
      this % bsubumnc = bsubumnc
      dims(1:2) = shape(bsubvmnc)
      ALLOCATE(this % bsubvmnc(dims(1),dims(2)))
      this % bsubvmnc = bsubvmnc

      dims(1:2) = shape(bsupumnc)
      ALLOCATE(this % bsupumnc(dims(1),dims(2)))
      this % bsupumnc = bsupumnc
      dims(1:2) = shape(bsupvmnc)
      ALLOCATE(this % bsupvmnc(dims(1),dims(2)))
      this % bsupvmnc = bsupvmnc

      dims(1:1) = shape(xm)
      ALLOCATE(this % xm(dims(1)))
      this % xm = xm
      dims(1:1) = shape(xn)
      ALLOCATE(this % xn(dims(1)))
      this % xn = xn
      dims(1:1) = shape(xm_nyq)
      ALLOCATE(this % xm_nyq(dims(1)))
      this % xm_nyq = xm_nyq
      dims(1:1) = shape(xn_nyq)
      ALLOCATE(this % xn_nyq(dims(1)))
      this % xn_nyq = xn_nyq

      dims(1:1) = shape(phi)
      ALLOCATE(this % phi(dims(1)))
      this % phi = phi
      dims(1:1) = shape(iotaf)
      ALLOCATE(this % iotaf(dims(1)))
      this % iotaf = iotaf

      IF (lasym) THEN
         dims(1:2) = shape(rmns)
         ALLOCATE(this % rmns(dims(1),dims(2)))
         this % rmns = rmns
         dims(1:2) = shape(zmnc)
         ALLOCATE(this % zmnc(dims(1),dims(2)))
         this % zmnc = zmnc
         dims(1:2) = shape(lmnc)
         ALLOCATE(this % lmnc(dims(1),dims(2)))
         this % lmnc = lmnc
         dims(1:2) = shape(gmns)
         ALLOCATE(this % gmns(dims(1),dims(2)))
         this % gmns = gmns
         dims(1:2) = shape(currumns)
         ALLOCATE(this % currumns(dims(1),dims(2)))
         this % currumns = currumns
         dims(1:2) = shape(currvmns)
         ALLOCATE(this % currvmns(dims(1),dims(2)))
         this % currvmns = currvmns
         dims(1:2) = shape(bsubumns)
         ALLOCATE(this % bsubumns(dims(1),dims(2)))
         this % bsubumns = bsubumns
         dims(1:2) = shape(bsubvmns)
         ALLOCATE(this % bsubvmns(dims(1),dims(2)))
         this % bsubvmns = bsubvmns

         dims(1:2) = shape(bsupumns)
         ALLOCATE(this % bsupumns(dims(1),dims(2)))
         this % bsupumns = bsupumns
         dims(1:2) = shape(bsupvmns)
         ALLOCATE(this % bsupvmns(dims(1),dims(2)))
         this % bsupvmns = bsupvmns
      END IF

!  Check variables
      CALL assert((lasym .EQV. state % fixp % lasym),sub_name //               &
     &  'lasym discrepancy')
      CALL assert_eq(ns,size(rmnc,2),                                          &
     &   state % fixp % ns_array(state % varp % ns_index),                     &
     &   sub_name // 'ns discrepancy')
      CALL assert_eq(mnmax,size(xm,1),sub_name // 'mnmax discrepancy')
      CALL assert_eq(mnmax_nyq,size(xm_nyq,1),sub_name //                      &
     &   'mnmax_nyq discrepancy')
      
      RETURN
      
      END SUBROUTINE eq_aux_1_construct
!-------------------------------------------------------------------------------
!  Construct an eq_aux_2
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_2_construct(this,fixp,aux1,kv)
!
!  11-09-2004
!  Note that kv information comes from the magnetic diagnostic mrf files. It
!  is communicated as an explicit argument because it is really not part of
!  the state information. It is needed to set up the v part of the suv grids.
!  
!  The coding right now ASSUMES that the grid in the toroidal direction is 
!  THE SAME for both the magnetic diagnostic response function, and for the
!  VMEC variables on the s-u-v grid. 
            
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_aux_2), INTENT (inout)       :: this
      TYPE (eq_param_fix), INTENT (in)      :: fixp
      TYPE (eq_aux_1), INTENT (in)          :: aux1
      INTEGER, INTENT(in)            :: kv

!  Declare local variables
      INTEGER :: ns, ju, n_field_periods
      INTEGER :: mnmax, mnmax_nyq
      LOGICAL :: lasym
      REAL(rprec)    :: fperiod
      REAL(rprec)    :: dels, delu, delv
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cosz, sinz
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosu, sinu,                  & 
     &                                            cosv, sinv
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sgrid, ugrid, vgrid

      INTEGER :: i, j, k, js             
      INTEGER, DIMENSION(3)   :: dims           

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_aux_2_construct: '

!  Start of executable code

!  Destroy this
      CALL eq_aux_2_destroy(this)
      
!----------------------------------------------------------------------
!-- get plasma contribution                                          --
!----------------------------------------------------------------------
!  local variable definitions
!  Consistency checks done in aux_1_construct
      lasym = fixp % lasym
      ns = size(aux1 % rmnc,2)
      ju = ns
      n_field_periods = fixp % nfp
      mnmax = size(aux1 % xm,1)
      mnmax_nyq = size(aux1 % xm_nyq,1)

      fperiod = twopi / n_field_periods
      dels = one / (ns - 1)
      delu = twopi / ju
      delv = fperiod / kv
!----------------------------------------------------------------------
!-- allocate arrays                                                  --
!----------------------------------------------------------------------
!  s, u, v, VMEC flux grids
      ALLOCATE (sgrid(ns),ugrid(ju),vgrid(kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc sgrid')

!  non-nyq variables
      ALLOCATE (this % rsuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc rsuv')

      ALLOCATE (this % zsuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc zsuv')

      ALLOCATE (this % gsuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc gsuv')

      ALLOCATE (this % currusuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc currusuv')

      ALLOCATE (this % currvsuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc currvsuv')

      ALLOCATE (this % rusuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc rusuv')

      ALLOCATE (this % zusuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc zusuv')

      ALLOCATE (this % rvsuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc rvsuv')

      ALLOCATE (this % zvsuv(ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc zvsuv')

!  Sin and cosine arrays
      ALLOCATE (cosz(mnmax), sinz(mnmax),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc cosz')

      ALLOCATE (cosu(mnmax, ju), sinu(mnmax, ju),
     &          cosv(mnmax, kv), sinv(mnmax, kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc cosu')

!  Surface formulation arrays
      ALLOCATE (this % bsubu(ns-1:ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc bsubu')

      ALLOCATE (this % bsubv(ns-1:ns,ju,kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc bsubv')

!----------------------------------------------------------------------
!-- set up (s,u,v) VMEC flux grids                                   --
!----------------------------------------------------------------------
      DO i = 1, ns
         sgrid(i) = (i - 1) * dels
      END DO
      DO i = 1, ju
         ugrid(i) = (i - 1) * delu
         cosu(:,i) = COS(aux1 % xm * ugrid(i))
         sinu(:,i) = SIN(aux1 % xm * ugrid(i))
      END DO
      DO i = 1, kv
         vgrid(i) = (i - 1) * delv
         cosv(:,i) = COS(aux1 % xn * vgrid(i))
         sinv(:,i) = SIN(aux1 % xn * vgrid(i))
      END DO
!
!       DO FOURIER MODE SUM IN INNERMOST LOOP
!       FIRST DO mnmax-sized ARRAYS (r, z)
!
      DO j = 1, ju
      DO k = 1, kv
!          zetamn = xm*ugrid(j) - xn*vgrid(k)
!          cosz = COS(zetamn)
!          sinz = SIN(zetamn)
         cosz = cosu(:,j) * cosv(:,k) + sinu(:,j) * sinv(:,k)
         sinz = sinu(:,j) * cosv(:,k) - cosu(:,j) * sinv(:,k)
         DO js = 1, ns
            this % rsuv(js,j,k) = SUM(aux1 % rmnc(:,js) * cosz)
            this % zsuv(js,j,k) = SUM(aux1 % zmns(:,js) * sinz)
            this % rusuv(js,j,k) =                                             &
     &         - SUM(aux1 % xm * aux1 % rmnc(:,js) * sinz)
            this % zusuv(js,j,k) =                                             &
     &         SUM(aux1 % xm * aux1 % zmns(:,js) * cosz)
            this % rvsuv(js,j,k) =                                             & 
     &         SUM(aux1 % xn * aux1 % rmnc(:,js) * sinz)
            this % zvsuv(js,j,k) =                                             &
     &         -SUM(aux1 % xn * aux1 % zmns(:,js) * cosz)
         END DO

!----------------------------------------------------------------------
!-- stellarator asymmetric terms                                     --
!----------------------------------------------------------------------
         IF (lasym) THEN
            DO js = 1, ns
               this % rsuv(js,j,k) = this % rsuv(js,j,k) +                     &
     &            SUM(aux1 % rmns(:,js)*sinz)
               this % zsuv(js,j,k) = this % zsuv(js,j,k) +                     &
     &            SUM(aux1 % zmnc(:,js)*cosz)
               this % rusuv(js,j,k) = this % rusuv(js,j,k) +                   &
     &            SUM(aux1 % xm * aux1 % rmns(:,js) * cosz)
               this % zusuv(js,j,k) = this % zusuv(js,j,k) -                   &
     &            SUM(aux1 % xm * aux1 % zmnc(:,js) * sinz)
               this % rvsuv(js,j,k) = this % rvsuv(js,j,k) -                   &
     &            SUM(aux1 % xn * aux1 % rmns(:,js) * cosz)
               this % zvsuv(js,j,k) = this % zvsuv(js,j,k) +                   & 
     &            SUM(aux1 % xn * aux1 % zmnc(:,js) * sinz)
            END DO
         ENDIF
      END DO
      END DO

!----------------------------------------------------------------------
!       NEXT DO mnmax_nyq-sized ARRAYS (g, curru,v, bsubu,v)
!----------------------------------------------------------------------
!  Deallocate and allocate space
      DEALLOCATE (cosz, sinz, cosu, cosv, sinu, sinv)

      ALLOCATE (cosz(mnmax_nyq), sinz(mnmax_nyq),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc cosz 2')

      ALLOCATE (cosu(mnmax_nyq, ju), sinu(mnmax_nyq, ju),
     &          cosv(mnmax_nyq, kv), sinv(mnmax_nyq, kv),stat=i)
      CALL assert_eq(0,i,sub_name // 'alloc cosu 2')

!  Set up cosine, sine terms
      DO i = 1, ju
         cosu(:,i) = COS(aux1 % xm_nyq(:) * ugrid(i))
         sinu(:,i) = SIN(aux1 % xm_nyq(:) * ugrid(i))
      END DO

      DO i = 1, kv
         vgrid(i) = (i - 1) * delv
         cosv(:,i) = COS(aux1 % xn_nyq(:) * vgrid(i))
         sinv(:,i) = SIN(aux1 % xn_nyq(:) * vgrid(i))
      END DO

!  Do FT summation
      DO j = 1, ju
      DO k = 1, kv
!          zetamn = xm_nyq*ugrid(j) - xn_nyq*vgrid(k)
!          cosz = COS(zetamn)
!          sinz = SIN(zetamn)
         cosz = cosu(:,j)*cosv(:,k) + sinu(:,j)*sinv(:,k)
         DO js = 1, ns
            this % gsuv(js,j,k) = SUM(aux1 % gmnc(:,js)*cosz)
            this % currusuv(js,j,k) = SUM(aux1 % currumnc(:,js) * cosz)
            this % currvsuv(js,j,k) = SUM(aux1 % currvmnc(:,js) * cosz)
         END DO

!  Surface integral terms
!         IF (lsurf) THEN
            DO js = ns-1, ns
               this % bsubu(js,j,k) = SUM(aux1 % bsubumnc(:,js) * cosz)
               this % bsubv(js,j,k) = SUM(aux1 % bsubvmnc(:,js) * cosz)
            END DO
!         END IF

!----------------------------------------------------------------------
!-- stellarator asymmetric terms                                     --
!----------------------------------------------------------------------
         IF (lasym) THEN
            sinz = sinu(:,j)*cosv(:,k) - cosu(:,j)*sinv(:,k)
            DO js = 1, ns
              this % gsuv(js,j,k) = this % gsuv(js,j,k) +                      & 
     &           SUM(aux1 % gmns(:,js) * sinz)
              this % currusuv(js,j,k) = this % currusuv(js,j,k) +              & 
     &           SUM(aux1 % currumns(:,js) * sinz)
              this % currvsuv(js,j,k) = this % currvsuv(js,j,k) +              & 
     &           SUM(aux1 % currvmns(:,js) * sinz)
            END DO
!  Surface integral terms
!            IF (lsurf) THEN
               DO js = ns-1, ns
                  this % bsubu(js,j,k) = this % bsubu(js,j,k) +                & 
     &               SUM(aux1 % bsubumns(:,js) * sinz)
                  this % bsubv(js,j,k) = this % bsubv(js,j,k) +                &
     &               SUM(aux1 % bsubvmns(:,js) * sinz)
               END DO
!            END IF
         ENDIF
      END DO
      END DO

      DEALLOCATE (cosz, sinz, cosu, cosv, sinu, sinv)

!----------------------------------------------------------------------
!-- Other variables                                                  --
!----------------------------------------------------------------------
      this % phitot = aux1 % phi(SIZE(aux1 % phi))
      this % mnyq = INT(MAXVAL(aux1 % xm_nyq))
      this % nnyq = INT(MAXVAL(ABS(aux1 % xn_nyq)))/ fixp % nfp
 
      RETURN
      
      END SUBROUTINE eq_aux_2_construct

!-------------------------------------------------------------------------------
!  Construct an eq_state
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_state_construct(this,code,version,s_id,l_id,               &
     &   fixp,varp,xc,xcdot,wout_filename,fsq_max)

      IMPLICIT NONE
!  Declare Arguments 
      TYPE (eq_state), INTENT (inout) :: this
      CHARACTER (len=*), INTENT(in)   :: code
      CHARACTER (len=*), INTENT(in)   :: version
      CHARACTER (len=*), INTENT(in)   :: s_id
      CHARACTER (len=*), INTENT(in)   :: l_id
      TYPE (eq_param_fix), INTENT(in) :: fixp
      TYPE (eq_param_var), INTENT(in) :: varp
!  Arguments specific to VMEC equilibrium
      REAL(rprec), DIMENSION(:) :: xc
      REAL(rprec), DIMENSION(:) :: xcdot
      CHARACTER(len=*) :: wout_filename
      REAL(rprec) :: fsq_max

!  Declare local variables
      INTEGER :: nd1, nd2, ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_state_construct: '

!  Start of executable code

!  Scalar assignments, common to all codes
      this % code = code
      this % version = version
      this % s_id = s_id
      this % l_id = l_id
      this % fixp = fixp
      this % varp = varp

!  Different coding, depending on equilibrium code
!  NB: JDH 10-21-04.There is a possible memory leakage here !!!
!        This should be on the list of issues !!!

      SELECT CASE (code) ! Different coding depending on code
      CASE ('VMEC2000')
         IF (ASSOCIATED(this % xc)) THEN
            IF (SIZE(this % xc) .ne. SIZE(xc)) THEN
               DEALLOCATE(this % xc, this % xcdot,STAT=ier1)
               CALL assert_eq(0,ier1,sub_name // 'dealloc xc')
               ALLOCATE(this % xc(1:SIZE(xc)),this % xcdot(1:SIZE(xc)),        &
     &            STAT=ier1)
               CALL assert_eq(0,ier1,sub_name // 'alloc xc 1')
            ENDIF
         ELSE
            ALLOCATE(this % xc(1:SIZE(xc)),this % xcdot(1:SIZE(xc)),           &
     &         STAT=ier1)
            CALL assert_eq(0,ier1,sub_name // 'alloc xc 2')
         ENDIF
         
         this % l_def_aux1 = .false.
         this % l_def_aux2 = .false.
         this % xc = xc
         this % xcdot = xcdot
         this % fsq_max = fsq_max
         this % wout_filename = wout_filename
            
      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized code: ',                     &
     &      char=code)
      END SELECT ! Different coding depending on code
      
      END SUBROUTINE eq_state_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy an eq_param_fix
!-------------------------------------------------------------------------------
      SUBROUTINE eq_param_fix_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_param_fix), INTENT(inout) :: this

!  Start of executable code
!  Get rid of all components
!  Variables for code='VMEC2000'
      this % lasym   = .false.
      this % lfreeb  = .false.
      this % lforbal = .false.
      this % lRFP = .false.
      this % mgrid_file = ' '
      this % nfp     = 0
      this % ncurr   = 0
      this % nvacskip = 0
      this % mpol    = 0
      this % ntor    = 0
      this % ntheta  = 0
      this % nzeta   = 0
      this % mns   = 0
      this % ns_array = 0
      this % nextcur = 0
      this % pcurr_type = ' '
      this % piota_type = ' '
      this % pmass_type = ' '
      this % niter   = 0
      this % nstep   = 0
      this % delt    = zero
      this % ftol_array = zero
      this % tcon0   = zero
      this % gamma   = zero

      END SUBROUTINE eq_param_fix_destroy

!-------------------------------------------------------------------------------
!  Destroy an eq_param_var
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_param_var_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_param_var), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_param_var_destroy: '

!  Start of executable code
!  Get rid of all components
!  Variables for code='VMEC2000'
      
      this % ns_index = 0

      this % ac = zero
      this % ac_aux_s = zero
      this % ac_aux_f = zero
      this % ai = zero
      this % ai_aux_s = zero
      this % ai_aux_f = zero
      this % am = zero
      this % am_aux_s = zero
      this % am_aux_f = zero

      this % bloat = 0
      
      IF (ASSOCIATED(this % rbc)) THEN
         DEALLOCATE(this % rbc,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc rbc')
      ENDIF

      IF (ASSOCIATED(this % zbs)) THEN
         DEALLOCATE(this % zbs,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc zbs')
      ENDIF

      IF (ASSOCIATED(this % zbc))  THEN
         DEALLOCATE(this % zbc,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc zbc')
      ENDIF

      IF (ASSOCIATED(this % rbs))  THEN
         DEALLOCATE(this % rbs,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc rbs')
      ENDIF

      this % extcur = zero         
      this % curtor = zero
      this % phiedge = zero
      this % pres_scale = zero

      END SUBROUTINE eq_param_var_destroy
!-------------------------------------------------------------------------------
!  Destroy an eq_aux_1
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_1_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_aux_1), INTENT (inout)       :: this

!  Start of executable code
!  Variables for code='VMEC'
      this % mgrid_mode = ' '
      this % ns         = 0
      this % mnmax      = 0
      this % mnmax_nyq  = 0

      IF (ASSOCIATED(this % rmnc)) THEN
         DEALLOCATE(this % rmnc, this % zmns, this % lmns,                     &
     &      this % gmnc, this % currumnc, this % currvmnc,                     &
     &      this % bsubumnc, this % bsubvmnc,                                  &
     &      this % bsupumnc, this % bsupvmnc)
      END IF
      IF (ASSOCIATED(this % rmns)) THEN
         DEALLOCATE(this % rmns, this % zmnc, this % lmnc,                     &
     &      this % gmns, this % currumns, this % currvmns,                     &
     &      this % bsubumns, this % bsubvmns,                                  &
     &      this % bsupumns, this % bsupvmns)
      END IF
      IF (ASSOCIATED(this % xm)) THEN
         DEALLOCATE(this % xm, this % xn, this % xm_nyq,                       &
     &      this % xn_nyq)
      END IF
      IF (ASSOCIATED(this % phi)) THEN
         DEALLOCATE(this % phi, this % iotaf)
      END IF

      END SUBROUTINE eq_aux_1_destroy

!-------------------------------------------------------------------------------
!  Destroy an eq_aux_2_destroy
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_2_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_aux_2), INTENT (inout)       :: this

!  Start of executable code

      IF (ASSOCIATED(this % rsuv)) THEN
         DEALLOCATE(this % rsuv, this % zsuv,                                  &
     &      this % gsuv, this % currusuv, this % currvsuv,                     &
     &      this % rusuv, this % zusuv, this % rvsuv, this % zvsuv,            &
     &      this % bsubu, this % bsubv)
      END IF
      this % phitot = zero
      this % mnyq   = 0
      this % nnyq   = 0
 
      END SUBROUTINE eq_aux_2_destroy

!-------------------------------------------------------------------------------
!  Destroy an eq_state
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_state_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_state), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_state_destroy: '

!  Start of executable code
!  Get rid of all components
      this % code = ' '
      this % version = ' '
      this % s_id = ' '
      this % l_id = ' '
!  Variables for code='VMEC2000'
!  NB: JDH 10-21-04. This deallocation could cause real problems if done incorrectly !!!
!        This should be on the list of issues !!!

      CALL eq_param_fix_destroy(this % fixp)
      CALL eq_param_var_destroy(this % varp)
      this % l_def_aux1 = .false.
      CALL eq_aux_1_destroy(this % aux1)
      this % l_def_aux2 = .false.
      CALL eq_aux_2_destroy(this % aux2)

      IF (ASSOCIATED(this % xc)) DEALLOCATE(this % xc,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'dealloc xc')
      IF (ASSOCIATED(this % xcdot)) DEALLOCATE(this % xcdot,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'dealloc xcdot')
      this % fsq_max = zero
      this % wout_filename = ' '

      END SUBROUTINE eq_state_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for eq_param_fix - Default OK as no Pointers
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Assignment for eq_param_var
!-------------------------------------------------------------------------------
      SUBROUTINE eq_param_var_assign(left,right)
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (eq_param_var), INTENT (out) :: left
      TYPE (eq_param_var), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_param_var_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right pointers are the same?. FIX IT'
      INTEGER :: ier1, ier2, ier3
      LOGICAL, DIMENSION(4) :: lassert
      INTEGER :: n, nm1
         
!  Start of executable code

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      lassert(1) = .not.ASSOCIATED(left % rbc,right % rbc)
      lassert(2) = .not.ASSOCIATED(left % zbs,right % zbs)
      lassert(3) = .not.ASSOCIATED(left % zbc,right % zbc)
      lassert(4) = .not.ASSOCIATED(left % rbs,right % rbs)
      CALL assert(lassert,sub_name // err_mess1)
      
!  Destroy left
      CALL eq_destroy(left)
      
!  Non-pointer variables.
!  Copy them over
      left % ns_index = right % ns_index
      left % ac = right % ac
      left % ac_aux_s = right % ac_aux_s
      left % ac_aux_f = right % ac_aux_f
      left % ai = right % ai
      left % ai_aux_s = right % ai_aux_s
      left % ai_aux_f = right % ai_aux_f
      left % am = right % am
      left % am_aux_s = right % am_aux_s
      left % am_aux_f = right % am_aux_f
      left % extcur = right % extcur
      left % bloat = right % bloat
      left % curtor = right % curtor
      left % phiedge = right % phiedge
      left % pres_scale = right % pres_scale

!  Allocate space for arrays (Were deallocated in _destroy)
      ALLOCATE(left % rbc(SIZE(right % rbc,1),                                 &
     &   SIZE(right % rbc,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rbc')
      left % rbc = right % rbc

      ALLOCATE(left % zbs(SIZE(right % zbs,1),                                 &
     &   SIZE(right % zbs,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc zbs')
      left % zbs = right % zbs

      ALLOCATE(left % zbc(SIZE(right % zbc,1),                                 &
     &   SIZE(right % zbc,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc zbc')
      left % zbc = right % zbc

      ALLOCATE(left % rbs(SIZE(right % rbs,1),                                 &
     &   SIZE(right % rbs,2)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rbs')
      left % rbs = right % rbs
         
      END SUBROUTINE eq_param_var_assign

!-------------------------------------------------------------------------------
!  Assignment for eq_aux_1
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_1_assign(left,right)
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (eq_aux_1), INTENT (inout) :: left
      TYPE (eq_aux_1), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_aux_1_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right pointers are the same?. FIX IT'
      INTEGER :: ier1, ier2, ier3
      LOGICAL :: lassert
      INTEGER :: n, nm1
      INTEGER :: dims(3)
         
!  Start of executable code

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      lassert = .not.ASSOCIATED(left % rmnc,right % rmnc)
      CALL assert(lassert,sub_name // err_mess1)

      CALL eq_aux_1_destroy(left)

!  Non-pointer variables.
!  Copy them over
      left % mgrid_mode = right % mgrid_mode
      left % mnmax = right % mnmax
      left % mnmax_nyq = right % mnmax_nyq
      left % ns = right % ns

!  Check that symmetric arrays are allocated
      IF (ASSOCIATED(right % rmnc)) THEN
      
         dims(1:2) = SHAPE(right % rmnc)
         ALLOCATE(left % rmnc(dims(1),dims(2)))
         left % rmnc = right % rmnc

         dims(1:2) = SHAPE(right % zmns)
         ALLOCATE(left % zmns(dims(1),dims(2)))
         left % zmns = right % zmns

         dims(1:2) = SHAPE(right % lmns)
         ALLOCATE(left % lmns(dims(1),dims(2)))
         left % lmns = right % lmns

         dims(1:2) = SHAPE(right % gmnc)
         ALLOCATE(left % gmnc(dims(1),dims(2)))
         left % gmnc = right % gmnc

         dims(1:2) = SHAPE(right % currumnc)
         ALLOCATE(left % currumnc(dims(1),dims(2)))
         left % currumnc = right % currumnc

         dims(1:2) = SHAPE(right % currvmnc)
         ALLOCATE(left % currvmnc(dims(1),dims(2)))
         left % currvmnc = right % currvmnc

         dims(1:2) = SHAPE(right % bsubumnc)
         ALLOCATE(left % bsubumnc(dims(1),dims(2)))
         left % bsubumnc = right % bsubumnc

         dims(1:2) = SHAPE(right % bsubvmnc)
         ALLOCATE(left % bsubvmnc(dims(1),dims(2)))
         left % bsubvmnc = right % bsubvmnc

         dims(1:2) = SHAPE(right % bsupumnc)
         ALLOCATE(left % bsupumnc(dims(1),dims(2)))
         left % bsupumnc = right % bsupumnc

         dims(1:2) = SHAPE(right % bsupvmnc)
         ALLOCATE(left % bsupvmnc(dims(1),dims(2)))
         left % bsupvmnc = right % bsupvmnc

         dims(1:1) = SHAPE(right % xm)
         ALLOCATE(left % xm(dims(1)))
         left % xm = right % xm

         dims(1:1) = SHAPE(right % xn)
         ALLOCATE(left % xn(dims(1)))
         left % xn = right % xn

         dims(1:1) = SHAPE(right % xm_nyq)
         ALLOCATE(left % xm_nyq(dims(1)))
         left % xm_nyq = right % xm_nyq

         dims(1:1) = SHAPE(right % xn_nyq)
         ALLOCATE(left % xn_nyq(dims(1)))
         left % xn_nyq = right % xn_nyq

         dims(1:1) = SHAPE(right % phi)
         ALLOCATE(left % phi(dims(1)))
         left % phi = right % phi

         dims(1:1) = SHAPE(right % iotaf)
         ALLOCATE(left % iotaf(dims(1)))
         left % iotaf = right % iotaf
      ELSE
         left % rmnc => null()
         left % zmns => null()
         left % lmns => null()
         left % gmnc => null()
         left % currumnc => null()
         left % currvmnc => null()
         left % bsubumnc => null()
         left % bsubvmnc => null()
         left % bsupumnc => null()
         left % bsupvmnc => null()
         left % xm => null()
         left % xn => null()
         left % xm_nyq => null()
         left % xn_nyq => null()
         left % phi => null()
         left % iotaf => null()
      ENDIF

!  Check that asymmetric arrays are associated
      IF (ASSOCIATED(right % rmns)) THEN

         dims(1:2) = SHAPE(right % rmns)
         ALLOCATE(left % rmns(dims(1),dims(2)))
         left % rmns = right % rmns

         dims(1:2) = SHAPE(right % zmnc)
         ALLOCATE(left % zmnc(dims(1),dims(2)))
         left % zmnc = right % zmnc

         dims(1:2) = SHAPE(right % lmnc)
         ALLOCATE(left % lmnc(dims(1),dims(2)))
         left % lmnc = right % lmnc

         dims(1:2) = SHAPE(right % gmns)
         ALLOCATE(left % gmns(dims(1),dims(2)))
         left % gmns = right % gmns

         dims(1:2) = SHAPE(right % currumns)
         ALLOCATE(left % currumns(dims(1),dims(2)))
         left % currumns = right % currumns

         dims(1:2) = SHAPE(right % currvmns)
         ALLOCATE(left % currvmns(dims(1),dims(2)))
         left % currvmns = right % currvmns

         dims(1:2) = SHAPE(right % bsubumns)
         ALLOCATE(left % bsubumns(dims(1),dims(2)))
         left % bsubumns = right % bsubumns

         dims(1:2) = SHAPE(right % bsubvmns)
         ALLOCATE(left % bsubvmns(dims(1),dims(2)))
         left % bsubvmns = right % bsubvmns

         dims(1:2) = SHAPE(right % bsupumns)
         ALLOCATE(left % bsupumns(dims(1),dims(2)))
         left % bsupumns = right % bsupumns

         dims(1:2) = SHAPE(right % bsupvmns)
         ALLOCATE(left % bsupvmns(dims(1),dims(2)))
         left % bsupvmns = right % bsupvmns
      ELSE
         left % rmns => null()
         left % zmnc => null()
         left % lmnc => null()
         left % gmns => null()
         left % currumns => null()
         left % currvmns => null()
         left % bsubumns => null()
         left % bsubvmns => null()
         left % bsupumns => null()
         left % bsupvmns => null()
      END IF
         
      END SUBROUTINE eq_aux_1_assign
!-------------------------------------------------------------------------------
!  Assignment for eq_aux_2
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_2_assign(left,right)
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (eq_aux_2), INTENT (inout) :: left
      TYPE (eq_aux_2), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_aux_2_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right pointers are the same?. FIX IT'
      INTEGER :: ier1, ier2, ier3
      LOGICAL :: lassert
      INTEGER :: n, nm1
      INTEGER :: dims(3)
         
!  Start of executable code

!  Non-pointer variables.
!  Copy them over
      left % phitot = right % phitot
      left % mnyq = right % mnyq
      left % nnyq = right % nnyq

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      lassert = .not.ASSOCIATED(left % bsubu,right % bsubu)
      CALL assert(lassert,sub_name // err_mess1)
      
!  Check that arrays are allocated
      IF (ASSOCIATED(right % bsubu)) THEN
      
         dims(1:3) = SHAPE(right % bsubu)
         ALLOCATE(left % bsubu(dims(1),dims(2),dims(3)))
         left % bsubu = right % bsubu
      
         dims(1:3) = SHAPE(right % bsubv)
         ALLOCATE(left % bsubv(dims(1),dims(2),dims(3)))
         left % bsubv = right % bsubv
      
         dims(1:3) = SHAPE(right % rsuv)
         ALLOCATE(left % rsuv(dims(1),dims(2),dims(3)))
         left % rsuv = right % rsuv
      
         dims(1:3) = SHAPE(right % zsuv)
         ALLOCATE(left % zsuv(dims(1),dims(2),dims(3)))
         left % zsuv = right % zsuv
      
         dims(1:3) = SHAPE(right % gsuv)
         ALLOCATE(left % gsuv(dims(1),dims(2),dims(3)))
         left % gsuv = right % gsuv
      
         dims(1:3) = SHAPE(right % currusuv)
         ALLOCATE(left % currusuv(dims(1),dims(2),dims(3)))
         left % currusuv = right % currusuv
      
         dims(1:3) = SHAPE(right % currvsuv)
         ALLOCATE(left % currvsuv(dims(1),dims(2),dims(3)))
         left % currvsuv = right % currvsuv
      
         dims(1:3) = SHAPE(right % rusuv)
         ALLOCATE(left % rusuv(dims(1),dims(2),dims(3)))
         left % rusuv = right % rusuv
      
         dims(1:3) = SHAPE(right % zusuv)
         ALLOCATE(left % zusuv(dims(1),dims(2),dims(3)))
         left % zusuv = right % zusuv
      
         dims(1:3) = SHAPE(right % rvsuv)
         ALLOCATE(left % rvsuv(dims(1),dims(2),dims(3)))
         left % rvsuv = right % rvsuv
      
         dims(1:3) = SHAPE(right % zvsuv)
         ALLOCATE(left % zvsuv(dims(1),dims(2),dims(3)))
         left % zvsuv = right % zvsuv

      ELSE
         left % bsubu => null()
         left % bsubv => null()
         left % rsuv => null()
         left % zsuv => null()
         left % gsuv => null()
         left % currusuv => null()
         left % currvsuv => null()
         left % rusuv => null()
         left % zusuv => null()
         left % rvsuv => null()
         left % zvsuv => null()
      ENDIF
         
      END SUBROUTINE eq_aux_2_assign
!-------------------------------------------------------------------------------
!  Assignment for eq_state
!-------------------------------------------------------------------------------
      SUBROUTINE eq_state_assign(left,right)
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (eq_state), INTENT (out) :: left
      TYPE (eq_state), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_state_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right pointers are the same?. FIX IT'
      INTEGER :: ier1, ier2, ier3
      LOGICAL :: lassert
      INTEGER :: n, nm1
      INTEGER :: dims(3)
         
!  Start of executable code

!  Non-pointer variables.
!  Copy them over
      left % code = right % code
      left % version = right % version
      left % s_id = right % s_id
      left % l_id = right % l_id
      left % l_def_aux1 = right % l_def_aux1
      left % l_def_aux2 = right % l_def_aux2
      left % wout_filename = right % wout_filename
      left % fsq_max = right % fsq_max

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      lassert = .not.ASSOCIATED(left % xc,right % xc)
      CALL assert(lassert,sub_name // err_mess1)
      
!  Assignments for sub-structures
      left % fixp = right % fixp
      left % varp = right % varp
      IF (right % l_def_aux1) left % aux1 = right % aux1
      IF (right % l_def_aux2) left % aux2 = right % aux2

      IF (ASSOCIATED(right % xc)) THEN
      
         dims(1:1) = SHAPE(right % xc)
         ALLOCATE(left % xc(dims(1)))
         left % xc = right % xc
      
         dims(1:1) = SHAPE(right % xcdot)
         ALLOCATE(left % xcdot(dims(1)))
         left % xcdot = right % xcdot

      ELSE
         left % xc => null()
         left % xcdot => null()
      ENDIF
         
      END SUBROUTINE eq_state_assign
          
!*******************************************************************************
! SECTION VII.    AUX DEFINITION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  The construction subroutines for the auxiliary are written so that the 
!  aux structures are independent of the state structure. However, the only use
!  contemplated at this time is for the aux structures to be components of
!  eq_states. 
!
!  In this section are subroutines that DEFINE the aux structure within a state.
!  (I am not sure that DEFINE is the best word, but I want to distinguish it 
!  from CONSTRUCT)
!  JDH 11-28-2004
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Define an eq_aux_1
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_1_define(this)
     
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_state), INTENT (inout)          :: this

!  Start of executable code
      CALL eq_aux_1_construct(this % aux1,this)
      this % l_def_aux1 = .true.
      
      RETURN
      
      END SUBROUTINE eq_aux_1_define
!-------------------------------------------------------------------------------
!  Define an eq_aux_2
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_aux_2_define(this,kv)

!  11-09-2004
!  Note that kv information comes from the magnetic diagnostic mrf files. It
!  is communicated as an explicit argument because it is really not part of
!  the state information.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_state), INTENT (inout)          :: this
      INTEGER, INTENT(in)               :: kv

!  Start of executable code
      IF (.not. this % l_def_aux1) THEN
         CALL eq_aux_1_define(this)
      END IF
      CALL eq_aux_2_construct(this % aux2,this % fixp,this % aux1,kv)
      this % l_def_aux2 = .true.
      
      RETURN
      
      END SUBROUTINE eq_aux_2_define

!*******************************************************************************
! SECTION VIII.   OTHER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Increase the ns_index
!    For control of the multi-grid capabilities of VMEC
!
!-------------------------------------------------------------------------------
      SUBROUTINE eq_ns_index_inc(this,lsuccess)

      IMPLICIT NONE

!  Declare Arguments
!  this       eq_state that will have the ns_index increased
!  lsuccess   logical, .T. if incremented correctly, 
!                      .F. if have now exceeded ns_index_max
      TYPE (eq_state), INTENT (inout)          :: this
      LOGICAL, INTENT(out)                     :: lsuccess

!  Declare local variables
      INTEGER :: ns_index_max, i
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_ns_index_inc: '

!  Start of executable code
!  First, find ns_index_max
      ns_index_max = SIZE(this % fixp % ns_array) 
      DO i = 1,SIZE(this % fixp % ns_array)
         IF (this % fixp % ns_array(i) .eq. 0) THEN
            ns_index_max = i - 1
            EXIT
         END IF
      END DO

      CALL assert(ns_index_max .lt. SIZE(this % fixp % ns_array),               &
     &   sub_name // 'No zeros in ns_array', 'Warn')

      CALL assert(this % varp % ns_index .le. ns_index_max,                    &    
     &   sub_name // 'ns_index too big')

      this % varp % ns_index = this % varp % ns_index + 1
      
      lsuccess = .TRUE.
      IF (this % varp % ns_index .gt. ns_index_max) THEN
         this % varp % ns_index = this % varp % ns_index - 1
         lsuccess = .FALSE.
      ENDIF

      RETURN
      
      END SUBROUTINE eq_ns_index_inc
      
!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 06.16.04
!     Start on module
!
!  JDH 09-21-2004
!     Steve Hirshman pared down the lists of variables in the param structures.
!     Fix up the construct and destroy subroutines. 
!     Added subroutine eq_state_point.
!
!  JDH 09-22-2004
!     Compiles.

!  JDH 10-27-2004
!    Added xcdot, fsq and wout_filename to _state. Now treating xc and xcdot as
!    allocatable arrays. Extra work in copying in and out, but cleaner logic, so
!    go with the cleaner logic, until timing becomes an issue.

!  JDH 11-04-2004
!    Changed fsq to fsq_max

!  JDH 11-26-2004
!    Changed so that eq_state has other derived types as components
!    (_aux_ pieces now in this module)

!  JDH 11-28-2004
!    Added nextcur to eq_param_fix. Added l_der aux1 and l_def_aux2 to eq_state
!    Added subroutines eq_aux_1_define and eq_aux_2_define
!    Cleaned up a little bit

!  JDH 11-30-2004
!    Moved variable consistency checks in to aux_1_construct.
!    Eliminated 'USE read_wout_mod ...' in aux_2_construct

!  JDH 12-06-2004
!    Moved ns_index from _fix to _var. Added mgrid_mode to _aux_1.
!    Added eq_ns_index_inc.

!  JDH 07-06-2005
!    Added lmns, lmnc, and phi arrays to eq_aux_1. Added some comments.
!    Added phitot to eq_aux_2.
!
!  JDH 07-13-2005. Corrected dimensions of am, ac, and ai to (0:10)
!    (In fixp)
!
!  JDH 07-21-2005
!    Added iotaf array to eq_aux_1.

!  JMS 08-03-2007
!    Added  ns, mnmax, mnmax_nyq  bsupumnc, bsupumns,  bsupvmnc, bsupvmns to eq_aux_1
!       and  mnyq, nnyq to eq_aux_2
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-08-04
!    Added subroutines eq_aux_1_assign, eq_aux_2_assign, and eq_state_assign
!
!  JDH 2010-12-06
!    Added variables to derived types:
!       eq_param_fix – pcurr_type, piota_type, pmass_type, lRFP
!       eq-param_var – bloat, ac_aux_s, ac_aux_f, ai_aux_s, ai_aux_f, 
!                      am_aux_s, am_aux_f
!    Did NOT update eq_cdf, which is out of date, and unused.
           
      END MODULE eq_T
