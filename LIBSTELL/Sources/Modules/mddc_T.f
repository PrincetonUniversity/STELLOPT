
!*******************************************************************************
!  File mddc_T.f
!  Contains module mddc_T
!  Defines derived-types: mddc_desc, mddc_mrf
!  A type of Diagnostic - Magnetic Diagnostic-Dot Coil
!   (A magnetic diagnostic coil, with original data from a diagnostic-dot file)
!
!*******************************************************************************
!  MODULE mddc_T
!    (MDDC Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE mddc_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!-------------------------------------------------------------------------------
      USE bsc_T
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: type_len=10      
      INTEGER, PARAMETER, PRIVATE :: sn_len=30      
      INTEGER, PARAMETER, PRIVATE :: ln_len=80      
      INTEGER, PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)  MDDC Description:
!       mddc_desc  
!     Type of diagnostic specified by  % d_type = 'mddc'.
!
! 2)   MDDC Magnetic Response Functions:
!          mddc_mrf
!
!  A mddc_mrf is a specialized structure, particular to magnetic diagnostics. The
!  model signal computation is an integration over the plasma volume, and
!  much of the integrand can be pre-computed, knowing only information about the
!  magnetic diagnostic. A mddc_mrf contains this pre-computed information, so 
!  that the model signal computation will be faster. The model signal
!  also contains a contribution due to the external field-coil groups. This 
!  coil-response information is also contained in a mddc_mrf.
!
!  Note that mddc_mrf is declared first, so that the type is known to a 
!  mddc_desc.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Declare type mddc_mrf
!         ---- Identification Variables ----
!    code_name        Character: name of the code which computed the responses
!    code_version     Character: version of the code
!    date_run         Character: date and time when the code was run
!    field_coils_id   Character: identifier of field-coils 
!
!         ----  Coil Response Function Variables----
!    n_field_cg         number of field-coil groups (external currents)
!                         (size of the mresponse_extcur array)
!    rdiag_coilg_1      array of diagnostic - field-coil-group responses
!                         (Single row of the rdiag_coilg array) 
!                         ^ USE AS ALLOCATABLE ARRAY!    
!    extcur_mg          array of external currents - 'MGRID'
!                         Used for normalization with 'raw' or 'scaled' mode 
!                         ^ USE AS ALLOCATABLE ARRAY!    
!
!         ---- Plasma Response Grid Variables ----
!    ir                 number of grid points in R, plasma grid
!    jz                 number of grid points in z, plasma grid
!    kp                 number of phi planes per field period in plasma grid
!    kp_store           number of phi planes actually store in plasma grid
!                          (With lstell_sym = .true., don't need to store all planes)
!    rmin               Minimum R for plasma grid
!    rmax               Maximum R for plasma grid
!    zmin               Minimum z for plasma grid
!    zmax               Maximum z for plasma grid
!    n_field_periods    Number of field periods
!    lstell_sym         Logical - True for stellarator symmetry
!
!          ---- Plasma Response Arrays ---- 
!    a_r                R component of plasma response function
!    a_f                phi component of plasma response function
!    a_z                Z component of plasma response function
!                         ^ USE AS ALLOCATABLE ARRAYS!    
!-------------------------------------------------------------------------------
      TYPE mddc_mrf
         CHARACTER(len=80)                  :: code_name
         CHARACTER(len=80)                  :: code_version
         CHARACTER(len=80)                  :: date_run
         CHARACTER(len=80)                  :: field_coils_id

         INTEGER                     :: n_field_cg
         REAL(rprec), DIMENSION(:), POINTER :: rdiag_coilg_1 => null()
         REAL(rprec), DIMENSION(:), POINTER :: extcur_mg => null()

         INTEGER                     :: ir
         INTEGER                     :: jz
         INTEGER                     :: kp    
         INTEGER                     :: kp_store    
         REAL(rprec)                        :: rmin
         REAL(rprec)                        :: rmax       
         REAL(rprec)                        :: zmin
         REAL(rprec)                        :: zmax
         INTEGER                     :: n_field_periods
         LOGICAL                            :: lstell_sym

         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_f => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_z => null()
      END TYPE mddc_mrf
!-------------------------------------------------------------------------------
!  Declare type mddc_desc
!       s_name          character, short name of diagnostic
!       l_name          character, long name of diagnostic
!       units           character, physical units that the data is measured in
!       sigma_default   real, default value of the uncertainty in the data
!       mddc_type       character, keyword from the diagnostic_dot file
!       l_mdcoil_def    logical, definition status of the mdcoil component     
!       flux_factor     real, factor to convert from flux to appropriate units
!       mdcoil          type bsc_coil, description of coil
!       mrf             type mddc_mrf, magnetic response functions
!!-------------------------------------------------------------------------------
      TYPE mddc_desc
         CHARACTER (len=sn_len)         :: s_name                                 
         CHARACTER (len=ln_len)         :: l_name
         CHARACTER (len=units_len)      :: units                                 
         CHARACTER (len=30)             :: mddc_type                                 
         LOGICAL                        :: l_mdcoil_def
         REAL(rprec)                    :: sigma_default
         REAL(rprec)                    :: flux_factor
         TYPE (bsc_coil)                :: mdcoil 
         TYPE (mddc_mrf)                :: mrf 
      END TYPE mddc_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE mddc_desc_assign,                                    &
     &                    mddc_mrf_assign  
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE mddc_construct
         MODULE PROCEDURE mddc_desc_construct,                                 &
     &                    mddc_mrf_construct
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE mddc_destroy
         MODULE PROCEDURE mddc_desc_destroy,                                   &
     &                    mddc_mrf_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE mddc_write
         MODULE PROCEDURE mddc_desc_write,                                     &
     &                    mddc_mrf_write
      END INTERFACE

!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a mddc_desc
!
!  For d_type = 'mddc' (magnetic mddc-dot coil)
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_desc_construct(this,s_name,l_name,units,                 &
     &   sigma_default,mddc_type,mdcoil,mrf,flux_factor)

!  NB. 
!  The mdcoil argument is assigned to the 'this' component. Do NOT call this
!  subroutine with this % mdcoil as the mdcoil argument.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (mddc_desc), INTENT(inout)            :: this
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      CHARACTER (len=*), INTENT(in)              :: mddc_type
      REAL(rprec), INTENT(in)                    :: sigma_default
      TYPE (bsc_coil), INTENT(in), TARGET        :: mdcoil  ! Why a TARGET ? 2007-06-11
      TYPE (mddc_mrf), INTENT(in), OPTIONAL      :: mrf
      REAL(rprec), INTENT(in), OPTIONAL          :: flux_factor

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_desc_construct: '

!  Start of executable code

!  Destroy the mdcoil component
      CALL bsc_destroy(this % mdcoil)

!  Destroy the mrf component
      CALL mddc_mrf_destroy(this % mrf)

!  Scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      this % mddc_type =  TRIM(ADJUSTL(mddc_type))
      this % l_mdcoil_def = .TRUE.
      IF (PRESENT(flux_factor)) THEN
         this % flux_factor = flux_factor
      ELSE
         this % flux_factor = one
      ENDIF

!  Derived Type Assignments
      this % mdcoil = mdcoil
      IF (PRESENT(mrf)) THEN
         this % mrf = mrf
      ENDIF
      
      END SUBROUTINE mddc_desc_construct

!-------------------------------------------------------------------------------
!  Construct a mddc_mrf
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_mrf_construct(this,code_name,code_version,               &
     &  date_run,field_coils_id,rdiag_coilg_1,extcur_mg,kp,                    &
     &  rmin,rmax,zmin,zmax,n_field_periods,lstell_sym,a_r,a_f,a_z)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (mddc_mrf), INTENT(inout)             :: this
      CHARACTER(len=*), INTENT(in)               :: code_name
      CHARACTER(len=*), INTENT(in)               :: code_version
      CHARACTER(len=*), INTENT(in)               :: date_run
      CHARACTER(len=*), INTENT(in)               :: field_coils_id
      REAL(rprec), DIMENSION(:), INTENT(in)      :: rdiag_coilg_1
      REAL(rprec), DIMENSION(:), INTENT(in)      :: extcur_mg
      INTEGER, INTENT(in)                 :: kp    
      REAL(rprec), INTENT(in)                    :: rmin
      REAL(rprec), INTENT(in)                    :: rmax       
      REAL(rprec), INTENT(in)                    :: zmin
      REAL(rprec), INTENT(in)                    :: zmax
      INTEGER, INTENT(in)                 :: n_field_periods
      LOGICAL, INTENT(in)                        :: lstell_sym
      REAL(rprec), DIMENSION(:,:,:), INTENT(in)  :: a_r
      REAL(rprec), DIMENSION(:,:,:), INTENT(in)  :: a_f
      REAL(rprec), DIMENSION(:,:,:), INTENT(in)  :: a_z

!  Declare local variables
      INTEGER           :: ir1, ir2, ir3, if1, if2, if3, iz1,           &
     &    iz2, iz3
      INTEGER           :: ier1, ier2, ier3
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_mrf_construct: '

!  Start of executable code

!  Destroy existing arrays
      CALL mddc_mrf_destroy(this)
            
!  Scalar variables
      this % code_name = ADJUSTL(code_name)
      this % code_version = ADJUSTL(code_version)
      this % date_run = ADJUSTL(date_run)
      this % field_coils_id = ADJUSTL(field_coils_id)
      this % kp = kp
      this % rmin = rmin
      this % rmax = rmax
      this % zmin = zmin
      this % zmax = zmax
      this % n_field_periods = n_field_periods
      this % lstell_sym = lstell_sym

!  Array Sizes
      this % n_field_cg = SIZE(rdiag_coilg_1)
      ir1 = SIZE(a_r,1)
      ir2 = SIZE(a_r,2)
      ir3 = SIZE(a_r,3)
      if1 = SIZE(a_f,1)
      if2 = SIZE(a_f,2)
      if3 = SIZE(a_f,3)
      iz1 = SIZE(a_z,1)
      iz2 = SIZE(a_z,2)
      iz3 = SIZE(a_z,3)
      CALL assert_eq(ir1,if1,iz1,sub_name // 'a_ first dims different')
      CALL assert_eq(ir2,if2,iz2,sub_name // 'a_ 2nd dims different')
      CALL assert_eq(ir3,if3,iz3,sub_name // 'a_ 3rd dims different')
      this % ir = ir1
      this % jz = ir2
      this % kp_store = ir3
      CALL assert_eq(this % n_field_cg,SIZE(extcur_mg),                        &
     &   sub_name // 'rd - extcur dims different')

!  Allocate space for arrays
      ALLOCATE(this % rdiag_coilg_1(this % n_field_cg),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rdiag_coilg_1')
      ALLOCATE(this % extcur_mg(this % n_field_cg),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc extcur_mg')

      ALLOCATE(this % a_r(ir1,ir2,ir3),STAT=ier1)
      ALLOCATE(this % a_f(ir1,ir2,ir3),STAT=ier2)
      ALLOCATE(this % a_z(ir1,ir2,ir3),STAT=ier3)
      CALL assert_eq(0,ier1,ier2,ier3,sub_name // 'alloc a_')

!  Move arrays
      this % rdiag_coilg_1 = rdiag_coilg_1
      this % extcur_mg = extcur_mg
      this % a_r = a_r
      this % a_f = a_f
      this % a_z = a_z
     
      END SUBROUTINE mddc_mrf_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a mddc_desc
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (mddc_desc), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_desc_destroy: '

!  Start of executable code

!  Destroy scalar components
      this % s_name = ' '
      this % l_name = ' '
      this % units = ' '
      this % mddc_type = ' '
      this % sigma_default = zero
      this % flux_factor = zero
      this % l_mdcoil_def = .FALSE.

!  Destroy Derived Types
      CALL bsc_destroy(this % mdcoil)
      CALL mddc_destroy(this % mrf)

      END SUBROUTINE mddc_desc_destroy

!-------------------------------------------------------------------------------
!  Destroy a mddc_mrf
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_mrf_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (mddc_mrf), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_mrf_destroy: '
      INTEGER :: ier1

!  Start of executable code

!  Get rid of all components
 
!  Scalar variables
      this % code_name = ' '
      this % code_version = ' '
      this % date_run = ' '
      this % field_coils_id = ' '
      this % kp = 0
      this % rmin = zero
      this % rmax = zero
      this % zmin = zero
      this % zmax = zero
      this % n_field_periods = 0
      this % lstell_sym = .false.

!  Array Sizes
      this % n_field_cg = 0
      this % ir = 0
      this % jz = 0
      this % kp_store = 0

!  Deallocate space for arrays
      IF (ASSOCIATED(this % rdiag_coilg_1)) THEN
         DEALLOCATE(this % rdiag_coilg_1,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc rdiag_coilg_1')
      ENDIF
      IF (ASSOCIATED(this % extcur_mg)) THEN
         DEALLOCATE(this % extcur_mg,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc extcur_mg')
      ENDIF
      IF (ASSOCIATED(this % a_r)) THEN
         DEALLOCATE(this % a_r,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc a_r')
      ENDIF
      IF (ASSOCIATED(this % a_f)) THEN
         DEALLOCATE(this % a_f,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc a_f')
      ENDIF
      IF (ASSOCIATED(this % a_z)) THEN
         DEALLOCATE(this % a_z,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc a_z')
      ENDIF

      END SUBROUTINE mddc_mrf_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for mddc_desc
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_desc_assign(left,right)

!  12-11-04. Can't get by with intrinsic assignment, because intrinsic 
!  assignment for the mdcoil component would  give incorrect results.

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (mddc_desc), INTENT (inout) :: left
      TYPE (mddc_desc), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_desc_assign: '
         
!  Start of executable code
      left % s_name = right % s_name
      left % l_name = right % l_name
      left % units = right % units
      left % mddc_type = right % mddc_type
      left % l_mdcoil_def = right % l_mdcoil_def
      left % sigma_default = right % sigma_default
      left % flux_factor = right % flux_factor
      left % mdcoil = right % mdcoil
      left % mrf = right % mrf
         
      END SUBROUTINE mddc_desc_assign

!-------------------------------------------------------------------------------
!  Assignment for mddc_mrf
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_mrf_assign(left,right)
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (mddc_mrf), INTENT (inout) :: left
      TYPE (mddc_mrf), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_mrf_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right pointers are the same?. FIX IT'
      INTEGER :: ier1, ier2, ier3
      LOGICAL, DIMENSION(5) :: lassert
         
!  Start of executable code

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      lassert(1) = .not.ASSOCIATED(left % a_r,right % a_r)
      lassert(2) = .not.ASSOCIATED(left % a_f,right % a_f)
      lassert(3) = .not.ASSOCIATED(left % a_z,right % a_z)
      lassert(4) = .not.ASSOCIATED(left % rdiag_coilg_1,                       &
     &                            right % rdiag_coilg_1)
      lassert(5) = .not.ASSOCIATED(left % extcur_mg,right % extcur_mg)
      CALL assert(lassert,sub_name // err_mess1)
      
!  Destroy left
      CALL mddc_mrf_destroy(left)

!  Scalar variables
      left % code_name = right % code_name
      left % code_version = right % code_version
      left % date_run = right % date_run
      left % field_coils_id = right % field_coils_id
      left % kp = right % kp
      left % rmin = right % rmin
      left % rmax = right % rmax
      left % zmin = right % zmin
      left % zmax = right % zmax
      left % n_field_periods = right % n_field_periods
      left % lstell_sym = right % lstell_sym
      left % n_field_cg = right % n_field_cg
      left % ir = right % ir
      left % jz = right % jz
      left % kp_store = right % kp_store

!  Allocate space for arrays (Were deallocated in _destroy)
      ALLOCATE(left % rdiag_coilg_1(left % n_field_cg),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rdiag_coilg_1')
      ALLOCATE(left % extcur_mg(left % n_field_cg),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc extcur_mg')

      ALLOCATE(left % a_r(left % ir,left % jz,left % kp_store),                &
     &   STAT=ier1)
      ALLOCATE(left % a_f(left % ir,left % jz,left % kp_store),                &
     &   STAT=ier2)
      ALLOCATE(left % a_z(left % ir,left % jz,left % kp_store),                &
     &   STAT=ier3)
      CALL assert_eq(0,ier1,ier2,ier3,sub_name // 'alloc a_')

!  Move arrays
!  JDH 2010-07-20. IF statements to avoid segmentation fault with
!  gfortran compiler.
      IF ( ASSOCIATED(right % rdiag_coilg_1)) THEN
         left % rdiag_coilg_1 = right % rdiag_coilg_1
      ENDIF
      IF ( ASSOCIATED(right % extcur_mg)) THEN
        left % extcur_mg = right % extcur_mg
      ENDIF
      IF ( ASSOCIATED(right % a_r)) THEN
        left % a_r = right % a_r
      ENDIF
      IF ( ASSOCIATED(right % a_f)) THEN
        left % a_f = right % a_f
      ENDIF
      IF ( ASSOCIATED(right % a_z)) THEN
        left % a_z = right % a_z
      ENDIF
         
      END SUBROUTINE mddc_mrf_assign
          
!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a mddc_desc
!-------------------------------------------------------------------------------

      SUBROUTINE mddc_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (mddc_desc), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      CHARACTER (len=60) :: id

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(10) :: fmt1 = (/                  &
     & '(" start mddc_desc write, called with id = ",a)      ',                &
     & '(" s_name = ",a)                                     ',                &
     & '(" l_name = ",a)                                     ',                &
     & '(" units = ",a)                                      ',                &
     & '(" l_mdcoil_def = ",L1)                              ',                &
     & '(" mddc_type = ",a)                                  ',                &
     & '(" bsc_coil s_name = ",a)                            ',                &
     & '(" flux_factor = ",es12.5)                           ',                &
     & '(" sigma_default = ",es12.5)                         ',                &
     & '(" end mddc_desc write, called with id = ",a)        '                 &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

!  Select Case of Verbosity Level
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,*) this % s_name
         WRITE(iou,*) this % l_name
         WRITE(iou,*) this % units
         WRITE(iou,*) this % l_mdcoil_def
         WRITE(iou,*) this % mddc_type
         WRITE(iou,*) this % mdcoil % s_name
         WRITE(iou,*) this % flux_factor
         WRITE(iou,*) this % sigma_default
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % s_name
         WRITE(iou,fmt1(3)) this % l_name
         WRITE(iou,fmt1(4)) this % units
         WRITE(iou,fmt1(5)) this % l_mdcoil_def
         WRITE(iou,fmt1(6)) this % mddc_type
         WRITE(iou,fmt1(7)) this % mdcoil % s_name
         WRITE(iou,fmt1(8)) this % flux_factor
         WRITE(iou,fmt1(9)) this % sigma_default
         WRITE(iou,fmt1(10)) id
      
      END SELECT

      END SUBROUTINE mddc_desc_write

!-------------------------------------------------------------------------------
!  Write out the contents of a mddc_mrf
!-------------------------------------------------------------------------------
!
      SUBROUTINE mddc_mrf_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (mddc_mrf), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, n_data
      INTEGER      :: i1, i2, i3

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(21) :: fmt1 = (/                  &
     & '(" start mddc_mrf write, called with id = ",a)       ',                &
     & '(" code_name  = ",a)                                 ',                &
     & '(" code_version = ",a)                               ',                &
     & '(" date_run = ",a)                                   ',                &
     & '(" field_coils_id = ",a)                             ',                &
     & '(" number of field-coil groups (n_field_cg) = ",i4)  ',                &
     & '(" index   rdiag_coilg_1: ",/,(1x,i4,3x,es12.5))     ',                &
     & '(" index   extcur_mg: ",/,(1x,i4,3x,es12.5))         ',                &
     & '(" number of grid points in R (ir) = ",i4)           ',                &
     & '(" number of grid points in z (jz) = ",i4)           ',                &
     & '(" number of grid points in phi (kp) = ",i4)         ',                &
     & '(" number of g. p. in phi stored (kp_store) = ",i4)  ',                &
     & '(" minimum R in grid (rmin) = ",es12.5)              ',                &
     & '(" maximum R in grid (rmax) = ",es12.5)              ',                &
     & '(" minimum Z in grid (zmin) = ",es12.5)              ',                &
     & '(" maximum Z in grid (zmax) = ",es12.5)              ',                &
     & '(" number of field periods (n_field_periods) = ",i4) ',                &
     & '(" Stellarator symmetry logical (lstell_sym) = ",l1) ',                &
     & '(" Three indices for a_ are ",i4,2x,i4,2x,i4)        ',                &
     & '(" a_r, a_f, a_z = ",3(3x,es12.5))                   ',                &
     & '(" end mddc_mrf write, called with id = ",a)         '                 &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF
      
!  Index values for single array value print out
      i1 = this % ir / 2
      i2 = this % jz / 2
      i3 = this % kp_store / 2

!  Select Case of Verbosity Level
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,*) this % code_name
         WRITE(iou,*) this % code_version
         WRITE(iou,*) this % date_run
         WRITE(iou,*) this % field_coils_id
         WRITE(iou,*) this % n_field_cg
         WRITE(iou,*) (i,this % rdiag_coilg_1(i),i=1,this % n_field_cg)
         WRITE(iou,*) (i,this % extcur_mg(i),i=1,this % n_field_cg)
         WRITE(iou,*) this % ir
         WRITE(iou,*) this % jz
         WRITE(iou,*) this % kp
         WRITE(iou,*) this % kp_store
         WRITE(iou,*) this % rmin
         WRITE(iou,*) this % rmax
         WRITE(iou,*) this % zmin
         WRITE(iou,*) this % zmax
         WRITE(iou,*) this % n_field_periods
         WRITE(iou,*) this % lstell_sym
         WRITE(iou,*) i1, i2, i3
         WRITE(iou,*) this % a_r(i1,i2,i3), this % a_f(i1,i2,i3),              &
     &      this % a_z(i1,i2,i3)
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % code_name
         WRITE(iou,fmt1(3)) this % code_version
         WRITE(iou,fmt1(4)) this % date_run
         WRITE(iou,fmt1(5)) this % field_coils_id
         WRITE(iou,fmt1(6)) this % n_field_cg
         WRITE(iou,fmt1(7)) (i,this % rdiag_coilg_1(i),                        &
     &      i=1,this % n_field_cg)
         WRITE(iou,fmt1(8)) (i,this % extcur_mg(i),                            &
     &      i=1,this % n_field_cg)
         WRITE(iou,fmt1(9)) this % ir
         WRITE(iou,fmt1(10)) this % jz
         WRITE(iou,fmt1(11)) this % kp
         WRITE(iou,fmt1(12)) this % kp_store
         WRITE(iou,fmt1(13)) this % rmin
         WRITE(iou,fmt1(14)) this % rmax
         WRITE(iou,fmt1(15)) this % zmin
         WRITE(iou,fmt1(16)) this % zmax
         WRITE(iou,fmt1(17)) this % n_field_periods
         WRITE(iou,fmt1(18)) this % lstell_sym
         WRITE(iou,fmt1(19)) i1, i2, i3
         WRITE(iou,fmt1(20)) this % a_r(i1,i2,i3), this % a_f(i1,i2,i3),       &
     &      this % a_z(i1,i2,i3)
         WRITE(iou,fmt1(21)) id

      END SELECT

      END SUBROUTINE mddc_mrf_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2007-06-11. First version of mddc_T. Copied and edited from diagnostic_T
!
!  
!-------------------------------- (diagnostic_T comments below) ----------------
! JDH 07-16-04. Modifying bsc.f to get diagnostic_mod.f
! JDH 08-11-04. More modifications. File diagnostic_T.f
!
!  JDH 08-16-2004
!     Add comments for collection, _ptr
!     Think about n_data, and where it belongs
!  JDH 08-19-2004
!     Eliminated n_data as a component
!  JDH 08-23-2004
!     Cleaned up sigma logic a bit.
!  JDH 09-10-2004
!     Added mddc_type component
!  JDH 12-11-2004
!     Removed 'pointer' attribute from mdcoil component of diagnostic_desc. 
!     Added l_mdcoil_def component to diagnostic_desc. Added subroutine 
!     diagnostic_desc_assign.
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-01-21
!    SPH Jan 2008 eliminated iprec. Completed elimination.
!    Initialized STAT variables in mddc_data_destroy
!       
!  JDH 2009-06-15
!    Eliminated mddc_data derived type - not needed.
!
!  JDH 2010-07-20
!    Added IF test for association in mddc_mrf_assign, to avoid
!    segmentation fault with gfortran compiler. Thanks to J Geiger.

       
      END MODULE mddc_T
