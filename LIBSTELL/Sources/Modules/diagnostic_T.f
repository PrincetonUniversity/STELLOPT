!*******************************************************************************
!  File diagnostic_T.f
!  Contains module diagnostic_T
!  Defines derived-types: diagnostic_desc

!*******************************************************************************
!  MODULE diagnostic_T
!    (Diagnostic Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE diagnostic_T

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
      USE mddc_T
      USE sxrch_T             ! GJH 2009-01-20
      USE ipch_T              ! JDH 2012-03-15
      USE thscte_T            ! JDH 2011-10-23
      USE extcurz_T           ! GLT 10-sep-2012
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
!   Diagnostic Description:
!       diagnostic_desc  
!     Type of diagnostic specified by  % d_type.
!     Allowable values of d_type:
!       mddc  - Magnetic 'Diagnostic-Dot' Coil
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type diagnostic_desc
!   Common to all d_types                                                     
!       d_type          character, type of diagnostic
!       s_name          character, short name of diagnostic
!       l_name          character, long name of diagnostic
!       units           character, physical units that the data is measured in
!       sigma_default   real, default value of the uncertainty in the data
!
!  Derived Types for various diagnostic types
!      mddc       Magnetic Diagnostic Dot Coil
!      sxrch      Soft X-Ray CHord
!      ipch       Interferometry-Polarimetry CHord
!      thscte     THomson SCattering TE
!      extcurz    EXTernal CURrent along Z
!-------------------------------------------------------------------------------
      TYPE diagnostic_desc
         CHARACTER (len=type_len)  :: d_type
         CHARACTER (len=sn_len)    :: s_name                                 
         CHARACTER (len=ln_len)    :: l_name
         CHARACTER (len=units_len) :: units                                 
         REAL(rprec)               :: sigma_default
         TYPE (mddc_desc)          :: mddc
         TYPE (sxrch_desc)         :: sxrch      !GJH 2010-01-20
         TYPE (ipch_desc)          :: ipch       ! JDH 2012-03-15
         TYPE (thscte_desc)        :: thscte     ! JDH 2011-10-23
         TYPE (extcurz_desc)       :: extcurz    ! GLT 10-sep-2012
      END TYPE diagnostic_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE diagnostic_desc_assign
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_construct
         MODULE PROCEDURE diagnostic_desc_construct_mddc,                      &
     &                    diagnostic_desc_construct_sxrch,                     &     
     &                    diagnostic_desc_construct_ipch,                      &     
     &                    diagnostic_desc_cnstrct_thscte,                      &
     &                    diagnostic_desc_construct_extcurz
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_destroy
         MODULE PROCEDURE diagnostic_desc_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_write
         MODULE PROCEDURE diagnostic_desc_write
      END INTERFACE

!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a diagnostic_desc with mddc diagnostic
!
!  For d_type = 'mddc' (magnetic diagnostic-dot coil)
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_desc_construct_mddc(this,d_type,s_name,            &
     &   l_name, units,sigma_default,mddc)

!  NB. 
!  The mddc argument is assigned to the 'this' component. Do NOT call this
!  subroutine with this % mddc as the mddc argument.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (diagnostic_desc), INTENT(inout)      :: this
      CHARACTER (len=*), INTENT(in)              :: d_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      REAL(rprec), INTENT(in)                    :: sigma_default
      TYPE (mddc_desc), INTENT(in)               :: mddc

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_construct_mddc: '

!  Start of executable code

!  Destroy the mddc component
      CALL mddc_destroy(this % mddc)

!  Scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default

!  Different coding, depending on d_type
      SELECT CASE (TRIM(ADJUSTL(d_type)))
      CASE ('mddc')
         this % d_type = 'mddc'
         this % mddc =  mddc

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char=d_type)
      END SELECT ! Different coding depending on d_type
      

      END SUBROUTINE diagnostic_desc_construct_mddc

!======================================================================
      SUBROUTINE diagnostic_desc_construct_sxrch(this,d_type,s_name,           &
     &   l_name, units,sigma_default,sxrch)
!-------------------------------------------------------------------------------
!  Construct a diagnostic_desc with sxrch diagnostic
!
!  For d_type = 'sxrch' (Soft X-Ray Chord diagnostic)
! 
!-------------------------------------------------------------------------------
     
!  The sxrch argument is assigned to the 'this' component. Do NOT call this
!  subroutine with this % sxrch as the mddc argument.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (diagnostic_desc), INTENT(inout)      :: this
      CHARACTER (len=*), INTENT(in)              :: d_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      REAL(rprec), INTENT(in)                    :: sigma_default
      TYPE (sxrch_desc), INTENT(in)               :: sxrch
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_construct_sxrch: '

!  Start of executable code

!  Destroy the sxrch component
      CALL sxrch_desc_destroy(this % sxrch)

!  Scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default
      
!  Different coding, depending on d_type
      SELECT CASE (TRIM(ADJUSTL(d_type)))
      CASE ('sxrch')
         this % d_type = 'sxrch'
         this % sxrch =  sxrch

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char=d_type)
      END SELECT ! Different coding depending on d_type
     
     
      END SUBROUTINE diagnostic_desc_construct_sxrch

!======================================================================

      SUBROUTINE diagnostic_desc_construct_extcurz(this,d_type,                &
     &   s_name,l_name,units,sigma_default,extcurz,s0,u0)

      IMPLICIT NONE

      ! arguments
      TYPE (diagnostic_desc), INTENT(inout) :: this
      CHARACTER (len=*), INTENT(in)         :: d_type
      CHARACTER (len=*), INTENT(in)         :: s_name
      CHARACTER (len=*), INTENT(in)         :: l_name
      CHARACTER (len=*), INTENT(in)         :: units
      REAL(rprec), INTENT(in)               :: sigma_default
      TYPE (extcurz_desc), INTENT(in)       :: extcurz
      REAL(rprec), INTENT(in)               :: s0, u0
      
      ! local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_construct_extcurz: '

      ! destroy the extcurz component
      CALL extcurz_desc_destroy(this % extcurz)

      ! scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default
      
      ! extcurz specific coding
      SELECT CASE (TRIM(ADJUSTL(d_type)))
      CASE ('extcurz')
         this % d_type = 'extcurz'
         this % extcurz = extcurz
         this % extcurz % s0 = s0
         this % extcurz % u0 = u0
      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char=d_type)
      END SELECT

      END SUBROUTINE diagnostic_desc_construct_extcurz

!======================================================================
      SUBROUTINE diagnostic_desc_construct_ipch(this,d_type,s_name,            &
     &   l_name, units,sigma_default,ipch)
!-------------------------------------------------------------------------------
!  Construct a diagnostic_desc with ipch diagnostic
!
!  For d_type = 'ipch' (Interferometry-Polarimetry Chord diagnostic)
! 
!-------------------------------------------------------------------------------
     
!  The ipch argument is assigned to the 'this' component. Do NOT call this
!  subroutine with this % ipch as the ipch argument.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (diagnostic_desc), INTENT(inout)      :: this
      CHARACTER (len=*), INTENT(in)              :: d_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      REAL(rprec), INTENT(in)                    :: sigma_default
      TYPE (ipch_desc), INTENT(in)               :: ipch
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_construct_ipch: '

!  Start of executable code

!  Destroy the ipch component
      CALL ipch_desc_destroy(this % ipch)

!  Scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default
      
!  Different coding, depending on d_type
      SELECT CASE (TRIM(ADJUSTL(d_type)))
      CASE ('ipch')
         this % d_type = 'ipch'
         this % ipch =  ipch

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char=d_type)
      END SELECT ! Different coding depending on d_type
     
     
      END SUBROUTINE diagnostic_desc_construct_ipch

!======================================================================
      SUBROUTINE diagnostic_desc_cnstrct_thscte(this,d_type,s_name,          &
     &   l_name, units,sigma_default,thscte)
!-------------------------------------------------------------------------------
!  Construct a diagnostic_desc with thscte diagnostic
!
!  For d_type = 'thscte' (Thomson Scattering Te)
! 
!-------------------------------------------------------------------------------
     
!  The thscte argument is assigned to the 'this' component. Do NOT call this
!  subroutine with this % thscte as the mddc argument.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (diagnostic_desc), INTENT(inout)      :: this
      CHARACTER (len=*), INTENT(in)              :: d_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      REAL(rprec), INTENT(in)                    :: sigma_default
      TYPE (thscte_desc), INTENT(in)             :: thscte
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_cnstrct_thscte: '

!  Start of executable code

!  Destroy the thscte component
      CALL thscte_desc_destroy(this % thscte)

!  Scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default
      
!  Different coding, depending on d_type
      SELECT CASE (TRIM(ADJUSTL(d_type)))
      CASE ('thscte')
         this % d_type = 'thscte'
         this % thscte =  thscte

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char=d_type)
      END SELECT ! Different coding depending on d_type
     
     
      END SUBROUTINE diagnostic_desc_cnstrct_thscte

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a diagnostic_desc
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (diagnostic_desc), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_destroy: '

!  Start of executable code

!  Get rid of all components
      this % s_name = ' '
      this % l_name = ' '
      this % units = ' '
      
!  Different coding, depending on d_type
      SELECT CASE (TRIM(ADJUSTL(this % d_type)))
      CASE ('mddc')
         this % d_type = ' '
         CALL mddc_destroy(this % mddc)
      
      CASE ('sxrch')                              ! GJH 2010-01-20
         this % d_type = ' '
         CALL sxrch_desc_destroy(this % sxrch)
         
      CASE ('ipch')                              ! JDH 2012-03-16
         this % d_type = ' '
         CALL ipch_desc_destroy(this % ipch)
         
      CASE ('thscte')
         this % d_type = ' '
         CALL thscte_desc_destroy(this % thscte)
         
      CASE ('extcurz')
         this % d_type = ' '
         CALL extcurz_desc_destroy(this % extcurz)
         
      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char=this % d_type)
      END SELECT ! Different coding depending on d_type

      END SUBROUTINE diagnostic_desc_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for diagnostic_desc
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_desc_assign(left,right)

!  12-11-04. Can't get by with intrinsic assignment, because intrinsic 
!  assignment for the mdcoil component would  give incorrect results.

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (diagnostic_desc), INTENT (inout) :: left
      TYPE (diagnostic_desc), INTENT (in) :: right
      
!  Declare temporary variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_desc_assign: '
         
!  Start of executable code
      left % d_type = right % d_type
      left % s_name = right % s_name
      left % l_name = right % l_name
      left % units = right % units
      left % sigma_default = right % sigma_default
      left % mddc = right % mddc
      left % sxrch = right % sxrch          ! GJH 2010-01-20
      left % ipch = right % ipch            ! JDH 2012-03-16
      left % thscte = right % thscte        ! JDH 2011-10-23
      left % extcurz = right % extcurz      ! GLT 10-sep-2012
         
      END SUBROUTINE diagnostic_desc_assign
          
!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a diagnostic_desc
!-------------------------------------------------------------------------------

      SUBROUTINE diagnostic_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_desc), INTENT (in) :: this
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

!  Declare Format array  modified GJH 2010-01-15
      CHARACTER(len=*), PARAMETER, DIMENSION(10) :: fmt1 = (/                  &
     & '(" start diagnostic_desc write, called with id = ",a)',                &
     & '(" d_type = ",a)                                     ',                &
     & '(" s_name = ",a)                                     ',                &
     & '(" l_name = ",a)                                     ',                &
     & '(" units = ",a)                                      ',                &
     & '(" mddc s_name = ",a)                                ',                &
     & '(" sxrch s_name = ",a)                               ',                & 
     & '(" ipch s_name = ",a)                                ',                & 
     & '(" thscte s_name = ",a)                              ',                & 
     & '(" end diagnostic_desc write, called with id = ",a)  '                 &
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
         WRITE(iou,*) this % d_type
         WRITE(iou,*) this % s_name
         WRITE(iou,*) this % l_name
         WRITE(iou,*) this % units
         WRITE(iou,*) this % mddc % s_name
         WRITE(iou,*) this % sxrch % chord_name           ! GJH 2010-01-20
         WRITE(iou,*) this % ipch % chord_name           ! JDH 2012-03-16
         WRITE(iou,*) this % thscte % chord_name          ! JDH 2011-10-23
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % d_type
         WRITE(iou,fmt1(3)) this % s_name
         WRITE(iou,fmt1(4)) this % l_name
         WRITE(iou,fmt1(5)) this % units
         WRITE(iou,fmt1(6)) this % mddc % s_name
         WRITE(iou,fmt1(7)) this % sxrch % chord_name    ! GJH 2010-01-20
         WRITE(iou,fmt1(8)) this % ipch % chord_name           ! JDH 2012-03-16
         WRITE(iou,fmt1(9)) this % thscte % chord_name          ! JDH 2011-10-23
         WRITE(iou,fmt1(10)) id
      
      END SELECT

      END SUBROUTINE diagnostic_desc_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 2007-06-11. Modified so that mddc_desc derived type is in diagnostic_desc
!    In preparation for different types of diagnostics
!
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
!    Initialized STAT variables in diagnostic_data_destroy
!
!  JDH 2009-06-15
!    Eliminated diagnostic_data derived type - not needed.
!
!  JDH 2011-08-01
!    Refactor sxrc -> sxrch
!
!  JDH 2011-10-23
!    Add thscte
!
!  JDH 2012-03-16
!    Add ipch
!
!  JDH 2012-06-08
!    Change diagnostic_desc_construct_thscte to diagnostic_desc_cnstrct_thscte
!    as first form had 32 characters.
!
!  GLT 10-sep-2012
!    added coding for extcurz
!

      END MODULE diagnostic_T
