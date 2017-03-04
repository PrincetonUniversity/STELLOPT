!*******************************************************************************
!  File geometric_T.f
!  Contains module geometric_T
!  Defines derived-types: geometric_desc

!*******************************************************************************
!  MODULE geometric_T
!    (geometric Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE geometric_T

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
      USE edge_limit_T
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
      INTEGER, PARAMETER, PRIVATE :: name_len=30      
      INTEGER, PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
!   geometric Description:
!       geometric_desc  
!     Type of geometric specified by  % g_type.
!     Allowable values of g_type:
!       'edge_limit'  - edge of the plasma (s = 1 surface)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type geometric_desc
!   Common to all g_types                                                     
!       g_type          character, type of geometric
!       name            character, name of geometric
!       units           character, physical units that the data is measured in
!       sigma_default   real, default value of the uncertainty in the data
!
!  Derived Types for various geometric types
!  Now, only edge_limit_desc - edge of the plasma
!      el_desc
!-------------------------------------------------------------------------------
      TYPE geometric_desc
         CHARACTER (len=type_len)  :: g_type
         CHARACTER (len=name_len)  :: name                                 
         CHARACTER (len=units_len) :: units                                 
         REAL(rprec)               :: sigma_default
         TYPE (edge_limit_desc)    :: el_desc 
      END TYPE geometric_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE geometric_desc_assign
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE geometric_construct
         MODULE PROCEDURE geometric_desc_construct
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE geometric_destroy
         MODULE PROCEDURE geometric_desc_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE geometric_write
         MODULE PROCEDURE geometric_desc_write
      END INTERFACE

!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a geometric_desc
!
!  For g_type = 'edge_limit' (edge - limiter)
!-------------------------------------------------------------------------------
      SUBROUTINE geometric_desc_construct(this,g_type,name,                    &
     &   units,sigma_default,el_desc)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (geometric_desc), INTENT(inout)       :: this
      CHARACTER (len=*), INTENT(in)              :: g_type
      CHARACTER (len=*), INTENT(in)              :: name
      CHARACTER (len=*), INTENT(in)              :: units
      REAL(rprec), INTENT(in)                    :: sigma_default
      TYPE (edge_limit_desc), INTENT(in)         :: el_desc

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'geometric_desc_construct_mddc: '

!  Start of executable code

!  Destroy the mddc component
      CALL edge_limit_destroy(this % el_desc)

!  Scalar assignments
      this % name = TRIM(ADJUSTL(name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default

!  Different coding, depending on g_type
      SELECT CASE (TRIM(ADJUSTL(g_type)))
      CASE ('edge_limit')
         this % g_type = 'edge_limit'
         this % el_desc =  el_desc

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized g_type: ',                   &
     &      char=g_type)
      END SELECT ! Different coding depending on g_type
      

      END SUBROUTINE geometric_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a geometric_desc
!-------------------------------------------------------------------------------
      SUBROUTINE geometric_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (geometric_desc), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'geometric_desc_destroy: '

!  Start of executable code

!  Get rid of all components
      this % name = ' '
      this % units = ' '
      
!  Different coding, depending on d_type
      SELECT CASE (TRIM(ADJUSTL(this % g_type)))
      CASE ('edge_limit')
         this % g_type = ' '

!  edge_limit component
         CALL edge_limit_destroy(this % el_desc)

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized g_type: ',                   &
     &      char=this % g_type)
      END SELECT ! Different coding depending on g_type

      END SUBROUTINE geometric_desc_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for geometric_desc
!-------------------------------------------------------------------------------
      SUBROUTINE geometric_desc_assign(left,right)

!  12-11-04. Can't get by with intrinsic assignment, because intrinsic 
!  assignment for the edge_limit_desc component might give incorrect results.

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (geometric_desc), INTENT (inout) :: left
      TYPE (geometric_desc), INTENT (in) :: right
      
!  Declare temporary variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'geometric_desc_assign: '
         
!  Start of executable code
      left % g_type = right % g_type
      left % name = right % name
      left % units = right % units
      left % sigma_default = right % sigma_default
      left % el_desc = right % el_desc
         
      END SUBROUTINE geometric_desc_assign
          
!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents of a geometric_desc
!-------------------------------------------------------------------------------

      SUBROUTINE geometric_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (geometric_desc), INTENT (in) :: this
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
      CHARACTER(len=*), PARAMETER, DIMENSION(6) :: fmt1 = (/                   &
     & '(" start geometric_desc write, called with id = ",a)',                 &
     & '(" g_type = ",a)                                    ',                 &
     & '(" name = ",a)                                      ',                 &
     & '(" units = ",a)                                     ',                 &
     & '(" el_desc name = ",a)                              ',                 &
     & '(" end geometric_desc write, called with id = ",a)  '                  &
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
         WRITE(iou,*) this % g_type
         WRITE(iou,*) this % name
         WRITE(iou,*) this % units
         WRITE(iou,*) this % el_desc % name
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % g_type
         WRITE(iou,fmt1(3)) this % name
         WRITE(iou,fmt1(4)) this % units
         WRITE(iou,fmt1(5)) this % el_desc % name
         WRITE(iou,fmt1(6)) id
      
      END SELECT

      END SUBROUTINE geometric_desc_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-01-21.  First version of geometric_T. Copied and edited from
!      diagnostic_T
!
!  ---------- Below - comments from diagnostic_T  ------------------------------
!  
!  JDH 2007-06-11. Modified so that mddc_desc derived type is in geometric_desc
!    In preparation for different types of geometrics
!
! JDH 07-16-04. Modifying bsc.f to get geometric_mod.f
! JDH 08-11-04. More modifications. File geometric_T.f
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
!     Removed 'pointer' attribute from mdcoil component of geometric_desc. 
!     Added l_mdcoil_def component to geometric_desc. Added subroutine 
!     geometric_desc_assign.
!  JMS 6-22-2007
!     Added ipsl as a subTYPE of geometric_desc. 
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-01-21
!    SPH Jan 2008 eliminated iprec. Completed elimination.
!    Initialized STAT variables in geometric_data_destroy
!
!  JDH 2009-06-15
!    Eliminated geometric_data derived type - not needed

      END MODULE geometric_T
