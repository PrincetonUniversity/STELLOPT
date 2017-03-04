
!*******************************************************************************
!  File edge_limit_T.f
!  Contains module edge_limit_T
!  Defines derived-types: edge_limit_desc
!  A type of geometric signal
!
!*******************************************************************************
!  MODULE edge_limit_T
!    (edge_limit Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE edge_limit_T

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
      USE limiter_iso_T
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
! 1)  edge_limit Description:
!       edge_limit_desc  
!     Type of geometric signal specified by  % g_type = 'edge_limit'.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Declare type edge_limit_desc
!    name             character, name of edge_limit
!    units            character, physical units that the edge_limit is measured in
!    edge_limit_type  character, specifies type of edge_limit
!                        iso_fun      only edge_limit type at this time
!                        polygon      future edge_limit type
!    sigma_default    real, default value of the uncertainty in the data
!    l_on_edge        logical
!                       True - want edge of plasma to just touch limiter
!                       False - want edge of plasma to be inside limiter
!    lim_iso          type limiter_iso, constructed for edge_limit_type='iso_fun'
!-------------------------------------------------------------------------------
      TYPE edge_limit_desc
         CHARACTER (len=name_len)       :: name                                 
         CHARACTER (len=30)             :: edge_limit_type                                 
         CHARACTER (len=units_len)      :: units                                 
         REAL(rprec)                    :: sigma_default
         LOGICAL                        :: l_on_edge
         TYPE (limiter_iso)             :: lim_iso 
      END TYPE edge_limit_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE edge_limit_construct
         MODULE PROCEDURE edge_limit_desc_construct
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE edge_limit_destroy
         MODULE PROCEDURE edge_limit_desc_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE edge_limit_write
         MODULE PROCEDURE edge_limit_desc_write
      END INTERFACE

!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a edge_limit_desc
!
!  For g_type = 'edge_limit' (geometric edge_limit)
!-------------------------------------------------------------------------------
      SUBROUTINE edge_limit_desc_construct(this,name,edge_limit_type,          &
     &   units,sigma_default,l_on_edge,lim_iso)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (edge_limit_desc), INTENT(inout)          :: this
      CHARACTER (len=*), INTENT(in)                  :: name
      CHARACTER (len=*), INTENT(in)                  :: edge_limit_type
      CHARACTER (len=*), INTENT(in)                  :: units
      REAL(rprec), INTENT(in)                        :: sigma_default
      LOGICAL, INTENT(in)                            :: l_on_edge
      TYPE (limiter_iso), INTENT(in)                 :: lim_iso  !Test without TARGET, 2009-01-31

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'edge_limit_desc_construct: '

!  Start of executable code

!  Destroy the limiter_iso component
      CALL limiter_iso_destroy(this % lim_iso)

!  Scalar assignments
      this % name = TRIM(ADJUSTL(name))
      this % units = TRIM(ADJUSTL(units))
      this % l_on_edge = l_on_edge
      
      SELECT CASE(TRIM(ADJUSTL(edge_limit_type)))
      CASE ('iso_fun')
         this % edge_limit_type = 'iso_fun'

      CASE ('polygon')
         WRITE(*,*) 'FATAL: edge_limit_desc_construct:'
         WRITE(*,*) 'edge_limit_type = polygon is NYI'
         STOP

      CASE DEFAULT
         WRITE(*,*) 'FATAL: edge_limit_desc_construct:'
         WRITE(*,*) 'unrecognized edge_limit_type = ',                         &
     &      TRIM(ADJUSTL(edge_limit_type))
         STOP
      END SELECT 
      
      this % sigma_default = sigma_default

!  Derived Type Assignments
      this % lim_iso = lim_iso
      
      END SUBROUTINE edge_limit_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a edge_limit_desc
!-------------------------------------------------------------------------------
      SUBROUTINE edge_limit_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (edge_limit_desc), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'edge_limit_desc_destroy: '

!  Start of executable code

!  Destroy scalar components
      this % name = ' '
      this % units = ' '
      this % edge_limit_type = ' '
      this % sigma_default = zero
      this % l_on_edge = .false.

!  Destroy Derived Types
      CALL limiter_iso_destroy(this % lim_iso)

      END SUBROUTINE edge_limit_desc_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for edge_limit_desc
!    Intrinsic assignment adequate - no pointer components
!-------------------------------------------------------------------------------
          
!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a edge_limit_desc
!-------------------------------------------------------------------------------

      SUBROUTINE edge_limit_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (edge_limit_desc), INTENT (in) :: this
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
      CHARACTER(len=*), PARAMETER, DIMENSION(7) :: fmt1 = (/                   &
     & '(" start edge_limit_desc write, called with id = ",a)      ',          &
     & '(" name = ",a)                                             ',          &
     & '(" edge_limit_type = ",a)                                  ',          &
     & '(" units = ",a)                                            ',          &
     & '(" sigma_default = ",es12.5)                               ',          &
     & '(" l_on_edge = ",l1)                                       ',          &
     & '(" end edge_limit_desc write, called with id = ",a)        '           &
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
         WRITE(iou,*) this % name
         WRITE(iou,*) this % edge_limit_type
         WRITE(iou,*) this % units
         WRITE(iou,*) this % sigma_default
         WRITE(iou,*) this % l_on_edge
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % name
         WRITE(iou,fmt1(3)) this % edge_limit_type
         WRITE(iou,fmt1(4)) this % units
         WRITE(iou,fmt1(5)) this % sigma_default
         WRITE(iou,fmt1(6)) this % l_on_edge
         WRITE(iou,fmt1(7)) id
      
      END SELECT

      END SUBROUTINE edge_limit_desc_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-01-21. First version of edge_limit_T. Copied and edited from mddc_T
!
!  JDH 2009-01-31. Change limiter_function -> limiter_iso
!
!  JDH 2009-02-04. Add l_on_edge.
!       
!  JDH 2009-06-15
!    Eliminated edge_limit_data derived type - not needed.
!
      END MODULE edge_limit_T
