!*******************************************************************************
!  File boundary_T.f
!  Contains module boundary_T
!  Defines derived-types: boundary_desc
!  A type of Geometric - Boundary signal
!
!*******************************************************************************
!  MODULE ip_T
!    (ip Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES
! SECTION VIII. PRIVATE ROUTINES USED IN ip_T
!
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE boundary_T
      USE safe_open_mod

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      IMPLICIT NONE

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
!     boundary Description:
!       boundary_desc
!     Type of geometric specified by  % d_type = 'boundary'.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Declare type boundary_desc
!         n_index       - Toroidal coefficient index.
!         m_index       - Poloidal coefficient index.
!         coefficient   - Coefficient name. This should be either rbc or zbs for
!                         symmetric cases or rbs or zbc for asymmetric cases.
!-------------------------------------------------------------------------------
      TYPE boundary_desc
         INTEGER          :: n_index
         INTEGER          :: m_index
         CHARACTER(LEN=3) :: coefficientName
      END TYPE

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a boundary_desc
!
!  For d_type = 'boundary' (boundary signal)
!-------------------------------------------------------------------------------
      SUBROUTINE boundary_desc_construct(this, coefficientName,                &
     &                                   n_index, m_index)

      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Argument Declarations
!-------------------------------------------------------------------------------
      TYPE (boundary_desc), INTENT(inout) :: this
      INTEGER, INTENT(in)                 :: n_index
      INTEGER, INTENT(in)                 :: m_index
      CHARACTER(LEN=3), INTENT(in)        :: coefficientName

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      this % n_index = n_index
      this % m_index = m_index
      this % coefficientName = coefficientName

      END SUBROUTINE boundary_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy an boundary_desc
!
!  ARGUMENT
!  this     - an boundary_desc
!-------------------------------------------------------------------------------
      SUBROUTINE boundary_desc_destroy(this)

      TYPE (boundary_desc),INTENT(inout) :: this

      this % n_index = 0
      this % m_index = 0
      this % coefficientName = '   '

      END SUBROUTINE boundary_desc_destroy

!*******************************************************************************
! SECTION VI.   ASSIGNMENT SUBROUTINES
!
!   These are not needed because the intrinsic assignments work
!*******************************************************************************

!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a boundary_desc
!  if iou and filaname are present - write to file
!  if iou and filename are not present - write to stdout (screen)
!
! THIS NEEDS MODIFYING TO BE ABLE TO APPEND RECORDS AND NOT OVERWRITE FILES
!-------------------------------------------------------------------------------

      SUBROUTINE boundary_desc_write(this,iounit,filename)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments
! this     - boundary
! iou      - output io unit number
! filename - output file name
!-------------------------------------------------------------------------------

      TYPE (boundary_desc),INTENT(in)   :: this
      INTEGER, OPTIONAL,INTENT(in)      :: iounit
      CHARACTER*300,OPTIONAL,INTENT(in) :: filename
!-------------------------------------------------------------------------------
! Local Variables
! iou  - iounit to use
! istat - status of file opening
!-------------------------------------------------------------------------------
      INTEGER :: iou = 6
      INTEGER :: istat = 0      !status of safe_open call

      IF (PRESENT(iounit).AND.PRESENT(filename)) THEN
         iou=iounit
         CALL safe_open(iou,istat,filename,'replace','formatted')
         WRITE(iou,*) 'n_index         - ', this % n_index
         WRITE(iou,*) 'm_index         - ', this % m_index
         WRITE(iou,*) 'coefficientName - ', this % coefficientName
      ELSE
         WRITE(*,*) 'n_index         - ', this % n_index
         WRITE(*,*) 'm_index         - ', this % m_index
         WRITE(*,*) 'coefficientName - ', this % coefficientName
      END IF

      END SUBROUTINE boundary_desc_write
!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      END MODULE boundary_T
