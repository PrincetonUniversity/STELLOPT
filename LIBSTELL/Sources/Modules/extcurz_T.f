!*******************************************************************************
!  File extcurz_T.f
!  Contains module extcurz_T
!  Defines derived-types: extcurz_desc
!
!*******************************************************************************
!  MODULE extcurz_T
!    (extcurz Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES
! SECTION VIII. PRIVATE ROUTINES USED IN SXR_T
!
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      MODULE extcurz_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds, ONLY: rprec
      USE stel_constants, ONLY: pi, zero
      USE safe_open_mod ! from LIBSTELL/MODULES

      IMPLICIT NONE

      PRIVATE rprec, pi, zero

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
!     extcurz Description:
!       extcurz_desc  
!     Type of diagnostic specified by  % d_type = 'extcurz'.
!
!*******************************************************************************

      TYPE extcurz_desc
        REAL(rprec) :: s0, u0 ! s,u values describing the circuit
      END TYPE extcurz_desc

      CONTAINS

!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************

      SUBROUTINE extcurz_desc_construct(this,s0,u0)

         IMPLICIT NONE

         TYPE (extcurz_desc), INTENT(inout) :: this
         REAL(rprec), INTENT(in) :: s0, u0

         this % s0 = s0
         this % u0 = u0

      END SUBROUTINE extcurz_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************

      SUBROUTINE extcurz_desc_destroy(this)

         TYPE (extcurz_desc), INTENT(inout) :: this

         this % s0 = zero
         this % u0 = zero

      END SUBROUTINE extcurz_desc_destroy

!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a extcurz_desc
!  if iou and filaname are present - write to file
!  if iou and filename are not present - write to stdout (screen)
!
! THIS NEEDS MODIFYING TO BE ABLE TO APPEND RECORDS AND NOT OVERWRITE FILES
!-------------------------------------------------------------------------------

      SUBROUTINE extcurz_desc_write(this,iounit,filename)

      IMPLICIT NONE

      TYPE (extcurz_desc),   INTENT(in) :: this
      INTEGER, OPTIONAL,       INTENT(in) :: iounit
      CHARACTER*300, OPTIONAL, INTENT(in) :: filename

      INTEGER :: iou = 6
      INTEGER :: istat = 0      !status of safe_open call
      
      IF (PRESENT(iounit).AND.PRESENT(filename)) THEN
         iou=iounit
         CALL safe_open(iou,istat,filename,'replace','formatted')
         WRITE(iou,*) 's0 - ', this % s0
         WRITE(iou,*) 'u0 - ', this % u0
      ELSE
         WRITE(*,*) 's0 - ', this % s0
         WRITE(*,*) 'u0 - ', this % u0
      END IF
      
      END SUBROUTINE extcurz_desc_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  GLT 10-sep-2012
!     first draft
!

      END MODULE extcurz_T

