!============================================================================
      MODULE ptr_links
!-------------------------------------------------------------------
! Module to define linked list for int/pol parsing routines
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!-------------------------------------------------------------------
      USE ip_beamline  ! module containing int_pol TYPE definition
      IMPLICIT NONE

      TYPE :: beamlink
        REAL :: value1
        TYPE(ip_beam) :: beam
        TYPE(beamlink) , POINTER :: p => null()
      END TYPE

      TYPE :: iplink
        REAL :: value1
        TYPE(int_pol) :: ip
        TYPE(iplink) , POINTER :: p => null()
      END TYPE

      TYPE :: ip_ptr
        TYPE(int_pol), POINTER ::  p => null()
      END TYPE

      END MODULE ptr_links
