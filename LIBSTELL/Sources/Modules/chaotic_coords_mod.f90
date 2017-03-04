!-----------------------------------------------------------------------
!     Module:        chaotic_coords_mod
!     Authors:       B. Israli (byi2000@columbia.edu)
!     Date:          06/16/2015
!     Description:   This module blah blah blah.
!-----------------------------------------------------------------------
      MODULE chaotic_coords_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      
!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER            :: integer_test
      DOUBLE PRECISION   :: RAXIS_CHAOS, ZAXIS_CHAOS
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0
      
!-----------------------------------------------------------------------
!     Subroutines
!         chaotic_coords_dump:  Dump some stuff to a file.
!         ccoords_find_axis:    Finds the magnetic axis.
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE chaotic_coords_dump(filename,istat)
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER :: iunit, ik
      CALL safe_open(iunit,istat,'chaotic_coords_dump.'//TRIM(filename),'unknown','formatted')

      WRITE(iunit,*) '#  V0(ik,1), V0(ik,2), V0(ik,3), FN(ik,1), FN(ik,2), FN(ik,3),',&
                        'V(ik,1), V(ik,2), V(ik,3), W(ik,1), W(ik,2), W(ik,3), d(ik)'
      CLOSE(iunit)
      RETURN
      END SUBROUTINE chaotic_coords_dump

      SUBROUTINE ccoords_find_axis(Raxis,Zaxis,iota,BFIELD_CC)
      DOUBLE PRECISION, INTENT(inout) :: Raxis
      DOUBLE PRECISION, INTENT(inout) :: Zaxis
      DOUBLE PRECISION, INTENT(out)   :: iota
      EXTERNAL                        :: BFIELD_CC
      INTEGER :: iunit, ik
      RETURN
      END SUBROUTINE ccoords_find_axis

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE chaotic_coords_mod
