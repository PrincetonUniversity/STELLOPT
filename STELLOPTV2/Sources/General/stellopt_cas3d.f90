!-----------------------------------------------------------------------
!     Subroutine:    stellopt_cas3d
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine calculates the boostrap current
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_cas3d(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_utils
      ! Boozer Coordinate transformation
      use safe_open_mod
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      integer, parameter :: nfax = 13
      INTEGER ::  ier, iunit_rzuv
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  STABILITY (CAS3D) CALCULATION  -------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','vmec2000_oneeq')
            CALL safe_open(iunit_rzuv, ier, 'for_cas3d_rzuv.dat','replace','formatted')
            WRITE(iunit_rzuv,*) 

         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  STABILITY (CAS3D) CALCULATION DONE  ----------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_cas3d
