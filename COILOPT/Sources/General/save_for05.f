      SUBROUTINE save_for05(for05file)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE coilsnamin
      USE safe_open_mod
      USE gade_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER for05file*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iunit=15, ierr
!-----------------------------------------------

      CALL safe_open(iunit, ierr, for05file, 'replace', 'formatted')
      CALL write_coilsin (iunit, ierr)
      IF (nopt_alg .GT. 0) CALL write_gade_nml(iunit)
      CLOSE(iunit)

      END SUBROUTINE save_for05
