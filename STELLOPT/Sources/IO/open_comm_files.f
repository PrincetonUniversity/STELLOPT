      SUBROUTINE open_comm_files (mboz, nboz, ns_booz, ns_booz_max,
     1  ns_surf, ns_surf_max, lbootsj_opt, extension, iflag)
      USE stel_kinds
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: mboz, nboz, ns_booz(*), ns_booz_max,
     1   ns_surf(*), ns_surf_max
      INTEGER :: iflag
      LOGICAL :: lbootsj_opt
      CHARACTER(LEN=*), INTENT(in)  :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_booz=24, unit_boot=26
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit, k

      iunit = unit_booz

      k=22
      IF (ns_booz_max .ge. 1) THEN
         CALL safe_open(iunit, k, 'in_booz.' // TRIM(extension), 
     1                  'replace', 'formatted')
         IF (k .ne. 0) THEN
            PRINT *,'ERROR opening in_booz'
            iflag = -14
            RETURN
         END IF

         WRITE(iunit, *) mboz, nboz           !!User may have input these through &OPTIMUM NAMELIST
         WRITE(iunit, *) TRIM(extension)
         WRITE(iunit, *) ns_booz(1:ns_booz_max)
         CLOSE(iunit)
      END IF

      IF (lbootsj_opt) THEN
         iunit = unit_boot
         CALL safe_open(iunit, k, 'in_bootsj.'//TRIM(extension),
     1             'replace', 'formatted')
         IF (k .ne. 0) THEN
            iflag = -14
            RETURN
         END IF
         WRITE(iunit, *) TRIM(extension)
         WRITE(iunit, *) ns_surf(1:ns_surf_max)
         CLOSE(iunit)
      END IF

      END SUBROUTINE open_comm_files
