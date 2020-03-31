      SUBROUTINE close_all_files
      USE vparams, ONLY: nmac, nthreed
      IMPLICIT NONE
C-----------------------------------------------

      IF (nthreed .gt. 0) CLOSE (nthreed)
      IF (nmac .gt. 0) CLOSE (nmac)

      END SUBROUTINE close_all_files
