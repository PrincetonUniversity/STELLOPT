      SUBROUTINE init_saddle_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE saddle_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, nsv, nc
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

      nc = nsad_coils_per_period
      nsad_currents = 0
      nvariables = 0

      IF (nc .le. 0) RETURN

      num_cursad = MAXVAL(nsad_group(1:nsmid))
        
!     initialize the variables to values of unique coil parameters
!     and count the number of variables

      nsv = 0

      IF (lsadcur) THEN
         DO i=1, num_cursad
            IF (ls_cur(i)) THEN
               nsv = nsv + 1
               xvariables(nsv) = cursad(i)
            END IF
         END DO
      END IF

      nsad_currents = nsv
      nvariables = nsv

      END SUBROUTINE init_saddle_currents
