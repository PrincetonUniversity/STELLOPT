      SUBROUTINE load_saddle_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE saddle_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, nsv
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

!     Load the coil currents with values from variables in
!     optimization

      nsv = 0

      IF (lsadcur) THEN
!     Vary saddle currents
         DO i = 1, num_cursad
            IF (ls_cur(i)) THEN
               nsv = nsv + 1
               cursad(i) = xvariables(nsv)
            END IF
         END DO
      END IF

      nvariables = nsv

      END SUBROUTINE load_saddle_currents
