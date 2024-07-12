      SUBROUTINE load_modular_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)

!-----------------------------------------------

!     Load the coil currents with values from optimization variables

      n = 0

      IF (lmodcur) THEN
         DO i = 1, nmid
            n = n + 1
            curmod(i) = xvariables(n)
         END DO
      END IF

      nvariables = n

      END SUBROUTINE load_modular_currents
