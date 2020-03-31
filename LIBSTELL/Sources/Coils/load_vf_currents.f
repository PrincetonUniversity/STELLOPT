      SUBROUTINE load_vf_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------
!     Load the coil currents with values from optimization variables

      nvariables = 0

      IF (lvfc) THEN
        DO i=1, num_vf
          IF (lcc_vf(i)) THEN
             nvariables = nvariables + 1
             cc_vf(i) = xvariables(nvariables)
          END IF
        END DO
      END IF

      END SUBROUTINE load_vf_currents
