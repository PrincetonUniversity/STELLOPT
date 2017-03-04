      SUBROUTINE init_modular_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nc, i, n
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

      nvariables = 0
      nmod_currents = 0
      nc = nmod_coils_per_period

      IF (.not.lmodcur .or. nc.le.0) RETURN


!     Initialize the variables to values of coil currents
!     and count the number of variables

      n = 0
      DO i = 1, nmid
         n = n + 1
         xvariables(n) = curmod(i)
      END DO

      nmod_currents = n
      nvariables = n

      END SUBROUTINE init_modular_currents
