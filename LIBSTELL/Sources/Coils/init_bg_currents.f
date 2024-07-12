      SUBROUTINE init_bg_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bcoils_mod
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, m, iunit, istat
      INTEGER, INTENT(out) :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------
      mbcoils = 0
      iunit = 77
      CALL safe_open(iunit, istat, TRIM(bcoil_file), 'old',
     1                       'formatted')
      IF (istat .eq. 0) THEN
         READ (iunit, *) mbcoils
         CLOSE (iunit)
      ELSE
         PRINT *, 'Background file - ', TRIM(bcoil_file),
     1            ' - could not be opened'
      END IF

!     Set variable background currents based on index mc_bg

      mc_max = 0
      DO i = 1, mbcoils
         m = mc_bg(i)
         IF (m .gt. mc_max) THEN
            cc_bg(m) = bcoil_cur(i)
            mc_max = m
         END IF
      END DO

!
!     lbcoil_cur == .false. is equivalent to mc_max = 0
!
      IF (mc_max == 0) lbcoil_cur = .false.

      IF (lbcoil_cur) THEN
         DO m = 1, mc_max
            xvariables(m) = cc_bg(m)
         END DO
      END IF

      nvariables = mc_max

      END SUBROUTINE init_bg_currents
