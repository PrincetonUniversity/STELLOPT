      SUBROUTINE load_bg_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bcoils_mod
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, m
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

!     Load variable background currents cc_bg into bcoil_cur
!     based on INDEX mc_bg

      IF (lbcoil_cur) THEN

         DO m = 1, mc_max
            cc_bg(m) = xvariables(m)
         END DO

         DO i = 1, mbcoils
            m = mc_bg(i)
            IF (m .gt. 0) THEN
               bcoil_cur(i) = cc_bg(m)
            END IF
         END DO

         nvariables = mc_max

      ELSE

         nvariables = 0

      END IF

      END SUBROUTINE load_bg_currents
