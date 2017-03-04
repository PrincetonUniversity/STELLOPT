      SUBROUTINE load_bg_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bcoils_mod
      USE biotsavart
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, m
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
      REAL(rprec)  :: current, current_first
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

      ! ADDED by SAL
         DO i = 1, mbcoils
            DO m = 1, coil_group(i) % ncoil
               current = coil_group(i) % coils(m) % current
               IF (m .eq. 1) current_first = current
               IF (current_first .ne. 0.0) coil_group(i) % 
     1             coils(m) % current = (current/current_first)
     1               *bcoil_cur(i)
            END DO
         END DO

      ELSE

         nvariables = 0

      END IF

      END SUBROUTINE load_bg_currents
