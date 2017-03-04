      SUBROUTINE init_bg_currents (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bcoils_mod
      USE safe_open_mod
      USE biotsavart
      USE mpi_params      
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, m, j, iunit, istat
      INTEGER, INTENT(out) :: nvariables
      REAL(rprec)  :: current, current_first
      REAL(rprec) :: xvariables(*)
      LOGICAL :: lcoil_open
!-----------------------------------------------
      mbcoils = 0
      iunit = 77

      CALL cleanup_biotsavart
      INQUIRE(FILE=TRIM(bcoil_file),NUMBER=iunit,
     1        OPENED=lcoil_open)
      IF (lcoil_open) CLOSE(iunit)
      CALL parse_coils_file(TRIM(bcoil_file))
      mbcoils = SIZE(coil_group) !SAL
      DO i = 1, mbcoils
         DO j = 1, coil_group(i) % ncoil
            current = coil_group(i) % coils(j) % current
            IF (j .eq. 1) current_first = current
            IF (current_first .ne. 0.0) coil_group(i) 
     1             % coils(j) % current = (current/current_first)
     2                                   *bcoil_cur(i)
         END DO
      END DO
      IF (myid .eq.master) THEN
         DO i = 1, mbcoils
            IF (ABS(bcoil_cur(i)).ge. 1.0E6) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',
     1               coil_group(i)%ncoil,'  EXTCUR = ',
     2               bcoil_cur(i)/1.0E6,' [MA]'
            ELSE IF (ABS(bcoil_cur(i)).ge. 1.0E3) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',
     1               coil_group(i)%ncoil,'  EXTCUR = ',
     2               bcoil_cur(i)/1.0E3,' [kA]'
            ELSE
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',
     1               coil_group(i)%ncoil,'  EXTCUR = ',
     2               bcoil_cur(i),' [A]'
            END IF
         END DO
      END IF
      CALL FLUSH(6)

!      CALL safe_open(iunit, istat, TRIM(bcoil_file), 'old',
!     1                       'formatted')
!      IF (istat .eq. 0) THEN
!         READ (iunit, *) mbcoils
!         CLOSE (iunit)
!      ELSE
!         PRINT *, 'Background file - ', TRIM(bcoil_file),
!     1            ' - could not be opened'
!      END IF

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
