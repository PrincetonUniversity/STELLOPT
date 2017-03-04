      SUBROUTINE chisq_jinvar (sigma, nrad, num, nopt, nu,
     1    sqrt_nsurf, iflag, extension, command, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, home_dir, exe_suffix
      USE optim_params, ONLY: NumJinvariant, lj_invariant
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt, nu, nrad
      INTEGER :: num, iflag
      CHARACTER(LEN=*) :: extension, command
      REAL(rprec) :: sigma(1,NumJinvariant), sqrt_nsurf
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER :: unit_jinvar = 30
      CHARACTER(LEN=20) :: temp
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, zero = 0, one = 1
      REAL(rprec), PARAMETER :: max_variation = 2.e-2_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit, ieps, n, istat, NJinvar
      CHARACTER(len=LEN_TRIM(home_dir)+20) :: version
      LOGICAL :: ex
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ljinvar
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: J_Invariant
      REAL(rprec) :: avg_Jinvar
      REAL(rprec) :: sqrt_nusurf
C-----------------------------------------------
      IF (ALL(ABS(sigma(1,1:NumJinvariant)) .ge. bigno)) RETURN
!
!        COMPUTE J-INVARIANT AT NUMJINVARIANT VALUES OF ep/mu RANGING FROM SLIGHTLY ABOVE
!        THE TRAPPED-PASSING BOUNDARY TO SLIGHTLY BELOW THE
!        DEEPLY TRAPPED-FORBIDDEN BOUNDARY.  THE PARAMETERS epl AND epu
!        DETERMINE DISTANCE TO THESE BOUNDARIES.
!
      version = TRIM(home_dir) // 'xj_invariant' // TRIM(exe_suffix)

      IF (nopt .gt. 0) THEN
!
!     RUN J-INVARIANT CODE TO CREATE OUTPUT FILE
!
         iunit = unit_jinvar

         WRITE (temp,'(1x,i3,1x,i3,1x,i3,a2)')
     1       nrad, NumJinvariant, nu, command

         CALL load_physics_codes (version, 'boozmn', temp,
     1         'j_invar_out', extension, iunit, iflag)
         IF (iflag .ne. 0) RETURN

         IF (nu .le. 0) RETURN
         sqrt_nusurf = SQRT(REAL(nu, rprec))
         ALLOCATE (J_Invariant(nu), ljinvar(nu))
!
!        Target ALL non-zero Jinvars to (d Jinvar/du) = 0
!
         DO ieps = 1, NumJinvariant
            READ (iunit, *, iostat=istat) (J_Invariant(n), n = 1, nu)
            IF (istat .ne. 0) THEN
               iflag = -18
               RETURN
            END IF

            ljinvar = J_Invariant .gt. zero
            avg_Jinvar = SUM(J_Invariant, mask=ljinvar)
            nJinvar = count(ljinvar)
            IF (nJInvar .gt. 0) avg_Jinvar = avg_Jinvar/nJinvar
            WHERE (.not.ljinvar) J_Invariant = avg_Jinvar

            IF (ABS(sigma(1,ieps)) .ge. bigno) CYCLE
            DO n = 1, nu                                 !!Read in nu from file...
              num = num+1
              index_array(num) = ivar_jinvar
              wegt(num) = max_variation * sqrt_nsurf * avg_Jinvar
     1                   * sigma(1,ieps)
              chisq_target(num) = avg_Jinvar
              chisq_match(num) = J_Invariant(n)
            END DO
         END DO

         CLOSE (iunit)                                   !!Opened in CALL to load_physics....
         DEALLOCATE (J_Invariant, ljinvar)
      ELSE
         INQUIRE(file=TRIM(version), exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            IF (lscreen) PRINT *,
     1          'xj_invariant file not found in ' // TRIM(home_dir)
            lj_invariant = .false.
         ELSE
            IF (nu .eq. 0) RETURN
            DO ieps = 1, NumJinvariant
               IF (ABS(sigma(1,ieps)) .ge. bigno) CYCLE
               DO n = 1, nu
                  num = num+1
                  IF (nopt .eq. -2) chisq_descript(num) = 
     1                              descript(ivar_jinvar)
               END DO
            END DO
         END IF
      END IF

      END SUBROUTINE chisq_jinvar
