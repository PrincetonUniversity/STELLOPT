      SUBROUTINE chisq_ballooning (sigma, target_bal, num, ns_surf, 
     1     nrad, nopt, iflag, extension, command)
      USE stel_kinds
      USE chisq_mod
      USE optim_params, ONLY: bal_theta0, bal_zeta0,
     1    lpres_prof_opt, lballoon_opt, nini_theta, nini_zeta, nini_tot
      USE optim, ONLY: home_dir, version_opt, exe_suffix
      USE safe_open_mod
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt, ns_surf(*)
      INTEGER, INTENT(inout) :: num, iflag
      REAL(rprec) :: sigma(*), target_bal(*)
      CHARACTER(LEN=*) :: extension, command
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ij, jk, jrad, k, iunit, jcount, ns_ball, istat
      REAL(rprec) :: grate_ball(nrad), init_zeta, init_theta
      CHARACTER(len=LEN_TRIM(home_dir)+20) :: version
      LOGICAL :: ex
C-----------------------------------------------
!
!     BALLOONING STABILITY CRITERION (GRATE>0, unstable)
!
      version = TRIM(home_dir) // 'xcobravmec' // TRIM(exe_suffix)
      IF (nopt. gt. 0) THEN
         iunit = unit_cobra
         CALL safe_open(iunit, k, 'in_cobra.'//TRIM(extension),
     1          'replace', 'formatted')
         IF (k .ne. 0) THEN
            iflag = -12
            RETURN
         ENDIF
         WRITE(iunit, *,iostat=k) nini_zeta, nini_theta, TRIM(extension)
         WRITE(iunit, *,iostat=k) (bal_zeta0(jk), jk = 1, nini_zeta)
         WRITE(iunit, *,iostat=k) (bal_theta0(ij), ij = 1, nini_theta)
         WRITE(iunit, *,iostat=k) nrad
         WRITE(iunit, *,iostat=k) ns_surf(1:nrad)
         CLOSE (iunit)
         IF (k .ne. 0) THEN
            PRINT *,' chisq_balloon error for file: ', TRIM(extension)
            iflag = -12
            RETURN
         END IF

!        Read in COBRA results:
!        ns_ball    =  INDEX of VMEC radial surface at which growth rate is computed
!                      (should be same as ns_balloon_max=nrad, so no need to SAVE it ...)
!        grate_ball(j) = growth rate at initial position in ij, lk loop at surface "j"
!                      (note that INDEX-j in list <=> ns_ball in VMEC indexing)

         CALL load_physics_codes (version, 'in_cobra', command,
     1             'cobra_grate', extension, iunit, iflag)
         IF (iflag .ne. 0) RETURN

         DO ij = 1, nini_theta
            DO jk = 1, nini_zeta
               READ (iunit, *, iostat=k) init_zeta, init_theta, jcount
               IF (jcount.ne.nrad .and. myid.eq.master)                 !MPI
     1            PRINT *,' JCOUNT = ', jcount,' != NS_SURF_max (= ',
     2                    nrad,' READING COBRA_GRATE FILE'
               IF (k .eq. 0) READ (iunit, *, iostat=k) (ns_ball,
     1              grate_ball(jrad), jrad = 1, jcount)
               IF (k .ne. 0) THEN
                  IF (myid .eq. master)                                   !MPI
     1               WRITE (6, *) 'Error reading cobra_grate in ',
     2               ' load_targe for file: ', TRIM(extension)
                  iflag = -12
                  RETURN
               END IF

               DO jrad = 1, nrad
                  num = num + 1
                  index_array(num) = ivar_balloon
                  wegt(num) = 0.1_dp*sigma(jrad)
                  chisq_target(num) = target_bal(jrad)
                  chisq_match(num) = MAX(grate_ball(jrad),
     1                                   target_bal(jrad))
               END DO
            END DO
         END DO

         CLOSE(unit=iunit)                                             !Keep grate file

      ELSE
         INQUIRE(file=TRIM(version), exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            IF (myid .eq. master)                                       !MPI
     1          PRINT *,'xcobra file not found in ' // TRIM(home_dir)
            lpres_prof_opt = .false.
            lballoon_opt = .false.
         ELSE IF (nrad .ge. 1) THEN
               DO ij = 1, nini_theta
                  DO jk = 1, nini_zeta
                     DO jrad = 1, nrad
                        num = num + 1
                        IF (nopt .eq. -2) 
     1                     chisq_descript(num) = descript(ivar_balloon)
                     END DO
                  END DO
               END DO
         ENDIF
         iflag = 0
      ENDIF

      END SUBROUTINE chisq_ballooning
