      SUBROUTINE chisq_neo (sigma, ivar, num, ns_neo, ns_max,
     1                      nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: home_dir, nrad, exe_suffix
      USE optim_params, ONLY: lneo_opt
      USE system_mod, ONLY: system
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ivar, nopt, iflag, num, ns_max
      REAL(rprec), DIMENSION(nrad) :: sigma
      INTEGER, DIMENSION(ns_max) :: ns_neo
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_neo = 27
      REAL(rprec), PARAMETER :: zero = 0
      INTEGER :: iunit, istat, jcount, nsval, isurf, write_prog
      REAL(rprec) :: ripple_eff, dummy_neo
      CHARACTER(LEN=200) :: version
      LOGICAL :: ex
C-----------------------------------------------
      version = TRIM(home_dir) // 'xneo' // TRIM(exe_suffix)

      IF (nopt > 0) THEN

         iunit = unit_neo
         write_prog = 0

         IF (lscreen) THEN
           write_prog = 1
           WRITE(6,*)
     1           'Running NEO for effective ripple calculations'
         ENDIF

!        WRITE out input control filed neo_param.ext for surface control
         CALL write_neoparam (ns_max, ns_neo, extension,
     1        write_prog, istat)
         IF (istat .ne. 0) THEN
            iflag = -27
            RETURN
         ENDIF

         CALL load_physics_codes(version, " ", TRIM(extension),
     1         'neo_out', extension, iunit, iflag)
         IF (iflag .ne. 0) RETURN

         IF (lscreen) WRITE(6,*) '  ns    ripple diffusion (eps_h**1.5)'

         DO jcount = 1, ns_max

              nsval = ns_neo(jcount)
              READ(unit_neo,*,iostat=istat) isurf, ripple_eff,
     1            dummy_neo, dummy_neo, dummy_neo, dummy_neo
              IF (istat .ne. 0) THEN
                 iflag = -27
                 RETURN
              END IF
              IF( lscreen) WRITE (6, 1000) nsval, ripple_eff
              num = num + 1
              index_array(num)  = ivar
              chisq_target(num) = zero
              chisq_match(num)  = ripple_eff
              wegt(num) = sigma(nsval)

          ENDDO
          CLOSE(unit_neo)

      ELSE
         INQUIRE(file=version, exist=ex, iostat=istat)
         IF (istat.ne.0 .or. .not.ex) THEN
            PRINT *, 'xneo file not found in ' // TRIM(home_dir)
            lneo_opt = .false.
         ELSE
            DO jcount = 1, ns_max
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
            ENDDO
         END IF
      ENDIF

 1000 FORMAT(2x,i3,3x,1p,e12.3)

      END SUBROUTINE chisq_neo
