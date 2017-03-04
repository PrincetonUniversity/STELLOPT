      SUBROUTINE chisq_bootsj (sigma, num, nsurf, nsurf_max,
     1    nrad, nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim_params, ONLY: lbootsj_opt, jboot, lseedcur
      USE optim, ONLY: home_dir, jdotb_opt, exe_suffix, remove, cat
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nrad, nopt, nsurf_max, nsurf(*)
      INTEGER :: num, iflag
      REAL(rprec) :: sigma(*)
      CHARACTER(LEN=*) :: extension
      LOGICAL lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_boot=26
      REAL(rprec), PARAMETER :: zero=0, p5=0.5_dp, one=1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jrad, iunit, k, nsval, nsin, iwrite, istat
      REAL(rprec) :: bsnorm, jdotb_norm, jdotb_targ, sjmax, frac_nustar
      CHARACTER(len=LEN_TRIM(home_dir)+20+LEN_TRIM(extension)) ::
     1         version
      LOGICAL :: ex
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: pseed_t, pfboot
C-----------------------------------------------

!
!     BOOTSTRAP CURRENT, <J dot B>, READ IN FROM BOOTSJ FILE
!
      version = TRIM(home_dir) // 'xbootsj' // TRIM(exe_suffix)
      IF (nopt .gt. 0) THEN
!
!     RUN BOOTSTRAP CODE TO CREATE OUTPUT FILE
!
         iunit = unit_boot

         CALL load_physics_codes (version, 'in_bootsj', comm1,
     1         'jBbs', extension, iunit, iflag)
         IF (iflag .ne. 0) RETURN

!     READ IN jdotb from xbootsj - produced data file, jBbs
         READ(iunit,'(a)', iostat=k)                                     !Discard header

         IF (nrad .gt. 0) THEN
            jdotb_norm = SUM(ABS(jdotb_opt(:nrad)))/nrad
         ELSE
            jdotb_norm = one
         END IF
         IF (jdotb_norm .eq. zero) jdotb_norm = one

         iwrite = iunit+1
         version = 'jboot.'//TRIM(extension)
         CALL safe_open (iwrite, k, version, 'replace', 'formatted')
         IF (k .ne. 0) THEN
            PRINT *,' Error opening file: ', TRIM(version)
            iflag = -13
            RETURN
         END IF
         WRITE (iwrite, "(2a,/,14x,a,42x,11('-'),/,61x,a)")
     1     '      S    <J*B> (xbootsj)   <J*B> (vmec)',
     2     '    FRAC(NU-STAR)   <J*B>(boot)', '*FRAC','<J*B>(tok)'

         DO jrad = 1, nsurf_max
            READ (iunit, *, iostat=k) nsin, jdotb_targ, bsnorm
            IF (k .ne. 0) THEN
               iflag = -13
               RETURN
            END IF
            nsval = nsurf(jrad)
            IF (nsin .ne. nsval) STOP 'nsin != nsval in chisq_bootsj'
            num = num + 1
            index_array(num) = ivar_bootsj
            wegt(num) = sigma(nsval)*jdotb_norm

            sjmax = REAL(nsval - 1,rprec)/(nrad - 1)
            frac_nustar = pfboot(sjmax)
            jdotb_targ = jdotb_targ * frac_nustar                ! mimic nu-star effect: pfboot = 1/(1+nu-star)

            IF (jboot.eq.1 .and. ABS(bsnorm).gt.zero)
     1         jdotb_targ = jdotb_targ / bsnorm                  ! Use axisymmetric value

            IF (lseedcur)
     1         jdotb_targ = jdotb_targ + pseed_t(sjmax)          ! Seed current (MCZ: Sept 98)

            chisq_target(num) = jdotb_targ                       ! Half radial grid
            chisq_match(num) = p5*(jdotb_opt(nsval)
     1                       +     jdotb_opt(nsval-1))           ! VMEC on full grid
            WRITE(iwrite,'(f8.2,1p,e15.2,2x,3e15.2,6x,3e15.2,
     1                     6x,3e15.2)') sjmax, chisq_target(num),
     2                     chisq_match(num), frac_nustar, bsnorm
         END DO

         CLOSE (iwrite)
         IF (lscreen) CALL system (cat // TRIM(version))
         CLOSE(iunit, status='delete')                           ! Remove jBbs.ext file
         version = "answers_plot" // "." // TRIM(extension)
         INQUIRE (file=version, exist=ex, iostat=istat)
         IF (ex) CALL system (remove // TRIM(version))

      ELSE
        INQUIRE(file=TRIM(version), exist=ex, iostat=istat)
        IF ((.not.ex) .or. (istat .ne. 0)) THEN
           lbootsj_opt = .false.
           IF (lscreen)
     1        PRINT *,'xbootsj file not found in ' // TRIM(home_dir)
           RETURN
        END IF
        IF (nsurf_max .ge. 1) THEN
           DO jrad = 1, nsurf_max
              num = num + 1
              IF (nopt .eq. -2) 
     1           chisq_descript(num) = descript(ivar_bootsj)
           END DO
        END IF
      ENDIF

      END SUBROUTINE chisq_bootsj
