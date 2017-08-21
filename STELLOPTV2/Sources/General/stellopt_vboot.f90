!-----------------------------------------------------------------------
!     Subroutine:    stellopt_vboot
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/24/2017
!     Description:   This subroutine calculates a (par)VMEC equilibrium
!                    with self-consistent bootstrap current profile, by
!                    iterating between VMEC and a bootstrap current
!                    code, either BOOTSJ or SFINCS. This subroutine is
!                    used when EQUIL_OPTION='VBOOT' in the OPTIMUM
!                    namelist.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_vboot(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_vars
      USE equil_utils
      USE safe_open_mod
      USE read_wout_mod, ONLY: ns
      use read_wout_mod, only: jdotb_vmec => jdotb
      USE vmec_input, ONLY: curtor_vmec => curtor
      USE parambs, ONLY: rhoar, dibs, irup, aibs, d_rho
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL :: lscreen_local, lfirst_pass
      INTEGER :: ik, nc, ibootlog
      REAL(rprec) :: val, jboot, val_last
      REAL(rprec), PARAMETER :: smooth_fac=0.05
      LOGICAL, DIMENSION(nsd) :: lbooz_sav
      REAL(rprec), DIMENSION(21) :: coefs
      REAL(rprec), ALLOCATABLE :: sfarr(:),sarr(:),farr(:)
      REAL(rprec), DIMENSION(5) :: f_out
      REAL(rprec), PARAMETER :: smooth_frac=0.75
      REAL(rprec), DIMENSION(5), PARAMETER :: s_out=(/0.0,0.25,0.50,0.75,1.0/)
      INTEGER :: vboot_iteration ! MJL
      CHARACTER(len=4) :: iteration_string ! MJL
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      lscreen_local = .FALSE.
      lfirst_pass = .TRUE.
      IF (lscreen) lscreen_local = .TRUE.
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  VBOOT CALCULATION  -------------------------'

      ! Handle boozer flags
      lbooz_sav = lbooz

      ! Open log file
      ibootlog = 12
      CALL safe_open(ibootlog, iflag, 'boot_fit.'//trim(proc_string), 'replace','formatted')
      val_last = 0
      vboot_iteration = -1 !MJL
      DO
         vboot_iteration = vboot_iteration + 1 ! MJL
         ! Run VMEC
         iflag = 0
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen_local)
         iflag = ier_paraexe
         IF (iflag .ne.0) RETURN
         !PRINT *,'-1-',iflag
         ! Begin MJL additions
         print *,"Here comes proc_string:",proc_string
         write (iteration_string,fmt="(i4.4)") vboot_iteration
         !call copy_txtfile('wout_'//trim(proc_string)//".nc", 'wout_'//trim(proc_string)//"_vboot"//trim(iteration_string)//".nc")
         !if (myworkid==master) call system('cp wout_'//trim(proc_string)//".nc wout_"//trim(proc_string)//"_vboot"//trim(iteration_string)//".nc")
         call system('cp wout_'//trim(proc_string)//".nc wout_"//trim(proc_string)//"_vboot"//trim(iteration_string)//".nc")
         ! End MJL additions

         ! Load Equilibrium
         CALL stellopt_load_equil(lscreen_local,iflag)

         print *,"Here comes jdotb from vmec:" ! MJL
         print *,jdotb_vmec ! MJL

         ! Don't do anything if pressure is zero
         IF (wp <= 0 .or. beta<=0) EXIT

         ! BOOTSJ requires Boozer coordinates, but SFINCS does not.
         !IF (TRIM(bootcalc_type) == 'bootsj') THEN
         IF (.true.) THEN
            ! Call BOOZER Transformation
            lbooz(1:ns) = .TRUE.
            lbooz(1)    = .FALSE.
            CALL stellopt_paraexe('booz_xform',proc_string,lscreen_local); iflag = ier_paraexe
            IF (iflag .ne.0) RETURN
         END IF

         ! Run the bootstrap current code
         CALL tolower(bootcalc_type)
         SELECT CASE (TRIM(bootcalc_type))
         CASE ('bootsj')
            CALL stellopt_paraexe('bootsj',proc_string,lscreen_local); iflag = ier_paraexe
         CASE ('sfincs')
            CALL stellopt_paraexe('sfincs',proc_string,lscreen_local); iflag = ier_paraexe
         CASE DEFAULT
            PRINT *,"Error! Invalid bootcalc_type:",bootcalc_type
            STOP
         END SELECT
         IF (iflag .ne.0) RETURN
         dibs = dibs * 1D6 ! Convert megaAmperes to Amperes.
         aibs = aibs * 1D6 ! Convert megaAmperes to Amperes.
         !PRINT *,'-3-',iflag

!!$         print *,"Here comes dibs before smoothing:" ! MJL
!!$         print *,dibs !MJL

!!$         ! Smooth Current Profile
!!$         ALLOCATE(sfarr(irup))
!!$         sfarr = 0
!!$         CALL smoothg(dibs,irup,smooth_fac,sfarr)
!!$         dibs = sfarr
!!$         DEALLOCATE(sfarr)

!!$         print *,"Here comes dibs after smoothing:" !MJL
!!$         print *,dibs  !MJL

         ! Calculate a scalar difference between the stellopt j_bootstrap profile and the
         ! new profile from the bootstrap current code.
         val = 0
         DO ik = 2, irup-1
            iflag = 0
            CALL get_equil_bootj(rhoar(ik),jboot,iflag)
            val = val + ABS(dibs(ik)-jboot)
         END DO
         val = val/(irup-2)

         ! Set up fitting arrays
         ALLOCATE(sarr(irup+2), farr(irup+2))
         sarr(1) = 0; sarr(irup+2) = 1;
         sarr(2:irup+1) = rhoar
         farr(2:irup+1) = dibs
         farr(1) = farr(2) + (sarr(1)-sarr(2))*(farr(3)-farr(2))/(sarr(3)-sarr(2))
         farr(irup+2) = farr(irup+1) + (sarr(irup+2)-sarr(irup+1))*(farr(irup+1)-farr(irup))/(sarr(irup+1)-sarr(irup))

         ! Print Bootstrap to Log
         ALLOCATE(sfarr(irup+2))
         DO ik = 1, irup+2
            CALL get_equil_bootj(sarr(ik),sfarr(ik),iflag)
         END DO
         IF (lfirst_pass) WRITE(ibootlog,'(512(1X,E20.10))')  (sarr(ik),ik=1,irup+2)
         !WRITE(ibootlog,'(512(1X,E20.10))')  0.0,(dibs(ik),  ik=1,irup)  ! Bootstrap
         WRITE(ibootlog,'(512(1X,E20.10))')  (farr(ik),  ik=1,irup+2)  ! Bootstrap
         WRITE(ibootlog,'(512(1X,E20.10))')  (sfarr(ik), ik=1,irup+2) ! Current Fit
         CALL FLUSH(ibootlog)
         DEALLOCATE(sfarr)

         ! Print to screen
         IF (lscreen) THEN ! Only on first pass through
            DO ik = 1, 5
               jboot = s_out(ik) ! First argument declared in/out
               CALL get_equil_bootj(jboot,f_out(ik),iflag)
            END DO
            IF (lfirst_pass) WRITE(6,'(A,5f10.2,A20,A20)') 'S:   ',s_out,'Total [MA]','Error'
            WRITE(6,'(5X,5f10.2,2E20.10)') f_out/1E3,aibs(irup)/1E6, val
            CALL FLUSH(6)
         END IF

         !Test for exit
         IF (ABS(val-val_last)/ABS(val) < 0.01) EXIT
         val_last = val
         !IF (val < 0.5E3) EXIT

         ! MJL: this next step does not make sense because farr has dimensions (Amps) and is ~ 10^6, whereas jboot (returned by get_equil_bootj)
         ! on the 1st iteration is ~1 since it has not yet been scaled by curtor. And if you wanted to blend the old and new current profiles, there should be a (1-smooth_frac) factor too.
!!$         ! Modify the bootstrap current for the fit
!!$         DO ik = 1, irup+2
!!$            CALL get_equil_bootj(sarr(ik),jboot,iflag)
!!$            !farr(ik) = jboot + (farr(ik) - jboot)*smooth_frac
!!$            farr(ik) = jboot*(1-smooth_frac) + (farr(ik) - jboot)*smooth_frac
!!$         END DO
         

!!$         if (lfirst_pass) then
!!$            print *,"First pass, so let's call fit_profile an extra time."
!!$            print *,"bootj_aux_f before fitting:",bootj_aux_f ! MJL
!!$            nc = 21
!!$            coefs(1:nc) = bootj_aux_f(1:nc)
!!$            ik = irup+2
!!$            CALL fit_profile(bootj_type,ik,sarr,farr,nc,coefs(1:nc))
!!$            print *,"bootj_aux_f after fitting:",coefs(1:nc) ! MJL
!!$         end if

         print *,"bootj_aux_f before fitting:",bootj_aux_f ! MJL
         !  Fit profile to bootstrap
         nc = 21
         coefs(1:nc) = bootj_aux_f(1:nc)
         ik = irup+2
         CALL fit_profile(bootj_type,ik,sarr,farr,nc,coefs(1:nc))
         DEALLOCATE(sarr,farr)

         ! Update Coefficients
         bootj_aux_f(1:21) = coefs(1:21)
         print *,"bootj_aux_f after fitting:",bootj_aux_f ! MJL

         ! Update curtor
         ! This guarantees that curtor is consistent with the bootstrap profile.
         ! Note that in the future we can adjust for the fraction here.
         val = 0
         DO ik = 1, irup
            CALL get_equil_bootj(rhoar(ik),jboot,iflag)
            val = val + jboot*d_rho(ik)
         END DO
         curtor_vmec = val 
         !PRINT *,val

         ! Setup for next pass
         lscreen_local = .FALSE.
         lfirst_pass = .FALSE.
      END DO
      
      ! Finish up
      CLOSE(ibootlog)
      lbooz = lbooz_sav
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  VBOOT CALCULATION DONE  ----------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_vboot
