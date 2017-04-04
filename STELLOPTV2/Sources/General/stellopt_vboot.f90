!-----------------------------------------------------------------------
!     Subroutine:    stellopt_vboot
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/24/2017
!     Description:   This subroutine calculates a self-consistent
!                    VMEC Bootstrap equilibrium
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
      DO
         ! Run VMEC
         iflag = 0
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen_local)
         iflag = ier_paraexe
         !PRINT *,'-1-',iflag

         ! Load Equilibrium
         CALL stellopt_load_equil(lscreen_local,iflag)
         lbooz(1:ns) = .TRUE.
         lbooz(1)    = .FALSE.
         CALL stellopt_paraexe('booz_xform',proc_string,lscreen_local); iflag = ier_paraexe
         !PRINT *,'-2-',iflag

         ! Run BOOTSTRAP
         CALL stellopt_paraexe('bootsj',proc_string,lscreen_local); iflag = ier_paraexe
         dibs = dibs * 1D6 ! Get in A
         aibs = aibs * 1D6 ! Get in A
         !PRINT *,'-3-',iflag

         ! Smooth Current Profile
         ALLOCATE(sfarr(irup))
         sfarr = 0
         CALL smoothg(dibs,irup,smooth_fac,sfarr)
         dibs = sfarr
         DEALLOCATE(sfarr)

         ! Calculated bootstrap error
         val = 0
         DO ik = 2, irup-1
            iflag = 0
            CALL get_equil_bootj(rhoar(ik),jboot,iflag)
            val = val + ABS(dibs(ik)-jboot)
         END DO
         val = val/(irup-2)

         ! Setup fitting arrays
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

         ! Modify the bootstrap current for the fit
         DO ik = 1, irup+2
            CALL get_equil_bootj(sarr(ik),jboot,iflag)
            farr(ik) = jboot + (farr(ik) - jboot)*smooth_frac
         END DO
         

         !  Fit profile to bootstrap
         nc = 21
         coefs(1:nc) = bootj_aux_f(1:nc)
         ik = irup+2
         CALL fit_profile(bootj_type,ik,sarr,farr,nc,coefs(1:nc))
         DEALLOCATE(sarr,farr)

         ! Update Coefficients
         bootj_aux_f(1:21) = coefs(1:21)

         ! Update curtor
         ! This garuntees that curtor is consistent with the bootstrap profile
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
