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
      REAL(rprec) :: val, jboot
      REAL(rprec), PARAMETER :: smooth_fac=0.7
      LOGICAL, DIMENSION(nsd) :: lbooz_sav
      REAL(rprec), DIMENSION(21) :: coefs
      REAL(rprec), ALLOCATABLE :: sfarr(:),sarr(:),farr(:)
      REAL(rprec), DIMENSION(5) :: f_out
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

      DO
         ! Run VMEC
         iflag = 0
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen_local)
         iflag = ier_paraexe

         ! Load Equilibrium
         CALL stellopt_load_equil(lscreen_local,iflag)
         lbooz(1:ns) = .TRUE.
         lbooz(1)    = .FALSE.
         CALL stellopt_paraexe('booz_xform',proc_string,lscreen_local); iflag = ier_paraexe

         ! Run BOOTSTRAP
         CALL stellopt_paraexe('bootsj',proc_string,lscreen_local); iflag = ier_paraexe
         dibs = dibs * 1D6 ! Get in A
         aibs = aibs * 1D6 ! Get in A

         ! Smooth Current Profile
         !ALLOCATE(sfarr(irup))
         !sfarr = 0
         !CALL smoothg(dibs,irup,smooth_fac,sfarr)
         !dibs = sfarr
         !DEALLOCATE(sfarr)

         ! Calculated bootstrap fraction
         DO ik = 1, irup
            iflag = 0
            CALL get_equil_bootj(rhoar(ik),jboot,iflag)
            val = MAX(val, ABS(dibs(ik)-jboot))
         END DO
         PRINT *,val

         !Test for exit
         IF (val < 1.0E3) EXIT

         ! Setup fitting arrays
         ALLOCATE(sarr(irup+1), farr(irup+1))
         sarr(1) = 0; sarr(irup+2) = 1;
         sarr(2:irup+1) = rhoar
         farr(2:irup+1) = dibs
         farr(1) = farr(2) + (sarr(1)-sarr(2))*(farr(3)-farr(2))/(sarr(3)-sarr(2))
         farr(irup+2) = farr(irup+1) + (sarr(irup+2)-sarr(irup+1))*(farr(irup+1)-farr(irup))/(sarr(irup+1)-sarr(irup))
         coefs = bootj_aux_f(1:21)

         !DO ik = 1, irup
         !   PRINT *,dibs(ik),farr(ik+1)
         !END DO

         ! Print Bootstrap to Log
         ALLOCATE(sfarr(irup+1))
         DO ik = 1, irup+1
            CALL get_equil_bootj(sarr(ik),sfarr(ik),iflag)
         END DO
         IF (lfirst_pass) WRITE(ibootlog,'(512(1X,E20.10))')  (sarr(ik),ik=1,irup+1)
         !WRITE(ibootlog,'(512(1X,E20.10))')  0.0,(dibs(ik),  ik=1,irup)  ! Bootstrap
         WRITE(ibootlog,'(512(1X,E20.10))')  (farr(ik),  ik=1,irup+1)  ! Bootstrap
         WRITE(ibootlog,'(512(1X,E20.10))')  (sfarr(ik), ik=1,irup+1) ! Current Fit
         CALL FLUSH(ibootlog)
         DEALLOCATE(sfarr)

         ! Print to screen
         IF (lscreen) THEN ! Only on first pass through
            DO ik = 1, 5
               jboot = s_out(ik) ! First argument declared in/out
               CALL get_equil_bootj(jboot,f_out(ik),iflag)
            END DO
            IF (lfirst_pass) WRITE(6,'(5f10.2)') s_out
            WRITE(6,'(5f10.2,E20.10)') f_out/1E3,aibs(irup)/1E6
            CALL FLUSH(6)
         END IF

         !  Fit profile to bootstrap
         CALL fit_profile(bootj_type,irup+2,sarr,farr,11,coefs(1:11))
         DEALLOCATE(sarr,farr)

         ! Update Coefficients
         bootj_aux_f(1:21) = bootj_aux_f(1:21) - (bootj_aux_f(1:21) - coefs)*0.25
         curtor_vmec = curtor_vmec - (curtor_vmec - aibs(irup))*0.25
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
