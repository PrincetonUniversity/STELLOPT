!-----------------------------------------------------------------------
!     Subroutine:    stellopt_vboot
!     Authors:       S. Lazerson (lazerson@pppl.gov) and Matt Landreman
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
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      LOGICAL :: lscreen_local, lfirst_pass
      INTEGER :: ik, nc, ibootlog, ier, j, radius_index, unit_out=12, Nradii
      REAL(rprec) :: vboot_convergence_factor, jboot, val_last
      REAL(rprec), PARAMETER :: smooth_fac=0.05
      LOGICAL, DIMENSION(nsd) :: lbooz_sav
      REAL(rprec), DIMENSION(21) :: coefs
      REAL(rprec), ALLOCATABLE :: sfarr(:),sarr(:),farr(:)
      REAL(rprec), DIMENSION(5) :: f_out
      REAL(rprec), PARAMETER :: smooth_frac=0.75
      REAL(rprec), DIMENSION(5), PARAMETER :: s_out=(/0.0,0.25,0.50,0.75,1.0/)
      INTEGER :: vboot_iteration
      CHARACTER(len=4) :: iteration_string
      REAL(rprec) :: temp1, temp2, ds_fine, scale_factor, curtor_beam, curtor_bootstrap
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sfincs_s_with_0
      INTEGER, PARAMETER :: Ns_fine = 1000
      REAL(rprec), DIMENSION(Ns_fine) :: s_fine_full, s_fine_half, J_dot_B_flux_surface_average_fine, d_p_d_s_fine_half
      REAL(rprec), DIMENSION(Ns_fine) :: B_squared_flux_surface_average_fine_half, B_squared_flux_surface_average_fine_full
      REAL(rprec), DIMENSION(Ns_fine) :: integrating_factor_full, integrating_factor_half, d_p_d_s_fine, integrand, isigng_I_F_full, isigng_I_full
      REAL(rprec), DIMENSION(Ns_fine) :: integrating_factor_half_approximate, pressure_fine_half
      REAL(rprec), DIMENSION(21) :: J_dot_B_flux_surface_average_fit, B_squared_flux_surface_average_fit, sfincs_AC_fit
      REAL(rprec), DIMENSION(Ns_fine) :: sfincs_ac_half, sfincs_ac_low_beta_limit, AC_profile_fine, AC_fit_results
      CHARACTER(LEN=32) :: B_squared_flux_surface_average_profile_type
      LOGICAL :: exit_after_next_vmec_run

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

      ! Here we use the same convention as in VMEC: half-mesh quantities use arrays with the same size as full-mesh quantities,
      ! but the first array element is 0.
      s_fine_full = [( (j-1.0d+0)/(Ns_fine-1), j=1, Ns_fine )]
      ds_fine = s_fine_full(2) - s_fine_full(1)
      s_fine_half(1) = 0
      s_fine_half(2:Ns_fine) = s_fine_full(1:(Ns_fine-1)) + ds_fine/2

      ! Evaluate the total beam current, which we will need later for computing the new CURTOR for vmec:
      DO radius_index = 1,Ns_fine
         CALL eval_profile(s_fine_full(radius_index), beamj_type, AC_profile_fine(radius_index), beamj_aux_s, beamj_aux_f, ier)
      END DO
      curtor_beam = (SUM(AC_profile_fine) - 0.5*AC_profile_fine(1) - 0.5*AC_profile_fine(Ns_fine)) * ds_fine ! Integration using the trapezoid rule.

      Nradii = MINLOC(sfincs_s(2:),DIM=1) ! Number of radii at which sfincs is run, if bootcalc_type='sfincs'.

      ! Open log file
      ibootlog = 12
      CALL safe_open(ibootlog, iflag, 'boot_fit.'//trim(proc_string), 'replace','formatted')

      ! Beginning of the main iteration:
      val_last = 0
      vboot_iteration = -1
      ier = 0
      exit_after_next_vmec_run = .false.
      AC_profile_fine = 0
      DO
         vboot_iteration = vboot_iteration + 1

         ! Run VMEC
         iflag = 0
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen_local)
         iflag = ier_paraexe
         IF (iflag .ne.0) THEN
            PRINT *,"WARNING: paravmec returned with an error flag: iflag =",iflag
            RETURN
         END IF
         ! Save vmec wout files from each iteration:
         WRITE (iteration_string,fmt="(i4.4)") vboot_iteration
         CALL system('cp wout_'//trim(proc_string)//".nc wout_"//trim(proc_string)//"_vboot"//trim(iteration_string)//".nc")
         !CALL stellopt_paraexe('paravmec_write',trim(proc_string)//"_vboot"//trim(iteration_string),.true.) 

         ! Load Equilibrium
         CALL stellopt_load_equil(lscreen_local,iflag)

         ! Don't do anything if pressure is zero
         IF (wp <= 0 .or. beta<=0) EXIT

         ! Run the bootstrap current code
         CALL tolower(bootcalc_type)
         SELECT CASE (TRIM(bootcalc_type))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Beginning of code specific to bootsj.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CASE ('bootsj')

            IF (irup > Ns_fine) STOP "irup must be <= Ns_fine."

            IF (vboot_iteration>0) THEN
               ! Log results in the boot_log file:
               WRITE(ibootlog,'(a,512(1X,E20.10))')  "AC_profile_fine: ",(AC_profile_fine(ik), ik=1,irup)
               CALL FLUSH(ibootlog)
            END IF

            IF (exit_after_next_vmec_run) EXIT

            lbooz(1:ns) = .TRUE.
            lbooz(1)    = .FALSE.
            CALL stellopt_paraexe('booz_xform',proc_string,lscreen_local); iflag = ier_paraexe
            IF (iflag .ne.0) RETURN

            CALL stellopt_paraexe('bootsj',proc_string,lscreen_local); iflag = ier_paraexe

            IF (iflag .ne.0) RETURN
            dibs = dibs * 1D6 ! Convert megaAmperes to Amperes.
            aibs = aibs * 1D6 ! Convert megaAmperes to Amperes.
            
            ! Log results in the boot_log file:
            IF (vboot_iteration == 0) THEN
               WRITE(ibootlog,'(a,512(1X,E20.10))')  "s=rhoar: ", rhoar
               ! Evaluate the initial stellopt bootstrap AC profile on the bootsj s (a.k.a. rhoar) grid.
               AC_profile_fine = 0
               DO radius_index = 1,irup
                  CALL eval_profile(rhoar(radius_index), bootj_type, AC_profile_fine(radius_index), bootj_aux_s, bootj_aux_f, ier)
               END DO
               WRITE(ibootlog,'(a,512(1X,E20.10))')  "AC_profile_fine: ",(AC_profile_fine(ik), ik=1,irup)
            END IF
            WRITE(ibootlog,'(a,512(1X,E20.10))')  "dibs: ",dibs
            CALL FLUSH(ibootlog)

            ! Fit the resulting AC profile
            CALL fit_profile(bootj_type, irup, rhoar, dibs, 21, bootj_aux_f)
            ! Evaluate the fit:
            DO radius_index = 1, irup
               CALL eval_profile(rhoar(radius_index), bootj_type, AC_fit_results(radius_index), bootj_aux_s, bootj_aux_f, ier)
            END DO

            !Test if the vboot iterations have converged, by comparing "L_1 norm of the difference between the last two AC profiles" to "L_1 norm of the latest AC profile":
            vboot_convergence_factor = SUM(ABS(AC_fit_results(1:irup) - AC_profile_fine(1:irup))) / SUM(ABS(AC_fit_results(1:irup)))
            AC_profile_fine = AC_fit_results
                        
            ! In this next line, perhaps the first and last point should be treated differently for greater accuracy?
            ds_fine = rhoar(2)-rhoar(1)
            curtor_bootstrap = SUM(AC_fit_results(1:irup)) * ds_fine
            curtor_vmec = curtor_bootstrap + curtor_beam

            IF (vboot_convergence_factor < vboot_tolerance) THEN
               WRITE(6,"(a,i4,a,es10.3,a,es10.3,a)") "Vboot iteration",vboot_iteration,": ctor=",curtor_vmec,", vboot convergence factor=",vboot_convergence_factor,". Tolerance achieved."
               exit_after_next_vmec_run = .true. ! VMEC is cheap, so always finish the vboot iteration with 1 last vmec run.
            ELSE
               WRITE(6,"(a,i4,a,es10.3,a,es10.3)")   "Vboot iteration",vboot_iteration,": ctor=",curtor_vmec,", vboot convergence factor=",vboot_convergence_factor
            END IF

            WRITE(ibootlog,"(4(a,es22.15))") "curtor_bootstrap = ",curtor_bootstrap," curtor_beam = ",curtor_beam," curtor_total = ",curtor_vmec," vboot convergence factor = ",vboot_convergence_factor
            CALL FLUSH(ibootlog)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of code specific to bootsj. Beginning of code specific to sfincs.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CASE ('sfincs')

!DEC$ IF DEFINED (SFINCS)
            IF (vboot_iteration==0) THEN
               ! Evaluate the initial AC profile on the fine grid.
               DO radius_index = 2,Ns_fine
                  CALL eval_profile(s_fine_half(radius_index), bootj_type, AC_profile_fine(radius_index), bootj_aux_s, bootj_aux_f, ier)
               END DO
               AC_profile_fine(1) = 0
            END IF

            ! Log results in the boot_log file:
            IF (lfirst_pass) WRITE(ibootlog,'(a,512(1X,E20.10))')  "s: ", (s_fine_half(ik),ik=2,Ns_fine,10)
            ! Note that when bootcalc_type='sfincs', the order of the lines in the boot_log file is different from the choice made previously for bootsj!
            WRITE(ibootlog,'(a,512(1X,E20.10))')  "AC_profile_fine: ",(AC_profile_fine(ik), ik=2,Ns_fine,10)
            CALL FLUSH(ibootlog)

            IF (exit_after_next_vmec_run) EXIT

            CALL stellopt_paraexe('sfincs',proc_string,lscreen_local); iflag = ier_paraexe

            ! Save sfincs files from each iteration
            CALL system('cp -r sfincs_'//trim(proc_string)//" sfincs_"//trim(proc_string)//"_vboot"//trim(iteration_string))
            CALL copy_txtfile('sfincs_results_before_profile_fitting.'//TRIM(proc_string), 'sfincs_results_before_profile_fitting.'//TRIM(proc_string)//"_vboot"//trim(iteration_string))

            ! The next section of code implements the formulae in the note computing_vmec_AC_profile_from_a_bootstrap_current_code.tex,
            ! for converting <j dot B> to VMEC's curtor and AC profile.

            ! Fit a function to <j dot B>. We will use the same functional type as bootj_type so that
            ! bootj_aux_f can provide a good initial guess for <j dot B>.
            ALLOCATE(sfincs_s_with_0(Nradii+1))
            sfincs_s_with_0 = 0
            sfincs_s_with_0(2:) = sfincs_s(1:Nradii)
            ! For an initial guess at the fit coefficients for <j dot B>, use the previous dI/ds profile:
            ! AC = -d I / d s = -2 pi (d psi / d s) <j dot B> / <B^2> + (small d p / d s term)
            J_dot_B_flux_surface_average_fit = bootj_aux_f(1:21)
            scale_factor = -SUM(sfincs_B_squared_flux_surface_average(2:(Nradii+1)))/(Nradii*(-phiedge))
            CALL scale_profile(bootj_type, J_dot_B_flux_surface_average_fit, scale_factor)
            CALL fit_profile(bootj_type, Nradii+1, sfincs_s_with_0, sfincs_J_dot_B_flux_surface_average, 21, &
                 J_dot_B_flux_surface_average_fit)
            DEALLOCATE(sfincs_s_with_0)

            ! Fit a profile to <B^2>:
            B_squared_flux_surface_average_fit = 0
            B_squared_flux_surface_average_profile_type = 'power_series'
            B_squared_flux_surface_average_fit(1) = SUM(sfincs_B_squared_flux_surface_average(2:(Nradii+1))) / Nradii ! Initial guess for constant term in the polynomial fit.
            B_squared_flux_surface_average_fit(2:3) = B_squared_flux_surface_average_fit(1) / 10 ! Set to a small nonzero value, so 'fit_profile' detects the correct polynomial order.
            CALL fit_profile(B_squared_flux_surface_average_profile_type, Nradii, sfincs_s, sfincs_B_squared_flux_surface_average(2:(Nradii+1)), 21, &
                 B_squared_flux_surface_average_fit)
            DO radius_index = 1,Ns_fine
               CALL eval_profile(s_fine_half(radius_index), B_squared_flux_surface_average_profile_type, &
                    B_squared_flux_surface_average_fine_half(radius_index), bootj_aux_s, B_squared_flux_surface_average_fit, ier)
               CALL eval_profile(s_fine_full(radius_index), B_squared_flux_surface_average_profile_type, &
                    B_squared_flux_surface_average_fine_full(radius_index), bootj_aux_s, B_squared_flux_surface_average_fit, ier)
            END DO

            ! Evaluate dp/ds:
            DO radius_index = 2,Ns_fine
               CALL get_equil_p(s_fine_half(radius_index), pressure_fine_half(radius_index), ier, d_p_d_s_fine_half(radius_index))
            END DO
            CALL get_equil_p(s_fine_full(1), temp1, ier)
            d_p_d_s_fine_half(1) = 0
            ! Sanity test: compute the integrating factor in the approximation that <B^2> is constant:
            integrating_factor_half_approximate = exp(mu0 * (pressure_fine_half-temp1) / (SUM(sfincs_B_squared_flux_surface_average(2:(Nradii+1))) / Nradii))

            ! Compute the integrating factor, eq (19):
            integrand = 0
            DO radius_index = 2,Ns_fine
               integrand(radius_index) = integrand(radius_index-1) + ds_fine * d_p_d_s_fine_half(radius_index) / B_squared_flux_surface_average_fine_half(radius_index)
            END DO
            integrating_factor_full = exp(mu0 * integrand)
            integrating_factor_half(1) = 0
            integrating_factor_half(2:Ns_fine) = 0.5d+0 * (integrating_factor_full(1:(Ns_fine-1)) + integrating_factor_full(2:Ns_fine))

            ! Form the integrand to eq (20) on the half mesh:
            DO radius_index = 2,Ns_fine
               CALL eval_profile(s_fine_half(radius_index), bootj_type, j_dot_B_flux_surface_average_fine(radius_index), bootj_aux_s, J_dot_B_flux_surface_average_fit, ier) ! <J dot B>
            END DO
            j_dot_B_flux_surface_average_fine(1) = 0
            B_squared_flux_surface_average_fine_half(1) = 0
            ! Here comes the integrand of eq (20). There are 2 factors of isigng that multiply together to give +1:
            !   * The fact that phi = 2 pi psi isigng.
            !   * The fact that AC = isigng * dI/ds
            sfincs_AC_low_beta_limit = phiedge * j_dot_B_flux_surface_average_fine / B_squared_flux_surface_average_fine_half
            integrand = sfincs_AC_low_beta_limit * integrating_factor_half
            integrand(1) = 0

            ! Integrate to get isigng*I(s)*F(s) on the full mesh:
            isigng_I_F_full = 0
            DO j=2,Ns_fine
               isigng_I_F_full(j) = isigng_I_F_full(j-1) + integrand(j) * ds_fine
            END DO

            isigng_I_full = isigng_I_F_full / integrating_factor_full
            ! Differentiate isigng * I(s) to get AC(s):
            sfincs_AC_half = 0
            DO j=2,Ns_fine
               sfincs_AC_half(j) = (isigng_I_full(j) - isigng_I_full(j-1)) / ds_fine
            END DO

            ! Log results in the boot_log file:
            WRITE(ibootlog,'(a,512(1X,E20.10))')  "sfincs_AC_half: ",(sfincs_AC_half(ik),  ik=2,Ns_fine,10)  ! Output of the bootstrap current code
            CALL FLUSH(ibootlog)

            ! Fit the resulting AC profile
            CALL fit_profile(bootj_type, Ns_fine-1, s_fine_half(2:Ns_fine), sfincs_AC_half(2:Ns_fine), 21, bootj_aux_f)
            ! Evaluate the fit:
            DO radius_index = 2,Ns_fine
               CALL eval_profile(s_fine_half(radius_index), bootj_type, AC_fit_results(radius_index), bootj_aux_s, bootj_aux_f, ier)
            END DO
            AC_fit_results(1) = 0

            !Test if the vboot iterations have converged, by comparing "L_1 norm of the difference between the last two AC profiles" to "L_1 norm of the latest AC profile":
            vboot_convergence_factor = SUM(ABS(AC_fit_results - AC_profile_fine)) / SUM(ABS(AC_fit_results))
            AC_profile_fine = AC_fit_results

            !IF (myworkid == master) THEN
            ier=0
            CALL safe_open(unit_out,ier,'sfincs_results_after_profile_fitting.'//TRIM(proc_string),'REPLACE','formatted')
            WRITE (unit_out,"(a)") &
                 "s (normalized toroidal flux), " // &
                 "Toroidal current derivative AC (Amperes), " // &
                 "Low beta approximation of toroidal current derivative AC (Amperes), " // &
                 "Fit to toroidal current derivative AC (Amperes), " // &
                 "Fit to <j dot B> (Tesla Amperes / meters^2), " // &
                 "Fit to <B^2> (Tesla^2), " // &
                 "d pressure / d s (Pascals), " // &
                 "Integrating factor (eq (19), dimensionless), " // &
                 "Constant-<B^2> approximation for the integrating factor (dimensionless), " // &
                 "Integrand of eq (20)"
            DO radius_index = 2, Ns_fine
               WRITE (unit_out,"(10(es24.15))") s_fine_half(radius_index), sfincs_AC_half(radius_index), sfincs_AC_low_beta_limit(radius_index), &
                    AC_fit_results(radius_index), j_dot_B_flux_surface_average_fine(radius_index), B_squared_flux_surface_average_fine_half(radius_index), &
                    d_p_d_s_fine_half(radius_index), integrating_factor_half(radius_index), integrating_factor_half_approximate(radius_index), integrand(radius_index)
            END DO
            CLOSE (unit_out)
            !END IF

            CALL copy_txtfile('sfincs_results_after_profile_fitting.'//TRIM(proc_string), 'sfincs_results_after_profile_fitting.'//TRIM(proc_string)//"_vboot"//trim(iteration_string))

            ! In this next line, perhaps the first and last point should be treated differently for greater accuracy?
            curtor_bootstrap = SUM(AC_fit_results) * ds_fine
            curtor_vmec = curtor_bootstrap + curtor_beam

            IF (vboot_convergence_factor < vboot_tolerance) THEN
               WRITE(6,"(a,i4,a,es10.3,a,es10.3,a)") "Vboot iteration",vboot_iteration,": ctor=",curtor_vmec,", vboot convergence factor=",vboot_convergence_factor,". Tolerance achieved."
               exit_after_next_vmec_run = .true. ! VMEC is computationally very cheap compared to sfincs, so always finish the vboot iteration with 1 last vmec run.
            ELSE
               WRITE(6,"(a,i4,a,es10.3,a,es10.3)")   "Vboot iteration",vboot_iteration,": ctor=",curtor_vmec,", vboot convergence factor=",vboot_convergence_factor
            END IF

            WRITE(ibootlog,"(4(a,es22.15))") "curtor_bootstrap = ",curtor_bootstrap," curtor_beam = ",curtor_beam," curtor_total = ",curtor_vmec," vboot convergence factor = ",vboot_convergence_factor
            CALL FLUSH(ibootlog)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of code specific to SFINCS.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEC$ ELSE
            STOP "Error! STELLOPT was compiled without SFINCS support, but bootcalc_type was set to sfincs."
!DEC$ ENDIF
         CASE DEFAULT
            PRINT *,"Error! Invalid bootcalc_type:",bootcalc_type
            STOP
         END SELECT

         ! Setup for next pass
         !lscreen_local = .FALSE.
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
