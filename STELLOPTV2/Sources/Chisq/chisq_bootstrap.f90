!-----------------------------------------------------------------------
!     Subroutine:    chisq_bootstrap
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between <J*B> in equilibria
!                    vs <J*B> as calculated via a bootstrap code.
!                    Note target_bootstrap should be set to zero so
!                    that chisq is the difference between VMEC and
!                    bootstrap calcualtion.
!                    Note the bootsj namelist needs teti set to <= 0.0
!                    to use the STELLOPT Ti profile.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_bootstrap(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: equil_type, beamj_aux_s, bootj_aux_s,&
                               beamj_aux_f, bootj_aux_f
      USE equil_vals, ONLY: curtor
      USE equil_utils, ONLY: get_equil_jdotb, get_equil_nustar, &
                             eval_prof_spline, get_equil_bootj, &
                             get_equil_beamj
      USE parambs, ONLY: ajBbs, bsnorm, rhoar, l_boot_all, aibs
      USE safe_open_mod, ONLY: safe_open
      USE bootsj_input
      USE mpi_params, ONLY: myid, master
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in) ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: i, j, ik, iunit, ier, dex1, dex2
      REAL(rprec) :: avg_jdotb, nustar, facnu, s, j_beam, j_boot, &
                     boot_jdotb, fac_boot, beam_jdotb, sig_temp
      CHARACTER(len = 256)   :: temp_str
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'BOOTSTRAP ',ik,10
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  RHO  AVG_JDOTB  BEAM_JDOTB  BOOT_JDOTB  AJBBS  FACNU  BSNORM'
      IF (niter >= 0) THEN
         dex1 = MINLOC(beamj_aux_s(2:),DIM=1)
         dex2 = MINLOC(bootj_aux_s(2:),DIM=1)
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               CALL get_equil_nustar(rhoar(ik),nustar,ier)
               IF (ier == -1) nustar=0
               facnu = 1./(1.+nustar)
               ier   = 0
               CALL get_equil_jdotb(rhoar(ik),avg_jdotb,ier)
               ! Scale avg_jdotb if net toroidal current is present
               beam_jdotb = 0.0
               boot_jdotb = 0.0
               IF ((dex1>4) .and. (dex2>4)) THEN
                  s = rhoar(ik)
                  CALL get_equil_beamj(s,j_beam,ier)
                  CALL get_equil_bootj(s,j_boot,ier)
                  !CALL eval_prof_spline(dex1,beamj_aux_s(1:dex1),beamj_aux_f(1:dex1),s,j_beam,ier)
                  !CALL eval_prof_spline(dex2,bootj_aux_s(1:dex2),bootj_aux_f(1:dex2),s,j_boot,ier)
                  beam_jdotb = avg_jdotb * j_beam/(j_beam+j_boot) * SIGN(1.0_rprec,curtor)
                  boot_jdotb = avg_jdotb * j_boot/(j_beam+j_boot) * SIGN(1.0_rprec,curtor)
               ELSE IF ((ABS(curtor) > 0) .and. (l_boot_all)) THEN
                  boot_jdotb = avg_jdotb * aibs(ik)/curtor
               END IF
               IF (ABS(ajBbs(ik)) < 1.0E-30) ajBbs(ik) = 0.0
               IF (ABS(bsnorm(ik)) < 1.0E-30) bsnorm(ik) = 0.0
               vals(mtargets)    = ABS(boot_jdotb-ajBbs(ik)*facnu)
               sig_temp          = ABS(ajBbs(ik)*facnu)
               !IF (sig_temp > 1.0E3) sigmas(mtargets) = sig_temp*10.
               IF (sig_temp > 1.0E4) sigmas(mtargets) = sig_temp*10.
               IF (sig_temp > 1.0E5) sigmas(mtargets) = sig_temp*1.0
               IF (sig_temp > 1.0E6) sigmas(mtargets) = sig_temp*0.2
               IF (iflag == 1) WRITE(iunit_out,'(10ES22.12E3)') target(ik),sigmas(mtargets),vals(mtargets),rhoar(ik),avg_jdotb,beam_jdotb,boot_jdotb,ajBbs(ik),facnu,bsnorm(ik)
            END IF
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_bootstrap
            END IF
         END DO
         iunit=12
         CALL safe_open(iunit,iflag,'input.'//TRIM(id_string),'old','formatted')
         IF (iflag < 0) RETURN
         CALL read_namelist (iunit, iflag, 'bootin')
         !IF (iflag < 0 .and. niter == -2 .and. myid == master) THEN
         IF (iflag < 0 .and. myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
            WRITE(6,*) '  BOOTIN Namelist not found     '
            WRITE(6,*) ' '
            ik=0
            CALL write_bootsj_input(6,ik)
            WRITE(6,*) ' '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         END IF
         CLOSE(iunit)
         IF (iflag < 0) RETURN
      END IF
      ! Save the output files
!      IF ((iflag == 1) .and. (lkeep_mins)) THEN
!         WRITE(temp_str,'(i5.5)') niter
!         temp_str = '.'//TRIM(temp_str)
!         SELECT CASE(TRIM(equil_type))
!            CASE('vmec2000')
!               CALL RENAME('answers_plot.'//TRIM(proc_string_old),'answers_plot.'//TRIM(id_string)//TRIM(temp_str),ier)
!               CALL RENAME('answers.'//TRIM(proc_string_old),'answers.'//TRIM(id_string)//TRIM(temp_str),ier)
!               CALL RENAME('jBbs.'//TRIM(proc_string_old),'jBbs.'//TRIM(id_string)//TRIM(temp_str),ier)
!         END SELECT
!      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_bootstrap
