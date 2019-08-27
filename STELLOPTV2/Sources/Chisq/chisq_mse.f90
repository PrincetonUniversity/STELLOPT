!-----------------------------------------------------------------------
!     Subroutine:    chisq_mse
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate the difference between equilibrium and
!                    experimentally measured MSE signals.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_mse(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars
      USE equil_utils
      USE biotsavart
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nprof)
      REAL(rprec), INTENT(in)    ::  sigma(nprof)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        dex         Length of ne_aux_s array
!        ne_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      LOGICAL ::  lcoil_open = .false.
      INTEGER ::  ik, j, ier, dex, iunit
      REAL(rprec) :: Bx, By, Bz, Br, Bphi, Er, Ez, gradphi, mse_val
      REAL(rprec) :: Bx_vac, By_vac, Bz_vac
      REAL(rprec) :: current, current_first, num_test
      REAL(rprec) :: xvec(3),bvec(3)
      REAL(rprec) :: a_temp(6)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'MSE ',ik,9
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  ER  EZ  MSE_VAL  '
      dex = MINLOC(phi_aux_s(2:),DIM=1)
      IF (niter >= 0) THEN
         IF (ANY(s_mse > 0)) lreset_s = .false.
         ! Load Coil
         IF (ANY(lextcur_opt) .and. ANY(lmse_extcur)) THEN
            CALL cleanup_biotsavart                  !Do this in case the coil has already been read
            INQUIRE(FILE=TRIM(magdiag_coil),NUMBER=iunit,OPENED=lcoil_open)
            IF (lcoil_open) CLOSE(iunit)
            CALL parse_coils_file(TRIM(magdiag_coil))
            DO ik = 1, SIZE(extcur)
               DO j = 1, coil_group(ik) % ncoil
                  current = coil_group(ik) % coils(j) % current
                  IF (j .eq. 1) current_first = current
                  IF (current_first .ne. zero) coil_group(ik) % coils(j) % current = (current/current_first)*extcur(ik)
               END DO
            END DO
         END IF
         ! Now calculate values
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ! For now we disable the phi contribution to the MSE signal
            ! Get E
            !IF (dex > 1) THEN
            !   CALL get_equil_E(r_mse(ik),phi_mse(ik),z_mse(ik),Er,Ez,ier)
            !ELSE
               Er = 0.0
               Ez = 0.0
            !END IF
            ! Get S (for plotting purposes only)
            ier = 0
            !PRINT *,ik,r_mse(ik),phi_mse(ik),z_mse(ik)
            CALL get_equil_s(r_mse(ik),phi_mse(ik),z_mse(ik),s_mse(ik),ier)
            ! Get B
            Br = 0.0; Bphi = 0.0; Bz = 0.0
            CALL get_equil_Bcyl(r_mse(ik),phi_mse(ik),z_mse(ik),Br,Bphi,Bz,ier)
            ! Calc VMEC MSE response
            IF (ier == 0) THEN
               num_test  = (  a2_mse(ik) * Bphi + a3_mse(ik) * Br &
                            + a4_mse(ik) * Bz   + a6_mse(ik) * Ez &
                            + a7_mse(ik) * Er)
               IF (num_test == 0) THEN
                  mse_val = 0.0
               ELSE
                  mse_val = ATAN((  a1_mse(ik) * Bz   + a5_mse(ik) * Er)/num_test)
               END IF
            ELSE
               mse_val = 0.0
            END IF
            ! Calc vac_mse
            IF (ANY(lextcur_opt) .and. ANY(lmse_extcur)) THEN
               Br   = 0.0
               Bphi = 0.0
               Bz   = 0.0
               DO j = 1, SIZE(coil_group)
                  IF (.not.lmse_extcur(j)) CYCLE ! Skip if we don't turn the current on
                  CALL bfield(r_mse(ik),phi_mse(ik),z_mse(ik),bvec(1),bvec(2),bvec(3),IG=j)
                  Br   = Br + bvec(1)
                  Bphi = Bphi + bvec(2)
                  Bz   = Bz + bvec(3)
               END DO
               num_test  = (  a2_mse(ik) * Bphi + a3_mse(ik) * Br + a4_mse(ik) * Bz )
               IF (num_test == 0.0) THEN
                  vac_mse(ik) = 0.0
               ELSE
                  vac_mse(ik) = ATAN((a1_mse(ik) * Bz)/ num_test)
               END IF
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)+vac_mse(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = mse_val
            IF (iflag == 1) WRITE(iunit_out,'(9ES22.12E3)') r_mse(ik),phi_mse(ik),z_mse(ik),s_mse(ik),targets(mtargets),sigma(ik),Er,Ez,vals(mtargets)
         END DO
         IF (lreset_s) s_mse(:) = -1.0
         IF (ANY(lextcur_opt)) CALL cleanup_biotsavart
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               lneed_output = .true.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_mse
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_mse
