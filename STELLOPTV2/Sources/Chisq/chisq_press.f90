!-----------------------------------------------------------------------
!     Subroutine:    chisq_press
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between equilbrium total
!                    plasma pressure and target total plams pressure
!-----------------------------------------------------------------------
      SUBROUTINE chisq_press(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      
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
      INTEGER ::  ik, ier, dex
      REAL(rprec) :: p_val, u_val, p_norm, s0
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'PRESS ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         s0 = 0.0; ier = 0; p_norm = 0.0; p_val = 0.0; u_val = 0.0
         CALL get_equil_p(s0,p_norm,ier)
         IF (ANY(s_ne > 0)) lreset_s = .false.
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ! GET s if necessary
            IF (lreset_s) THEN
               CALL get_equil_s(r_press(ik),phi_press(ik),z_press(ik),s_press(ik),ier,u_val)
               IF (ier < 0) THEN; iflag = ier; RETURN; END IF
            END IF
            IF (s_press(ik) <= 1.0 .and. s_press(ik) >= 0.0) THEN
               CALL get_equil_p(s_press(ik),p_val,ier)
               IF (ier < 0) THEN; iflag = ier; RETURN; END IF
            ELSE
               p_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)/norm_press
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = p_val/p_norm
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_press(ik),phi_press(ik),z_press(ik),s_press(ik),target(ik),sigma(ik),p_val/p_norm
         END DO
         IF (lreset_s) s_press(:) = -1.0
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_press
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_press
