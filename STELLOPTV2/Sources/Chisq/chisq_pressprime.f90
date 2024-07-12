!-----------------------------------------------------------------------
!     Subroutine:    chisq_pressprime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/27/20
!     Description:   Calculate difference between equilibrium pressure
!                    gradient and a target value
!-----------------------------------------------------------------------
      SUBROUTINE chisq_pressprime(target,sigma,niter,iflag)
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
      REAL(rprec) :: p_val, u_val, s0, pp_val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'PRESSPRIME ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         s0 = 0.0; ier = 0; p_val = 0.0; u_val = 0.0; pp_val =0.0
         IF (ANY(s_pressprime > 0)) lreset_s = .false.
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ! GET s if necessary
            IF (lreset_s) THEN
               CALL get_equil_s(r_pressprime(ik),phi_pressprime(ik),z_pressprime(ik),s_pressprime(ik),ier,u_val)
               IF (ier < 0) THEN; iflag = ier; RETURN; END IF
            END IF
            IF (s_pressprime(ik) <= 1.0 .and. s_pressprime(ik) >= 0.0) THEN
               CALL get_equil_p(s_press(ik),p_val,ier,pp_val)
               IF (ier < 0) THEN; iflag = ier; RETURN; END IF
            ELSE
               pp_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = pp_val
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_pressprime(ik),phi_pressprime(ik),z_pressprime(ik),s_pressprime(ik),target(ik),sigma(ik),pp_val
         END DO
         IF (lreset_s) s_pressprime(:) = -1.0
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_pressprime
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_pressprime
