!-----------------------------------------------------------------------
!     Subroutine:    chisq_te
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between equilbrium electron
!                    temperature and target electron temperature
!-----------------------------------------------------------------------
      SUBROUTINE chisq_te(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: ndatafmax, lte_f_opt, te_aux_s, te_aux_f
      USE equil_utils
      USE EZspline_obj
      USE EZspline
      
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
!        te_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  ik, ier, dex, flag, u2
      REAL(rprec) :: s, xu, xv, R1, Z1, R2, Z2, dist, temp
      REAL(rprec) :: te_val, targ_temp, sig_temp
      REAL(rprec) :: d1(nu_max)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'TE ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  VAL'
      dex = MINLOC(te_aux_s(2:),DIM=1)
      IF (niter >= 0) THEN
         IF (ANY(s_te > 0)) lreset_s = .false.
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ! GET s if necessary
            IF (lreset_s) THEN
               ier = 0
               CALL get_equil_s(r_te(ik),phi_te(ik),z_te(ik),s_te(ik),ier)
            END IF
            IF (s_te(ik) <= 1.0 .and. s_te(ik) >= 0.0 .and. ier == 0) THEN
               ier = 0
               CALL get_equil_te(s_te(ik),TRIM(te_type),te_val,ier)
               targ_temp = target(ik)
               sig_temp  = sigma(ik)
            ELSE
               te_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            !targets(mtargets) = targ_temp
            !sigmas(mtargets)  = sig_temp
            vals(mtargets)    = te_val
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_te(ik),phi_te(ik),z_te(ik),s_te(ik),target(ik),sigma(ik),te_val
         END DO
         IF (lreset_s) s_te(:) = -1.0
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_te
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_te
