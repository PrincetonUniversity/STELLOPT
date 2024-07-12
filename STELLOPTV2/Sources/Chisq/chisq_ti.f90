!-----------------------------------------------------------------------
!     Subroutine:    chisq_ti
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between equilbrium ion
!                    temperature and target ion temperature
!-----------------------------------------------------------------------
      SUBROUTINE chisq_ti(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: ti_aux_s, ti_aux_f
      USE equil_utils
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nprof)
      REAL(rprec), INTENT(in)    ::  sigma(nprof)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        ti_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  ik, ier, dex
      REAL(rprec) :: ti_val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'TI ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  TI'
      dex = MINLOC(ti_aux_s(2:),DIM=1)
      IF (niter >= 0) THEN
         IF (ANY(s_ti > 0)) lreset_s = .false.
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ! GET s if necessary
            IF (lreset_s) THEN
               ier = 0
               CALL get_equil_s(r_ti(ik),phi_ti(ik),z_ti(ik),s_ti(ik),ier)
            END IF
            IF (s_ti(ik) >= 0.0) THEN
               ier = 0
               CALL get_equil_ti(s_ti(ik),TRIM(ti_type),ti_val,ier)
            ELSE
               ti_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = ti_val
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_ti(ik),phi_ti(ik),z_ti(ik),s_ti(ik),target(ik),sigma(ik),ti_val
         END DO
         IF (lreset_s) s_ti(:) = -1.
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_ti
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_ti
