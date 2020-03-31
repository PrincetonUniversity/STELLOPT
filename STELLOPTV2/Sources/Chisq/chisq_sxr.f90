!-----------------------------------------------------------------------
!     Subroutine:    chisq_sxr
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/19/2013
!     Description:   Calculate difference measured and equilibrium
!                    soft X-Ray signals
!-----------------------------------------------------------------------
      SUBROUTINE chisq_sxr(target,sigma,niter,iflag)
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
      INTEGER ::  ik, ier, dex
      REAL(rprec) :: u_val, r_try, targ_scale, val_scale
      REAL(rprec) :: sxr_val(nprof), targ_val(nprof)
      REAL(rprec) :: x0(3), x1(3)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'SXR',ik,9
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  R0  PHI0  Z0  R1  PHI1  Z1'
      val_scale = 1.0
      targ_scale = 1.0
      IF (niter >= 0) THEN
         sxr_val = 0.0
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            x0(1)=r0_sxr(ik); x1(1)=r1_sxr(ik)
            x0(2)=phi0_sxr(ik); x1(2)=phi1_sxr(ik)
            x0(3)=z0_sxr(ik); x1(3)=z1_sxr(ik)
            CALL line_int(fcn_sxr,x0,x1,sxr_val(ik))
         END DO
         val_scale  = MAXVAL(sxr_val)
         targ_scale = MAXVAL(target)
         IF (val_scale == 0) val_scale = 1.0
         IF (targ_scale == 0) targ_scale = 1.0
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            x0(1)=r0_sxr(ik); x1(1)=r1_sxr(ik)
            x0(2)=phi0_sxr(ik); x1(2)=phi1_sxr(ik)
            x0(3)=z0_sxr(ik); x1(3)=z1_sxr(ik)
            mtargets = mtargets + 1
            targets(mtargets) = target(ik) / targ_scale
            sigmas(mtargets)  = sigma(ik)  / targ_scale
            vals(mtargets)    = sxr_val(ik) / val_scale
            IF (iflag == 1) WRITE(iunit_out,'(9ES22.12E3)') target(ik)/targ_scale,sigma(ik)/targ_scale,sxr_val(ik)/val_scale,x0(1),x0(2),x0(3),x1(1),x1(2),x1(3)
         END DO
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_sxr
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_sxr
