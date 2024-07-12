!-----------------------------------------------------------------------
!     Subroutine:    chisq_line_te
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference measured and equilibrium
!                    line integrated electron temperature.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_line_te(target,sigma,niter,iflag)
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
!        ne_val      Holds profile evaulation
!-----------------------------------------------------------------------
      INTEGER ::  ik
      REAL(rprec) :: te_val, te_length
      REAL(rprec) :: x0(3), x1(3)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'TELINE',ik,9
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  LINE_TE  R0  PHI0  Z0  R1  PHI1  Z1'
      IF (niter >= 0) THEN
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            x0(1)=r0_te_line(ik); x1(1)=r1_te_line(ik)
            x0(2)=phi0_te_line(ik); x1(2)=phi1_te_line(ik)
            x0(3)=z0_te_line(ik); x1(3)=z1_te_line(ik)
            te_val = 0.0
            CALL line_int(fcn_linete,x0,x1,te_val,LENGTH=te_length)
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = te_val/te_length
            IF (iflag == 1) WRITE(iunit_out,'(9ES22.12E3)') target(ik),sigma(ik),te_val/te_length,x0(1),x0(2),x0(3),x1(1),x1(2),x1(3)
         END DO
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_line_te
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_line_te
