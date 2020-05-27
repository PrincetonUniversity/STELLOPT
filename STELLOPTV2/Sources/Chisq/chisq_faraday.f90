!-----------------------------------------------------------------------
!     Subroutine:    chisq_faraday
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/09/2012
!     Description:   Calculate difference measured and equilibrium
!                    faraday rotation measurements
!-----------------------------------------------------------------------
      SUBROUTINE chisq_faraday(target,sigma,niter,iflag)
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
      REAL(rprec) :: faraday_val, u_val, r_try
      REAL(rprec) :: x0(3), x1(3)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno_ne)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'FARADAY',ik,9
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  R0  PHI0  Z0  R1  PHI1  Z1'
      IF (niter >= 0) THEN
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno_ne) CYCLE
            x0(1)=r0_faraday(ik); x1(1)=r1_faraday(ik)
            x0(2)=phi0_faraday(ik); x1(2)=phi1_faraday(ik)
            x0(3)=z0_faraday(ik); x1(3)=z1_faraday(ik)
            CALL line_int(fcn_faraday,x0,x1,faraday_val)
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = faraday_val
            IF (iflag == 1) WRITE(iunit_out,'(9ES22.12E3)') target(ik),sigma(ik),faraday_val,x0(1),x0(2),x0(3),x1(1),x1(2),x1(3)
         END DO
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno_ne) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_faraday
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_faraday
