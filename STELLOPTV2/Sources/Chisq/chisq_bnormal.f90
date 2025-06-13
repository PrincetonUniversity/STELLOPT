!-----------------------------------------------------------------------
!     Subroutine:    chisq_bnormal
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   Calculate the difference between the normal field
!                    on the equilibrium surface and zero.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_bnormal(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: bnormal_total
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(in)    ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'BNORMAL_TOTAL ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  BNORMAL_TOTAL'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = SUM(bnormal_total)
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,SUM(bnormal_total)
      ELSE
         IF (sigma < bigno) THEN
            lneed_bnormal = .TRUE.
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_bnormal
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_bnormal
