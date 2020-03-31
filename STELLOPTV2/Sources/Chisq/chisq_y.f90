!-----------------------------------------------------------------------
!     Subroutine:    chisq_y
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          10/28/2018
!     Description:   This is a chisq subroutine for a test variable
!                    called y.  It does not depend on any equilibrium
!                    test and is simply for testing the optimizers on
!                    a known problem.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_y(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'TEST_Y ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  YVAL'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         targets(mtargets) = target_y
         sigmas(mtargets)  = sigma_y
         vals(mtargets)    = yval
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,vals(mtargets)
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_y
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_y
