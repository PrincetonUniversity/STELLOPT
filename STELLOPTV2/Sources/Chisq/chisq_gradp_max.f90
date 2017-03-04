!-----------------------------------------------------------------------
!     Subroutine:    chisq_gradp_max
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Constrain pressure gradient to be below target.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_gradp_max(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: aspect
      
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
      INTEGER :: ik 
      REAL(rprec) :: gradp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'GRADP_MAX ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  GRADP  CHI'
      IF (niter >= 0) THEN
         gradp = 0.0
         
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         ! Impose a limit on the calculated aspect ratio, below aspec_max
         ! it should match target but above it should match grow
         vals(mtargets)    = target + 1 + TANH((gradp-target)/(width_gradp_max))
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,aspect,vals(mtargets)
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_gradp_max
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_gradp_max
