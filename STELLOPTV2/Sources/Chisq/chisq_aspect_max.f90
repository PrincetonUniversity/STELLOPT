!-----------------------------------------------------------------------
!     Subroutine:    chisq_aspect_max
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Constrain aspect ratio to be below target.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_aspect_max(target,sigma,niter,iflag)
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
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'ASPECT_MAX ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET SIGMA ASPECT_MAX CHI'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         ! Impose a limit on the calculated aspect ratio, below aspec_max
         ! it should match target but above it should match grow
         vals(mtargets)    = target + 1 + TANH((aspect-target)/(width_aspect_max))
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,aspect,vals(mtargets)
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_aspect_max
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_aspect_max
