!-----------------------------------------------------------------------
!     Subroutine:    chisq_b0
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/17/2017
!     Description:   Calculate difference between equilibrium magnetic
!                    axis field and the target magnetic axis field.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_b0(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: Baxis
      
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
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'B0 ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  B0'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = Baxis
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,Baxis
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_b0
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_b0
