!-----------------------------------------------------------------------
!     Subroutine:    chisq_wp
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between equilibrium stored
!                    energy and target stored energy.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_wp(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: wp
      
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
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'WP ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         targets(mtargets) = target_wp
         sigmas(mtargets)  = sigma_wp
         vals(mtargets)     = wp
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,wp
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_wp
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_wp
