!-----------------------------------------------------------------------
!     Subroutine:    chisq_mercier
!     Authors:       J. Schmitt (jcschmitt@auburn.edu)
!     Date:          2018
!     Description:   Calculate difference between desired and actual
!                    mercier stability growth rate.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_mercier(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: mercier_criterion
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(in)    ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ik 
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I5.5))') 'MERCIER_CRITERION ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  CHI  MERCIER_CRITERION  #'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               IF (.not. ALLOCATED(mercier_criterion)) PRINT *, "<----ERROR::mercier_criterion not allocated"
               vals(mtargets)    = -1 * MIN(target(ik),mercier_criterion(ik)) ! Assume positive to be stable
               IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3,2X,I3.3)') target(ik),sigma(ik),vals(mtargets),mercier_criterion(ik),ik
            END IF
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_mercier_criterion
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_mercier

