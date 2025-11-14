!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_baxis
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/13/2025
!     Description:   Calculate the angle between the vacuum field on
!                    axis and the VMEC mangetic axis. Note that this
!                    only works for beta=0 VMEC equilibria.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_baxis(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: baxis_total
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
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
      ik = nv_bnormal
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8.8))') 'BAXIS ',ik,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  BAXIS'
      IF (niter >= 0) THEN
         DO ik = 1, nv_bnormal
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)     = baxis_total(ik)
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,baxis_total(ik)
         END DO
      ELSE
         IF (sigma < bigno) THEN
            lneed_bnormal = .TRUE.
            DO ik = 1, nv_bnormal
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_baxis
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_baxis
