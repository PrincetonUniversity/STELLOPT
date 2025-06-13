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
      ik = nu_bnormal*nv_bnormal
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I8.8))') 'BNORMAL ',ik,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  BNORMAL'
      IF (niter >= 0) THEN
         DO ik = 1, nu_bnormal*nv_bnormal
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)     = bnormal_total(ik)
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,bnormal_total(ik)
         END DO
      ELSE
         IF (sigma < bigno) THEN
            lneed_bnormal = .TRUE.
            DO ik = 1, nu_bnormal*nv_bnormal
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_bnormal
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_bnormal
