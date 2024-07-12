!-----------------------------------------------------------------------
!     Subroutine:    chisq_jtor
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/05/2013
!     Description:   Calculates the difference between Jcurv in an
!                    equilibrium and a desired Jcurv profile.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_jcurv(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils, ONLY: get_equil_jcurv, shat
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER     :: ik,ier
      REAL(rprec) :: local_jcurv
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'JTOR ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  S'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            ier = 0
            CALL get_equil_jcurv(shat(ik),local_jcurv,ier)
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = local_jcurv
            IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ik),sigma(ik),local_jcurv,shat(ik)
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_jcurv
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_jcurv
