!-----------------------------------------------------------------------
!     Subroutine:    chisq_magwell
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          07/05/2017
!     Description:   This subroutine calculates the magnetic well
!                    parameter as outline at:
!                    https://fusion.gat.com/pubs-ext/ComPlasmaPhys/A22135.pdf
!                    Here W > 0 implies stability
!-----------------------------------------------------------------------
      SUBROUTINE chisq_magwell(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE equil_utils
      USE stellopt_targets
      
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
      INTEGER     :: ik, ier
      REAL(rprec) :: modb, Bsqav, dBsqav, p, pp, W, Vp,temp1, temp2, &
                     rhosqav, V, Bav
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'MAGWELL ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MAGWELL  <B**2>  P  dP/ds'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            CALL get_equil_Bav(rho(ik),Bav,Bsqav,ier,dBsqav)
            CALL get_equil_p(rho(ik),p,ier,pp)
            CALL get_equil_volume(rho(ik),V,ier,Vp)
            W = V*(2*mu0*pp/Vp+dBsqav)/Bsqav
            ! Output value
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = W
            IF (iflag == 1) WRITE(iunit_out,'(6ES22.12E3,2X,I3.3)') target(ik),sigma(ik),W, Bsqav, p, pp, ik
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_magwell
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_magwell
