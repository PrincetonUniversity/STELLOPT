!-----------------------------------------------------------------------
!     Subroutine:    chisq_gs2_ptsm3d
!     Authors:       A. Bader (abader@engr.wisc.edu)
!     Date:          07/2018
!     Description:   This subroutine targets the energy transfer
!           between stable and unstable modes using the ptsm3d module
!           Geometric quantities are calculated in the gs2 framework
!-----------------------------------------------------------------------
      SUBROUTINE chisq_gs2_ptsm3d(target,sigma,niter,iflag)
      !     Libraries
      !-----------------------------------------------------------------------

      USE stellopt_runtime
      USE stellopt_targets  

      !-----------------------------------------------------------------------
      !     Input/Output Variables
      !
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag


      !----------------BEGIN SUBROUTINE --------------

      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'P2     ',1,3
      IF (niter >= 0) THEN

        write (*,*) 'GS2 PTSM3D TEST FUNCTION'
        mtargets = mtargets + 1
        targets(mtargets) = target
        sigmas(mtargets) = sigma
        vals(mtargets) = 1.0
        IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,1.0

      ELSE
        IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_gs2_ptsm3d
        END IF
      END IF
      RETURN
      END SUBROUTINE chisq_gs2_ptsm3d
