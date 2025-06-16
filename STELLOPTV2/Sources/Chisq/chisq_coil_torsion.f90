!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_torsion
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   Calculate the mean local coil torision.
!                    https://en.wikipedia.org/wiki/Torsion_of_a_curve
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_torsion(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE spline_coils_mod, ONLY: get_coil_torsion
      
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
      INTEGER :: ncoilgroups, ntotal_coils,i,j,nc,nc1
      REAL(rprec) :: torsion_mean, torsion_min, torsion_max, hs
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc,yc,zc,xcp,ycp,zcp, &
         xcpp,ycpp,zcpp,xcppp,ycppp,zcppp,torsion
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_TORSION ',1,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MEAN  MIN  MAX'
      IF (niter >= 0) THEN
         torsion_min = 0.0; torsion_max = 0.0; torsion_mean = 0.0
         CALL get_coil_torsion(torsion_mean,torsion_max,torsion_min)
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = torsion_mean
         IF (iflag == 1) WRITE(iunit_out,'(5ES22.12E3)') target,sigma,torsion_mean,torsion_min,torsion_max
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_coil_torsion
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_torsion
