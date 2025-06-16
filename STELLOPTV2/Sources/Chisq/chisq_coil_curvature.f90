!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_curvature
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   Calculate the mean local coil curvature.
!                    https://en.wikipedia.org/wiki/Curvature
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_curvature(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE spline_coils_mod, ONLY: get_coil_curvature
      
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
      REAL(rprec) :: curve_mean, curve_min, curve_max, hs
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc,yc,zc,xcp,ycp,zcp,xcpp,ycpp,zcpp,curve
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_CURVATURE ',1,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MEAN  MIN  MAX'
      IF (niter >= 0) THEN
         curve_min = 0.0; curve_max = 0.0; curve_mean = 0.0
         CALL get_coil_curvature(curve_mean,curve_max,curve_min)
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = curve_mean
         IF (iflag == 1) WRITE(iunit_out,'(5ES22.12E3)') target,sigma,curve_mean,curve_min,curve_max
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_coil_curvature
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_curvature
