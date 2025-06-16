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
      USE stellopt_vars, ONLY: nw_coil, nh_coil
      USE biotsavart, ONLY: coil_group
      
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
         ntotal_coils = 0
         ncoilgroups = SIZE(coil_group)
         DO i = 1, ncoilgroups
            !DO j = 1, coil_group(i)%ncoil
            ! We only do the first coil in each group
            DO j = 1, nw_coil*nh_coil
               nc = SIZE(coil_group(i)%coils(j)%xnod,2)
               nc1 = nc - 1
               hs  = 1.0D+00/nc1
               ALLOCATE(xc(nc),yc(nc),zc(nc))
               ALLOCATE(xcp(nc),ycp(nc),zcp(nc))
               ALLOCATE(xcpp(nc),ycpp(nc),zcpp(nc))
               ALLOCATE(curve(nc))
               xc = coil_group(i)%coils(j)%xnod(1,:)
               yc = coil_group(i)%coils(j)%xnod(2,:)
               zc = coil_group(i)%coils(j)%xnod(3,:)
               xcp(1:nc1) = xc(2:nc) - xc(1:nc1)
               ycp(1:nc1) = yc(2:nc) - yc(1:nc1)
               zcp(1:nc1) = zc(2:nc) - zc(1:nc1)
               ! Note that x(1)=x(nc) so we need
               xcp(nc) = xcp(1)
               ycp(nc) = ycp(1)
               zcp(nc) = zcp(1)
               xcp = xcp * hs; ycp = ycp * hs; zcp = zcp * hs;
               xcpp(1:nc1) = xcp(2:nc) - xcp(1:nc1)
               ycpp(1:nc1) = ycp(2:nc) - ycp(1:nc1)
               zcpp(1:nc1) = zcp(2:nc) - zcp(1:nc1)
               xcpp(nc) = xcpp(1)
               ycpp(nc) = ycpp(1)
               zcpp(nc) = zcpp(1)
               xcpp = xcpp * hs; ycpp = ycpp * hs; zcpp = zcpp * hs;
               curve = SQRT((zcpp*ycp-ycpp*zcp)**2 &
                     + (xcpp*zcp-zcpp*xcp)**2 &
                     + (ycpp*xcp-xcpp*ycp)**2) &
                     / (xcp*xcp+ycp*ycp+zcp*zcp)**(3.0/2.0)
               curve_min = curve_min + MINVAL(curve)
               curve_max = curve_max + MAXVAL(curve)
               curve_mean = curve_mean + SUM(curve)/nc
               ntotal_coils = ntotal_coils + 1
               DEALLOCATE(xc,yc,zc)
               DEALLOCATE(xcp,ycp,zcp)
               DEALLOCATE(xcpp,ycpp,zcpp,curve)
            END DO
         END DO
         curve_min = curve_min/ntotal_coils
         curve_max = curve_max/ntotal_coils
         curve_mean = curve_mean/ntotal_coils
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
