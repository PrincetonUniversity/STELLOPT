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
               ALLOCATE(xcppp(nc),ycppp(nc),zcppp(nc))
               ALLOCATE(torsion(nc))
               xc = coil_group(i)%coils(j)%xnod(1,:)
               yc = coil_group(i)%coils(j)%xnod(2,:)
               zc = coil_group(i)%coils(j)%xnod(3,:)
               xcp(1:nc1) = xc(2:nc) - xc(1:nc1)
               ycp(1:nc1) = yc(2:nc) - yc(1:nc1)
               zcp(1:nc1) = zc(2:nc) - zc(1:nc1)
               xcp = xcp * hs; ycp = ycp * hs; zcp = zcp * hs;
               xcp(nc) = xc(1) - xc(nc)
               ycp(nc) = yc(1) - yc(nc)
               zcp(nc) = zc(1) - zc(nc)
               xcpp(1:nc1) = xcp(2:nc) - xcp(1:nc1)
               ycpp(1:nc1) = ycp(2:nc) - ycp(1:nc1)
               zcpp(1:nc1) = zcp(2:nc) - zcp(1:nc1)
               xcpp(nc) = xcp(1) - xcp(nc)
               ycpp(nc) = ycp(1) - ycp(nc)
               zcpp(nc) = zcp(1) - zcp(nc)
               xcpp = xcpp * hs; ycpp = ycpp * hs; zcpp = zcpp * hs;
               xcppp(1:nc1) = xcpp(2:nc) - xcpp(1:nc1)
               ycppp(1:nc1) = ycpp(2:nc) - ycpp(1:nc1)
               zcppp(1:nc1) = zcpp(2:nc) - zcpp(1:nc1)
               xcppp(nc) = xcpp(1) - xcpp(nc)
               ycppp(nc) = ycpp(1) - ycpp(nc)
               zcppp(nc) = zcpp(1) - zcpp(nc)
               xcppp = xcppp * hs; ycppp = ycppp * hs; zcppp = zcppp * hs;
               torsion = ((zcpp*ycp-ycpp*zcp)*xcppp &
                       +  (xcpp*zcp-zcpp*xcp)*ycppp &
                       +  (ycpp*xcp-xcpp*ycp)*zcppp) &
                       / ((zcpp*ycp-ycpp*zcp)**2 &
                       +  (xcpp*zcp-zcpp*xcp)**2 &
                       +  (ycpp*xcp-xcpp*ycp)**2)
               torsion_min = torsion_min + MINVAL(torsion)
               torsion_max = torsion_max + MAXVAL(torsion)
               torsion_mean = torsion_mean + SUM(torsion)/nc
               ntotal_coils = ntotal_coils + 1
               DEALLOCATE(xc,yc,zc)
               DEALLOCATE(xcp,ycp,zcp)
               DEALLOCATE(xcpp,ycpp,zcpp)
               DEALLOCATE(xcppp,ycppp,zcppp,torsion)
            END DO
         END DO
         torsion_min = torsion_min/ntotal_coils
         torsion_max = torsion_max/ntotal_coils
         torsion_mean = torsion_mean/ntotal_coils
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
