!-----------------------------------------------------------------------
!     Subroutine:    chisq_coilcoil_distance
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   Calculate the minimum coil_coil distance
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coilcoil_distance(target,sigma,niter,iflag)
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
      INTEGER :: ncoilgroups, i1, i2, j1, j2, n1, k, nc1, nc2
      REAL(rprec) :: dist_min
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc1,yc1,zc1,xc2,yc2,zc2
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: x2d,y2d,z2d,d2d
      REAL(rprec), PARAMETER :: DCC_EXP = -5.0
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COILCOIL_DISTANCE ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MINDIST_INV  MINDIST'
      IF (niter >= 0) THEN
         dist_min = 1.0D+30
         ncoilgroups = SIZE(coil_group)
         !--------------------------------------------------------------
         !    Self-distance loop
         !    This loop compares a given coil group to all other coils
         !    in that same group. This catches the self-coil
         !    intersections.
         !--------------------------------------------------------------
         DO i1 = 1, ncoilgroups
            DO j1 = 1, nw_coil*nh_coil
               n1 = j1+nw_coil*nh_coil
               DO j2 = n1, coil_group(i1)%ncoil
                  nc1 = SIZE(coil_group(i1)%coils(j1)%xnod,2)
                  nc2 = SIZE(coil_group(i1)%coils(j2)%xnod,2)
                  ALLOCATE(x2d(nc1,nc2),y2d(nc1,nc2),z2d(nc1,nc2),d2d(nc1,nc2))
                  FORALL(k=1:nc1) x2d(k,:) = coil_group(i1)%coils(j1)%xnod(1,k)
                  FORALL(k=1:nc1) y2d(k,:) = coil_group(i1)%coils(j1)%xnod(2,k)
                  FORALL(k=1:nc1) z2d(k,:) = coil_group(i1)%coils(j1)%xnod(3,k)
                  FORALL(k=1:nc2) x2d(:,k) = x2d(:,k) - coil_group(i1)%coils(j2)%xnod(1,k)
                  FORALL(k=1:nc2) y2d(:,k) = y2d(:,k) - coil_group(i1)%coils(j2)%xnod(2,k)
                  FORALL(k=1:nc2) z2d(:,k) = z2d(:,k) - coil_group(i1)%coils(j2)%xnod(3,k)
                  d2d = x2d*x2d+y2d*y2d+z2d*z2d
                  !WHERE(d2d < 1.0E-6) d2d = 1.0E6
                  !dist_min = MIN(SQRT(MINVAL(d2d)),dist_min)
                  !WHERE(d2d > target) d2d = target*100.0
                  dist_min = MIN((MINVAL(d2d)),dist_min)
                  DEALLOCATE(x2d,y2d,z2d,d2d)
               END DO
            END DO
         END DO
         !--------------------------------------------------------------
         !    In this loop we compare coils in different groups. 
         !--------------------------------------------------------------
         DO i1 = 1, ncoilgroups
            n1 = i1+1
            DO i2 = n1,ncoilgroups
               DO j1 = 1, coil_group(i1)%ncoil
                  DO j2 = 1, coil_group(i2)%ncoil
                     nc1 = SIZE(coil_group(i1)%coils(j1)%xnod,2)
                     nc2 = SIZE(coil_group(i2)%coils(j2)%xnod,2)
                     ALLOCATE(x2d(nc1,nc2),y2d(nc1,nc2),z2d(nc1,nc2),d2d(nc1,nc2))
                     FORALL(k=1:nc1) x2d(k,:) = coil_group(i1)%coils(j1)%xnod(1,k)
                     FORALL(k=1:nc1) y2d(k,:) = coil_group(i1)%coils(j1)%xnod(2,k)
                     FORALL(k=1:nc1) z2d(k,:) = coil_group(i1)%coils(j1)%xnod(3,k)
                     FORALL(k=1:nc2) x2d(:,k) = x2d(:,k) - coil_group(i2)%coils(j2)%xnod(1,k)
                     FORALL(k=1:nc2) y2d(:,k) = y2d(:,k) - coil_group(i2)%coils(j2)%xnod(2,k)
                     FORALL(k=1:nc2) z2d(:,k) = z2d(:,k) - coil_group(i2)%coils(j2)%xnod(3,k)
                     d2d = x2d*x2d+y2d*y2d+z2d*z2d
                     !WHERE(d2d < 1.0E-6) d2d = 1.0E6
                     !dist_min = MIN(SQRT(MINVAL(d2d)),dist_min)
                     !WHERE(d2d > target) d2d = target*100.0
                     dist_min = MIN((MINVAL(d2d)),dist_min)
                     DEALLOCATE(x2d,y2d,z2d,d2d)
                  END DO
               END DO
            END DO
         END DO
         dist_min = SQRT(dist_min)
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = EXP(DCC_EXP*(dist_min-target))
         IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target,sigma,vals(mtargets),dist_min
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_coilcoil_distance
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilcoil_distance
