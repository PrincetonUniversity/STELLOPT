!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_distortion
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          02/15/2025
!     Description:   Computes the coil distortion using a similiar
!                    deffinition to the ONSET code.
!                    F = \frac{1}{2\pi}\int g(\kappa)\kappa(s)ds
!                    g(\kappa) = (\kappa-\kappa_min)/(\kappa_max-\kappa_min)
!                              = 0 \kappa \le \kappa_min
!                    Factor of 2*pi omitted since s goes from 0 to 1.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_distortion(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: nw_coil, nh_coil, rho_coil_kts
      USE spline_coils_mod, ONLY: get_coil_curvature, get_coil_ns
      
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
      INTEGER :: numcoilgroups, k, n, nc1
      REAL(rprec) :: curve, g, dcurve, val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_DISTORTION ',numcoilgroups*nw_coil*nh_coil,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DISTORTION  COILGROUP'
      IF (niter >= 0) THEN
         dcurve = distortion_curvature_max-distortion_curvature_min
         DO k = 1, numcoilgroups*nw_coil*nh_coil
            val = 0.0
            CALL get_coil_ns(nc1)
            DO n = 1, nc1
               CALL get_coil_curvature(k,n,curve)
               g = 0.0
               IF (curve > distortion_curvature_min) g = (curve - distortion_curvature_min) / dcurve
               val = val + g*curve
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)    = val/DBLE(nc1)
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') &
                      target,sigma,vals(mtargets),k
         END DO
      ELSE
         IF (sigma < bigno) THEN
            numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
            DO k = 1, numcoilgroups*nw_coil*nh_coil
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_distortion
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_distortion
