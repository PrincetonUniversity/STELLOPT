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
      USE stellopt_vars, ONLY: nw_coil, nh_coil, rho_coil_kts
      USE spline_coils_mod, ONLY: get_coil_curvature, get_coil_dl, &
                                  get_coil_ns
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
      INTEGER :: numcoilgroups, k, n, nc1
      REAL(rprec) :: curve, curve_hold, hypc, dl, L, val, curve_max, curve_min
      
      ! The following mimics the FOCUS code algorithm
      !   penfunc = 0 Minimize toward curve_k0 (from one side)
      !   penfunc = 1 Minimize toward curve_k0 (from both sides)
      INTEGER, PARAMETER :: penfun_curve = 0
      REAL(rprec), PARAMETER :: curve_k0 = 1.0 ! >= 0.0
      REAL(rprec), PARAMETER :: curve_k1 = 1.0 ! >= 0.0
      REAL(rprec), PARAMETER :: curve_alpha = 1.0 ! >= 0.0
      REAL(rprec), PARAMETER :: curve_beta = 2.0 ! >= 2.0
      REAL(rprec), PARAMETER :: curve_gamma = 1.0 ! >= 1.0
      REAL(rprec), PARAMETER :: curve_sigma = 1.0 ! >= 0.0
      ! Note if gamma == 1.0 then k1 = 0 
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_CURVATURE ',numcoilgroups*nw_coil*nh_coil,6
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MEAN  COILGROUP  MAX  MIN'
      IF (niter >= 0) THEN
         curve_max = 0.0; curve_min = bigno
         DO k = 1, numcoilgroups*nw_coil*nh_coil
            val = 0.0; L = 0.0
            CALL get_coil_ns(nc1)
            DO n = 1, nc1
               curve_hold = 0.0
               CALL get_coil_curvature(k,n,curve)
               CALL get_coil_dl(k,n,dl)
               curve_max = MAX(curve,curve_max)
               curve_min = MIN(curve,curve_min)
               IF (curve > curve_k0) THEN
                  IF (penfun_curve == 1) THEN
                     hypc = 0.5 * EXP( curve_alpha * ( curve - curve_k0 ) ) &
                          + 0.5 * EXP(-curve_alpha * ( curve - curve_k0 ) )
                     curve_hold = ( hypc - 1.0 )**2
                  ELSE
                     curve_hold = ( curve_alpha * ( curve - curve_k0 ) )**curve_beta
                  END IF
               END IF
               IF (curve > curve_k1) THEN
                  curve_hold = curve_hold + curve_sigma * ( (curve - curve_k1)**curve_gamma)
               END IF
               L   = L + dl
               val = val + curve_hold * dl
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)    = val/DBLE(nc1*L)
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3,2ES22.12E3)') &
                      target,sigma,vals(mtargets),k,curve_max,curve_min
         END DO
      ELSE
         IF (sigma < bigno) THEN
            numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
            DO k = 1, numcoilgroups*nw_coil*nh_coil
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_curvature
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_curvature
