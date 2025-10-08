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
      USE stellopt_vars, ONLY: nw_coil, nh_coil, rho_coil_kts
      USE spline_coils_mod, ONLY: get_coil_torsion, get_coil_dl
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
      REAL(rprec) :: tors, tors_hold, tors0, dl, L, val, hypc, tors_max, tors_min

      ! The following mimics the FOCUS code algorithm
      INTEGER, PARAMETER :: penfun_tors = 0
      REAL(rprec), PARAMETER :: tors_k0 = 1.0 ! >= 0.0
      REAL(rprec), PARAMETER :: tors_k1 = 1.0 ! >= 0.0
      REAL(rprec), PARAMETER :: tors_alpha = 1.0 ! >= 0.0
      REAL(rprec), PARAMETER :: tors_beta = 2.0 ! >= 2.0
      REAL(rprec), PARAMETER :: tors_gamma = 1.0 ! >= 1.0
      REAL(rprec), PARAMETER :: tors_sigma = 1.0 ! >= 0.0
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_TORSION ',numcoilgroups*nw_coil*nh_coil,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MEAN  COILGROUP  MAX  MIN'
      IF (niter >= 0) THEN
         numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
         tors_max = 0.0; tors_min = bigno
         DO k = 1, numcoilgroups*nw_coil*nh_coil
            val = 0.0; L = 0.0
            nc1 = SIZE(coil_group(k)%coils(1)%xnod,2)
            DO n = 1, nc1
               CALL get_coil_torsion(k,n,tors)
               CALL get_coil_dl(k,n,dl)
               tors_max = MAX(tors,tors_max)
               tors_min = MIN(tors,tors_min)
               IF (tors > tors_k0) THEN
                  IF (penfun_tors == 1) THEN
                     hypc = 0.5 * EXP( tors_alpha * ( tors - tors_k0 ) ) &
                          + 0.5 * EXP(-tors_alpha * ( tors - tors_k0 ) )
                     tors_hold = ( hypc - 1.0 )**2
                  ELSE
                     tors_hold = ( tors_alpha * ( tors - tors_k0 ) )**tors_beta
                  END IF
               END IF
               IF (tors > tors_k1) THEN
                  tors_hold = tors_hold + tors_sigma * ( (tors - tors_k1)**tors_gamma)
               END IF
               L   = L + dl
               val = val + tors_hold * dl
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)    = val/DBLE(nc1*L)
            IF (iflag == 1) WRITE(iunit_out,'(5ES22.12E3)') target,sigma,vals(mtargets),k,tors_max,tors_min
         END DO
      ELSE
         IF (sigma < bigno) THEN
            numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
            DO k = 1, numcoilgroups*nw_coil*nh_coil
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_torsion
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_torsion
