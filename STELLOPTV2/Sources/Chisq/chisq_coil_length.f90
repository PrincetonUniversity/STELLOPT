!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_length
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/29/2025
!     Description:   Calculate the coil length of each coil.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_length(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: nw_coil, nh_coil, rho_coil_kts
      USE spline_coils_mod, ONLY: get_coil_dl, get_coil_ns
      USE biotsavart, ONLY: coil_group
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(nsd), INTENT(in)    ::  target
      REAL(rprec), DIMENSION(nsd), INTENT(in)    ::  sigma
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
      REAL(rprec), PARAMETER :: curve_k1 = 0.25 ! >= 0.0
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
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_LENGTH ',numcoilgroups,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  COILGROUP  LENGTH'
      IF (niter >= 0) THEN
         curve_max = 0.0; curve_min = bigno
         DO k = 1, numcoilgroups
            L = 0.0
            CALL get_coil_ns(nc1)
            DO n = 1, nc1
               CALL get_coil_dl(k,n,dl)
               L   = L + dl
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = MAX(target(k),0.0_rprec)
            sigmas(mtargets)  = sigma(k)
            ! Pick form of functional based on target
            IF (target(k) <= 0) THEN ! Exponential
               vals(mtargets) = SQRT(EXP(L)*EXP(target(k)))
            ELSE
               vals(mtargets)    = L
            END IF
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3,1ES22.12E3)') &
                      targets(mtargets),sigma(k),vals(mtargets),k,L
         END DO
      ELSE
         DO k = 1, numcoilgroups
            IF (sigma(k) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_length
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_length
