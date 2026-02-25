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
      USE vsvd0, ONLY : nigroup
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(nigroup), INTENT(in)    ::  target
      REAL(rprec), DIMENSION(nigroup), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(in)    ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: numcoilgroups, k, n, nc1
      REAL(rprec) :: dl, L, val
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      k = COUNT(sigma<bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_LENGTH ',k,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  COILGROUP  LENGTH'
      IF (niter >= 0) THEN
         DO k = 1, nigroup
            IF (sigma(k)>=bigno) CYCLE
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
         DO k = 1, nigroup
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
