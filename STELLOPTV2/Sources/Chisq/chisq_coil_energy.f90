!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_energy
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          01/19/2026
!     Description:   Calculate the magnetic energy of a coil based on
!                    the formulas found in:
!                    https://arxiv.org/pdf/2310.12087
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_energy(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: nw_coil, nh_coil, width_coil, height_coil, &
      rho_coil_kts
      USE spline_coils_mod, ONLY: compute_coil_energy
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
      INTEGER :: k
      REAL(rprec) :: E
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      k   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_ENERGY ',k,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  COIL'
      IF (niter >= 0) THEN
         DO k = 1, nigroup
            IF (sigma(k)>=bigno) CYCLE
            CALL compute_coil_energy(nw_coil,nh_coil,width_coil,height_coil,k,E)
            mtargets = mtargets + 1
            targets(mtargets) = target(k)
            sigmas(mtargets)  = sigma(k)
            vals(mtargets)    = E
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') &
                      targets(mtargets),sigma(k),vals(mtargets),k
         END DO
      ELSE
         DO k = 1, nigroup
            IF (sigma(k) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_energy
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_energy
