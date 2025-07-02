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
      INTEGER :: numcoilgroups, k, n
      REAL(rprec) :: tors, tors_hold, tors0, dl, L, val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_TORSION ',numcoilgroups*nw_coil*nh_coil,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MEAN'
      IF (niter >= 0) THEN
         numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
         DO k = 1, numcoilgroups*nw_coil*nh_coil
            CALL get_coil_torsion(k,1,tors0)
            CALL get_coil_dl(k,1,dl)
            tors0 = tors0 * dl
            val = tors0; L = dl
            DO n = 2, 128
               CALL get_coil_torsion(k,n,tors)
               CALL get_coil_dl(k,n,dl)
               val = val + tors * dl
               L   = L + dl
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)    = (val-tors0)/DBLE(128*L)
            IF (iflag == 1) WRITE(iunit_out,'(5ES22.12E3)') target,sigma,vals(mtargets)
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
