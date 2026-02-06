!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_twist
!     Authors:       M.S. Rang (maxsrang@gmail.com)
!     Date:          02/02/2026
!     Description:   Calculate the mean local coil twist - the integral of the torsion
!                    https://en.wikipedia.org/wiki/Torsion_of_a_curve
!                    https://math.franklin.uga.edu/sites/default/files/users/user317/ShifrinDiffGeo.pdf
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_twist(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: nw_coil, nh_coil, rho_coil_kts
      USE spline_coils_mod, ONLY: get_coil_torsion, get_coil_dl, &
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
      REAL(rprec) :: tors, dl, L, val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'COIL_TWIST ',numcoilgroups*nw_coil*nh_coil,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  MEAN  COILGROUP  MAX  MIN'
      IF (niter >= 0) THEN
         DO k = 1, numcoilgroups*nw_coil*nh_coil
            val = 0.0; L = 0.0
            CALL get_coil_ns(nc1)
            DO n = 1, nc1
               CALL get_coil_torsion(k,n,tors)
               CALL get_coil_dl(k,n,dl)
               val = val + tors * dl
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)    = val
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3,2ES22.12E3)') & ! write format may be bugged (I deleted two printed variables in the next line)
                     target,sigma,vals(mtargets),k
         END DO
      ELSE
         IF (sigma < bigno) THEN
            numcoilgroups = COUNT(ANY(rho_coil_kts>0,DIM=2))
            DO k = 1, numcoilgroups*nw_coil*nh_coil
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coil_twist
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_twist
