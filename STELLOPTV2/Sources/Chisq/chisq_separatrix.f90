!-----------------------------------------------------------------------
!     Subroutine:    chisq_separatrix
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate the maximum distance between a target
!                    separatrix and the plasma equilibrium.  Assume
!                    That we define the toroidal angle of the
!                    separatrix in terms of a field period.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_separatrix(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nu_max,nv_max)
      REAL(rprec), INTENT(in)    ::  sigma(nu_max,nv_max)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: u, u2, v, nv2, ier, dex
      REAL(rprec) :: s, xu, xv, R1, Z1, R2, Z2, dist, temp
      REAL(rprec) :: d1(nu_max)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      dex = SUM(COUNT(sigma_separatrix < bigno,DIM=1))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I5.5))') 'SEPARATRIX ',dex,6
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DIST  R  PHI  Z'
      IF (niter >= 0) THEN
         DO v = 1, nv_max
            DO u = 1, nu_max
               IF (sigma_separatrix(u,v) >= bigno) CYCLE
               R1 = r_separatrix(u,v)
               Z1 = z_separatrix(u,v)
               s  = 1.0_rprec
               xv = phi_separatrix(u,v)
               IF (xv < 0.0) xv = xv + pi2
               !xv = MOD(nfp*xv,pi2)
               xv = xv/nfp
               dist = bigno
               d1 = bigno
               DO u2 = 1, nu_max
                  ier = 0
                  xu = pi2*REAL(u2-1)/REAL(nu_max)
                  CALL get_equil_RZ(s,xu,xv,R2,Z2,ier)
                  IF (ier == 0) d1(u2) = SQRT((R2 - R1)*(R2 - R1) + (Z2 - Z1)*(Z2 - Z1))
               END DO
               dist = MINVAL(d1,DIM=1)
               mtargets = mtargets + 1
               targets(mtargets) = target(u,v)
               sigmas(mtargets)  = sigma(u,v)
               vals(mtargets)    = dist
               IF (iflag == 1) WRITE(iunit_out,'(6ES22.12E3)') target(u,v),sigma(u,v),vals(mtargets),R1,phi_separatrix(u,v),Z1
            END DO
         END DO
      ELSE
         DO v = 1, nv_max
            DO u = 1, nu_max
               IF (sigma_separatrix(u,v) >= bigno) CYCLE
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_separatrix
            END DO
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_separatrix
