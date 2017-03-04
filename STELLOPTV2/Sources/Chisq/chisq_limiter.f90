!-----------------------------------------------------------------------
!     Subroutine:    chisq_limiter
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/25/2013
!     Description:   Calculates the distance from the equilibrium
!                    surface to a limiter surface.  If the equilibrium
!                    falls outside the surface, a large chisq is
!                    returned.  Target essentially sets the point and
!                    slope at which the effect of the limiter is felt.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_limiter(target,sigma,niter,iflag)
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
      INTEGER :: u, u2, v, nv2, ier, dex, nu_val
      REAL(rprec) :: s, xu, xv, R1, Z1, R2, Z2, dist, temp, lim_sgn, &
                     RA, ZA, rho_plasma, rho_limiter, RP, ZP, RN, ZN, &
                     RD, ZD, RX, ZX
      REAL(rprec) :: d1(nu_max)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      dex = SUM(COUNT(sigma_limiter < bigno,DIM=1))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I5.5))') 'LIMITER ',dex,6
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DIST  R  PHI  Z'
      IF (niter >= 0) THEN
         xu = 0.0
         xv = 0.0
         s  = 0.0
         CALL get_equil_RZ(s,xu,xv,RA,ZA,ier)
         R1 = r_limiter(1,1)
         Z1 = z_limiter(1,1)
         R2 = r_limiter(2,1)
         Z2 = z_limiter(2,1)
         xu = ATAN2(Z1-ZA,R1-RA)
         xv = ATAN2(Z2-ZA,R2-RA)
         lim_sgn = 1.0
         IF (xv < xu) lim_sgn = -1.0
         DO v = 1, nv_max
            nu_val = MAXLOC(sigma_limiter(:,v),DIM=1)
            !PRINT *,'nu_val',nu_val
            DO u = 1, nu_val-2
               IF (sigma_limiter(u,v) >= bigno) CYCLE
               s  = 1.0_rprec
               xv = phi_limiter(u,v)
               IF (xv < 0.0) xv = xv + pi2
               xv = MOD(nfp*xv,pi2)
               R1 = r_limiter(u,v)
               Z1 = z_limiter(u,v)
               R2 = r_limiter(u+1,v)
               Z2 = z_limiter(u+1,v)
               RN = lim_sgn*Z2-Z1
               ZN = -lim_sgn*(R2-R1)
               RD = R2-R1
               ZD = Z2-Z1
               temp = SQRT(RD*RD+ZD*ZD)
               RD = RD/temp
               ZD = ZD/temp
               d1 = 0.0
               !PRINT *,'R1,Z1',R1,Z1
               !PRINT *,'R2,Z2',R2,Z2
               !PRINT *,'RD,ZD',RD,ZD
               !PRINT *,'RD,ZD',RN,ZN
               DO u2 = 1, nu_max
                  ier = 0
                  xu = pi2*REAL(u2-1)/REAL(nu_max)
                  CALL get_equil_RZ(s,xu,xv,RP,ZP,ier)
                  ! Now X1 is the first point X2 is the second point an XP is the plasma point
                  dist =   ABS( (R2-R1)*(Z1-ZP)-(R1-RP)*(Z2-Z1) ) &
                         / SQRT( (R2-R1)*(R2-R1)+(Z2-Z1)*(Z2-Z1) )
                  ! Now (RN, ZN) points out of the vessel
                  ! and (RD, ZD) points along the vessel
                  ! THEN X = P - X1 + ((P-X1).D)D points from the line to the point p and is normal to the line segment
                  RX = RP - (R1 + ((RP-R1)*RD+(ZP-Z1)*ZD)*RD)
                  ZX = ZP - (Z1 + ((RP-R1)*RD+(ZP-Z1)*ZD)*ZD)
                  ! So if X.N is positive  the point is outside, if X.N is negative the point is inside
                  temp = RX*RN + ZX*ZN
                  d1(u2) = dist
                  IF (temp > 0) d1(u2) = -dist
                  !PRINT *,'dist',u2,d1(u2)
               END DO
               dist = MINVAL(d1,DIM=1)
               mtargets = mtargets + 1
               targets(mtargets) = target(u,v)
               sigmas(mtargets)  = sigma(u,v)
               IF (dist < target(u,v)) THEN
                  vals(mtargets) = dist
               ELSE
                  vals(mtargets) = target(u,v)
               END IF
               IF (iflag == 1) WRITE(iunit_out,'(6ES22.12E3)') target(u,v),sigma(u,v),vals(mtargets),R1,phi_limiter(u,v),Z1
            END DO
         END DO
      ELSE
         DO v = 1, nv_max
            DO u = 1, nu_max
               IF (sigma_limiter(u,v) >= bigno) CYCLE
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_limiter
            END DO
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_limiter
