!-----------------------------------------------------------------------
!     Subroutine:    chisq_LgradB
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion)
!     Date:          08/04/2025
!     Description:   Computes L_grad(B) used as a proxy for the
!                    coil plasma distance as defined in:
!                    https://dx.doi.org/10.1088/1361-6587/ad1a3e
!-----------------------------------------------------------------------
      SUBROUTINE chisq_lgradb(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stel_tools, ONLY: get_equil_LgradB
      USE equil_vals, ONLY: nfp
      
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
      INTEGER, PARAMETER :: nu_local = 128
      INTEGER, PARAMETER :: nv_local = 96
      INTEGER :: u, v, ier
      REAL(rprec) :: s,theta,zeta
      DOUBLE PRECISION :: lgradB, R, Z
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I6.6))') 'LGRADB ',nu_local*nv_local,8
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  LGRADB  S  THETA  PHI  R  Z'
      IF (niter >= 0) THEN
         s = 1.0_rprec
         lgradB = bigno
         DO u = 1, nu_local
            DO v = 1, nv_local
               theta = DBLE(u-1)/DBLE(nu_local)*pi2
               zeta  = DBLE(v-1)/DBLE(nv_local)*pi2
               ier = 0
               CALL get_equil_LgradB(s, theta, zeta, lgradB, ier, R_out=R, Z_out=Z)
               mtargets = mtargets + 1
               targets(mtargets) = target
               sigmas(mtargets)  = sigma
               vals(mtargets)    = lgradB
               IF (lgradB > target) sigmas(mtargets)  = sigma*10.0
               IF (iflag == 1) WRITE(iunit_out,'(8ES22.12E3)') target, sigmas(mtargets), lgradB, s, theta, zeta/nfp, R, Z
            END DO
         END DO
      ELSE
         IF (sigma < bigno) THEN
            DO u = 1, nu_local
               DO v = 1, nv_local
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets)=jtarget_lgradb
               END DO
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_lgradb
