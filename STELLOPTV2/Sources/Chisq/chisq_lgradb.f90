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
      USE stel_tools, ONLY: get_equil_RZ
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
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
      REAL(rprec), DIMENSION(nu_local,nv_local) :: Bxds
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'LGRADB ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  LGRADB'
      IF (niter >= 0) THEN
         s = 1.0_rprec
         DO u = 1, nu_local
            DO v = 1, nv_local
               theta = DBLE(u-1)/DBLE(nu_local)
               zeta  = DBLE(v-1)/DBLE(nv_local)
               ier = 0
               ! First compute the vector gradients
               CALL get_equil_LgradB(s,u,v,R,Z,ier,Rgrad,Zgrad)
         ! Compute the Frobenius norm
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = beta
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,beta
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_lgradb
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_lgradb
