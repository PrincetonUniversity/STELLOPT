!-----------------------------------------------------------------------
!     Subroutine:    chisq_kappa
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/18/2017
!     Description:   Calculate difference between equilibrium elongation
!                    and target elongation.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_kappa(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE stel_tools
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: nt = 359
      INTEGER     :: i, ier
      REAL(rprec) :: s,u,v,Rout,temp,Rin, kappa, phi
      REAL(rprec) :: Zarr(nt)
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'KAPPA ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET SIGMA KAPPA'
      IF (niter >= 0) THEN
         s=1; u=0; v=0; ier = 0
         CALL get_equil_RZ(s,u,v,Rout,temp,ier)
         temp = r0-0.1*(Rout-r0) ! 10% minor radius in negative direction
         ier = 0
         phi = v/nfp
         CALL get_equil_s(temp,phi,z0,s,ier,u)
         s=1; v=0; ier = 0
         CALL get_equil_RZ(s,u,v,Rin,temp,ier)
         DO i = 1, nt
            s=1; u = REAL((i-1))/REAL(nt); v= 0; ier = 0
            CALL get_equil_RZ(s,u,v,temp,Zarr(i),ier)
         END DO
         kappa = (MAXVAL(Zarr)-MINVAL(Zarr))/(Rout-Rin)
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = kappa
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,kappa
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_kappa
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_kappa
