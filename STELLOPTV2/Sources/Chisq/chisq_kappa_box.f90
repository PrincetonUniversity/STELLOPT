!-----------------------------------------------------------------------
!     Subroutine:    chisq_kappa_box
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/25/2017
!     Description:   Calculate difference between equilibrium elongation
!                    and target elongation. (Box Method)
!-----------------------------------------------------------------------
      SUBROUTINE chisq_kappa_box(target,sigma,niter,iflag)
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
!        i           Index over theta points
!        ier         Error flag
!        s           Radial coordiante [0,1]
!        u           Poloidal coordiante [0,2*pi]
!        v           Toroidal coordiante [0,2*pi] (over field period)
!        kappa       Elongation
!        Zarr(nt)    Array of Vertical points [m]
!        Rarr(nt)    Array of Radial points [m]
!
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: nt = 359
      INTEGER     :: i, ier
      REAL(rprec) :: s,u,v,kappa
      REAL(rprec) :: Zarr(nt), Rarr(nt)
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'KAPPA_BOX ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET SIGMA KAPPA PHI'
      IF (niter >= 0) THEN
         v = MOD(phi_kapa_box,pi2/nfp)*nfp
         DO i = 1, nt
            s=1; u = pi2*REAL((i-1))/REAL(nt); ier = 0
            CALL get_equil_RZ(s,u,v,Rarr(i),Zarr(i),ier)
         END DO
         kappa = (MAXVAL(Zarr)-MINVAL(Zarr))/(MAXVAL(Rarr)-MINVAL(Rarr))
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = kappa
         IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target,sigma,kappa,phi_kappa_box
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_kappa_box
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_kappa_box
