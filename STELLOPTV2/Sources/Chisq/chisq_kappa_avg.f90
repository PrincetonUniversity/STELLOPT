!-----------------------------------------------------------------------
!     Subroutine:    chisq_kappa_avg
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/25/2017
!     Description:   Calculate difference between equilibrium elongation
!                    and target elongation. (Avg. Method)
!-----------------------------------------------------------------------
      SUBROUTINE chisq_kappa_avg(target,sigma,niter,iflag)
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
!        j           Index over toroidal points
!        ier         Error flag
!        nz          Number of toroidal points
!        s           Radial coordiante [0,1]
!        u           Poloidal coordiante [0,2*pi]
!        v           Toroidal coordiante [0,2*pi] (over field period)
!        Rax         Radial Axis Location [m]
!        Zax         Vertical Axis Location [m]
!        phi         Toroidal Coordiante [0,2*pi] (over device)
!        temp        Temporary variable
!        Rout        Radial outboard point [m]
!        Rin         Raidal inboard point [m]
!        kappa       Elongation
!        Zarr(nt)    Array of Vertical points [m]
!
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: nt = 359
      INTEGER     :: i, j, ier, nz
      REAL(rprec) :: s, u, v, Rax, Zax, phi, Rout,temp,Rin, kappa
      REAL(rprec) :: Zarr(nt)
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'KAPPA_AVG ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET SIGMA KAPPA'
      IF (niter >= 0) THEN
         nz = 180/nfp
         kappa = 0
         DO j = 1, nz
            ! Need axis Point
            s = 0; u = 0; ier = 0
            v = pi2*REAL(j-1)/REAL(nz-1)
            CALL get_equil_RZ(s,u,v,Rax,Zax,ier)
            ! First calculated outboard point
            temp = Rax + 0.05*aminor
            ier = 0
            phi = v/nfp
            CALL get_equil_s(temp,phi,Zax,s,ier,u)
            s=1; ier = 0
            CALL get_equil_RZ(s,u,v,Rout,temp,ier) ! v on [0,2pi]
            ! Now Calculate Inboard point
            temp = Rax-0.1*(Rout-Rax) ! 10% minor radius in negative direction
            ier = 0
            phi = v/nfp
            CALL get_equil_s(temp,phi,Zax,s,ier,u)
            s=1; ier = 0
            CALL get_equil_RZ(s,u,v,Rin,temp,ier)
            ! Now Calculate Z
            DO i = 1, nt
               s=1; u = pi2*REAL((i-1))/REAL(nt-1); ier = 0
               CALL get_equil_RZ(s,u,v,temp,Zarr(i),ier)
            END DO
            kappa = kappa + (MAXVAL(Zarr)-MINVAL(Zarr))/(Rout-Rin)
         END DO
         kappa = kappa/nz
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)     = kappa
         IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target,sigma,kappa
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_kappa_avg
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_kappa_avg
