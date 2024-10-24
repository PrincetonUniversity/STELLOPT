!-----------------------------------------------------------------------
!     Subroutine:    chisq_quasiiso
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/08/2024
!     Description:   This chisquared functional replicates the QI
!                    functional described in:
!                    https://doi.org/10.1017/S002237782300065X
!                    and
!                    https://doi.org/10.1103/PRXEnergy.3.023010
!-----------------------------------------------------------------------
      SUBROUTINE chisq_quasiiso(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ik, mn, ier, l, lmin
      INTEGER, PARAMETER :: nalpha = 65 ! must be odd
      INTEGER, PARAMETER :: ntheta0 = 5
      INTEGER, PARAMETER :: nlambda = 4
      REAL(rprec) :: s_temp, iota0, phi, theta, deltaphi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: modb, modbs


!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'QUASIISO ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         ! Loop over surfaces
         DO ik = 1, nsd
            IF (ABS(sigma(ik)) >= bigno) CYCLE
            ! Get iota for the surface
            s_temp = DBLE(ik)/DBLE(ns_b)
            ier = 0
            CALL get_equil_iota(s_temp,iota0,ier)
            ! Calculate the distance to follow a field line
            deltaphi = pi2/iota0
            ALLOCATE(modb(nalpha),modbs(nalpha))
            ! Loop over inital values of theta0
            DO i = 1, ntheta0
               DO l = 1, nalpha
                  phi = deltaphi*(l-1)/(nalpha-1)
                  theta = theta0 + iota0*phi
                  DO mn = 1, mnmax_b
                     modb(l) = bmnc_(mn,ik)*COS(ixm_b(mn)*theta+ixn_b(mn)*phi/nfp_b)
                  END DO
               END DO
               ! Find min/max values
               lmin = MINLOC(modb)
               ! Sort the arrays
               lh = FLOOR(nalpha/2)+1
               modb = CSHIFT(modb,lmin-lh) !BOOSH
               ! Squash the array
               modbs = modb
               DO l = lh,nalpha
                  IF (modbs(l)<modbs(l-1)) modbs(l) = modbs(l-1)
               END DO
               DO l = 1, lh
                  IF (modbs(l)<modbs(l+1)) modbs(l) = modbs(l+1)
               END DO
               ! Stretch the array
               Bmin = MINVAL(modb)
               Bmax = MAXVAL(modb)
               DO l = 1, lh
                  modbs(l) = Bmin + (Bmax-Bmin)*(modbs(l)-Bmin)/(modbs(1)-Bmin)
               END DO
               DO l = lh, nalpha
                  modbs(l) = Bmin + (Bmax-Bmin)*(modbs(l)-Bmin)/(modbs(nalpha)-Bmin)
               END DO
               ! Now calculate the J's
               
            END DO

         mtargets = mtargets + 1
         targets(mtargets) = 0.0
         sigmas(mtargets)  = bigno
         vals(mtargets)    = 0.0
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,0.0,vals(mtargets)
         END DO
      ELSE
         ! Consistency check
         mboz = MAX(6*mpol, 2, mboz)             
         nboz = MAX(2*ntor-1, 0, nboz)     
         ! CALCULATE mnboz_b becasue we don't know it yet (setup_booz.f)
         mnboz_b = (2*nboz+1)*(mboz-1) + (nboz + 1)
         DO ik = 1, nsd
            IF (ABS(sigma(ik)) < bigno) THEN
               lbooz(ik) = .TRUE.
               DO mn = 1, mnboz_b
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets) = jtarget_qi
               END DO
            END IF
         END DO
      END IF
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE chisq_quasiiso
