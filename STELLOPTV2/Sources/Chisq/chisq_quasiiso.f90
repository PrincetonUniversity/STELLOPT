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
      USE equil_utils, ONLY: get_equil_iota
      USE read_boozer_mod
      USE vmec_input, ONLY: mpol, ntor
      
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
      INTEGER :: ik, mn, ier, i, l, m, lmin, lmax, i1, i2, lh, lr
      INTEGER, PARAMETER :: nalpha = 65 ! must be odd
      INTEGER, PARAMETER :: ntheta0 = 5
      INTEGER, PARAMETER :: nlambda = 4
      REAL(rprec) :: s_temp, iota0, phi, theta, deltaphi, ftemp, norm, &
         Bmax, Bmin, Bmir, lambda, theta0, phimin,phimax
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: modb, modbs, dl, &
         integral_A, integral_B
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: J_I, J_C


!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'QUASIISO ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  K'
      IF (niter >= 0) THEN
         ! Allocate helpers along field line
         ALLOCATE(modb(nalpha),modbs(nalpha),dl(nalpha),integral_A(nalpha),&
            integral_B(nalpha))
         ! Allocate J arrays
         ALLOCATE(J_C(nlambda,ntheta0),J_I(nlambda,ntheta0))
         ! Loop over surfaces
         DO ik = 1, nsd
            IF (ABS(sigma(ik)) >= bigno) CYCLE
            ! Get iota for the surface
            s_temp = DBLE(ik)/DBLE(ns_b)
            ier = 0
            CALL get_equil_iota(s_temp,iota0,ier)
            ! Loop over inital values of theta0
            J_C = 0.0; J_I = 0.0
            DO i = 1, ntheta0
               ! Calculate the distance to follow a field line
               deltaphi = pi2/iota0
               modb = 0.0
               modbs = 0.0
               theta0 = pi2*(i-1)/ntheta0
               !theta0 = 0.0
               DO l = 1, nalpha
                  phi = deltaphi*(l-1)/(nalpha-1)
                  theta = theta0 + iota0*phi
                  DO mn = 1, mnboz_b
                     !PRINT *,ik,mn,bmnc_b(mn,ik),ixm_b(mn),theta,ixn_b(mn),phi,nfp_b
                     modb(l) = modb(l)+bmnc_b(mn,ik)*COS(DBLE(ixm_b(mn))*theta+DBLE(ixn_b(mn))*phi/nfp_b)
                  END DO
               END DO
               !WRITE(290+i,*) modb; CALL FLUSH(290+i)
               ! Find min/max values
               lmin = MINLOC(modb,DIM=1)
               lmax = MAXLOC(modb,DIM=1)
               phimin = deltaphi*(lmin-1)/(nalpha-1)
               phimax = deltaphi*(lmax-1)/(nalpha-1)
               deltaphi = ABS(phimax-phimin)*2.0
               ! Now recalculate modb centered around lmin
               modb = 0
               DO l = 1, nalpha
                  phi = deltaphi*(l-lmin-1)/(nalpha-1) - deltaphi/2
                  theta = theta0 + iota0*phi
                  DO mn = 1, mnboz_b
                     !PRINT *,ik,mn,bmnc_b(mn,ik),ixm_b(mn),theta,ixn_b(mn),phi,nfp_b
                     modb(l) = modb(l)+bmnc_b(mn,ik)*COS(DBLE(ixm_b(mn))*theta+DBLE(ixn_b(mn))*phi/nfp_b)
                  END DO
               END DO
               ! Sort the arrays
               lh = FLOOR(nalpha*0.5)+1
               PRINT *,lh,lmin
               !WRITE(300+i,*) modb; CALL FLUSH(300+i)
               modb = CSHIFT(modb,lmin-lh) !BOOSH
               !WRITE(310+i,*) modb; CALL FLUSH(310+i)
               ! Squash the array
               modbs = modb
               DO l = lh,nalpha
                  IF (modbs(l)<modbs(l-1)) modbs(l) = modbs(l-1)
               END DO
               DO l = 1, lh
                  IF (modbs(l)<modbs(l+1)) modbs(l) = modbs(l+1)
               END DO
               !WRITE(320+i,*) modbs; CALL FLUSH(320+i)
               ! Stretch the array
               Bmin = MINVAL(modb)
               Bmax = MAXVAL(modb)
               !PRINT *,Bmin,Bmax
               DO l = 1, lh
                  modbs(l) = Bmin + (Bmax-Bmin)*(modbs(l)-Bmin)/(modbs(1)-Bmin)
               END DO
               DO l = lh, nalpha
                  modbs(l) = Bmin + (Bmax-Bmin)*(modbs(l)-Bmin)/(modbs(nalpha)-Bmin)
               END DO
               !WRITE(330+i,*) modbs; CALL FLUSH(330+i)
               ! We need dl (missing B.grad(phi) term)
               dl = modb*deltaphi/(nalpha-1)
               !WRITE(340+i,*) dl; CALL FLUSH(340+i)
               ! Now we evaluate integrals for different values of lambda
               ! Alpha is enclosing loop
               DO m = 1, nlambda
                  Bmir = Bmin+0.9*(Bmax-Bmin)*DBLE(m)/DBLE(nlambda)
                  lambda = 1.0/Bmir
                  i1 = COUNT(modbs(1:lh)>Bmir) + 1
                  i2 = COUNT(modbs(lh:nalpha)<Bmir) + lh - 1
                  print *, nlambda,i1,i2,Bmin,Bmir,Bmax
                  i1 = MAX(i1,1); i2 = MIN(i2,nalpha)
                  integral_A = 1.0 - lambda * modb
                  !WRITE(400+m,*) integral_A; CALL FLUSH(400+m)
                  integral_A = SIGN(1.0_rprec,integral_A)*SQRT(ABS(integral_A))*dl
                  !WRITE(410+m,*) integral_A; CALL FLUSH(410+m)
                  integral_B = 1.0 - lambda * modbs
                  !WRITE(420+m,*) integral_B; CALL FLUSH(420+m)
                  integral_B = SQRT(integral_B)*dl
                  !WRITE(430+m,*) integral_B; CALL FLUSH(430+m)
                  PRINT *,'A :',integral_A(i1:i2)
                  PRINT *,'B :',integral_B(i1:i2)
                  J_I(m,i) = SUM(integral_A(i1:i2))
                  J_C(m,i) = SUM(integral_B(i1:i2))
               END DO
            END DO
            WRITE(500,*) J_I
            WRITE(501,*) J_C
            ! Construct fuctional value
            DO m = 1, nlambda
               DO i1 = 1, ntheta0
                  DO i2 = 1, ntheta0
                     ftemp = ftemp + J_I(m,i1) - J_C(m,i2)
                  END DO
               END DO
            END DO
            norm = SUM(J_C+J_I)/(nlambda*ntheta0)
            mtargets = mtargets + 1
            targets(mtargets) = 0.0
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = ftemp/norm
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') target(ik),sigma(ik),vals(mtargets),ik
         END DO
         DEALLOCATE(modb,modbs,dl,integral_A,integral_B,J_C,J_I)
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
                  IF (niter == -2) target_dex(mtargets) = jtarget_quasiiso
               END DO
            END IF
         END DO
      END IF
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE chisq_quasiiso
