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

      LOGICAL, PARAMETER :: lwrite_out = .False.


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
               ! Calculate the initial distance to follow a field line
               deltaphi = pi2/iota0
               modb = 0.0
               modbs = 0.0
               theta0 = pi2*(i-1)/ntheta0
               DO l = 1, nalpha
                  phi = deltaphi*(l-1)/(nalpha-1)
                  theta = theta0 + iota0*phi
                  DO mn = 1, mnboz_b
                     modb(l) = modb(l)+bmnc_b(mn,ik)*COS(DBLE(ixm_b(mn))*theta+DBLE(ixn_b(mn))*phi/nfp_b)
                  END DO
               END DO
               IF (lwrite_out) WRITE(320,*) modb
               ! Find min/max modB values
               lmin = MINLOC(modb,DIM=1)
               lmax = MAXLOC(modb,DIM=1)
               ! Recalculate length of fieldline
               phimin = deltaphi*(lmin-1)/(nalpha-1)
               phimax = deltaphi*(lmax-1)/(nalpha-1)
               deltaphi = ABS(phimax-phimin)
               ! Now recalculate modb centered around lmin
               modb = 0
               DO l = 1, nalpha
                  phi = phimin-deltaphi + 2*deltaphi*(l-1)/(nalpha-1)
                  !phi = deltaphi*(l-lmin-1)/(nalpha-1) - deltaphi/2
                  theta = theta0 + iota0*phi
                  DO mn = 1, mnboz_b
                     modb(l) = modb(l)+bmnc_b(mn,ik)*COS(DBLE(ixm_b(mn))*theta+DBLE(ixn_b(mn))*phi/nfp_b)
                  END DO
               END DO
               IF (lwrite_out) WRITE(321,*) modb
               ! Sort the arrays
               !lh = FLOOR(nalpha*0.5)+1
               !modb = CSHIFT(modb,lmin-lh) !BOOSH
               lh = MINLOC(modb,DIM=1)
               lh = MIN(MAX(lh,1),nalpha)
               ! Squash the array
               modbs = modb
               DO l = lh, 1, -1
                  IF (modbs(l)<modbs(l+1)) modbs(l) = modbs(l+1)
               END DO
               DO l = lh,nalpha
                  IF (modbs(l)<modbs(l-1)) modbs(l) = modbs(l-1)
               END DO
               IF (lwrite_out) WRITE(322,*) modbs
               ! Stretch the array
               Bmin = MINVAL(modb)
               Bmax = MAXVAL(modb)
               DO l = 1, lh
                  modbs(l) = Bmin + (Bmax-Bmin)*(modbs(l)-Bmin)/(modbs(1)-Bmin)
               END DO
               DO l = lh, nalpha
                  modbs(l) = Bmin + (Bmax-Bmin)*(modbs(l)-Bmin)/(modbs(nalpha)-Bmin)
               END DO
               IF (lwrite_out) WRITE(323,*) modbs
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
                  IF (lwrite_out) print *, nlambda,i1,i2,Bmin,Bmir,Bmax
                  i1 = MAX(i1,1); i2 = MIN(i2,nalpha)
                  integral_A = 1.0 - lambda * modb
                  integral_A = SIGN(1.0_rprec,integral_A)*SQRT(ABS(integral_A))*dl
                  IF (lwrite_out) WRITE(324,*) integral_A
                  integral_B = 1.0 - lambda * modbs
                  integral_B = SQRT(integral_B)*dl
                  IF (lwrite_out) WRITE(325,*) integral_B
                  IF (lwrite_out) PRINT *,'A :',integral_A(i1:i2)
                  IF (lwrite_out) PRINT *,'B :',integral_B(i1:i2)
                  J_I(m,i) = SUM(integral_A(i1:i2))
                  J_C(m,i) = SUM(integral_B(i1:i2))
               END DO
            END DO
            IF (lwrite_out) WRITE(326,*) J_I
            IF (lwrite_out) WRITE(327,*) J_C
            ! Construct fuctional value
            ftemp = 0.0_rprec
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
