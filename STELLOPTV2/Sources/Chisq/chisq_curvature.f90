!-----------------------------------------------------------------------
!     Subroutine:    chisq_curvature
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/11/2013
!     Description:   This subroutine targets the curvature of the
!                    equilibrium edge to prevent cusps.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_curvature(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: i,j, ier
      INTEGER, PARAMETER :: nbdry_pts = 90
      REAL(rprec) :: s,u,v, xp,zp,xpp,zpp, denom
      REAL(rprec) :: R_grad(3), Z_grad(3)
      REAL(rprec), DIMENSION(nbdry_pts) :: kappa_avg, kappa_max, kurtosis,&
                                           kappa_sig, kappa_mu
      REAL(rprec) :: kappa(nbdry_pts,nbdry_pts)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'CURVATURE_KERT ',1,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  KURTOSIS  KURTOSIS_AVG  KURTOSIS_MAX  PHI'
      IF (niter >= 0) THEN
         s=1.0
         DO i = 1, nbdry_pts
            DO j = 1, nbdry_pts
               u = pi2*(i-1)/nbdry_pts
               v = pi2*(j-1)/nbdry_pts
               CALL get_equil_kappa(s,u,v,kappa(i,j),ier)
            END DO
         END DO
         kappa_avg = SUM(kappa,DIM=2)/nbdry_pts
         DO i = 1, nbdry_pts
            DO j = 1, nbdry_pts
              kappa_max(j) = MAX(kappa_max(j),kappa(i,j))
              kappa_sig(j) = kappa_sig(j) + (kappa(i,j)-kappa_avg(j)*kappa_avg(j))/nbdry_pts
              kappa_mu(j) = kappa_mu(j) + (kappa(i,j)-kappa_avg(j)*kappa_avg(j)*kappa_avg(j)*kappa_avg(j))/nbdry_pts
            END DO
         END DO
         kurtosis = kappa_mu/(kappa_sig*kappa_sig)
         DO i = 1, nbdry_pts
            mtargets = mtargets + 1
            v = pi2*(i-1)/nbdry_pts
            targets(mtargets) = 0.0
            sigmas(mtargets)  = bigno
            vals(mtargets)    = 0.44 + 0.5*TANH((kurtosis(i) -20)/15)
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') target,sigma,vals(mtargets),kurtosis(i),kappa_avg(i),kappa_max(i),v
         END DO
      ELSE
         IF (sigma < bigno) THEN
            DO i =1, nbdry_pts
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_curvature
            END DO
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_curvature
