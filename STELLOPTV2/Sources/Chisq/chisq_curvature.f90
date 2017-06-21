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
      INTEGER, PARAMETER :: ntheta_pts = 128
      INTEGER, PARAMETER :: nzeta_pts = 4
      REAL(rprec) :: s,u,v
      REAL(rprec), DIMENSION(nzeta_pts) :: kappa_avg, kappa_max, kurtosis,&
                                           kappa_sig, kappa_mu
      REAL(rprec) :: kappa(ntheta_pts,nzeta_pts)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'CURVATURE_KERT ',1,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  KURTOSIS  KURTOSIS_AVG  KURTOSIS_MAX  PHI'
      IF (niter >= 0) THEN
         s=1.0
         kappa = 0
         DO i = 1, ntheta_pts
            DO j = 1, nzeta_pts
               u = pi2*(i-1)/ntheta_pts
               v = pi2*(j-1)/nzeta_pts
               CALL get_equil_kappa(s,u,v,kappa(i,j),ier)
               !PRINT *,u,v,kappa(i,j),ier
            END DO
         END DO
         kappa_max = -1.0E30; kappa_sig = 0; kappa_mu = 0
         kappa_avg = SUM(kappa,DIM=1)/ntheta_pts
         DO j = 1, nzeta_pts
            DO i = 1, ntheta_pts
              kappa_max(j) = MAX(kappa_max(j),kappa(i,j))
              kappa_sig(j) = kappa_sig(j) + (kappa(i,j)-kappa_avg(j)*kappa_avg(j)) !/ntheta_pts
              kappa_mu(j) = kappa_mu(j) + (kappa(i,j)-kappa_avg(j)*kappa_avg(j)*kappa_avg(j)*kappa_avg(j)) !/ntheta_pts
            END DO
         END DO
         kappa_sig = kappa_sig / ntheta_pts
         kappa_mu   = kappa_mu / ntheta_pts
         kurtosis = kappa_mu/(kappa_sig*kappa_sig)
         DO j = 1, nzeta_pts
            mtargets = mtargets + 1
            v = pi2*(j-1)/nzeta_pts
            targets(mtargets) = target
            sigmas(mtargets)  = sigma
            vals(mtargets)    = 0.44 + 0.5*TANH((kurtosis(j) -20)/15)
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') target,sigma,vals(mtargets),kurtosis(j),kappa_avg(j),kappa_max(j),v
         END DO
      ELSE
         IF (sigma < bigno) THEN
            DO i =1, nzeta_pts
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
