!-----------------------------------------------------------------------
!     Subroutine:    chisq_jinvariant
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/11/2013
!     Description:   This subroutine calculates the difference between
!                    target and equilibrium J-Invariant for various
!                    values of ep/mu.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_jinvariant(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
!DEC$ IF DEFINED (JINV_OPT)
      USE BJdata
!DEC$ ENDIF
      
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
      LOGICAL :: fix_pitch
      LOGICAL, ALLOCATABLE :: ljinvar(:,:)
      INTEGER :: ik, j, k
      REAL(rprec) :: min_epmu, max_epmu, J_min, J_max, J_avg
      REAL(rprec), ALLOCATABLE :: epmu(:)
      REAL(rprec), ALLOCATABLE :: J_inv_all(:,:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!DEC$ IF DEFINED (JINV_OPT)
!      IF (iflag < 0) RETURN
!      ik   = COUNT(sigma < bigno)
!      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'J_INVARIANT ',1,4
!      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  J_INV_ALL  NS  NUM_EP_MU  NUMJSTAR'
!      IF (niter >= 0) THEN
!         device = "qas"
!         ALLOCATE(J_inv_all(NumJstar,num_ep_mu),epmu(num_ep_mu), &
!                  ljinvar(NumJstar,num_ep_mu))
!         DO ik = 1, nsd
!            IF (sigma(ik) >= bigno) CYCLE
!            CALL j_inv_calc(ik,num_ep_mu,NumJstar,fix_pitch,min_epmu,max_epmu,J_inv_all,epmu)
!            ljinvar = J_inv_all > 0.0
!            DO j = 1, num_ep_mu
!               J_min = MINVAL(J_inv_all(:,j))
!               J_max = MAXVAL(J_inv_all(:,j))
!               J_avg = SUM(J_inv_all(:,j),MASK=ljinvar(:,j))/COUNT(LJINVAR(:,j))
!               mtargets = mtargets + 1
!               DO k = 1, NumJstar
!                  targets(mtargets) = J_avg
!                  sigmas(mtargets)  = sigma(ik)
!                  vals(mtargets)    = J_inv_all(k,j)
!                  IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,3(2x,I4.4))') J_avg,sigma(ik),J_inv_all(k,j),ik,j,k
!               END DO
!            END DO
!         END DO
!         DEALLOCATE(J_inv_all,num_ep_mu,ljinvar)
!      ELSE
!         DO ik = 1, nsd
!            IF (sigma >= bigno) CYCLE
!            DO j =1, num_ep_mu
!               DO k = 1, NumJstar
!                  mtargets = mtargets + 1
!                  IF (niter == -2) target_dex(mtargets)=jtarget_jinvariant
!               END DO
!            END DO
!         END DO
!      END IF
!DEC$ ENDIF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_jinvariant
