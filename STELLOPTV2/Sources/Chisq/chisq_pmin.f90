!-----------------------------------------------------------------------
!     Subroutine:    chisq_pmin
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Constrain pressure to be above a minimum value.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_pmin(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: pres_spl, nrad, rho
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(in)    ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: pmin
      REAL(rprec), ALLOCATABLE :: p_temp(:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'PMIN ',1,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  PMIN  CHI'
      IF (niter >= 0) THEN
         pmin = 0.0
         ALLOCATE(p_temp(nrad))
         CALL EZspline_interp(pres_spl,nrad,rho,p_temp,ier)
         pmin = MINVAL(p_temp,DIM=1)
         DEALLOCATE(p_temp)
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         ! Impose a limit on the calculated pressure, below aspec_max
         ! it should match target but above it should match grow
         IF (pmin > target) THEN
            vals(mtargets) = target
         ELSE
            vals(mtargets) = target - (pmin-target)/width_pmin
         END IF
         !vals(mtargets)    = target + 1 - TANH((pmin-target)/(width_pmin))
         !PRINT *,target,sigma,pmin,vals(mtargets)
         IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target,sigma,pmin,vals(mtargets)
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_pmin
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_pmin
