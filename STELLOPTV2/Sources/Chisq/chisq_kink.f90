!-----------------------------------------------------------------------
!     Subroutine:    chisq_kink
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/10/2016
!     Description:   This subroutine handles loading the chi-squared
!                    values for kink stability as determined by
!                    TERPSICHORE.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_kink(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
!DEC$ IF DEFINED (TERPSICHORE)
      !USE tpr_param, ONLY: wp_global, wk_global, omega_global, growth_global
!DEC$ ENDIF
 
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsys)
      REAL(rprec), INTENT(in)    ::  sigma(nsys)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ik

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'KINK ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  WP/WK  WP  WK  OMEGA  GROWTH'
      IF (niter >= 0) THEN
         DO ik = 1, nsys
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            !vals(mtargets)    = growth_kink(ik)
            vals(mtargets)    = wp_kink(ik)/wk_kink(ik)  ! Rayleigh Quotient
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') target(ik),sigma(ik),vals(mtargets),wp_kink(ik),wk_kink(ik), omega_kink(ik), growth_kink(ik)
         END DO
      ELSE
         DO ik = 1, nsys
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_kink
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_kink
