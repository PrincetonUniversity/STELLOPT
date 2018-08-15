!-----------------------------------------------------------------------
!     Subroutine:    chisq_alpha_gc
!     Authors:       A. Bader, A. Ware, B. Faber
!     Date:          08/15/2018
!     Description:   This subroutine targets alpha particle confinement using 
!                    a guiding center calculation written by V.V. Nemov
!                    It is a more targeted metric than the generalized 
!                    BEAMS3D.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_alpha_gc(target,sigma,niter,iflag)
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

      INTEGER :: ik 

      !----------------BEGIN SUBROUTINE --------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'ALPHA CONFINEMENT ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DIFFERENCE  #'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               ! Here is where you call nemov's code
               ! The two variables you need are "ik" which is the flux surface label
               ! and "wout_"//TRIM(PROC_STRING)//".nc" which is the wout file
               ! the result for each flux surface will be set in 
               ! vals(mtargets) which right now holds a dummy value
               vals(mtargets)    = 0.0 !dummy value

               IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') target(ik),sigma(ik),vals(mtargets),ik
               IF (iflag == 1) CALL FLUSH(iunit_out)
            END IF
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_alpha_gc
            END IF
         END DO
      END IF
      END SUBROUTINE chisq_alpha_gc
