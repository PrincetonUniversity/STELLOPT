!-----------------------------------------------------------------------
!     Subroutine:    chisq_sfincs_bootstrap
!     Authors:       E. Paul (ejpaul@umd.edu)
!     Date:          05/04/2018
!     Description:
!-----------------------------------------------------------------------
      SUBROUTINE chisq_sfincs_bootstrap(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: sfincs_s, sfincs_J_dot_B_flux_surface_average, sfincs_B_squared_flux_surface_average
      USE safe_open_mod, ONLY: safe_open
      USE mpi_params, ONLY: myid, master
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in) ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ik, Nradii
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN

      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'SFINCS BOOTSTRAP ',ik,10
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL'
      ! Compute length of sfincs_s array
      Nradii = MINLOC(sfincs_s(2:),DIM=1)
      IF (niter >= 0) THEN
         ! Iterate over sfincs_s
         DO ik = 1, Nradii
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               ! Here +1 is needed - s = 0 corresponds to ik=1
               vals(mtargets) = sfincs_J_dot_B_flux_surface_average(ik+1)/sqrt(sfincs_B_squared_flux_surface_average(ik+1))
               IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ik), &
                  sigma(ik), vals(mtargets)
            END IF
         END DO
      ELSE
         IF (ANY(sigma < bigno)) THEN
           DO ik = 1, Nradii
              IF (sigma(ik) < bigno) THEN
                 mtargets = mtargets + 1
                 IF (niter == -2) THEN
                  target_dex(mtargets)=jtarget_sfincs_bootstrap
                 END IF
              END IF
           END DO
         END IF
      END IF
      !print *,"End of chisq_sfincs_bootstrap."
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_sfincs_bootstrap
