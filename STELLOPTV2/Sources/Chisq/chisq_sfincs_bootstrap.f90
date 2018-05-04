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
      USE stellopt_vars, ONLY: sfincs_s
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
      print *,"Hello world from chisq_bootstrap_sfincs.f90. niter=",niter
      IF (iflag < 0) RETURN

      ! This counts the number of elements of sigma that are < bigno
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'BOOTSTRAP ',ik,10
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  RHO  AVG_JDOTB  BEAM_JDOTB  BOOT_JDOTB  AJBBS  FACNU  BSNORM'
      IF (niter >= 0) THEN
         ! Compute length of sfincs_s array
         Nradii = MINLOC(sfincs_s(2:),DIM=1)
         print *,"In the chisq_sfincs_bootstrap (niter >= 0) block."
         ! Iterate over sfincs_s
         DO ik = 1, Nradii
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
         END DO
      ELSE
         print *,"In the chisq_sfincs_bootstrap (niter < 0) block."
         DO ik = 1, Nradii
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
            END IF
         END DO
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_sfincs_bootstrap
