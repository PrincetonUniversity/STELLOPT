!-----------------------------------------------------------------------
!     Subroutine:    chisq_penta_er
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!                    A. Coelho (antonio.coelho@gauss-fusion.com)
!     Date:          03/23/2026
!     Description:   Targeting of Er using PENTA code
!-----------------------------------------------------------------------
      SUBROUTINE chisq_penta_er(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: ER_PENTA
      USE penta_interface_mod, ONLY: read_penta_ion_params_namelist, &
         read_penta_run_params_namelist
      
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
      INTEGER :: ik, ij, ii, istat
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ! Print Header
      IF (iflag == 1) THEN
         ik   = COUNT(target_dex == jtarget_penta_er)
         WRITE(iunit_out,'(A,2(2X,I3.3))') 'PENTA_ER ', ik, 4
         WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  K'
      END IF
      IF (niter >= 0) THEN
         ii = 1
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = ER_PENTA(ii)
            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') target(ik),sigma(ik),vals(mtargets),ik
            IF (iflag == 1) CALL FLUSH(iunit_out)
            ii = ii + 1
         END DO
      ELSE
         DO ii = 1, nsd
            IF (sigma(ii) >= bigno) CYCLE
            lbooz(ii) = .TRUE.
            lneed_dkes(ii) = .TRUE.
            lneed_penta(ii) = .TRUE.
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_penta_er
            !DO ij = 1, nprof
            !   IF (E_dkes(ij) <= -bigno .or. nu_dkes(ij) <= -bigno) CYCLE
            !   nruns_dkes = nruns_dkes + 1
            !END DO
         END DO
         istat = 0
         CALL read_penta_ion_params_namelist("input."//TRIM(id_string),istat)
         CALL read_penta_run_params_namelist("input."//TRIM(id_string),istat)
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_penta_er
