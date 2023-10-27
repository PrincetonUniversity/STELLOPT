!-----------------------------------------------------------------------
!     Subroutine:    chisq_mercier
!     Authors:       E. Sanchez (edi.sanchez@ciemat.es))
!     Date:          08/05/2023
!     Description:   This subroutine calculates the mercier stability (DMercier) target contribution to Chisquare
!-----------------------------------------------------------------------
      SUBROUTINE chisq_dmercier(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE equil_utils
      USE stellopt_targets
      USE equil_vals, ONLY: Baxis
      USE read_wout_mod, ONLY: read_wout_file, DMerc

      USE safe_open_mod

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
      INTEGER     :: ik, ier, iCount, diagUnit, istat
      REAL(rprec) :: dmercier(nsd)
      character(256) :: fname
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      diagUnit = 223
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
        IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'DMER ', ik, 4
        IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DMER  k'
      IF (niter >= 0) THEN
        fname = 'dmerc_out' // trim(proc_string)
        istat=0
        ! Data are output to separated file for easier diagnostic
        CALL safe_open(diagUnit,istat,fname,'replace','formatted')
        IF (istat .ne. 0) THEN
         WRITE(6,*) 'Error opening output file:',TRIM(fname), istat
         iflag=-1
         RETURN
        END IF
        ! read Dmercier data
        CALL read_wout_file (proc_string, ier)
        DO ik = 1, nsd
          IF (sigma(ik) >= bigno) CYCLE
          mtargets = mtargets + 1
          targets(mtargets) = target(ik)
          sigmas(mtargets)  = sigma(ik)
          vals(mtargets)    = DMerc(ik)
            !write(*,*) 'rho(ik)' , rho(ik), 'mshear(ik): ',  DMerc(ik)
          IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') target(ik), sigma(ik), DMerc(ik), ik
          IF (iflag == 1) THEN
              WRITE(diagUnit,'(2ES22.12E3)') rho(ik), DMerc(ik)
          END IF
        END DO

        CALL FLUSH(6)
        CALL FLUSH(diagUnit)
        CLOSE(diagUnit)
      ELSE
         DO ik = 1, nsd
            ! define the radial positions at which the boozxform is run (lbooz = true))
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_dmer
            END IF
         END DO
      END IF
      !CALL DEALLOCATE_READ_WOUT
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_dmercier
