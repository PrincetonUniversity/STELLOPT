!-----------------------------------------------------------------------
!     Subroutine:    chisq_mbm
!     Authors:       J. L. Velasco (joseluis.velasco@ciemat.es) based on chisq_raderb00 by E. Sanchez
!     Date:          01/04/2022
!     Description:   This subroutine calculates the maximum value of B minus its minimum value at s=1
!-----------------------------------------------------------------------
      SUBROUTINE chisq_mbm(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE equil_utils
      USE stellopt_targets
      USE equil_vals, ONLY: Baxis

      USE safe_open_mod
      USE read_boozer_mod

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
      INTEGER     :: ik, ier, diagUnit, istat, iz, it, imn
      REAL(rprec) :: zeta, theta, dz, dt, Btemp, bmin, bmax(nsd), mbm(nsd)
      character(256) :: fname
      INTEGER, PARAMETER :: npoints=64
      REAL(rprec), PARAMETER :: twopi = 6.283185307179586476925286766559D0
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
      !----------------------------------------------------------------------

      diagUnit = 223
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
      IF (niter >= 0) THEN
        fname = 'mBm_out' // trim(proc_string)
        istat=0
        ! Data are output to separated file for easier diagnostic
        CALL safe_open(diagUnit,istat,fname,'replace','formatted')
        IF (istat .ne. 0) THEN
         WRITE(6,*) 'Error opening output file:',TRIM(fname), istat
         iflag=-1
         RETURN
        END IF
        CALL read_boozer_file (proc_string, ier)
        dz=twopi/npoints/nfp_b
        dt=twopi/npoints
        DO ik = 1, nsd
            ! positions at which the boozer transform was carried out
           IF (lbooz(ik)) THEN
              bmax(ik)=bmnc_b(1, ik)
              bmin=bmax(ik)
              DO iz=1,npoints
                 zeta=(iz-1)*dz
                 DO it=1,npoints
                    theta=(it-1)*dt
                    Btemp=0
                    DO imn=1,mnboz_b
                       Btemp=Btemp+bmnc_b(imn, ik)*cos(ixm_b(imn)*theta-ixn_b(imn)*zeta)
                    END DO
                    IF(Btemp.GT.bmax(ik)) THEN
                       bmax(ik)=Btemp
                    ELSE IF (Btemp.LT.bmin) THEN
                       bmin=Btemp
                    END IF
                 END DO
              END DO
           END IF

        END DO
        
        DO ik = 1, nsd
           IF (sigma(ik) >= bigno) CYCLE	
           !  normalization
           mbm(ik) = bmax(ik)/bmin
           mtargets = mtargets + 1
           targets(mtargets) = target(ik)
           sigmas(mtargets)  = sigma(ik)
           vals(mtargets)    = mbm(ik)
           
           IF (iflag == 1) THEN
              WRITE(diagUnit,'(2ES22.12E3)') rho(ik), mbm(ik)
           END IF

        END DO

        CALL FLUSH(6)
        CALL FLUSH(diagUnit)
        CLOSE(diagUnit)
      ELSE
         DO ik = 1, nsd
            ! define the radial positions at which the boozxform is run (lbooz = true))
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               IF(ik.GT.2) THEN
                  lbooz(ik-1) = .TRUE.
               END IF
               IF(ik.LT.nsd) THEN
                  lbooz(ik+1) = .TRUE.
               END IF
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_mbm
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_mbm
