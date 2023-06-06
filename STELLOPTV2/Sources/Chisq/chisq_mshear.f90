!-----------------------------------------------------------------------
!     Subroutine:    chisq_mshear
!     Authors:       E. Sanchez (edi.sanchez@ciemat.es))
!     Date:          15/11/2020
!     Description:   This subroutine calculates the magnetic shear target contribution to Chisquare
!-----------------------------------------------------------------------
      SUBROUTINE chisq_mshear(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE equil_utils
      USE stellopt_targets
      USE equil_vals, ONLY: Baxis

      USE safe_open_mod
      USE read_boozer_mod
      USE EZspline_obj
      USE EZspline

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
      REAL(rprec) :: msh(nsd)
      DOUBLE PRECISION :: iota_v, s_val
!      TYPE(EZspline1_r8) :: iot_spl
      character(256) :: fname
      REAL(rprec) :: iota, iotaprime
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      diagUnit = 223
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
        IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'MSHEAR ',ik, 4
        IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA   MSHEAR  k'
      IF (niter >= 0) THEN
        fname = 'msh_out' // trim(proc_string)
        istat=0
        ! Data are output to separated file for easier diagnostic
        CALL safe_open(diagUnit,istat,fname,'replace','formatted')
        IF (istat .ne. 0) THEN
         WRITE(6,*) 'Error opening output file:',TRIM(fname), istat
         iflag=-1
         RETURN
        END IF
        DO ik = 1, nsd
          IF (sigma(ik) >= bigno) CYCLE
          s_val = rho(ik) * rho(ik)
          CALL EZspline_isInDomain(iota_spl, s_val, ier)
          IF (ier .ne. 0) RETURN
          !write(*,*) 'mag.shear: '
          CALL EZspline_interp(iota_spl, s_val, iota_v, ier)
          !write(*,*) 'mag.shear: 2'
          CALL EZspline_derivative(iota_spl, 1, s_val, iotaprime, ier)
          msh(ik) =  2  * rho(ik) * iotaprime /iota_v
          !write(*,*) 'mag.shear: ', msh(ik)
          mtargets = mtargets + 1
          targets(mtargets) = target(ik)
          sigmas(mtargets)  = sigma(ik)
          vals(mtargets)    = msh(ik)
            !write(*,*) 'rho(ik)' , rho(ik), 'mshear(ik): ',  msh(ik)
          IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') target(ik), sigma(ik), msh(ik), ik
          IF (iflag == 1) THEN
              WRITE(diagUnit,'(2ES22.12E3)') rho(ik), msh(ik)
          END IF
        END DO
        !CALL EZspline_free(iota_spl,ier)

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
               IF (niter == -2) target_dex(mtargets)=jtarget_mshear
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_mshear
