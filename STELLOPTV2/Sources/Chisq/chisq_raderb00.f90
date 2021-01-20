!-----------------------------------------------------------------------
!     Subroutine:    chisq_raderb00
!     Authors:       E. Sanchez (edi.sanchez@ciemat.es))
!     Date:          15/11/2020
!     Description:   This subroutine calculates the radial derivative of B00 (in Boozer coordinates)
!-----------------------------------------------------------------------
      SUBROUTINE chisq_raderb00(target,sigma,niter,iflag)
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
      REAL(rprec) :: raderb00(nsd)
      DOUBLE PRECISION :: s_val(nsd), b00(nsd)
      TYPE(EZspline1_r8) :: b00prof_spl
      character(256) :: fname
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      diagUnit = 222
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
!        IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'RADERB0 ',ik, 3
!        IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA dB_00/ds'
      IF (niter >= 0) THEN
        fname = 'rdB00_out' // trim(proc_string)
        istat=0
        ! Data are output to separated file for easier diagnostic
        CALL safe_open(diagUnit,istat,fname,'replace','formatted')
        IF (istat .ne. 0) THEN
         WRITE(6,*) 'Error opening output file:',TRIM(fname), istat
         iflag=-1
         RETURN
        END IF
        CALL read_boozer_file (proc_string, ier)
        iCount =0
        DO ik = 1, nsd
            ! positions at which the boozer transform was carried out
            IF ( lbooz(ik) .eq. .true.) THEN
                iCount=iCount+1
                s_val(iCount) = rho(ik)
                b00(iCount) = bmnc_b(1, ik)
                !write(*,*) 's_val(ik): ', s_val(iCount)
            END IF
        END DO
!            write(*,*) 'lbooz(ik) .eq. .true.: ', iCount
        CALL setup_prof_spline(b00prof_spl, iCount, s_val, b00, ier)
        DO ik = 1, nsd
          IF (sigma(ik) >= bigno) CYCLE

          CALL EZspline_isInDomain(b00prof_spl, rho(ik), ier)
          IF (ier .ne. 0) RETURN
          CALL EZspline_derivative(b00prof_spl, 1, rho(ik), raderb00(ik), ier)
          !write(*,*) 'derivative: ', raderb00

        !  normalization
          raderb00(ik) = raderb00(ik)/ Baxis
          mtargets = mtargets + 1
          targets(mtargets) = target(ik)
          sigmas(mtargets)  = sigma(ik)
          vals(mtargets)    = raderb00(ik)
        !            write(*,*) 'rho(ik)' , rho(ik), 'raderb00(ik): ',  raderb00(ik), 'b00(ik)', b00(ik)
        !            IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target(ik), sigma(ik), raderb00(ik)
          IF (iflag == 1) THEN
              WRITE(diagUnit,'(2ES22.12E3)') rho(ik), raderb00(ik)
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
               IF (niter == -2) target_dex(mtargets)=jtarget_raderb00
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_raderb00
