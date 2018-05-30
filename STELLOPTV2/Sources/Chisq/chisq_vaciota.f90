!-----------------------------------------------------------------------
!     Subroutine:    chisq_iota
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/12/2012
!     Description:   Calculate difference between equilbrium rotational
!                    transform and target rotational transform.
!                    Although iota is not prescribed at real space
!                    points we include it for flexibility.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_vaciota(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nprof)
      REAL(rprec), INTENT(in)    ::  sigma(nprof)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        iota_val    Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  ik, ier
      REAL(rprec) :: iota_val, s11, s12, s21, s22
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'VACIOTA ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  VACIOTA'
      IF (niter >= 0) THEN
         IF (ANY(s_vaciota >= 0)) lreset_s = .false.
         ! GET s if necessary
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ier = 0
            IF (lreset_s) THEN
               CALL get_equil_s(r_vaciota(ik),phi_vaciota(ik),z_vaciota(ik),s_vaciota(ik),ier)
            END IF
            IF (s_vaciota(ik) <= 1.0 .and. s_vaciota(ik) >= 0.0 .and. ier == 0) THEN
               ier = 0
               CALL get_equil_sus(s_vaciota(ik),s11,s12,s21,s22,ier)
               iota_val = -s12/s11
            ELSE
               iota_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = iota_val
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_vaciota(ik),phi_vaciota(ik),z_vaciota(ik),s_vaciota(ik),target(ik),sigma(ik),iota_val
         END DO
         IF (lreset_s) s_vaciota(:) = -1.0
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_vaciota
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_vaciota
