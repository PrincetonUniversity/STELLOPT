!-----------------------------------------------------------------------
!     Subroutine:    chisq_coillen
!     Authors:       J. Breslau (jbreslau@pppl.gov)
!     Date:          8/22/2017
!     Description:   Calculates the difference between the length of each 
!                    coil and its target value.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coillen(target_l,sigma_l,target_v,sigma_v,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: lcoil_spline

!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target_l(nigroup), target_v(nigroup)
      REAL(rprec), INTENT(in)    ::  sigma_l(nigroup), sigma_v(nigroup)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag

!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER     :: ik
      REAL(rprec) :: local_len, rmssegvar
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT((sigma_l < bigno).OR.(sigma_v < bigno))
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'LEN_SEG_SD ',ik, 7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         DO ik = 1, nigroup
            IF (((sigma_l(ik) < bigno).OR.(sigma_v(ik) < bigno)) .AND. &
                 ANY(lcoil_spline(ik,:))) THEN
               CALL get_coil_length(ik, local_len, rmssegvar)
               IF (sigma_l(ik) < bigno) THEN
                  mtargets = mtargets + 1
                  targets(mtargets) = target_l(ik)
                  sigmas(mtargets)  = sigma_l(ik)
                  vals(mtargets)    = local_len
               END IF
               IF (sigma_v(ik) < bigno) THEN
                  mtargets = mtargets + 1
                  targets(mtargets) = target_v(ik)
                  sigmas(mtargets)  = sigma_v(ik)
                  vals(mtargets)    = rmssegvar
               END IF
               IF (iflag == 1) WRITE(iunit_out,'(I5,6ES22.12E3)') ik,&
                    target_l(ik), sigma_l(ik), local_len, target_v(ik), sigma_v(ik), rmssegvar
            END IF
         END DO
      ELSE
         DO ik = 1, nigroup
            IF (ANY(lcoil_spline(ik,:))) THEN
               IF (sigma_l(ik) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets)=jtarget_coillen
               END IF
               IF (sigma_v(ik) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets)=jtarget_coilsegvar
               END IF
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coillen
