!-----------------------------------------------------------------------
!     Subroutine:    chisq_coiltorvar
!     Authors:       J. Breslau (jbreslau@pppl.gov)
!     Date:          1/14/2019
!     Description:   Calculates the difference between RMS toroidal angle
!                    excursion for each coil and its target value.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coiltorvar(target,sigma,niter,iflag)
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
      REAL(rprec), INTENT(in)    ::  target(nigroup)
      REAL(rprec), INTENT(in)    ::  sigma(nigroup)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER     :: ik
      REAL(rprec) :: rmsdv
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'RMS COIL TOROIDAL EXCURSION ',ik,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               CALL get_coil_tor_excur(ik, thwt_coiltorvar(ik), rmsdv)
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = MAX(rmsdv, target(ik))  ! One-sided barrier
               IF (iflag == 1) WRITE(iunit_out,'(I5,3ES22.12E3)') ik,target(ik),sigma(ik),rmsdv
            END IF
         END DO
      ELSE
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coiltorvar
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coiltorvar
