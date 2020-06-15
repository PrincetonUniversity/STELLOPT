!-----------------------------------------------------------------------
!     Subroutine:    chisq_coilcrv
!     Authors:       J. Breslau (jbreslau@pppl.gov)
!     Date:          8/24/2017
!     Description:   Calculates the difference between maximum curvature
!                    for each coil and its target value.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coilcrv(target,sigma,niter,iflag)
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

      INTERFACE
         SUBROUTINE get_coil_maxcurv(icoil, maxcurv, s_max, u_max, v_max, uout)
           USE stel_kinds, ONLY: rprec
           INTEGER, INTENT(IN)           :: icoil
           REAL(rprec), INTENT(OUT)      :: maxcurv, s_max, u_max, v_max
           INTEGER, INTENT(IN), OPTIONAL :: uout
         END SUBROUTINE get_coil_maxcurv
      END INTERFACE
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER     :: ik
      REAL(rprec) :: maxcurv, smax, umax, vmax
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'MAX_COIL_CURVATURE ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  LOC  U  V'
      IF (niter >= 0) THEN
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               CALL get_coil_maxcurv(ik, maxcurv, smax, umax, vmax)
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = MAX(maxcurv, target(ik))  ! One-sided barrier
               IF (iflag == 1) WRITE(iunit_out,'(I5,6ES22.12E3)') &
                    ik,target(ik),sigma(ik),maxcurv,smax,umax,vmax
            END IF
         END DO
      ELSE
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coilcrv
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilcrv
