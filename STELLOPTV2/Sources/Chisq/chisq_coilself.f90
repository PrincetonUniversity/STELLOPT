!-----------------------------------------------------------------------
!     Subroutine:    chisq_coilself
!     Authors:       J. Breslau (jbreslau@pppl.gov)
!     Date:          8/28/2017
!     Description:   Calculates the number of places each coil intersects
!                    itself on a winding surface.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coilself(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: lcoil_spline, lwindsurf
      
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
      INTEGER     :: ik, nselfint

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (.NOT.lwindsurf) RETURN  !This routine only makes sense for embedded curves in 2D
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'SELF-INTS ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL'
      IF (niter >= 0) THEN
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               CALL get_coil_self_int(ik, nselfint)
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = REAL(nselfint, rprec)
               IF (iflag == 1) WRITE(iunit_out,'(I5,2ES22.12E3,I5)') ik,target(ik),sigma(ik),nselfint
            END IF
         END DO
      ELSE
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coilself
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilself
