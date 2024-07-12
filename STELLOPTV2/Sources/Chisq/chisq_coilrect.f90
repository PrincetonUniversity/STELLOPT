!-----------------------------------------------------------------------
!     Subroutine:    chisq_coilrect
!     Authors:       J. Breslau (jbreslau@pppl.gov)
!     Date:          6/4/2019
!     Description:   Calculates penalty based on excursion of coil from
!                    prescribed v bounds
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coilrect(target,sigma,niter,iflag)
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
      REAL(rprec) :: vmin, vmax
      INTEGER     :: ik

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (.NOT.ANY(lwindsurf)) RETURN  !This routine only makes sense for embedded curves in 2D
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'VBOUNDS ',ik,5
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VMIN  VMAX'
      IF (niter >= 0) THEN
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               CALL get_coil_vbounds(ik, coilrectduu(ik), coilrectdul(ik), vmin, vmax)
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = (MIN(vmin - coilrectvmin(ik), 0.0)/coilrectpfw)**2
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = (MAX(vmax - coilrectvmax(ik), 0.0)/coilrectpfw)**2
               IF (iflag == 1) WRITE(iunit_out,'(I5,4ES22.12E3)') ik,target(ik),sigma(ik),vmin,vmax
            END IF
         END DO
      ELSE
         DO ik = 1, nigroup
            IF ((sigma(ik) < bigno) .AND. ANY(lcoil_spline(ik,:))) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coilrect
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_coilrect
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coilrect
