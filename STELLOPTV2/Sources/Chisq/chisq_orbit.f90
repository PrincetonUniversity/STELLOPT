!-----------------------------------------------------------------------
!     Subroutine:    chisq_orbit
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          10/24/2014
!     Description:   This subroutine calculates chisqared for particle
!                    orbits.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_orbit(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
!DEC$ IF DEFINED (BEAMS3D_OPT)
      USE beams3d_input_mod, ONLY: read_beams3d_input, r_start_in
!DEC$ ENDIF
      
      
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
      LOGICAL :: mult_par = .true.
      INTEGER :: iunit, ik, i, j, nu, nv
      REAL(rprec) :: val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'ORBIT ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  LOST_FRAC  S '
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            IF (ALLOCATED(orbit_lost_frac)) THEN
               val = orbit_lost_frac(ik)
            ELSE
               val = bigno
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = val
            IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ik),sigma(ik),vals(mtargets),rho(ik)
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_orbit
         END DO
         iflag=0
!DEC$ IF DEFINED (BEAMS3D_OPT)
         r_start_in(1) = 1.0 ! So that the beams aren't read.
         CALL read_beams3d_input('input.'//TRIM(id_string), iflag)
!DEC$ ENDIF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_orbit
