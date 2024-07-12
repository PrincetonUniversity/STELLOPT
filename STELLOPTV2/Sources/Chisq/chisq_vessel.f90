!-----------------------------------------------------------------------
!     Subroutine:    chisq_vessel
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate the distance to the vessel wall.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_vessel(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils
      USE vessel_mod
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ier         Error Flag
!        dex         Length of ne_aux_s array
!        ne_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  loutside = .false.
      INTEGER ::  ik, ier, dex, u, v, iout
      REAL(rprec) :: dist, dist_max, dist_min, theta, zeta, R_eq, Z_eq
      INTEGER, PARAMETER :: nu_eq=30
      INTEGER, PARAMETER :: nv_eq=30
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'VESSEL ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  DIST'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         ! Default to chisq=0
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)    = target
         dist = 0.0
         dist_max = target
         dist_min = bigno
         DO u = 1, nu_eq
            DO v = 1, nv_eq
               theta = pi2*(u-1)/nu_eq
               zeta  = pi2*(v-1)/nv_eq
               CALL get_equil_RZ(1.0_rprec,theta,zeta,R_eq,Z_eq,ier)
               CALL vessel_dist_cyl(R_eq,zeta/nfp,Z_eq,dist,iout)
               IF (iout <= 0.0) THEN
                  loutside=.true.
                  dist_max = dist
               ELSE IF (dist < dist_min) THEN
                  dist_min=dist
               END IF
            END DO
         END DO
         IF (loutside) THEN
            vals(mtargets) = dist_max
         ELSE IF (dist_min < target) THEN
            vals(mtargets) = target-dist_min
         END IF
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,vals(mtargets)
      ELSE
         IF (sigma < bigno) THEN
            mtargets = mtargets + 1
            IF (niter == -1) CALL vessel_load_txt(TRIM(vessel_string),iflag)
            IF (niter == -2) target_dex(mtargets) = jtarget_vessel
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_vessel
