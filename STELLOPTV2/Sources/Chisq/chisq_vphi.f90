!-----------------------------------------------------------------------
!     Subroutine:    chisq_vphi
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference between equilbrium toroidal
!                    rotation and target toroidal rotation.
!                    Per Tony Cooper
!                    v_phi= omega(s) * R(s,u,v)*c_s*sqrt(machsq) / R_0
!                    c_s=sqrt(T*q/m) if T in [eV]
!-----------------------------------------------------------------------
      SUBROUTINE chisq_vphi(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars
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
!        ti_val      Holds profile evaulation
!-----------------------------------------------------------------------
      LOGICAL ::  lreset_s = .true.
      INTEGER ::  ik, ier, dex
      REAL(rprec) :: vphi_val,omega_val,ti_val,sound_speed,s_temp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'VPHI ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'R  PHI  Z  S  TARGET  SIGMA  VPHI'
      IF (niter >= 0) THEN
         IF (ANY(s_vphi > 0)) lreset_s = .false.
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            ! GET s if necessary
            IF (lreset_s) THEN
               ier = 0
               CALL get_equil_s(r_vphi(ik),phi_vphi(ik),z_vphi(ik),s_vphi(ik),ier)
            END IF
            IF (s_vphi(ik) <= 1.0 .and. s_vphi(ik) >= 0.0) THEN
               ier = 0
               ! WILL need to get R to calc VPHI from omega
               CALL get_equil_omega(s_vphi(ik),omega_val,ier)
               s_temp = 0.0
               CALL get_equil_ti(s_temp,TRIM(ti_type),ti_val,ier)
               sound_speed = SQRT(ti_val*qm_ratio)  ! Assume Ti in eV
               vphi_val = sound_speed*mach0*omega_val*r_vphi(ik)/rmajor  ! omega = vphi/r
            ELSE
               vphi_val = 0.0
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = vphi_val
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') r_vphi(ik),phi_vphi(ik),z_vphi(ik),s_vphi(ik),target(ik),sigma(ik),vphi_val
         END DO
         IF (lreset_s) s_vphi(:) = -1.
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_vphi
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_vphi
