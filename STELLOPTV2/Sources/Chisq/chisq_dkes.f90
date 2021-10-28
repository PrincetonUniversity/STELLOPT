!-----------------------------------------------------------------------
!     Subroutine:    chisq_dkes
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This is the chisquared routine which compares
!                    target to equilibrium DKES values.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_dkes(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils, ONLY: get_equil_phi,shat, phi_type
!DEC$ IF DEFINED (DKES_OPT)
      USE dkes_realspace, ONLY: DKES_L11p, DKES_L33p, DKES_L31p, &
                                DKES_L11m, DKES_L33m, DKES_L31m, &
                                DKES_scal11, DKES_scal33, DKES_scal31
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
      INTEGER :: ik, ier_phi
      REAL(rprec) :: dkes_efield, phi_temp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'DKES ',ik,15
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  S  NU  ER  L11p  L11m  L33p  L33m  L31p  L31m  SCAL11  SCAL33  SCAL31'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            targets(mtargets) = target_dkes(ik)
            sigmas(mtargets)  = sigma_dkes(ik)
!DEC$ IF DEFINED (DKES_OPT)
            CALL get_equil_phi(shat(ik),phi_type,phi_temp,ier_phi,dkes_efield)
            IF (ier_phi < 0) dkes_efield = 0.0
            vals(mtargets)    = 0.5*(DKES_L31p(ik)+DKES_L31m(ik))
            IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
               target(ik),sigma(ik),vals(mtargets),&
               shat(ik), nu_dkes(ik), dkes_efield, &
               DKES_L11p(ik),DKES_L11m(ik),DKES_L33p(ik),&
               DKES_L33m(ik),DKES_L31p(ik),DKES_L31m(ik),&
               DKES_scal11(ik),DKES_scal33(ik),DKES_scal31(ik)
!DEC$ ELSE
            vals(mtargets) = target_dkes(ik)
            IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
               target(ik),sigma(ik),vals(mtargets),&
               shat(ik), nu_dkes(ik), dkes_efield, &
               0.0, 0.0, 0.0,&
               0.0, 0.0, 0.0,&
               0.0, 0.0, 0.0
!DEC$ ENDIF
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_dkes
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_dkes
