!-----------------------------------------------------------------------
!     Subroutine:    chisq_dkes_31
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/04/2025
!     Description:   Chisquared for the DKES L31 coefficient.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_dkes_31(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals, ONLY: shat
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
      INTEGER :: ik, ij, ii
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ! Print Header
      IF (iflag == 1) THEN
         ik   = COUNT(target_dex == jtarget_dkes_31)
         WRITE(iunit_out,'(A,2(2X,I3.3))') 'DKES_31 ',ik,15
         WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  S  NU  ER  L11p  L11m  L33p  L33m  L31p  L31m  SCAL11  SCAL33  SCAL31'
      END IF
      IF (niter >= 0) THEN
         ik = 0 ! Always the first one
         DO ii = 1, nsd
            IF (sigma(ii) >= bigno) CYCLE
            DO ij = 1, nprof
               IF (E_dkes(ij) <= -bigno .or. nu_dkes(ij) <= -bigno) CYCLE
               ik = ik + 1
               mtargets = mtargets + 1
               targets(mtargets) = target_dkes_31(ii)
               sigmas(mtargets)  = sigma_dkes_31(ii)
!DEC$ IF DEFINED (DKES_OPT)
               vals(mtargets)    = 0.5*(DKES_L31p(ik)+DKES_L31m(ik))
               IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
                  targets(mtargets),sigmas(mtargets),vals(mtargets),&
                  shat(ii), nu_dkes(ij), E_dkes(ij), &
                  DKES_L11p(ik),DKES_L11m(ik),DKES_L33p(ik),&
                  DKES_L33m(ik),DKES_L31p(ik),DKES_L31m(ik),&
                  DKES_scal11(ik),DKES_scal33(ik),DKES_scal31(ik)
!DEC$ ELSE
               vals(mtargets) = target_dkes_31(ii)
               IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
                  targets(mtargets),sigmas(mtargets),vals(mtargets),&
                  shat(ii), nu_dkes(ij), E_dkes(ij), &
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0
!DEC$ ENDIF
            END DO
         END DO
      ELSE
         DO ii = 1, nsd
            IF (sigma(ii) >= bigno) CYCLE
            lbooz(ii) = .TRUE.
            lneed_dkes(ii) = .TRUE.
            DO ij = 1, nprof
               IF (E_dkes(ij) <= -bigno .or. nu_dkes(ij) <= -bigno) CYCLE
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_dkes_31
            END DO
         END DO
         ! SIGMA DKES may have already counted the nruns
         DO ii = 1, nsd
            IF ((sigma_dkes_11(ii) < bigno) .or. (sigma(ii) >= bigno)) CYCLE
               DO ij = 1, nprof
                  IF (E_dkes(ij) <= -bigno .or. nu_dkes(ij) <= -bigno) CYCLE
                  nruns_dkes = nruns_dkes + 1
               END DO
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_dkes_31
