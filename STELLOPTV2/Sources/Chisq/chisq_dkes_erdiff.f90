!-----------------------------------------------------------------------
!     Subroutine:    chisq_dkes_erdiff
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          12/19/2022
!     Description:   This is the chisquared which handles the proxy
!                    for temperature screening using DKES.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_dkes_erdiff(target,sigma,niter,iflag)
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
      INTEGER :: ik, ii
      REAL(rprec) :: fp, fm
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ! Count values
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'DKES_ERDIFF ',ik,25
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  S  NU  ERp ERm '&
                // ' L11p+  L11m+  L33p+  L33m+  L31p+  L31m+  SCAL11+  SCAL33+  SCAL31+'& 
                // ' L11p-  L11m-  L33p-  L33m-  L31p-  L31m-  SCAL11-  SCAL33-  SCAL31-'
      IF (niter >= 0) THEN
         ik = 0
         DO ii = 1, nsd
            IF (sigma(ii) >= bigno) CYCLE
            mtargets = mtargets + 1
            targets(mtargets) = target(ii)
            sigmas(mtargets)  = sigma(ii)
!DEC$ IF DEFINED (DKES_OPT)
            ik = ik + 1
            fp = 0.5*(DKES_L11p(ik)+DKES_L11m(ik))
            ik = ik + 1
            fm = 0.5*(DKES_L11p(ik)+DKES_L11m(ik))
            vals(mtargets) = fp - fm
            IF (iflag == 1) WRITE(iunit_out,'(25ES22.12E3)') &
                  targets(mtargets),sigmas(mtargets),vals(mtargets),&
                  shat(ii), nu_dkes_Erdiff, Ep_dkes_Erdiff, &
                  Em_dkes_erdiff,&
                  DKES_L11p(ik-1),DKES_L11m(ik-1),DKES_L33p(ik-1),&
                  DKES_L33m(ik-1),DKES_L31p(ik-1),DKES_L31m(ik-1),&
                  DKES_scal11(ik-1),DKES_scal33(ik-1),DKES_scal31(ik-1), &
                  DKES_L11p(ik),DKES_L11m(ik),DKES_L33p(ik),&
                  DKES_L33m(ik),DKES_L31p(ik),DKES_L31m(ik),&
                  DKES_scal11(ik),DKES_scal33(ik),DKES_scal31(ik)

!DEC$ ELSE
            vals(mtargets) = target(ii)
            IF (iflag == 1) WRITE(iunit_out,'(25ES22.12E3)') &
                  targets(mtargets),sigmas(mtargets),vals(mtargets),&
                  shat(ii), nu_dkes_Erdiff, Ep_dkes_Erdiff, &
                  Em_dkes_erdiff,&
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0,&
                  0.0, 0.0, 0.0
!DEC$ ENDIF
         END DO
      ELSE
         nruns_dkes = 0
         DO ii = 1, nsd
            IF (sigma(ii) >= bigno) CYCLE
            lbooz(ii) = .TRUE.
            mtargets = mtargets + 1
            nruns_dkes = nruns_dkes + 2
            E_dkes(1) = Ep_DKES_Erdiff
            E_dkes(2) = Em_DKES_Erdiff
            nu_dkes(1) = nu_dkes_Erdiff
            nu_dkes(2) = nu_dkes_Erdiff
            IF (niter == -2) target_dex(mtargets)=jtarget_dkes_erdiff
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_dkes_erdiff
