!-----------------------------------------------------------------------
!     Subroutine:    chisq_b10b11
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          02/15/2025
!     Description:   Calculate the Boozer B10/B11 ratio (Bmn) which is a
!                    proxy for bootstrap current.
!                    http://aip.scitation.org/doi/10.1063/1.860843
!-----------------------------------------------------------------------
      SUBROUTINE chisq_b10b11(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
      USE read_boozer_mod
      USE vmec_input, ONLY: mpol, ntor
      
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
      INTEGER :: ik, mn
      REAL(rprec) :: b00,b01,b11,b10
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'B10B11 ',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  B00  B01  B10  B11'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            b00 = 0.0; b10 = 0.0; b01 = 0.0; b11 = 1.0
            DO mn = 1, mnboz_b
               IF ((ixn_b(mn)/nfp_b == 0) .and. (ixm_b(mn) == 0)) b00 = bmnc_b(mn,ik)
               IF ((ixn_b(mn)/nfp_b == 1) .and. (ixm_b(mn) == 0)) b01 = bmnc_b(mn,ik)
               IF ((ixn_b(mn)/nfp_b == 0) .and. (ixm_b(mn) == 1)) b10 = bmnc_b(mn,ik)
               IF ((ixn_b(mn)/nfp_b == 1) .and. (ixm_b(mn) == 1)) b11 = bmnc_b(mn,ik)
            END DO
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets) = b10/b11
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') target(ik),sigmas(mtargets),vals(mtargets),b00,b01,b10,b11
         END DO
      ELSE
         ! Consistency check
         mboz = MAX(6*mpol, 2, mboz)             
         nboz = MAX(2*ntor-1, 0, nboz)  
         ! CALCULATE mnboz_b becasue we don't know it yet (setup_booz.f)
         mnboz_b = (2*nboz+1)*(mboz-1) + (nboz + 1)   
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_b10b11
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_b10b11
