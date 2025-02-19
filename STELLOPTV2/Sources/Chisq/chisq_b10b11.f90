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
!        lfound   Logical indicating that all harmonic indexes are found
!        ik       Radial index helper
!        mn       Harmonic index helper
!        bXX      B(m,n) values
!-----------------------------------------------------------------------
      LOGICAL :: lfound
      INTEGER :: ik, mn, mn00, mn01, mn11, mn10
      REAL(rprec) :: b00, b01, b11, b10
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'B10B11 ',ik,8
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  B00  B01  B10  B11  K'
      IF (niter >= 0) THEN
         ! First get the indices on the mn grid
         mn00 = -1; mn01 = -1; mn10 = -1; mn11 = -1; lfound = .FALSE.
         DO mn = 1, mnboz_b
            IF ((ixn_b(mn)/nfp_b == 0) .and. (ixm_b(mn) == 0)) mn00 = mn
            IF ((ixn_b(mn)/nfp_b == 1) .and. (ixm_b(mn) == 0)) mn01 = mn
            IF ((ixn_b(mn)/nfp_b == 0) .and. (ixm_b(mn) == 1)) mn10 = mn
            IF ((ixn_b(mn)/nfp_b == 1) .and. (ixm_b(mn) == 1)) mn11 = mn
         END DO
         IF ((mn00 > 0) .or. (mn01 > 0) .or. (mn10 > 0) .or. (mn11 > 0)) lfound = .TRUE.
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            b00 = 0.0; b10 = 0.0; b01 = 0.0; b11 = 1.0
            IF (lfound) THEN
               b00 = bmnc_b(mn00,ik)
               b01 = bmnc_b(mn01,ik)
               b10 = bmnc_b(mn10,ik)
               b11 = bmnc_b(mn11,ik)
            END IF
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets) = ABS(b10/b11) ! The ABS is becasue I don't care about the sign. (SAL)
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3,2X,I3.3)') target(ik),sigmas(mtargets),vals(mtargets),b00,b01,b10,b11,ik
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
