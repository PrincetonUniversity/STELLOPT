!-----------------------------------------------------------------------
!     Subroutine:    chisq_xmcts_bright
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference measured and equilibrium
!                    XMCTS signals.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_xmcts_bright(target,sigma,niter,iflag)
        !-----------------------------------------------------------------------
        !     Libraries
        !-----------------------------------------------------------------------
              USE stellopt_runtime
              USE stellopt_targets
              USE equil_utils
              
        !-----------------------------------------------------------------------
        !     Input/Output Variables
        !
        !-----------------------------------------------------------------------
              IMPLICIT NONE
              REAL(rprec), INTENT(in)    ::  target(nprof)
              REAL(rprec), INTENT(in)    ::  sigma(nprof)
              INTEGER,     INTENT(in)    ::  niter
              INTEGER,     INTENT(inout) ::  iflag
              
        !-----------------------------------------------------------------------
        !     Local Variables
        !        lreset_s    Gets set to true if using R,PHI,Z Specification
        !        ik          Dummy index
        !        ti_val      Holds profile evaulation
        !-----------------------------------------------------------------------
              INTEGER ::  ik
              REAL(rprec) :: xmcts_val, xmcts_length, etendu
              REAL(rprec) :: x0(3), x1(3)
        !----------------------------------------------------------------------
        !     BEGIN SUBROUTINE
        !----------------------------------------------------------------------
              IF (iflag < 0 ) RETURN
              ik = COUNT(sigma < bigno)
              IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'XMCTS_BRIGHT',ik,10
              IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  EQUIL  R0  PHI0  Z0  R1  PHI1  Z1  ETENDU'
              IF (niter >= 0) THEN
                 PRINT *,'integrating sxr emissivity with etendu'
                 DO ik = 1, nprof
                    IF (sigma(ik) >= bigno) CYCLE
                    x0(1)=r0_xmcts(ik); x1(1)=r1_xmcts(ik)
                    x0(2)=phi0_xmcts(ik); x1(2)=phi1_xmcts(ik)
                    x0(3)=z0_xmcts(ik); x1(3)=z1_xmcts(ik)
                    etendu = etendu_xmcts(ik)
                    xmcts_val = 0.0
                    CALL line_int(fcn_xmcts_bright,x0,x1,xmcts_val,LENGTH=xmcts_length)
                    xmcts_val = xmcts_val * etendu
                    !PRINT *,'XMCTS_LENGTH (',ik,')',xmcts_length
                    mtargets = mtargets + 1
                    targets(mtargets) = target(ik)
                    sigmas(mtargets)  = sigma(ik)
                    vals(mtargets)    = xmcts_val
                    IF (iflag == 1) WRITE(iunit_out,'(10ES22.12E3)') target(ik),sigma(ik),xmcts_val,x0(1),x0(2),x0(3),x1(1),x1(2),x1(3),etendu
                 END DO
              ELSE
                 DO ik = 1, nprof
                    IF (sigma(ik) < bigno) THEN
                       mtargets = mtargets + 1
                       IF (niter == -2) target_dex(mtargets) = jtarget_xmcts_bright
                    END IF
                 END DO
              END IF
              RETURN
        !----------------------------------------------------------------------
        !     END SUBROUTINE
        !----------------------------------------------------------------------
              END SUBROUTINE chisq_xmcts_bright
        