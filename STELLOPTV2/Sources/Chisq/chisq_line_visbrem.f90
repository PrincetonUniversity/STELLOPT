!-----------------------------------------------------------------------
!     Subroutine:    chisq_visbrem_zeff
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          04/18/2021
!     Description:   Calculate line integrated visual bremsstrahlung.
!                    For details of calculation see:
!                       Pavone, A., Hergenhahn, U., Krychowiak, 
!                       M., Hoefel, U., Kwak, S., Svensson, J., et al.
!                       (2019). Journal of Instrumentation, 14(10),
!                       C10003–C10003. 
!                       http://doi.org/10.1088/1748-0221/14/10/C10003
!                    Note that if signal is normalized to the
!                    calibration factor then calib_visbrem_line should
!                    be set to 1.0.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_line_visbrem(target,sigma,niter,iflag)
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
!        ne_val      Holds profile evaulation
!-----------------------------------------------------------------------
      INTEGER ::  ik
      REAL(rprec) :: visbrem_val, visbrem_length
      REAL(rprec) :: x0(3), x1(3)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'VISBREMLINE',ik,11
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VISBREMLINE LAMBDA CALIB  R0  PHI0  Z0  R1  PHI1  Z1'
      IF (niter >= 0) THEN
         DO ik = 1, nprof
            IF (sigma(ik) >= bigno) CYCLE
            x0(1)=r0_visbrem_line(ik); x1(1)=r1_visbrem_line(ik)
            x0(2)=phi0_visbrem_line(ik); x1(2)=phi1_visbrem_line(ik)
            x0(3)=z0_visbrem_line(ik); x1(3)=z1_visbrem_line(ik)
            visbrem_val = 0.0
            visbrem_lambda = lambda_visbrem_line(ik)
            CALL line_int(fcn_bremsstrahlung,x0,x1,visbrem_val,LENGTH=visbrem_length)
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
            vals(mtargets)    = visbrem_val*calib_visbrem_line(ik)
            IF (iflag == 1) WRITE(iunit_out,'(11ES22.12E3)') &
                                 target(ik),sigma(ik),visbrem_val, &
                                 lambda_visbrem_line(ik), calib_visbrem_line(ik), &
                                 x0(1),x0(2),x0(3),x1(1),x1(2),x1(3)
         END DO
      ELSE
         DO ik = 1, nprof
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets) = jtarget_visbrem_line
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_line_visbrem
