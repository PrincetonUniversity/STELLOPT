!-----------------------------------------------------------------------
!     Subroutine:    chisq_ece
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculate difference measured and equilibrium
!                    ECE reflected power.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_ece(target,sigma,niter,iflag)
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
      REAL(rprec), INTENT(in)    ::  target(nsys,nprof)
      REAL(rprec), INTENT(in)    ::  sigma(nsys,nprof)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        lreset_s    Gets set to true if using R,PHI,Z Specification
!        ik          Dummy index
!        ti_val      Holds profile evaulation
!-----------------------------------------------------------------------
      INTEGER ::  ii,ij,ik,ier,n
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0 ) RETURN
      ik = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'ECEREFLECT',ik,7
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  REFLECT  FREQ  RADTX  RADTO  MIX'
      IF (niter >= 0) THEN
      	DO ij = 1, nsys
      	   IF (ALL(sigma(ij,:) >= bigno,DIM=1)) CYCLE
            DO ik = 1, nprof
               IF (sigma(ij,ik) >= bigno) CYCLE
               mtargets = mtargets + 1
               targets(mtargets) = target(ij,ik)
               sigmas(mtargets)  = sigma(ij,ik)
               vals(mtargets)    = (radtx_ece(ij,ik)+mix_ece*radto_ece(ij,ik))*1000
               IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') target(ij,ik),sigma(ij,ik),vals(mtargets),freq_ece(ij,ik),radtx_ece(ij,ik)*1000,radto_ece(ij,ik)*1000,mix_ece
            END DO
         END DO
      ELSE
      	DO ij = 1, nsys
            DO ik = 1, nprof
               IF (sigma(ij,ik) < bigno) THEN
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets) = jtarget_ece
               END IF
            END DO
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_ece
