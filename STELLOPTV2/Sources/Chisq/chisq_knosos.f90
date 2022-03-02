!-----------------------------------------------------------------------
!     Subroutine:    chisq_knosos
!     Authors:       J.L. Velasco (joseluis.velasco@ciemat.es),
!                    building on chisq_dkes (by S. Lazerson)
!     Date:          last modified, 16/02/2020
!     Description:   This is the chisquared routine which compares
!                    target to equilibrium KNOSOS values.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_knosos(target,sigma,niter,iflag,jtarget_knosos)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils, ONLY: get_equil_phi,rho, phi_type
      USE knosos_stellopt_mod
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter,jtarget_knosos
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ik
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'KNOSOS',ik,15
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            targets(mtargets) = target(ik) 
            sigmas(mtargets)  = sigma(ik)
!            targets(mtargets) = target_knosos(ik) 
!            sigmas(mtargets)  = sigma_knosos(ik)  
!DEC$ IF DEFINED (KNOSOS_OPT)
            IF(jtarget_knosos.EQ.jtarget_knosos_1nu) THEN
               vals(mtargets)    = KNOSOS_1nu(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_snu) THEN
               vals(mtargets)    = KNOSOS_snu(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_sbp) THEN
               vals(mtargets)    = KNOSOS_sbp(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_gmc) THEN
               vals(mtargets)    = KNOSOS_gmc(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_gma) THEN
               vals(mtargets)    = KNOSOS_gma(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_qer) THEN
               vals(mtargets)    = KNOSOS_qer(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_vb0) THEN
               vals(mtargets)    = KNOSOS_vb0(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_vbm) THEN
               vals(mtargets)    = KNOSOS_vbm(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_vbb) THEN
               vals(mtargets)    = KNOSOS_vbb(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_wbw) THEN
               vals(mtargets)    = KNOSOS_wbw(ik)
            ELSE IF(jtarget_knosos.EQ.jtarget_knosos_dbo) THEN
               vals(mtargets)    = KNOSOS_dbo(ik)
            END IF
!            IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
!               target(ik),sigma(ik),vals(mtargets),&
!               rho(ik), KNOSOS(ik)
!DEC$ ELSE
            vals(mtargets) = target(ik)
!            IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
!               target(ik),sigma(ik),vals(mtargets),&
!               rho(ik), target_knosos(ik)
!DEC$ ENDIF
         END DO
      ELSE
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               IF(ik.GT.2) THEN
                  lbooz(ik-1) = .TRUE.
               ELSE IF(ik.LT.nsd) THEN
                  lbooz(ik+1) = .TRUE.
               END IF
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_knosos
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_knosos
