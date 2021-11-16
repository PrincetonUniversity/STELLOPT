!-----------------------------------------------------------------------
!     Subroutine:    chisq_stella
!     Authors:       J.L. Velasco (joseluis.velasco@ciemat.es),
!                    building on chisq_stella
!     Date:          last modified, 29/10/2020
!     Description:   This is the chisquared routine which compares
!                    target to equilibrium stella values.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_stella(target,sigma,niter,iflag,jtarget_stella)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_utils, ONLY: get_equil_phi,rho, phi_type
!DEC$ IF DEFINED (STELLA_OPT)
      USE stella_stellopt_mod
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter,jtarget_stella
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
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'STELLA',ik,15
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            mtargets = mtargets + 1
            targets(mtargets) = target(ik)
            sigmas(mtargets)  = sigma(ik)
!            targets(mtargets) = target_stella(ik)
!            sigmas(mtargets)  = sigma_stella(ik)
!DEC$ IF DEFINED (STELLA_OPT)
            IF(jtarget_stella.EQ.jtarget_stella_q1) THEN
               vals(mtargets)    = STELLA_q1(ik)
            ELSE IF(jtarget_stella.EQ.jtarget_stella_q2) THEN
               vals(mtargets)    = STELLA_q2(ik)
            END IF
!            IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
!               target(ik),sigma(ik),vals(mtargets),&
!               rho(ik), STELLA(ik)
!DEC$ ELSE
            vals(mtargets) = target(ik)
!            IF (iflag == 1) WRITE(iunit_out,'(15ES22.12E3)') &
!               target(ik),sigma(ik),vals(mtargets),&
!               rho(ik), target(ik)
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
               IF (niter == -2) target_dex(mtargets)=jtarget_stella
            END IF
         END DO
!!DEC$ IF DEFINED (STELLA_OPT)
!         CALL read_stellein_input(TRIM(id_string),iflag)
!         IF (iflag < 0 .and. niter == -2 .and. myid == master) THEN
!            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
!            WRITE(6,*) '  STELLA_IN Namelist not found     '
!            WRITE(6,*) ' '
!            k=0
!            CALL write_stellain_namelist(6,k)
!            WRITE(6,*) ' '
!            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!         END IF
!         IF (iflag /=0) RETURN
!!DEC$ ENDIF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_stella
