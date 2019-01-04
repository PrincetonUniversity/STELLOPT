!-----------------------------------------------------------------------
!     Subroutine:    chisq_neo
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/09/2012
!     Description:   Calculate difference between target and simulated
!                    effective ripple as calculated by NEO.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_neo(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars, ONLY: equil_type
      USE equil_vals, ONLY: eff_ripple
!DEC$ IF DEFINED (NEO_OPT)
      USE neo_input_mod, ONLY: read_neoin_input, write_neoin_namelist
!DEC$ ENDIF
      USE mpi_params, ONLY: myid, master
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: i, j, ik, ier
      CHARACTER(len = 256)   :: temp_str
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ik   = COUNT(sigma < bigno)
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'NEO ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  EFF_RIPPLE  #'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = eff_ripple(ik)
               IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3,2X,I3.3)') target(ik),sigma(ik),vals(mtargets),ik
               IF (iflag == 1) CALL FLUSH(iunit_out)
            END IF
         END DO
      ELSE            !12/28/18.so niter < 0 (stel_init > load_targ calls)
         DO ik = 1, nsd
            IF (sigma(ik) < bigno) THEN
               lbooz(ik) = .TRUE.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_neo
            END IF
         END DO
!DEC$ IF DEFINED (NEO_OPT)
!         CALL read_neoin_input(TRIM(id_string),iflag)  !12/28/18.(7m14b)mv to stel_init.
         IF (iflag < 0 .and. niter == -2 .and. myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
            WRITE(6,*) '  NEO_IN Namelist not found     '
            WRITE(6,*) ' '
            k=0
            CALL write_neoin_namelist(6,k)
            WRITE(6,*) ' '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         END IF
         IF (iflag /=0) RETURN
!DEC$ ENDIF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_neo
