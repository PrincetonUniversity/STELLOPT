!-----------------------------------------------------------------------
!     Subroutine:    chisq_totalbootstrap
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          03/04/2025
!     Description:   This subroutine calculates the total bootstrap
!                    current functional.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_totalbootstrap(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE safe_open_mod, ONLY: safe_open
      USE bootsj_input
      USE vmec_input, ONLY: ns_array
      USE mpi_params, ONLY: myid, master
      USE parambs, ONLY: aibs, irup
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: ns_max, iunit, ier

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'TOTALBOOTSTRAP ',1,3
      IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  TOTALBOOTSTRAP'
      IF (niter >= 0) THEN
         mtargets = mtargets + 1
         targets(mtargets) = target
         sigmas(mtargets)  = sigma
         vals(mtargets)    = aibs(irup)*1E6 ! aibs in MA
         IF (iflag == 1) WRITE(iunit_out,'(3ES22.12E3)') target,sigma,aibs(irup)*1E6
      ELSE
         IF (sigma < bigno) THEN
            ! Set LBOOZ to true over full radius
            ns_max = MAXVAL(ns_array)
            lbooz(2:ns_max) = .TRUE.
            mtargets = mtargets + 1
            IF (niter == -2) target_dex(mtargets)=jtarget_totalbootstrap
         END IF
         iunit=12
         CALL safe_open(iunit,iflag,'input.'//TRIM(id_string),'old','formatted')
         IF (iflag < 0) RETURN
         CALL read_namelist (iunit, iflag, 'bootin')
         IF (iflag < 0 .and. myid == master) THEN
            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
            WRITE(6,*) '  BOOTIN Namelist not found     '
            WRITE(6,*) ' '
            ier=0
            CALL write_bootsj_input(6,ier)
            WRITE(6,*) ' '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         END IF
         CLOSE(iunit)
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_totalbootstrap
