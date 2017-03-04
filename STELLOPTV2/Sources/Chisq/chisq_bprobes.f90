!-----------------------------------------------------------------------
!     Subroutine:    chisq_bprobes
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculates fit to B-Probe magnetic diagnostics.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_bprobes(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE diagno_input_mod, ONLY: read_diagno_input
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nprobes)
      REAL(rprec), INTENT(in)    ::  sigma(nprobes)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: iunit, dex, itmp, ik,ier
      REAL(rprec), ALLOCATABLE :: xp(:),yp(:),zp(:),modb(:),flux(:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (niter >= 0) THEN
         ! Read the bprobes file
         iunit = 12
         CALL safe_open(iunit,ier,'diagno_bth.'//TRIM(proc_string_old),'old','formatted')
         IF (ier < 0) THEN; iflag=ier; RETURN; END IF
         dex = max(COUNT(target /= 0.0),COUNT(sigma < bigno))
         ALLOCATE(xp(dex),yp(dex),zp(dex),modb(dex),flux(dex))
         DO ik = 1, dex
            READ(iunit,*) itmp,xp(ik),yp(ik),zp(ik),modb(ik),flux(ik)
         END DO
         CLOSE(iunit)
         IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'B_PROBES ',dex,7
         IF (iflag == 1) WRITE(iunit_out,'(A)') 'X  Y  Z  |B|  TARGET  SIGMA  FLUX'
         DO ik = 1, nprobes
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = flux(ik)
            END IF
            IF (iflag == 1 .and. ik <= dex) WRITE(iunit_out,'(7ES22.12E3)') xp(ik),yp(ik),zp(ik),modb(ik),target(ik),sigma(ik),flux(ik)
         END DO
         DEALLOCATE(xp,yp,zp,modb,flux)
      ELSE
         DO ik =1, nprobes
            IF (sigma(ik) < bigno) THEN
               lneed_magdiag = .true.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_bprobe
            END IF
         END DO
         CALL read_diagno_input(TRIM(id_string),iflag)
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_bprobes
