!-----------------------------------------------------------------------
!     Subroutine:    chisq_fluxloops
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   Calculates fit to fluxloop magnetic diagnostics.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_fluxloops(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE diagno_input_mod, ONLY: read_diagno_input, lskip_flux
      USE safe_open_mod
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
!
!-----------------------------------------------------------------------
      INTEGER :: iunit, dex, itmp, ik, ier
      REAL(rprec) :: flx_temp
      REAL(rprec), ALLOCATABLE :: flux(:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (niter >= 0) THEN
         ! Read the fluxloops file
         iunit = 12
         CALL safe_open(iunit,ier,'diagno_flux.'//TRIM(proc_string_old),'old','formatted')
         IF (ier < 0) THEN; iflag=ier; RETURN; END IF
         READ(iunit,'(I6)',IOSTAT=ier) dex
         IF (ier /= 0) THEN; iflag=ier; RETURN; END IF
         ALLOCATE(flux(dex))
         flux(:) = 0.0
         DO ik = 1, dex
            READ(iunit,'(1X,1E22.14)',IOSTAT=ier) flx_temp
            IF (ier == 0) flux(ik) = flx_temp
         END DO
         CLOSE(iunit)
         IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I3.3))') 'FLUXLOOPS ',dex,3
         IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  FLUX'
         DO ik = 1, nprof
            ! Reordered so we get a value for every flux loop
            IF (sigma(ik) < bigno) THEN
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = flux(ik)
            END IF
            IF (iflag == 1 .and. ik <= dex) WRITE(iunit_out,'(3ES22.12E3)') target(ik),sigma(ik),flux(ik)
         END DO
         DEALLOCATE(flux)
      ELSE
         CALL read_diagno_input(TRIM(id_string),iflag)
         lskip_flux = .true.
         DO ik =1, nprof
            IF (sigma(ik) < bigno) THEN
               lneed_magdiag = .true.
               lskip_flux(ik) = .false.
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_fluxloop
            END IF
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_fluxloops
