!-----------------------------------------------------------------------
!     Subroutine:    stellopt_run_vmec
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          12/11/2025
!     Description:   This subroutine runs the VMEC code twice, targeting
!                    the vacuum field on axis.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_run_vmec(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stellopt_runtime
      USE stellopt_globals, ONLY: b0_vac
      USE equil_utils, ONLY: eval_prof_spline, Baxis
      USE vmec_input, ONLY: curtor, pres_scale,phiedge
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        lscreen     Print to screen logical
!        iflag       Error flag
!-----------------------------------------------------------------------
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!        iflag       Error flag
!        iunit       File unit number
!-----------------------------------------------------------------------
      REAL(rprec) :: curtor_save, pres_scale_save, phiedge_new
      
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      iflag = 0
      curtor_save = curtor
      pres_scale_save = pres_scale
      IF (ABS(B0_VAC) > 0) THEN
         curtor = 0.0
         pres_scale = 0.0
      END IF
      CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
      iflag = ier_paraexe
      IF (ABS(B0_VAC) > 0) THEN
         CALL stellopt_load_equil(.FALSE.,iflag)
         phiedge_new = (b0_vac/Baxis)*phiedge
         IF (lscreen .and. lverb) THEN 
            WRITE(6,*)  '----- Recomputing PHIEDGE ------'
            WRITE(6,'(A,F7.3)') '      PHIEDGE(OLD): ', phiedge
            WRITE(6,'(A,F7.3)') '     BVACAXIS(OLD): ', Baxis
            WRITE(6,'(A,F7.3)') '  BVACAXIS(TARGET): ', b0_vac
            WRITE(6,'(A,F7.3)') '      PHIEDGE(NEW): ', phiedge_new
            WRITE(6,*)  '--------------------------------'
         END IF
         phiedge = phiedge_new
         curtor = curtor_save
         pres_scale = pres_scale_save
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
         iflag = ier_paraexe
      END IF
      IF (lscreen) WRITE(6,*)  '-------------------------  PARAVMEC CALCULATION DONE  -----------------------'
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE stellopt_run_vmec