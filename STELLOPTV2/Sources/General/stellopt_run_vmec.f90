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
      USE stellopt_vars, ONLY: lcreate_coils
      USE equil_utils, ONLY: eval_prof_spline, Baxis
      USE equil_vals, ONLY: volume
      USE vmec_input, ONLY: curtor, pres_scale,phiedge, lfreeb, extcur, &
         tvolume
      USE biotsavart, ONLY: parse_coils_file
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
         IF (lscreen .and. lverb) WRITE(6,*)  '----- Runing zero beta VMEC ------'
      END IF
      CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
      iflag = ier_paraexe
      ! First 
      IF (ABS(B0_VAC) > 0) THEN
         CALL stellopt_load_equil(.FALSE.,iflag)
         IF (lfreeb) THEN
            extcur = (b0_vac/Baxis)*extcur
            IF (lscreen .and. lverb) THEN 
               WRITE(6,*)  '----- Recomputing EXTCUR ------'
               WRITE(6,'(A,F7.3)') '     BVACAXIS(OLD): ', Baxis
               WRITE(6,'(A,F7.3)') '  BVACAXIS(TARGET): ', b0_vac
               WRITE(6,'(A,F7.3)') '     EXTCUR_FACTOR: ', b0_vac/Baxis
               WRITE(6,'(A,F7.3)') '      PHIEDGE(OLD): ', phiedge
               WRITE(6,*)  '--------------------------------'
            END IF
         ELSE
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
         END IF
         curtor = curtor_save
         pres_scale = pres_scale_save
         IF (lscreen .and. lverb) WRITE(6,*)  '----- Runing finite beta VMEC ------'
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
         iflag = ier_paraexe
      END IF
      ! Now adjust volume via PHIEDGE if in free boundary at full parameters
      IF (lfreeb .and. (tvolume .gt. 0.0)) THEN
         CALL stellopt_load_equil(.FALSE.,iflag)
         phiedge_new = phiedge * (tvolume/volume)
         IF (lscreen .and. lverb) THEN 
            WRITE(6,*)  '----- Recomputing PHIEDGE ------'
            WRITE(6,'(A,F7.3)') '      PHIEDGE(OLD): ', phiedge
            WRITE(6,'(A,F7.2)') '       VOLUME(OLD): ', volume
            WRITE(6,'(A,F7.2)') '    VOLUME(TARGET): ', tvolume
            WRITE(6,'(A,F7.3)') '      PHIEDGE(NEW): ', phiedge_new
            WRITE(6,*)  '--------------------------------'
         END IF
         phiedge = phiedge_new
         IF (lscreen .and. lverb) WRITE(6,*)  '----- Runing finite beta VMEC ------'
         CALL stellopt_paraexe('paravmec_run',proc_string,lscreen)
         iflag = ier_paraexe
      END IF
      ! VMEC deallocates the coil_group structure in STELLOPT
      IF (lfreeb .and. lcreate_coils) CALL parse_coils_file('coils.'//TRIM(proc_string))
      IF (lscreen) WRITE(6,*)  '-------------------------  PARAVMEC CALCULATION DONE  -----------------------'
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE stellopt_run_vmec