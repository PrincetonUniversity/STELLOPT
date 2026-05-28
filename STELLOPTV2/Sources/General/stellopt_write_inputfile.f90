!-----------------------------------------------------------------------
!     Subroutine:    stellopt_write_inputfile
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          11/12/2019
!     Description:   This subroutine handles writeing the full
!                    STELLOPT input file
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_write_inputfile(ncnt,lmin)
      USE stellopt_runtime
      USE stellopt_input_mod
      USE safe_open_mod, ONLY: safe_open
      USE vmec_input
      USE bootsj_input, ONLY: write_bootsj_input
      USE diagno_input_mod, ONLY: write_diagno_input
      USE fieldlines_input_mod, ONLY: write_fieldlines_namelist
!DEC$ IF DEFINED (NEO_OPT)
      USE neo_input_mod, ONLY: write_neoin_namelist
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
      USE beams3d_input_mod, ONLY: write_beams3d_namelist
!DEC$ ENDIF
!DEC$ IF DEFINED (AEOPT)
      USE trapped_avail_energy_mod, ONLY: write_avail_energy_nml
!DEC$ ENDIF
      USE penta_interface_mod, ONLY: write_ion_params_nml ,&
                                     write_run_params_nml
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        ncnt          Interation Identifier
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: ncnt
      LOGICAL, INTENT(in)    :: lmin

!----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER                :: ier
      CHARACTER(len = 256)   :: temp_str

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ier = 0
      iunit_out = 12
      WRITE(temp_str,'(i5.5)') ncnt
      IF (lmin) THEN
         proc_string = TRIM(id_string) // '_min'
      ELSE
         proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
      END IF
      CALL safe_open(iunit_out,ier,TRIM('input.'//TRIM(proc_string)),'unknown','formatted')
         SELECT CASE(TRIM(equil_type))
            CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot','vmec2000_oneeq','vmec_provided')
               CALL RESCALE_BOUNDARY ! Necssary for output file to have correct RBC/ZBS
               CALL write_indata_namelist(iunit_out,ier)
            CASE('test')
         END SELECT
      CALL write_optimum_namelist(iunit_out,ier)
      IF (lneed_magdiag) CALL write_diagno_input(iunit_out,ier)
      IF (ANY(sigma_bootstrap < bigno) .or. (sigma_totalbootstrap < bigno)) CALL write_bootsj_input(iunit_out,ier)
      IF (lpoincare) CALL write_fieldlines_namelist(iunit_out,ier)
!DEC$ IF DEFINED (NEO_OPT)
      IF (ANY(sigma_neo < bigno)) CALL write_neoin_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
      IF (ANY(sigma_orbit < bigno)) CALL write_beams3d_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (AEOPT)
      IF (ANY(sigma_txport < bigno)) CALL write_avail_energy_nml(iunit_out,ier)
!DEC$ ENDIF
      IF (ANY(lneed_penta)) THEN
         CALL write_ion_params_nml(iunit_out)
         CALL write_run_params_nml(iunit_out)
      ENDIF
      WRITE(iunit_out,'(A)') '&END'
      CLOSE(iunit_out)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_write_inputfile
