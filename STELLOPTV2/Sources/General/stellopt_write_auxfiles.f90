!-----------------------------------------------------------------------
!     Subroutine:    stellopt_write_auxfiles
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          11/12/2019
!     Description:   This subroutine handles writing and copying of
!                    all the auxilliary stellopt files.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_write_auxfiles
      USE stellopt_runtime
      USE stellopt_input_mod
      USE equil_utils, ONLY: move_txtfile, copy_txtfile, copy_boozer_file
      USE beams3d_runtime, ONLY: id_string_beams => id_string, lverb_beams => lverb
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        NONE
!----------------------------------------------------------------------
      IMPLICIT NONE

!----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        vcrtl_array VMEC Control Array
!----------------------------------------------------------------------
      LOGICAL                :: lfile_found
      INTEGER                :: ier, ik, ialpha
      INTEGER                ::  vctrl_array(5)
      CHARACTER(len = 256)   :: temp_str, reset_string

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ier = 0
      IF (ANY(sigma_txport < bigno)) THEN
         DO ik = 1, 256
            DO ialpha = 1, 256
               WRITE(temp_str,'(2(A,I3.3))') '_',ik,'_',ialpha
               CALL move_txtfile('gist_genet_'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'gist_genet_'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
               CALL move_txtfile('curv_stellopt_'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'curv_stellopt_'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
            END DO
            WRITE(temp_str,'(A,I3.3)') '_',ik
            CALL move_txtfile('txport_out.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'txport_out.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
            CALL move_txtfile('gist_geney_'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'gist_geney_'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
         END DO
      END IF
      lfile_found = .false.
      INQUIRE(FILE='parameters',EXIST=lfile_found)
      IF (lfile_found .AND. ANY(sigma_txport < bigno) .AND. (txport_proxy == 'gene_parallel')) THEN
         CALL move_txtfile('log_gene.'//TRIM(proc_string_old),'log_gene.'//TRIM(proc_string))
         DO ik = 1, 256
            DO ialpha = 1, 256
               CALL move_txtfile('gist_'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'gist_'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
               CALL move_txtfile('eigenvalues_'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'eigenvalues_'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
               CALL move_txtfile('parameters_'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'parameters_'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
            END DO
         END DO
      END IF
      CALL move_txtfile('mercier.'//TRIM(proc_string_old),'mercier.'//TRIM(proc_string))
      CALL copy_txtfile('diagno_bth.'//TRIM(proc_string_old),'diagno_bth.'//TRIM(proc_string))
      CALL copy_txtfile('diagno_flux.'//TRIM(proc_string_old),'diagno_flux.'//TRIM(proc_string))
      CALL copy_txtfile('diagno_seg.'//TRIM(proc_string_old),'diagno_seg.'//TRIM(proc_string))
      CALL move_txtfile('jBbs.'//TRIM(proc_string_old),'jBbs.'//TRIM(proc_string))
      CALL move_txtfile('answers_plot.'//TRIM(proc_string_old),'answers_plot.'//TRIM(proc_string))
      CALL move_txtfile('answers.'//TRIM(proc_string_old),'answers.'//TRIM(proc_string))
      CALL move_txtfile('neo_cur.'//TRIM(proc_string_old),'neo_cur.'//TRIM(proc_string))
      CALL move_txtfile('neolog.'//TRIM(proc_string_old),'neolog.'//TRIM(proc_string))
      CALL move_txtfile('neo_out.'//TRIM(proc_string_old),'neo_out.'//TRIM(proc_string))
      CALL move_txtfile('tprof.'//TRIM(proc_string_old),'tprof.'//TRIM(proc_string))
      CALL move_txtfile('jprof.'//TRIM(proc_string_old),'jprof.'//TRIM(proc_string))
      CALL move_txtfile('dprof.'//TRIM(proc_string_old),'dprof.'//TRIM(proc_string))
      CALL move_txtfile('boot_fit.'//TRIM(proc_string_old),'boot_fit.'//TRIM(proc_string))
      CALL move_txtfile('coils.'//TRIM(proc_string_old),'coils.'//TRIM(proc_string))
      CALL move_txtfile('bnorm.'//TRIM(proc_string_old),'bnorm.'//TRIM(proc_string))
      CALL move_txtfile('bnorm_real.'//TRIM(proc_string_old),'bnorm_real.'//TRIM(proc_string))
      CALL move_txtfile('coil_curvature.'//TRIM(proc_string_old),'coil_curvature.'//TRIM(proc_string))
      CALL copy_boozer_file(TRIM(proc_string_old),TRIM(proc_string))
      DO ik = 1, nsd
         WRITE(temp_str,'(A,I3.3)') '_s',ik
         CALL move_txtfile('input_dkes.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'input_dkes.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
         CALL move_txtfile('dkesout.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'dkesout.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
         CALL move_txtfile('opt_dkes.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'opt_dkes.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
         CALL move_txtfile('results.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'results.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
      END DO
      IF (ANY(sigma_orbit .lt. bigno)) THEN
         lverb_beams = .FALSE.
!         id_string_beams = TRIM(proc_string_old)
         CALL beams3d_read(TRIM(proc_string_old))
!         id_string_beams = TRIM(proc_string)
         CALL beams3d_write('GRID_INIT')
         CALL beams3d_write('TRAJECTORY_FULL')
         CALL beams3d_write('DIAG')
         CALL beams3d_free
         CALL move_txtfile('beams3d_diag_'//TRIM(proc_string_old)//'.txt',&
                           'beams3d_diag_'//TRIM(proc_string)//'.txt')
      END IF
!DEC$ IF DEFINED (TERPSICHORE)
      IF (ANY(sigma_kink < bigno)) THEN
         CALL move_txtfile('terpsichore_eq.'//TRIM(proc_string_old),&
                           'terpsichore_eq.'//TRIM(proc_string))
         DO ik = 1, nsys
            WRITE(temp_str,'(1(A,I2.2))') '_',ik-1
            CALL move_txtfile('terpsichore_16.'//TRIM(proc_string_old)//TRIM(temp_str),&
                              'terpsichore_16.'//TRIM(proc_string)//TRIM(temp_str))
            CALL move_txtfile('terpsichore_17.'//TRIM(proc_string_old)//TRIM(temp_str),&
                              'terpsichore_17.'//TRIM(proc_string)//TRIM(temp_str))
            CALL move_txtfile('terpsichore_19.'//TRIM(proc_string_old)//TRIM(temp_str),&
                              'terpsichore_19.'//TRIM(proc_string)//TRIM(temp_str))
            CALL move_txtfile('terpsichore_22.'//TRIM(proc_string_old)//TRIM(temp_str),&
                              'terpsichore_22.'//TRIM(proc_string)//TRIM(temp_str))
         END DO
      END IF
!DEC$ ENDIF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_write_auxfiles