!-----------------------------------------------------------------------
!     Subroutine:    stellopt_clean_up
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/26/2012
!     Description:   This subroutine handles cleaning up the function
!                    calls.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_clean_up(ncnt,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE stellopt_targets, ONLY: sigma_bootstrap, lbooz, numws
      USE safe_open_mod, ONLY: safe_open
      USE diagno_input_mod, ONLY: write_diagno_input
      USE gist_mod, ONLY: write_gist_namelist
      USE bootsj_input, ONLY: write_bootsj_input
      USE equil_utils, ONLY: move_txtfile, copy_txtfile, copy_boozer_file
!DEC$ IF DEFINED (NEO_OPT)
      USE neo_input_mod, ONLY: write_neoin_namelist
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
      USE beams3d_runtime, ONLY: id_string_beams => id_string, lverb_beams => lverb
      USE beams3d_input_mod, ONLY: write_beams3d_namelist
!DEC$ ENDIF
      USE vmec_input
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
!                             animec_flag, flow_flag
      USE fdjac_mod, ONLY: flag_singletask, flag_cleanup, &
                           JAC_CLEANUP => flag_cleanup_jac,&
                           LEV_CLEANUP => flag_cleanup_lev
      USE gade_mod, ONLY: GADE_CLEANUP, PSO_CLEANUP
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: ncnt
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL ::  lfile_found
      INTEGER ::  ier, ik, iunit, ctype, temp_max, ialpha, m, n
      INTEGER ::  vctrl_array(5)
      CHARACTER(len = 256)   :: temp_str, reset_string
      REAL(rprec), ALLOCATABLE :: fvec_temp(:)
      
!      INTEGER, PARAMETER :: JAC_CLEANUP = -100
!      INTEGER, PARAMETER :: LEV_CLEANUP = -101
!      INTEGER, PARAMETER :: PSO_CLEANUP = -300
      INTEGER, PARAMETER :: JUST_INPUT  = -110
      INTEGER, PARAMETER :: LAST_GO     = -500
      
      LOGICAL, PARAMETER :: lkeep_extra = .TRUE.  ! For debugging the code will keep the _optXXX files
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ctype = iflag
      iflag = 0
      iunit = 10
      iunit_out = 12
      ier = 0
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot')
            IF (ctype == PSO_CLEANUP) THEN
               IF (ncnt /= 1) THEN
                  WRITE(temp_str,'(i5.5)') ncnt
                  proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
                  CALL safe_open(iunit_out,iflag,TRIM('input.'//TRIM(proc_string)),'unknown','formatted')
                  CALL write_indata_namelist(iunit_out,ier)
                  CALL write_optimum_namelist(iunit_out,ier)
                  IF (lneed_magdiag) CALL write_diagno_input(iunit_out,ier)
!DEC$ IF DEFINED (TXPORT_OPT)
!                  IF (ANY(sigma_txport < bigno)) CALL write_gist_namelist(iunit_out,ier)
!DEC$ ENDIF
                  IF (ANY(sigma_bootstrap < bigno)) CALL write_bootsj_input(iunit_out,ier)
!DEC$ IF DEFINED (NEO_OPT)
                  IF (ANY(sigma_neo < bigno)) CALL write_neoin_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
                  IF (ANY(sigma_orbit < bigno)) CALL write_beams3d_namelist(iunit_out,ier)
!DEC$ ENDIF
                  WRITE(iunit_out,'(A)') '&END'
                  CLOSE(iunit_out)
               END IF
               IF (lkeep_mins) THEN
                  WRITE(temp_str,'(i5.5)') ncnt
                  proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
                  SELECT CASE(TRIM(equil_type))
                     CASE ('vmec2000_old','animec','flow','satire')
                          vctrl_array(1) = output_flag ! Output to file
                          vctrl_array(2) = 0     ! vmec error flag  
                          vctrl_array(3) = 0    ! Use multigrid
                          vctrl_array(4) = 0     ! Iterative 
                          vctrl_array(5) = myid ! Output file sequence number
                          !IF (TRIM(equil_type)=='animec') vctrl_array(1) = vctrl_array(1) + animec_flag
                          !IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') vctrl_array(1) = vctrl_array(1) + flow_flag
                         CALL runvmec(vctrl_array,proc_string,.false.,MPI_COMM_SELF,'')
                     CASE('parvmec','paravmec','vmec2000','vboot')
                         CALL stellopt_paraexe('paravmec_write',proc_string,.false.)
                  END SELECT
                  ier=vctrl_array(2)
                  iflag = ier
                  IF (ier == successful_term_flag) iflag = 0
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
                  IF (lcoil_geom) THEN
                     CALL move_txtfile('coils.'//TRIM(proc_string_old),'coils.'//TRIM(proc_string))
                  END IF
                  DO ik = 1, nsd
                     WRITE(temp_str,'(A,I3.3)') '_s',ik
                     CALL move_txtfile('dkesout.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'dkesout.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
                     CALL move_txtfile('opt_dkes.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'opt_dkes.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
                     CALL move_txtfile('results.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'results.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
                  END DO
!DEC$ IF DEFINED (BEAMS3D_OPT)
                  IF (ANY(sigma_orbit .lt. bigno)) THEN
                     lverb_beams = .FALSE.
                     id_string_beams = TRIM(proc_string_old)
                     CALL beams3d_read
                     id_string_beams = TRIM(proc_string)
                     CALL beams3d_write('GRID_INIT')
                     CALL beams3d_write('TRAJECTORY_FULL')
                     CALL beams3d_write('DIAG')
                     CALL beams3d_free
                     CALL move_txtfile('beams3d_diag_'//TRIM(proc_string_old)//'.txt',&
                                       'beams3d_diag_'//TRIM(proc_string)//'.txt')
                  END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (COILOPTPP)
                  IF (sigma_coil_bnorm < bigno) THEN
                     CALL move_txtfile('bnorm.'//TRIM(proc_string_old),&
                                       'bnorm.'//TRIM(proc_string))
                     CALL move_txtfile('coilopt_params.'//TRIM(proc_string_old),&
                                       'coilopt_params.'//TRIM(proc_string))
                     CALL copy_txtfile('b_norm_eq_'//TRIM(proc_string_old)//'.dat',&
                                       'b_norm_eq_'//TRIM(proc_string)//'.dat')
                     CALL copy_txtfile('b_norm_final_'//TRIM(proc_string_old)//'.dat',&
                                       'b_norm_final_'//TRIM(proc_string)//'.dat')
                     CALL move_txtfile('b_norm_init_'//TRIM(proc_string_old)//'.dat',&
                                       'b_norm_init_'//TRIM(proc_string)//'.dat')
                     DO ik = 0, numws-1
                        WRITE(temp_str,'(I3.3)') ik
                        CALL copy_txtfile('coil_spline'//TRIM(temp_str)//'_'//TRIM(proc_string_old)//'.out',&
                                          'coil_spline'//TRIM(temp_str)//'_'//TRIM(proc_string)//'.out')
                     END DO
                  END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (REGCOIL)
                 ! OUTPUT FILES SHOULD BE WRITTEN HERE - Use the regcoil
                 ! functions to write the hdf5 output file
                 ! This is inside of the PSO loop. Should be
                 ! duplicated, or broken out to a subroutine
                 ! WRITE *, '<----- REGCOIL Output files missing -----'
!DEC$ ENDIF
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
                        !CALL move_txtfile('terpsichore_23'//TRIM(proc_string_old)//TRIM(temp_str),&
                        !                  'terpsichore_23'//TRIM(proc_string)//TRIM(temp_str))
                     END DO
                  END IF
!DEC$ ENDIF
               END IF
               ! Now open the Output file
               ALLOCATE(fvec_temp(mtargets))
               CALL safe_open(iunit_out,iflag,TRIM('stellopt.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
               iflag = 1
               IF (ncnt == 0) WRITE(iunit_out,'(A,1X,F5.2)') 'VERSION',STELLOPT_VERSION
               WRITE(iunit_out,'(A,1X,I5.5)') 'ITER',ncnt
               CALL stellopt_load_targets(mtargets,fvec_temp,iflag,ncnt)          ! Count
               WRITE(iunit_out,'(A,2(2X,I8))') 'TARGETS ',mtargets,1
               WRITE(iunit_out,'(A)') 'TARGETS'
               WRITE(iunit_out,'(ES22.12E3)') targets(1:mtargets)
               WRITE(iunit_out,'(A,2(2X,I8))') 'SIGMAS ',mtargets,1
               WRITE(iunit_out,'(A)') 'SIGMAS'
               WRITE(iunit_out,'(ES22.12E3)') sigmas(1:mtargets)
               WRITE(iunit_out,'(A,2(2X,I8))') 'VALS ',mtargets,1
               WRITE(iunit_out,'(A)') 'VALUES'
               WRITE(iunit_out,'(ES22.12E3)') vals(1:mtargets)
               CLOSE(iunit_out)
               DEALLOCATE(fvec_temp)
            ELSE IF ((ctype == LEV_CLEANUP) .or. (ctype == GADE_CLEANUP)) THEN
               IF (ncnt /= 1 .or. ctype == GADE_CLEANUP) THEN
                  ! Write the input file
                  WRITE(temp_str,'(i5.5)') ncnt
                  proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
                  CALL safe_open(iunit_out,iflag,TRIM('input.'//TRIM(proc_string)),'unknown','formatted')
                  CALL write_indata_namelist(iunit_out,ier)
                  CALL write_optimum_namelist(iunit_out,ier)
                  IF (lneed_magdiag) CALL write_diagno_input(iunit_out,ier)
!DEC$ IF DEFINED (TXPORT_OPT)
!                  IF (ANY(sigma_txport < bigno)) CALL write_gist_namelist(iunit_out,ier)
!DEC$ ENDIF
                  IF (ANY(sigma_bootstrap < bigno)) CALL write_bootsj_input(iunit_out,ier)
!DEC$ IF DEFINED (NEO_OPT)
                  IF (ANY(sigma_neo < bigno)) CALL write_neoin_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
                  IF (ANY(sigma_orbit < bigno)) CALL write_beams3d_namelist(iunit_out,ier)
!DEC$ ENDIF
                  WRITE(iunit_out,'(A)') '&END'
                  CLOSE(iunit_out)
               END IF
               ! Overwrite the restart file
               proc_string = 'reset_file'
               SELECT CASE(TRIM(equil_type))
                  CASE ('vmec2000_old','animec','flow','satire')
                     vctrl_array(1) = output_flag ! Output to file
                     vctrl_array(2) = 0     ! vmec error flag  
                     vctrl_array(3) = 0    ! Use multigrid
                     vctrl_array(4) = 0     ! Iterative 
                     vctrl_array(5) = myid ! Output file sequence number
                     !IF (TRIM(equil_type)=='animec') vctrl_array(1) = vctrl_array(1) + animec_flag
                     !IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') vctrl_array(1) = vctrl_array(1) + flow_flag
                     CALL runvmec(vctrl_array,proc_string,.false.,MPI_COMM_SELF,'')
                  CASE('parvmec','paravmec','vmec2000','vboot')
                     CALL stellopt_paraexe('paravmec_write',proc_string,.false.)
               END SELECT
               iflag = ier_paraexe
               iflag = ier
               IF (ier_paraexe == successful_term_flag) iflag = 0
!DEC$ IF DEFINED (COILOPTPP)
               IF (sigma_coil_bnorm < bigno .and. (proc_string.ne.proc_string_old) ) THEN
               	  DO ik = 0, numws-1
               	  WRITE(temp_str,'(I3.3)') ik
                     CALL copy_txtfile('coil_spline'//TRIM(temp_str)//'_'//TRIM(proc_string_old)//'.out',&
                                       'coil_spline'//TRIM(temp_str)//'_'//TRIM(proc_string)//'.out')
                  END DO
               END IF
!DEC$ ENDIF
               ! Keep minimum states
               IF (lkeep_mins) THEN
                  WRITE(temp_str,'(i5.5)') ncnt
                  proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
                  SELECT CASE(TRIM(equil_type))
                     CASE('parvmec','paravmec','vmec2000','vboot')
                         CALL stellopt_paraexe('paravmec_write',proc_string,.false.)
                  END SELECT
                  ier=vctrl_array(2)
                  iflag = ier
                  IF (ier == successful_term_flag) iflag = 0
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
                  CALL copy_boozer_file(TRIM(proc_string_old),TRIM(proc_string))
                  IF (lcoil_geom) THEN
                     CALL SYSTEM('mv coils.'//TRIM(proc_string_old)//' coils.'//TRIM(proc_string))
                  END IF
                  IF (ANY(sigma_dkes < bigno)) THEN
                     DO ik = 1, nsd
                        WRITE(temp_str,'(A,I3.3)') '_s',ik
                        CALL move_txtfile('dkesout.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'dkesout.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
                        CALL move_txtfile('opt_dkes.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'opt_dkes.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
                        CALL move_txtfile('results.'//TRIM(proc_string_old)//TRIM(ADJUSTL(temp_str)),'results.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str)))
                     END DO
                  END IF
!DEC$ IF DEFINED (BEAMS3D_OPT)
                  IF (ANY(sigma_orbit .lt. bigno)) THEN
                     lverb_beams = .FALSE.
                     id_string_beams = TRIM(proc_string_old)
                     CALL beams3d_read
                     id_string_beams = TRIM(proc_string)
                     CALL beams3d_write('GRID_INIT')
                     CALL beams3d_write('TRAJECTORY_FULL')
                     CALL beams3d_write('DIAG')
                     CALL beams3d_free
                     CALL move_txtfile('beams3d_diag_'//TRIM(proc_string_old)//'.txt',&
                                       'beams3d_diag_'//TRIM(proc_string)//'.txt')
                  END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (COILOPTPP)
                  IF (sigma_coil_bnorm < bigno) THEN
                     CALL move_txtfile('bnorm.'//TRIM(proc_string_old),&
                                       'bnorm.'//TRIM(proc_string))
                     CALL move_txtfile('coilopt_params.'//TRIM(proc_string_old),&
                                       'coilopt_params.'//TRIM(proc_string))
                     CALL copy_txtfile('b_norm_eq_'//TRIM(proc_string_old)//'.dat',&
                                       'b_norm_eq_'//TRIM(proc_string)//'.dat')
                     CALL copy_txtfile('b_norm_final_'//TRIM(proc_string_old)//'.dat',&
                                       'b_norm_final_'//TRIM(proc_string)//'.dat')
                     CALL move_txtfile('b_norm_init_'//TRIM(proc_string_old)//'.dat',&
                                       'b_norm_init_'//TRIM(proc_string)//'.dat')
                     DO ik = 0, numws-1
                        WRITE(temp_str,'(I3.3)') ik
                        CALL copy_txtfile('coil_spline'//TRIM(temp_str)//'_'//TRIM(proc_string_old)//'.out',&
                                          'coil_spline'//TRIM(temp_str)//'_'//TRIM(proc_string)//'.out')
                     END DO
                  END IF
!DEC$ ENDIF
!DEC$ IF DEFINED (REGCOIL)
                  ! Currently inside of LEV and GADE cleanup loop, and 
                  ! 'Keeping the mins' section
                  IF ( ANY(sigma_regcoil_chi2_b < bigno) .and. &
                     ( ANY(lregcoil_rcws_rbound_c_opt) .or. ANY(lregcoil_rcws_rbound_s_opt) .or. &
                       ANY(lregcoil_rcws_zbound_c_opt) .or. ANY(lregcoil_rcws_zbound_s_opt) ) ) THEN
                     !print *, '<---In LEV/GADE cleanup.'
                     !print *, '<---proc_string_old = ', proc_string_old
                     !print *, '<---proc_string = ', proc_string
                     CALL copy_txtfile('regcoil_nescout.'//TRIM(proc_string_old),&
                                       'regcoil_nescout.'//TRIM(proc_string))
                  END IF
!DEC$ ENDIF
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
                        !CALL move_txtfile('terpsichore_23'//TRIM(proc_string_old)//TRIM(temp_str),&
                        !                  'terpsichore_23'//TRIM(proc_string)//TRIM(temp_str))
                     END DO
                  END IF
!DEC$ ENDIF
               END IF
               ! Now open the Output file
               ALLOCATE(fvec_temp(mtargets))

               !WRITE COIL KNOT FILE
               IF (ANY(lcoil_spline)) THEN
                  CALL safe_open(iunit_out,iflag,TRIM('knots.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
                  IF (ncnt == 0) WRITE(iunit_out,'(A)') 'COIL KNOTS'
                  WRITE(iunit_out,'(A,1X,I5.5)') 'ITER',ncnt
                  DO n = LBOUND(lcoil_spline,DIM=1), UBOUND(lcoil_spline,DIM=1)
                     IF (ANY(lcoil_spline(n,:))) THEN
                        WRITE(iunit_out,'(2X,A,2X,I5.5)') 'COIL', n
                        ik = MINLOC(coil_splinesx(n,:),DIM=1) - 1
                        IF (lwindsurf) THEN
                           WRITE(iunit_out,"(4X,'u =',4(2X,ES22.12E3))") (coil_splinefx(n,m), m = 1, ik)
                           WRITE(iunit_out,"(4X,'v =',4(2X,ES22.12E3))") (coil_splinefy(n,m), m = 1, ik)
                        ELSE
                           WRITE(iunit_out,"(4X,'x =',4(2X,ES22.12E3))") (coil_splinefx(n,m), m = 1, ik)
                           WRITE(iunit_out,"(4X,'y =',4(2X,ES22.12E3))") (coil_splinefy(n,m), m = 1, ik)
                           WRITE(iunit_out,"(4X,'z =',4(2X,ES22.12E3))") (coil_splinefz(n,m), m = 1, ik)
                        END IF !lwindsurf
                     END IF
                  END DO !n
                  CLOSE(iunit_out)
               END IF

               CALL safe_open(iunit_out,iflag,TRIM('stellopt.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
               iflag = 1
               IF (ncnt == 0) WRITE(iunit_out,'(A,1X,F5.2)') 'VERSION',STELLOPT_VERSION
               WRITE(iunit_out,'(A,1X,I5.5)') 'ITER',ncnt
               CALL stellopt_load_targets(mtargets,fvec_temp,iflag,ncnt)          ! Count
               WRITE(iunit_out,'(A,2(2X,I8))') 'TARGETS ',mtargets,1
               WRITE(iunit_out,'(A)') 'TARGETS'
               WRITE(iunit_out,'(ES22.12E3)') targets(1:mtargets)
               WRITE(iunit_out,'(A,2(2X,I8))') 'SIGMAS ',mtargets,1
               WRITE(iunit_out,'(A)') 'SIGMAS'
               WRITE(iunit_out,'(ES22.12E3)') sigmas(1:mtargets)
               WRITE(iunit_out,'(A,2(2X,I8))') 'VALS ',mtargets,1
               WRITE(iunit_out,'(A)') 'VALUES'
               WRITE(iunit_out,'(ES22.12E3)') vals(1:mtargets)
               CLOSE(iunit_out)
               DEALLOCATE(fvec_temp)
            ELSE IF (ctype == JAC_CLEANUP) THEN
            ELSE IF (ctype == JUST_INPUT) THEN
               ! Write the input file
               WRITE(temp_str,'(i5.5)') ncnt
               proc_string = TRIM(id_string) // '.' // TRIM(ADJUSTL(temp_str))
               CALL safe_open(iunit_out,iflag,TRIM('input.'//TRIM(proc_string)),'unknown','formatted')
               CALL write_indata_namelist(iunit_out,ier)
               CALL write_optimum_namelist(iunit_out,ier)
               IF (lneed_magdiag) CALL write_diagno_input(iunit_out,ier)
!DEC$ IF DEFINED (TXPORT_OPT)
               IF (ANY(sigma_txport < bigno)) CALL write_gist_namelist(iunit_out,ier)
!DEC$ ENDIF
               IF (ANY(sigma_bootstrap < bigno)) CALL write_bootsj_input(iunit_out,ier)
!DEC$ IF DEFINED (NEO_OPT)
               IF (ANY(sigma_neo < bigno)) CALL write_neoin_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
                  IF (ANY(sigma_orbit < bigno)) CALL write_beams3d_namelist(iunit_out,ier)
!DEC$ ENDIF
               WRITE(iunit_out,'(A)') '&END'
               CLOSE(iunit_out)
            ELSE IF (ctype == LAST_GO) THEN
               CALL safe_open(iunit_out,iflag,TRIM('input.'//TRIM(id_string)//'_min'),'unknown','formatted')
               CALL write_indata_namelist(iunit_out,ier)
               CALL write_optimum_namelist(iunit_out,ier)
               IF (lneed_magdiag) CALL write_diagno_input(iunit_out,ier)
!DEC$ IF DEFINED (TXPORT_OPT)
!                  IF (ANY(sigma_txport < bigno)) CALL write_gist_namelist(iunit_out,ier)
!DEC$ ENDIF
                  IF (ANY(sigma_bootstrap < bigno)) CALL write_bootsj_input(iunit_out,ier)
!DEC$ IF DEFINED (NEO_OPT)
                  IF (ANY(sigma_neo < bigno)) CALL write_neoin_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
                  IF (ANY(sigma_orbit < bigno)) CALL write_beams3d_namelist(iunit_out,ier)
!DEC$ ENDIF
!DEC$ IF DEFINED (REGCOIL)
                 ! JCS to do: If needed put regcoil items here.
!DEC$ ENDIF

               WRITE(iunit_out,'(A)') '&END'
               CLOSE(iunit_out)
               ! Overwrite the restart file
               proc_string = 'reset_file'
               vctrl_array(1) = output_flag ! Output to file
               vctrl_array(2) = 0     ! vmec error flag  
               vctrl_array(3) = 0    ! Use multigrid
               vctrl_array(4) = 0     ! Iterative 
               vctrl_array(5) = myid ! Output file sequence number
               !IF (TRIM(equil_type)=='animec') vctrl_array(1) = vctrl_array(1) + animec_flag
               !IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') vctrl_array(1) = vctrl_array(1) + flow_flag
               CALL runvmec(vctrl_array,proc_string,.false.,'')
               ier=vctrl_array(2)
               iflag = ier
               IF (ier == successful_term_flag) iflag = 0
               ! Now open the Output file
               ALLOCATE(fvec_temp(mtargets))
               CALL safe_open(iunit_out,iflag,TRIM('stellopt.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
               iflag = 1
               WRITE(iunit_out,'(A)') 'ITER MIN'
               ik = 999
               CALL stellopt_load_targets(mtargets,fvec_temp,iflag,ik)    ! Count
               WRITE(iunit_out,'(A,2(2X,I8))') 'TARGETS ',mtargets,1
               WRITE(iunit_out,'(A)') 'TARGETS'
               WRITE(iunit_out,'(ES22.12E3)') targets(1:mtargets)
               WRITE(iunit_out,'(A,2(2X,I8))') 'SIGMAS ',mtargets,1
               WRITE(iunit_out,'(A)') 'SIGMAS'
               WRITE(iunit_out,'(ES22.12E3)') sigmas(1:mtargets)
               WRITE(iunit_out,'(A,2(2X,I8))') 'VALS ',mtargets,1
               WRITE(iunit_out,'(A)') 'VALUES'
               WRITE(iunit_out,'(ES22.12E3)') vals(1:mtargets)
               CLOSE(iunit_out)
               DEALLOCATE(fvec_temp)
               IF (.not.lkeep_extra) THEN
                  temp_max = max(nvars,numprocs)
                  DO ik = 0, temp_max, numprocs-1
                     WRITE(temp_str,'(i5)') ik
                     proc_string = TRIM(id_string) // '_opt' // TRIM(ADJUSTL(temp_str))
                     OPEN(iunit,FILE=TRIM('answers_plot.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('answers.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('boozmn_'//TRIM(proc_string)//'.nc'),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('diagno_bth.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('diagno_flux.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('jBbs.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('mercier.'//TRIM(proc_string)),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('jxbout_'//TRIM(proc_string)//'.nc'),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                     OPEN(iunit,FILE=TRIM('wout_'//TRIM(proc_string)//'.nc'),STATUS='unknown',IOSTAT=ier)
                     IF (ier == 0) CLOSE(iunit,STATUS='delete')
                  END DO
               END IF
            END IF
         CASE('spec')
      END SELECT
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_clean_up
