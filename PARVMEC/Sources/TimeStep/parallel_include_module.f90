      MODULE parallel_include_module

        USE stel_kinds
        USE vmec_input, ONLY:lfreeb
        USE mpi_inc

        USE parallel_vmec_module, ONLY: PARVMEC
        USE parallel_vmec_module, ONLY: LV3FITCALL
        USE parallel_vmec_module, ONLY: SKS_ALLGATHER
        USE parallel_vmec_module, ONLY: BCYCLIC
        USE parallel_vmec_module, ONLY: THOMAS
        USE parallel_vmec_module, ONLY: TOFU

        USE parallel_vmec_module, ONLY: num_grids 
        USE parallel_vmec_module, ONLY: grid_procs
        USE parallel_vmec_module, ONLY: grid_size
        USE parallel_vmec_module, ONLY: grid_time
        USE parallel_vmec_module, ONLY: vgrid_time
        USE parallel_vmec_module, ONLY: f3d_time
        USE parallel_vmec_module, ONLY: f3d_num
        USE parallel_vmec_module, ONLY: RUNVMEC_COMM_WORLD
        USE parallel_vmec_module, ONLY: NS_COMM
        USE parallel_vmec_module, ONLY: VAC_COMM
        USE parallel_vmec_module, ONLY: TWODCOMM
        USE parallel_vmec_module, ONLY: px, py
        USE parallel_vmec_module, ONLY: NS_RESLTN
        USE parallel_vmec_module, ONLY: MPI_ERR
        USE parallel_vmec_module, ONLY: rank, nranks
        USE parallel_vmec_module, ONLY: lactive, vlactive
        USE parallel_vmec_module, ONLY: grank, gnranks
        USE parallel_vmec_module, ONLY: vrank, vnranks
        USE parallel_vmec_module, ONLY: nuvmin, nuvmax
        USE parallel_vmec_module, ONLY: nuv3min, nuv3max
        USE parallel_vmec_module, ONLY: tlglob, trglob
        USE parallel_vmec_module, ONLY: t1lglob, t1rglob
        USE parallel_vmec_module, ONLY: t2lglob, t2rglob
        USE parallel_vmec_module, ONLY: tlglob_arr, trglob_arr
        USE parallel_vmec_module, ONLY: nuv3min_arr,nuv3max_arr 
        USE parallel_vmec_module, ONLY: trow, brow
        USE parallel_vmec_module, ONLY: lcol, rcol
        USE parallel_vmec_module, ONLY: par_ntmax
        USE parallel_vmec_module, ONLY: par_mpol1
        USE parallel_vmec_module, ONLY: par_ntor 
        USE parallel_vmec_module, ONLY: par_ns 
        USE parallel_vmec_module, ONLY: par_nuv 
        USE parallel_vmec_module, ONLY: par_nuv3 

        USE parallel_vmec_module, ONLY: STOPMPI

        USE parallel_vmec_module, ONLY: SetVacuumPartitions
        USE parallel_vmec_module, ONLY: SetVacuumCommunicator


        USE parallel_vmec_module, ONLY: Parallel2Serial2X
        USE parallel_vmec_module, ONLY: Parallel2Serial4X

        USE parallel_vmec_module, ONLY: Serial2Parallel4X

        USE parallel_vmec_module, ONLY: Gather1XArray
        USE parallel_vmec_module, ONLY: Gather2XArray
        USE parallel_vmec_module, ONLY: Gather4XArray

        USE parallel_vmec_module, ONLY: PrintParallelIJSPArray
        USE parallel_vmec_module, ONLY: PrintSerialIJSPArray
        USE parallel_vmec_module, ONLY: PrintParallelMNSPArray
        USE parallel_vmec_module, ONLY: PrintSerialMNSPArray
        USE parallel_vmec_module, ONLY: PrintNSArray
        USE parallel_vmec_module, ONLY: PrintOutLinearArray

        USE parallel_vmec_module, ONLY: bcyclic_comp_time
        USE parallel_vmec_module, ONLY: bcyclic_comm_time
        USE parallel_vmec_module, ONLY: waitall_time
        USE parallel_vmec_module, ONLY: dgemm_time
        USE parallel_vmec_module, ONLY: dgemv_time
        USE parallel_vmec_module, ONLY: dgetrf_time
        USE parallel_vmec_module, ONLY: dgetrs_time
        USE parallel_vmec_module, ONLY: ForwardSolveLoop_time
        USE parallel_vmec_module, ONLY: t

        USE parallel_vmec_module, ONLY: allgather_time
        USE parallel_vmec_module, ONLY: allreduce_time
        USE parallel_vmec_module, ONLY: broadcast_time
        USE parallel_vmec_module, ONLY: sendrecv_time
        USE parallel_vmec_module, ONLY: scatter_time

        LOGICAL :: INFILEOUT=.FALSE.

        !Common timers
        REAL(dp) :: mgrid_file_read_time=0
        REAL(dp) :: total_time=0
        REAL(dp) :: evolve_time=0
        REAL(dp) :: restart_time=0
        REAL(dp) :: eqsolve_time=0
        REAL(dp) :: fileout_time=0
        REAL(dp) :: runvmec_time=0
        REAL(dp) :: init_parallel_time=0
        REAL(dp) :: init_radial_time=0
        REAL(dp) :: allocate_ns_time=0
        REAL(dp) :: profile1d_time=0
        REAL(dp) :: profile3d_time=0
        REAL(dp) :: before_main_loop_time=0
        REAL(dp) :: reset_params_time=0
        REAL(dp) :: vsetup_time=0
        REAL(dp) :: readin_time=0
        REAL(dp) :: fixarray_time=0

        REAL(dp) :: myenvvar_time=0
        REAL(dp) :: init_MPI_time=0
        REAL(dp) :: read_namelist_time=0
        REAL(dp) :: get_args_time=0 
        REAL(dp) :: safe_open_time=0

        REAL(dp) :: jacob1=0, jacob2=0

        REAL(dp) :: res_getfsq=0
        REAL(dp) :: res_scalfor=0

        INTEGER :: nfunct3d=0
        REAL(dp) :: old_vacuum_time=0

        REAL(dp) :: fo_funct3d_time=0
        REAL(dp) :: fo_eqfor_time=0
        REAL(dp) :: fo_wrout_time=0
        REAL(dp) :: fo_prepare_time=0
        REAL(dp) :: fo_par_call_time=0

        !Vacuum, Scalpot variables
        INTEGER :: blksize_scp, numjs_vac
        INTEGER, ALLOCATABLE, DIMENSION(:) :: counts_scp, disps_scp,    &
                                              counts_vac, disps_vac, lindx_scp

        !Parallel timers
        REAL(dp) :: totzsps_time=0
        REAL(dp) :: totzspa_time=0
        REAL(dp) :: jacobian_time=0
        REAL(dp) :: bcovar_time=0
        REAL(dp) :: alias_time=0
        REAL(dp) :: forces_time=0
        REAL(dp) :: tomnsps_time=0
        REAL(dp) :: tomnspa_time=0
        REAL(dp) :: symrzl_time=0
        REAL(dp) :: symforces_time=0
        REAL(dp) :: residue_time=0
        REAL(dp) :: tridslv_time=0
        REAL(dp) :: scalfor_time=0
        REAL(dp) :: funct3d_time=0
        REAL(dp) :: gmres_time=0
        REAL(dp) :: guess_axis_time=0
        REAL(dp) :: vacuum_time=0
        REAL(dp) :: precal_time=0
        REAL(dp) :: surface_time=0
        REAL(dp) :: bextern_time=0
        REAL(dp) :: scalpot_time=0
        REAL(dp) :: solver_time=0
        REAL(dp) :: analyt_time=0
        REAL(dp) :: tolicu_time=0
        REAL(dp) :: belicu_time=0
        REAL(dp) :: becoil_time=0
        REAL(dp) :: greenf_time=0
        REAL(dp) :: fourp_time=0
        REAL(dp) :: fouri_time=0

        REAL(dp) :: init_time=0
        REAL(dp) :: setup_time=0
        REAL(dp) :: forwardsolve_time=0
        REAL(dp) :: backwardsolve_time=0
        REAL(dp) :: finalize_time=0

        !Serial timers
        REAL(dp) :: s_totzsps_time=0
        REAL(dp) :: s_totzspa_time=0
        REAL(dp) :: s_jacobian_time=0
        REAL(dp) :: s_bcovar_time=0
        REAL(dp) :: s_alias_time=0
        REAL(dp) :: s_forces_time=0
        REAL(dp) :: s_tomnsps_time=0
        REAL(dp) :: s_tomnspa_time=0
        REAL(dp) :: s_symrzl_time=0
        REAL(dp) :: s_symforces_time=0
        REAL(dp) :: s_residue_time=0
        REAL(dp) :: s_tridslv_time=0
        REAL(dp) :: s_scalfor_time=0
        REAL(dp) :: s_gmres_time=0
        REAL(dp) :: s_guess_axis_time=0
        REAL(dp) :: s_vacuum_time=0
        REAL(dp) :: s_precal_time=0
        REAL(dp) :: s_surface_time=0
        REAL(dp) :: s_bextern_time=0
        REAL(dp) :: s_scalpot_time=0
        REAL(dp) :: s_solver_time=0
        REAL(dp) :: s_analyt_time=0
        REAL(dp) :: s_tolicu_time=0
        REAL(dp) :: s_belicu_time=0
        REAL(dp) :: s_becoil_time=0
        REAL(dp) :: s_greenf_time=0
        REAL(dp) :: s_fourp_time=0
        REAL(dp) :: s_fouri_time=0

        REAL(rprec) :: maxvalue, minvalue
        INTEGER :: maxrank, minrank

        !V3FIT 
        INTEGER :: RUNVMEC_PASS=0

!        INTEGER :: v3fgcomm, v3freccomm, v3feqcomm
!        INTEGER :: v3fgnranks, v3frecnranks, v3feqnranks
!        INTEGER :: v3fgrank, v3frecrank, v3feqrank 
!        LOGICAL :: lreconactive=.FALSE.

CONTAINS

        !------------------------------------------------
        ! Print out parallel timing information 
        !------------------------------------------------
        SUBROUTINE PrintTimes
          IMPLICIT NONE
          DOUBLE PRECISION :: ind(2,1), outd(2,1)
          DOUBLE PRECISION :: mt(nranks)
          LOGICAL :: LMINMAXFILE=.FALSE.

          IF (PARVMEC) THEN
#if defined(MPI_OPT)
            ind(1,1) = total_time; ind(2,1) = DBLE(rank)
            CALL MPI_Reduce(ind,outd,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,0,NS_COMM,MPI_ERR)
            IF(grank.EQ.0) THEN
              minvalue = outd(1,1); minrank = outd(2,1)
            END IF
            CALL MPI_Reduce(ind,outd,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,0,NS_COMM,MPI_ERR)
            IF(grank.EQ.0) THEN
              maxvalue = outd(1,1); maxrank = outd(2,1)
            END IF
            CALL MPI_Bcast(minvalue,1,MPI_DOUBLE_PRECISION,0,NS_COMM,MPI_ERR)
            CALL MPI_Bcast(maxvalue,1,MPI_DOUBLE_PRECISION,0,NS_COMM,MPI_ERR)
            CALL MPI_Bcast(minrank,1,MPI_INTEGER,0,NS_COMM,MPI_ERR)
            CALL MPI_Bcast(maxrank,1,MPI_INTEGER,0,NS_COMM,MPI_ERR)
            
            IF (.FALSE..AND.ABS(maxvalue-minvalue).GE.0.0) THEN
              LMINMAXFILE=.TRUE.
              IF (grank.EQ.0) THEN
                WRITE(*,*)
                WRITE(*,*)'---------------------------------------------------'
                WRITE(*,*) 'Min and max runtimes exceed tolerance....'
                WRITE(*,'(A12,F15.6,A12,I6)')'Max. time: ',maxvalue,'Rank :',maxrank
                WRITE(*,'(A12,F15.6,A12,I6)')'Min. time: ',minvalue,'Rank :',minrank
                WRITE(*,*)'---------------------------------------------------'
                WRITE(*,*)
              END IF
            END IF
#endif
          ELSE
            totzsps_time = s_totzsps_time
            totzspa_time = s_totzspa_time
            jacobian_time = s_jacobian_time
            bcovar_time = s_bcovar_time
            alias_time = s_alias_time
            forces_time = s_forces_time
            tomnsps_time = s_tomnsps_time
            tomnspa_time = s_tomnspa_time
            symrzl_time = s_symrzl_time
            symforces_time = s_symforces_time
            residue_time = s_residue_time
            tridslv_time = s_tridslv_time
            scalfor_time = s_scalfor_time
            gmres_time = s_gmres_time
            guess_axis_time = s_guess_axis_time
            vacuum_time = s_vacuum_time
            precal_time = s_precal_time
            surface_time = s_surface_time
            bextern_time = s_bextern_time
            scalpot_time = s_scalpot_time
            solver_time = s_solver_time
            analyt_time = s_analyt_time
            tolicu_time = s_tolicu_time
            belicu_time = s_belicu_time
            becoil_time = s_becoil_time
            greenf_time = s_greenf_time
            fourp_time = s_fourp_time
            fouri_time = s_fouri_time
            maxvalue=total_time; minvalue=total_time
          END IF

          IF (PARVMEC.AND.LMINMAXFILE) THEN
            IF(grank.EQ.minrank) THEN
              CALL WriteTimes('min-timings.txt')
            END IF
            IF(grank.EQ.maxrank) THEN
              CALL WriteTimes('max-timings.txt')
            END IF
          END IF

          IF(grank.EQ.0) THEN
            CALL WriteTimes('timings.txt')
          END IF

        END SUBROUTINE PrintTimes 
!------------------------------------------------

!------------------------------------------------
        SUBROUTINE WriteTimes(fname)
          IMPLICIT NONE
          CHARACTER(*), INTENT(IN) :: fname
          CHARACTER(100) :: cfname
          CHARACTER(50) :: ciam, cnprocs
          INTEGER :: TFILE, i, istat, igrd

          WRITE(ciam,*) rank; WRITE(cnprocs,*) nranks
          ciam=ADJUSTL(ciam); cnprocs=ADJUSTL(cnprocs)
          TFILE = 9*nranks+grank+1000
          OPEN(UNIT=TFILE, FILE=fname, STATUS="REPLACE", ACTION="WRITE",&
            &FORM="FORMATTED",POSITION="APPEND", IOSTAT=istat)

          IF (PARVMEC) THEN
            WRITE(TFILE,*) '====================== PARALLEL TIMINGS ====================' 
          ELSE
            WRITE(TFILE,*) '======================= SERIAL TIMINGS =====================' 
          END IF

          WRITE(TFILE,'(A20,A4,F15.6)') 'total', ': ', total_time
          WRITE(TFILE,'(A20,A4,I15)') 'rank', ': ', rank
          WRITE(TFILE,'(A20,A4,F15.6)') 'mgrid file read time', ': ', mgrid_file_read_time
          WRITE(TFILE,'(A20,A4,I15)') 'No. of procs', ': ', gnranks
          WRITE(TFILE,*) 

          IF (lfreeb) THEN
            DO igrd=1, num_grids
              WRITE(TFILE,'(A20,A4,3I15,F15.6)') '--- non-vacuum', ': ',&
                f3d_num(igrd), grid_size(igrd), grid_procs(igrd), f3d_time(igrd)-vgrid_time(igrd)
            END DO
            WRITE(TFILE,*) 

            WRITE(TFILE,'(A20,A4,I15)') 'VNRANKS ',': ', vnranks 
            DO igrd=1, num_grids
              WRITE(TFILE,'(A20,A4,I15,F15.6)') '--- vacuum ', ': ',&
                grid_size(igrd),vgrid_time(igrd)
            END DO
            WRITE(TFILE,*) 
          ELSE
            DO igrd=1, num_grids
              WRITE(TFILE,'(A20,A4,3I15,F15.6)') '--- non-vacuum', ': ',&
                f3d_num(igrd), grid_size(igrd), grid_procs(igrd), f3d_time(igrd)
            END DO
            WRITE(TFILE,*) 
          END IF
        
          WRITE(TFILE,'(A20,A4,F15.6)') 'runvmec', ': ', runvmec_time
          WRITE(TFILE,*) 

          WRITE(TFILE,'(A20,A4,F15.6)') 'init radial', ': ',init_radial_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'eqsolve', ': ', eqsolve_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'fileout', ': ', fileout_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') 'evolve', ': ', evolve_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'funct3d', ': ', funct3d_time
          WRITE(TFILE,'(A20,A4,I15)') 'nfunct3d', ': ', nfunct3d
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') 'totzsps', ': ', totzsps_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'totzspa', ': ', totzspa_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'symrzl', ': ', symrzl_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'jacobian', ': ', jacobian_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'bcovar', ': ', bcovar_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'vacuum', ': ', vacuum_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') '- precal', ': ', precal_time
          WRITE(TFILE,'(A20,A4,F15.6)') '- surface', ': ', surface_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') '- bextern', ': ', bextern_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') '-- becoil', ': ', becoil_time
          WRITE(TFILE,'(A20,A4,F15.6)') '-- tolicu', ': ', tolicu_time
          WRITE(TFILE,'(A20,A4,F15.6)') '-- belicu', ': ', belicu_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') '- scalpot', ': ', scalpot_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') '-- analyt', ': ', analyt_time
          WRITE(TFILE,'(A20,A4,F15.6)') '-- greenf', ': ', greenf_time
          WRITE(TFILE,'(A20,A4,F15.6)') '-- fourp', ': ', fourp_time
          WRITE(TFILE,'(A20,A4,F15.6)') '-- fouri', ': ', fouri_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') '- solver', ': ', solver_time
          WRITE(TFILE,*) 
          WRITE(TFILE,'(A20,A4,F15.6)') 'alias', ': ', alias_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'forces', ': ', forces_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'symforces', ': ', symforces_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'tomnsps', ': ', tomnsps_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'tomnspa', ': ', tomnspa_time
          WRITE(TFILE,'(A20,A4,F15.6)') 'residue', ': ', residue_time
          WRITE(TFILE,'(A20,A4,F15.6)') '-- tridslv', ': ', tridslv_time
          WRITE(TFILE,*) 
          WRITE(TFILE,*) '============================================================' 
          IF (PARVMEC) THEN
            WRITE(TFILE,*) 
            WRITE(TFILE,'(A20,A4,F15.6)') 'allgather', ': ', allgather_time
            WRITE(TFILE,'(A20,A4,F15.6)') 'allreduce', ': ', allreduce_time
            WRITE(TFILE,'(A20,A4,F15.6)') 'broadcast', ': ', broadcast_time
            WRITE(TFILE,'(A20,A4,F15.6)') 'sendrecv', ': ', sendrecv_time
            WRITE(TFILE,*) 
            WRITE(TFILE,*) '============================================================' 
          END IF
          CALL FLUSH(TFILE)
          CLOSE(TFILE)

        END SUBROUTINE WriteTimes
!------------------------------------------------

      END MODULE parallel_include_module
!------------------------------------------------

