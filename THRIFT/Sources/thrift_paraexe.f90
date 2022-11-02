!-----------------------------------------------------------------------
!     Subroutine:    thrift_paraexe
!     Authors:       L. van Ham, S. Lazerson
!     Date:          11/XX/22
!     Description:   This subroutine is only called by those
!                    processes which will execute parallel sub
!                    codes.  These processes do not execute
!                    the main blocks of the thrift code but
!                    rather, they sit around waiting for work from
!                    their master.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_paraexe(in_parameter_1,in_parameter_2,lscreen_local)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_input_mod
      USE mpi_params
      USE mpi_inc
!DEC$ IF DEFINED (SKS2)
      USE parallel_vmec_module, ONLY: &
            InitializeParallel, FinalizeParallel, grank, &
            PARVMEC, LV3FITCALL,&
            InitRunVmec,FinalizeRunVmec,RUNVMEC_COMM_WORLD,NS_RESLTN,&
            FinalizeSurfaceComm, NS_COMM, LIFFREEB
!DEC$ ENDIF
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
      USE vmec_input, ONLY:  ns_array, lfreeb, write_indata_namelist, &
                             lnyquist, bcast_indata_namelist
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(inout)    :: in_parameter_1
      CHARACTER(LEN=*), INTENT(inout)    :: in_parameter_2
      LOGICAL, INTENT(inout)        :: lscreen_local
      
!-----------------------------------------------------------------------
!     Local Variables
!        lscreen_local   Controls printing to the screen
!        pass            Counts how many times this code is called.
!        gene_dir_string GENE parameter file directory
!        gene_str        GENE parameter file extension
!        ch_in           GENE Checkpoint string
!        code_str         Code to run
!        file_str         Filename to pass to the code
!----------------------------------------------------------------------
      CHARACTER(256) :: code_str,file_str
      CHARACTER(128) :: gene_dir_str, gene_str, ch_in
      CHARACTER(len = 128)    :: reset_string
      INTEGER        :: dex, ier, iunit
      INTEGER :: ictrl(5)
      LOGICAL :: lhit
      INTEGER :: i, myseq
      DOUBLE PRECISION :: x0,y0,z0,x1,y1,z1,xw,yw,zw, rr0,zz0,phi0
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (TRIM(in_parameter_1) /= 'exit' .and. ier_paraexe /= 0) RETURN
      code_str = TRIM(in_parameter_1)
      file_str = TRIM(in_parameter_2)
      ierr_mpi = 0
      DO
         ! First get the name of the code blah
         ier_paraexe = 0; ierr_mpi = 0; ier = 0
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BARRIER1',ierr_mpi)
         CALL MPI_BCAST(code_str,256,MPI_CHARACTER,master,MPI_COMM_MYWORLD,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BCAST1',ierr_mpi)
         CALL MPI_BCAST(file_str,256,MPI_CHARACTER,master,MPI_COMM_MYWORLD,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BCAST2',ierr_mpi)
!DEC$ ENDIF

         ! Now run the proper code
         CALL tolower(code_str)
         SELECT CASE (TRIM(code_str))
            CASE('parvmec_init')
               myseq=myid
               CALL MPI_BCAST(myseq,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
              ! Now make initializing VMEC call which preforms allocations
               ictrl(1) = restart_flag + readin_flag + reset_jacdt_flag
               ictrl(2) = 0
               ictrl(3) = 50
               ictrl(4) = 0
               ictrl(5) = myid
               PARVMEC = .TRUE.
               NS_RESLTN = 0 ! Need to do this otherwise situations arrise which cause problems.
               reset_string =''
               !IF (TRIM(equil_type)=='animec') ictrl(1) = ictrl(1) + animec_flag
               !IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') ictrl(1) = ictrl(1) + flow_flag
               IF (myworkid==master) THEN
                  CALL safe_open(iunit,ier,'threed1.'//TRIM(file_str),'unknown','formatted')
                  CLOSE(iunit)
               END IF
               CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BARRIER',ierr_mpi)
               CALL runvmec(ictrl,file_str,.false.,MPI_COMM_MYWORLD,reset_string)
               CALL FinalizeSurfaceComm(NS_COMM)
               CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD)
               ier_paraexe=ictrl(2)
               in_parameter_2 = TRIM(file_str)
            CASE('paravmec_run')
!DEC$ IF DEFINED (SKS2)
               ! Broadcast the sequence number
               myseq=myid          ! MPI
               CALL MPI_BCAST(myseq,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BCAST2b',ierr_mpi)
               ! Now update the namelists
               !CALL bcast_indata_namelist(master,MPI_COMM_MYWORLD,ier)
               IF (ier .eq. 0) THEN
                  !CALL thrift_reinit_vmec
                  IF (myworkid == master) THEN
                     iunit = 37; ier = 0
                     CALL safe_open(iunit,ier,TRIM('temp_input.'//TRIM(file_str)),'unknown','formatted')
                     CALL write_indata_namelist(iunit,ier)
                     CALL FLUSH(iunit)
                  END IF
                  ! Setup ICTRL Array
                  ictrl(1) = restart_flag+timestep_flag+output_flag+reset_jacdt_flag
                  ictrl(2) = 0; ictrl(3) = -1; ictrl(4) = -1; ictrl(5) = myseq
                  PARVMEC = .TRUE.
                  ! Setup reset_string
                  reset_string =''
                  lhit = .FALSE.
                  NS_RESLTN = 0 ! Need to do this otherwise situations arrise which cause problems.
                  CALL runvmec(ictrl,file_str,lscreen_local,MPI_COMM_MYWORLD,reset_string)
                  CALL FinalizeSurfaceComm(NS_COMM)
                  CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD)
                  ier=ictrl(2)
                  CALL MPI_BCAST(ier,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BCAST2d',ierr_mpi)
                  IF (  ier == successful_term_flag  .or. &
                        ier == norm_term_flag) THEN
                     IF (myworkid == master) CLOSE(UNIT=iunit,STATUS='delete')
                     ier = 0
                  ELSE
                     IF (myworkid == master) CLOSE(UNIT=iunit)
                     ier = -1
                  END IF
               END IF
               in_parameter_2 = TRIM(file_str)
               ier_paraexe = ier
!DEC$ ENDIF
            CASE('paravmec_write')
!DEC$ IF DEFINED (SKS2)
               CALL MPI_BCAST(myseq,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BCAST2b',ierr_mpi)
               ictrl(1) = output_flag
               ictrl(2) = 0     ! vmec error flag  
               ictrl(3) = 0    ! Use multigrid
               ictrl(4) = 0
               ictrl(5) = myseq ! Output file sequence number
               reset_string =''
               NS_RESLTN = 0 ! Need to do this otherwise situations arrise which cause problems.
               CALL runvmec(ictrl,file_str,lscreen_local,MPI_COMM_MYWORLD,reset_string)
               LIFFREEB  = .FALSE. ! Already deallocated from before and we need to reset stuff
               CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD) ! We don't allocate the vacuum communicator when we write
               ier=ictrl(2)
               CALL MPI_BCAST(ier,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: BCAST2d',ierr_mpi)
               ier_paraexe = ier
!DEC$ ENDIF
!DEC$ IF DEFINED (TRAVIS)
!            CASE('travis')
!               proc_string = file_str
!               CALL thrift_travis(lscreen_local,ier)
!DEC$ ENDIF            
            CASE('booz_xform')
               proc_string = file_str
               ier = 0
               !CALL thrift_toboozer(lscreen_local,ier)
               ier_paraexe = ier
            CASE('bootsj')
               proc_string = file_str
               ier = 0
               !CALL thrift_bootsj(lscreen_local,ier)
               ier_paraexe = ier
            CASE('sfincs')
!               proc_string = file_str
!               ier = 0
!               CALL thrift_sfincs(lscreen_local,ier)
!               ier_paraexe = ier
            CASE('exit')  ! we send this when we want to terminate the code (everyone leaves)
               CALL MPI_COMM_FREE(MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'thrift_paraexe: FREE',ierr_mpi)
               RETURN
            CASE DEFAULT
               PRINT *,"Error! thrift_paraexe called with unknown argument: ",TRIM(code_str)
               STOP
         END SELECT
         IF (myworkid == master) RETURN ! The master process of the Communicator can leave
      END DO
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_paraexe
