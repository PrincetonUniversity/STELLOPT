!-----------------------------------------------------------------------
!     Subroutine:    stellopt_paraexe
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/18/14
!     Description:   This subroutine is only called by those
!                    processes which will execute parallel sub
!                    codes.  These processes do not execute
!                    the main blocks of the stellopt code but
!                    rather, they sit around waiting for work from
!                    their master.
!     GENE NOTE:     Must compile GENE with -assume noold_unit_star 
!                    under INTEL (hydra).
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_paraexe(in_parameter_1,in_parameter_2,lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_vals, ONLY: kx_gene
      USE wall_mod, ONLY: wall_free
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
      USE vmec_input, ONLY:  ns_array, lfreeb, write_indata_namelist
!DEC$ IF DEFINED (GENE)
      USE gene_subroutine, ONLY: rungene
      USE par_in, ONLY: diagdir,file_extension, beta_gene => beta
      USE par_other, ONLY: print_ini_msg
      USE parameters_IO, ONLY: read_parameters, write_parameters, PAR_OUT_FILE, &
                               PARFILE
      USE stellopt_interface, ONLY: bcast_parameters
      USE geometry, ONLY: geomdir, magn_geometry, geomfile
      USE file_io, ONLY: erase_stop_file
      USE discretization, ONLY: mype_gl, n_procs_sim, n_procs_s, n_procs_w, &
                                n_procs_v, n_procs_x, n_procs_y, n_procs_z
      USE communications, ONLY: omp_level, my_mpi_comm_world
      USE coordinates, ONLY: kx_center
!DEC$ ENDIF
!DEC$ IF DEFINED (BEAMS3D_OPT)
      ! BEAMS3D Libraries
      USE beams3d_runtime, ONLY: id_string_beams => id_string, &
            lread_input_beams => lread_input, lvmec_beams => lvmec, &
            lverb_beams => lverb, lbeam_beams => lbeam, &
            lpies_beams => lpies, lspec_beams => lspec, &
            lmgrid_beams => lmgrid, lascot_beams => lascot, &
            lvessel_beams => lvessel, lcoil_beams => lcoil, &
            lrestart_grid_beams => lrestart_grid, lrestart_particles_beams => lrestart_particles, &
            lbeam_simple_beams => lbeam_simple, &
            lplasma_only_beams => lplasma_only, lascot4_beams => lascot4, &
            lbbnbi_beams => lbbnbi, lascotfl_beams => lascotfl, &
            lcollision_beams => lcollision, lw7x_beams => lw7x, &
            coil_string_beams => coil_string, mgrid_string_beams => mgrid_string,&
            vessel_string_beams => vessel_string, restart_string_beams => restart_string, &
            lraw_beams => lraw, nbeams_beams => nbeams, &
            lvac_beams => lvac, lhitonly, nparticles_start, &
            vll_start_in, R_start_in, Z_start_in, PHI_start_in, mu_start_in, &
            mu_start_in, charge_in, mass_in, t_end_in, Zatom_in, &
            TE_AUX_S_BEAMS => TE_AUX_S, TE_AUX_F_BEAMS => TE_AUX_F, &
            NE_AUX_S_BEAMS => NE_AUX_S, NE_AUX_F_BEAMS => NE_AUX_F, &
            TI_AUX_S_BEAMS => TI_AUX_S, TI_AUX_F_BEAMS => TI_AUX_F, nprocs_beams, &
            ZEFF_AUX_S_BEAMS => ZEFF_AUX_S, ZEFF_AUX_F_BEAMS => ZEFF_AUX_F
      USE beams3d_lines, ONLY: nparticles_beams => nparticles, R_lines, Z_lines,&
            PHI_lines, vll_lines, moment_lines, neut_lines
      USE beams3d_grid, ONLY: nte, nne, nti, B_R, B_PHI, B_Z, raxis, zaxis, phiaxis,&
                              BR_spl, BZ_spl, BPHI_spl, MODB_spl, rmin, rmax, zmin, &
                              zmax, phimin, phimax, nzeff
      USE wall_mod, ONLY: wall_free
      USE beams3d_input_mod, ONLY: BCAST_BEAMS3D_INPUT
!DEC$ ENDIF
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(inout)    :: in_parameter_1
      CHARACTER(LEN=*), INTENT(inout)    :: in_parameter_2
      LOGICAL, INTENT(inout)        :: lscreen
      
!-----------------------------------------------------------------------
!     Local Variables
!        lscreen         Controls printing to the screen
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
      IF (ier_paraexe /= 0) RETURN
      code_str = TRIM(in_parameter_1)
      file_str = TRIM(in_parameter_2)
      ierr_mpi = 0
      print *,"Hello from stellopt_paraexe. code_str=",trim(code_str)
      DO
         ! First get the name of the code blah
         ier_paraexe = 0; ierr_mpi = 0; ier = 0
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BARRIER1',ierr_mpi)
         CALL MPI_BCAST(code_str,256,MPI_CHARACTER,master,MPI_COMM_MYWORLD,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BCAST1',ierr_mpi)
         CALL MPI_BCAST(file_str,256,MPI_CHARACTER,master,MPI_COMM_MYWORLD,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BCAST2',ierr_mpi)
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
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BARRIER',ierr_mpi)
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
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BCAST2b',ierr_mpi)
               ! Now update the namelists
               IF (myworkid == master) CALL stellopt_prof_to_vmec(file_str,ier)
               CALL stellopt_bcast_vmec(master,MPI_COMM_MYWORLD,ier)
               IF (ier .eq. 0) THEN
                  CALL stellopt_reinit_vmec
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
                  CALL runvmec(ictrl,file_str,lscreen,MPI_COMM_MYWORLD,reset_string)
                  CALL FinalizeSurfaceComm(NS_COMM)
                  CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD)
                  ier=ictrl(2)
                  CALL MPI_BCAST(ier,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BCAST2d',ierr_mpi)
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
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BCAST2b',ierr_mpi)
               ictrl(1) = output_flag
               ictrl(2) = 0     ! vmec error flag  
               ictrl(3) = 0    ! Use multigrid
               ictrl(4) = 0
               ictrl(5) = myseq ! Output file sequence number
               reset_string =''
               NS_RESLTN = 0 ! Need to do this otherwise situations arrise which cause problems.
               CALL runvmec(ictrl,file_str,lscreen,MPI_COMM_MYWORLD,reset_string)
               LIFFREEB  = .FALSE. ! Already deallocated from before and we need to reset stuff
               CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD) ! We don't allocate the vacuum communicator when we write
               ier=ictrl(2)
               CALL MPI_BCAST(ier,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: BCAST2d',ierr_mpi)
               ier_paraexe = ier
!DEC$ ENDIF
            CASE('gene_parallel')  ! Parallel Gene
!DEC$ IF DEFINED (MPI_OPT) .AND. DEFINED (GENE)
               ! First read the parameters file
               file_extension = ''
               geomdir = '.'
               magn_geometry = 'gist'
               geomfile = '.'
               beta_gene = 0.0
               print_ini_msg = .false.
               PARFILE = 30
               CALL read_parameters('')

               ! Setup the communicator
               my_mpi_comm_world = MPI_COMM_MYWORLD
               omp_level = MPI_THREAD_SINGLE
               mype_gl = myworkid
               CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, n_procs_sim, ierr_mpi )
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe',ierr_mpi)

               ! Initialize the code
               file_extension = TRIM(file_str)
               diagdir = './'
               geomdir = '.'
               magn_geometry = 'gist'
               gene_dir_str = '.'
               gene_str = TRIM(file_str)
               geomfile = 'gist_genet'//TRIM(file_str)
               ch_in = 'no'
               kx_center = kx_gene

               ! Do the initalization call
               IF (myworkid == master) THEN
                  CALL write_parameters
                  CALL read_parameters(gene_dir_str)
               END IF
               CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe',ierr_mpi)

               ! Update everyone
               CALL bcast_parameters
               gene_dir_str = 'skip_parfile'  ! This causes rungene to skip reading the parameters file.
               ! We need to write the parameters file so we can read it
               !CALL write_parameters
               !PARFILE = 30
               !IF (myworkid == master) CALL write_parameters

               ! Run the Code
               IF (lscreen) print_ini_msg = .true.
!!!!!!!!!!!!!!  The following line crashes the code
               IF (lscreen) WRITE(6,'(a)') ' ::::::::::::::::::::        LINEAR GENE RUN     ::::::::::::::::::::::'
               IF (.not.lscreen .and. myworkid == master) THEN
                  CLOSE(UNIT=6)
                  OPEN(UNIT=6,FILE='log_gene.'//TRIM(gene_str(2:128)))
               END IF
               CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe',ierr_mpi)
               CALL rungene(MPI_COMM_MYWORLD,gene_dir_str,gene_str,ch_in)
               CALL erase_stop_file
               IF (.not.lscreen .and. myworkid == master) THEN 
                  CLOSE(UNIT=6)
                  OPEN(UNIT=6,FILE=TRIM(screen_str),POSITION='APPEND')
               END IF
               IF (lscreen) WRITE(6,'(a)') ' ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
               CALL FLUSH(6)
!DEC$ ENDIF
            CASE('beams3d')
!DEC$ IF DEFINED (BEAMS3D_OPT)
               CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, nprocs_beams, ierr_mpi )
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe',ierr_mpi)

               ! Set vars so BEAMS3D knows it's being called from stellopt
               CALL MPI_COMM_DUP(MPI_COMM_MYWORLD, MPI_COMM_BEAMS, ierr_mpi)
               CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_BEAMS, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
               CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
               CALL BCAST_BEAMS3D_INPUT(master,MPI_COMM_MYWORLD,ierr_mpi)

               lverb_beams        = .FALSE.
               lvmec_beams        = .TRUE.  ! Use VMEC Equilibria
               lpies_beams        = .FALSE.
               lspec_beams        = .FALSE.
               lcoil_beams        = .FALSE.
               lmgrid_beams       = .FALSE.
               lascot_beams       = .FALSE.
               lascotfl_beams     = .FALSE.
               lascot4_beams      = .FALSE.
               lbbnbi_beams       = .FALSE.
               lraw_beams         = .FALSE.
               lvessel_beams      = .FALSE.
               lvac_beams         = .FALSE.
               lrestart_grid_beams     = .FALSE.
               lrestart_particles_beams     = .FALSE.
               lbeam_simple_beams = .FALSE.
               lhitonly           = .TRUE. ! Set to true to smaller files.
               IF (lscreen) lhitonly = .FALSE.
               lplasma_only_beams = .TRUE.
               lbeam_beams        = .FALSE.
               lread_input_beams  = .FALSE.
               lcollision_beams   = .FALSE.
               lw7x_beams   = .FALSE.
               id_string_beams    = TRIM(file_str)
               coil_string_beams  = ''
               mgrid_string_beams = ''
               vessel_string_beams = ''
               restart_string_beams = ''
               IF (myworkid .eq. master) lverb_beams = lscreen


               ! Broadcast some variables
               ierr_mpi = 0
               CALL MPI_BCAST(lhitonly,1,MPI_LOGICAL, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(nne,1,MPI_INTEGER, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(nte,1,MPI_INTEGER, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(nti,1,MPI_INTEGER, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(nzeff,1,MPI_INTEGER, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(nparticles_start,1,MPI_INTEGER, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(rmin,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(rmax,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(zmin,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(zmax,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(phimin,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(phimax,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(TE_AUX_S_BEAMS,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(TE_AUX_F_BEAMS,nte,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(NE_AUX_S_BEAMS,nne,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(NE_AUX_F_BEAMS,nne,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(TI_AUX_S_BEAMS,nti,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(TI_AUX_F_BEAMS,nti,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(ZEFF_AUX_S_BEAMS,nzeff,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(ZEFF_AUX_F_BEAMS,nzeff,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(R_start_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(Z_start_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(PHI_start_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(vll_start_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(charge_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(mass_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(Zatom_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(mu_start_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(t_end_in,nparticles_start,MPI_REAL8, master, MPI_COMM_MYWORLD,ierr_mpi)
               nparticles_beams = nparticles_start

               ! Print some stuff to screen
               IF (lscreen) THEN
                  WRITE(6, '(a,f5.2)') 'BEAMS3D Version ', BEAMS3D_VERSION
               END IF

               ! Initialize the code
               CALL beams3d_init
               CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe',ierr_mpi)

               ! Follow particles
               CALL beams3d_follow
               nbeams_beams = 1  ! Do this so the read in cleanup doesn't fail
               CALL beams3d_write('TRAJECTORY_PARTIAL')

               ! Calculated particle locations in flux space
               CALL beams3d_diagnostics

               ! Deallocate Arrays
               CALL beams3d_free(MPI_COMM_SHARMEM)
               CALL wall_free(ier,MPI_COMM_SHARMEM)
              
               ! Free the Shared Memory region
               CALL MPI_COMM_FREE(MPI_COMM_SHARMEM,ierr_mpi)
               CALL MPI_COMM_FREE(MPI_COMM_BEAMS,ierr_mpi)

               IF (lverb_beams) WRITE(6, '(A)') '----- BEAMS3D DONE -----'

!DEC$ ENDIF
!DEC$ IF DEFINED (TRAVIS)
            CASE('travis')
               proc_string = file_str
               CALL stellopt_travis(lscreen,ier)
!DEC$ ENDIF
            CASE('coilopt++')
               CALL stellopt_coiloptpp(file_str,lscreen)
            CASE('regcoil_chi2_b')
               CALL stellopt_regcoil_chi2_b(lscreen,ier)
            CASE('terpsichore')
               proc_string = file_str
               ier = 0
               CALL stellopt_kink(lscreen,ier)
            CASE('booz_xform')
               proc_string = file_str
               ier = 0
               CALL stellopt_toboozer(lscreen,ier)
               ier_paraexe = ier
            CASE('bootsj')
               proc_string = file_str
               ier = 0
               CALL stellopt_bootsj(lscreen,ier)
               ier_paraexe = ier
            CASE('sfincs')
               proc_string = file_str
               ier = 0
               CALL stellopt_sfincs(lscreen,ier)
               ier_paraexe = ier
            CASE('cobra')
               proc_string = file_str
               ier = 0
               CALL stellopt_balloon(lscreen,ier)
               ier_paraexe = ier
            CASE('diagno')
               proc_string = file_str
               ier = 0
               CALL stellopt_magdiag(lscreen,ier)
               ier_paraexe = ier
            CASE('neo')
               proc_string = file_str
               ier = 0
               CALL stellopt_neo(lscreen,ier)
               ier_paraexe = ier
            CASE('write_mgrid')
               CALL stellopt_write_mgrid(MPI_COMM_MYWORLD,file_str,lscreen)
            CASE('mango_init')
               CALL stellopt_mango_init
            CASE('mango_finalize')
               CALL stellopt_mango_finalize
            CASE('exit')  ! we send this when we want to terminate the code (everyone leaves)
               !PRINT *,'myid: ',myid,' exiting stellopt_paraexe'
               CALL MPI_COMM_FREE(MPI_COMM_MYWORLD,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_paraexe: FREE',ierr_mpi)
               RETURN
            CASE DEFAULT
               PRINT *,"Error! stellopt_paraexe called with unknown argument: ",TRIM(code_str)
               STOP
         END SELECT
         !lscreen = .false.
         IF (myworkid == master) RETURN ! The master process of the Communicator can leave
      END DO
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_paraexe
