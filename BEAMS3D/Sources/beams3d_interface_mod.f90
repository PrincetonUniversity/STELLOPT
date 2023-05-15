!-----------------------------------------------------------------------
!     Module:        BEAMS3D_INTERFACE_MOD
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/15/2022
!     Description:   Module provides routines for handling various
!                    initialization tasks.
!-----------------------------------------------------------------------
MODULE BEAMS3D_INTERFACE_MOD
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_runtime
      USE mpi_params
      USE mpi_inc
      USE wall_mod, ONLY: wall_free
#if defined(LHDF5)
      USE hdf5
#endif

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: vmajor, vminor, liblen
      INTEGER :: h5major, h5minor, h5rel, h5par
      INTEGER :: mpi_info_beams3d
      CHARACTER(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_lib_name
      
!-----------------------------------------------------------------------
!     Subroutines
!         beams3d_init_mpi:         MPI initialization
!         beams3d_init_HDF5:        HDF5 initialization
!         beams3d_init_constants:   Constants initialization
!         beams3d_init_commandline: Handle the command line
!-----------------------------------------------------------------------
CONTAINS

      SUBROUTINE beams3d_init_mpi
      IMPLICIT NONE
      myworkid = master
#if defined(MPI_OPT)
      CALL MPI_INIT(ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_INIT_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_BEAMS, myworkid, ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_BEAMS, nprocs_beams, ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_SIZE_ERR, 'beams3d_main', ierr_mpi)
      CALL beams3d_init_mpi_split(MPI_COMM_BEAMS)
      CALL MPI_GET_VERSION(vmajor,vminor,ierr_mpi)
      CALL MPI_GET_LIBRARY_VERSION(mpi_lib_name,liblen,ierr_mpi)
      ! Now we set some info
      CALL MPI_INFO_CREATE(mpi_info_beams3d, ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_beams3d, "IBM_largeblock_io", "true",    ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_beams3d, "stripping_unit",    "1048576", ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_beams3d, "romio_ds_read",     "disable", ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_beams3d, "romio_ds_write",    "disable", ierr_mpi)
#endif
      RETURN
      END SUBROUTINE beams3d_init_mpi

      SUBROUTINE beams3d_init_mpi_split(comm)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: comm
#if defined(MPI_OPT)
      CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
#endif
      END SUBROUTINE beams3d_init_mpi_split

      SUBROUTINE beams3d_cleanup
      IMPLICIT NONE
      INTEGER :: ier
      ! Clean up
      ier = 0
      CALL beams3d_free(MPI_COMM_SHARMEM)
      CALL wall_free(ier,MPI_COMM_BEAMS)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_main', ierr_mpi)
      ierr_mpi=0
      CALL MPI_INFO_FREE(mpi_info_beams3d, ierr_mpi)
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR, 'beams3d_main', ierr_mpi)
#endif
      IF (lverb) WRITE(6, '(A)') '----- BEAMS3D DONE -----'
      RETURN
      END SUBROUTINE

      SUBROUTINE beams3d_init_hdf5
      IMPLICIT NONE
      INTEGER :: ier
#if defined(LHDF5)
      CALL H5GET_LIBVERSION_F(h5major, h5minor, h5rel, ier)
      h5par = 0
#endif
#if defined(HDF5_PAR)
      h5par = 1
#endif
      RETURN
      END SUBROUTINE beams3d_init_hdf5

      SUBROUTINE beams3d_init_constants
      IMPLICIT NONE
      pi = 4.0 * ATAN(1.0)
      pi2 = 8.0 * ATAN(1.0)
      invpi2 = 1./pi2
      mu0 = (16.0E-7) * ATAN(1.0)
      to3 = REAL(2)/REAL(3)
      lverb = .true.
      lread_input = .true.
      RETURN
      END SUBROUTINE beams3d_init_constants

      SUBROUTINE beams3d_init_commandline
      IMPLICIT NONE
      INTEGER :: numargs, i, ier
      INTEGER, parameter :: arg_len = 256
      CHARACTER*(arg_len) :: arg1
      CHARACTER*(arg_len), allocatable, dimension(:) :: args
      lverb = .false.
      IF (myworkid == master) THEN
        numargs = 0
        i = 0
        arg1 = ''
        limas = .false.
        lverb = .true.
        lvmec = .false.
        lpies = .false.
        lspec = .false.
        leqdsk = .false.
        lcoil = .false.
        lmgrid = .false.
        lvessel = .false.
        lvac = .false.
        lrestart_grid = .false.
        lrestart_particles = .false.
        lhitonly  = .false.
        lplasma_only = .false.
        lraw = .false.
        ldepo = .false.
        lbeam_simple = .false.
        lcollision = .false.
        lw7x = .false.
        lascot = .false.
        lfidasim = .false.
        lfidasim2 = .false.
        lascot4 = .false.
        lbbnbi = .false.
        lascotfl = .false.
        lrandomize = .false.
        lsuzuki = .false.
        lfusion = .false.
        lfusion_alpha = .false.
        id_string = ''
        coil_string = ''
        mgrid_string = ''
        vessel_string = ''
        restart_string = ''
        bbnbi_string = ''
        eqdsk_string = ''

        ! First Handle the input arguments
        CALL GETCARG(1, arg1, numargs)
        ALLOCATE(args(numargs))
        ! Cycle through Arguments
        i = 1
        DO WHILE (i <= numargs)
            call GETCARG(i, args(i), numargs)
            select case (args(i))
            case ("-noverb") ! No Verbose Output
                lverb = .false.
            case ("-vac") ! Vacuum Fields Only
                lvac = .true.
            case ("-ascot","-ascot5")
                lascot = .true.
            case ("-ascot_fl","-ascot5_fl")
                lascot = .true.
                lascotfl = .true.
            case ("-ascot4")
                lascot4 = .true.
            case ("-fidasim")
                lfidasim = .true.
            case ("-fidasim2")
                lfidasim2 = .true.
            case ("-vmec")
                i = i + 1
                lvmec = .true.
                CALL GETCARG(i, id_string, numargs)
            case ("-pies")
                i = i + 1
                lpies = .true.
                CALL GETCARG(i, id_string, numargs)
            case ("-spec")
                i = i + 1
                lspec = .true.
                CALL GETCARG(i, id_string, numargs)
            case ("-hint")
                i = i + 1
                lhint = .true.
                CALL GETCARG(i, id_string, numargs)
            case ("-eqdsk")
                i = i + 1
                leqdsk = .true.
                CALL GETCARG(i, id_string, numargs)
                i = i + 1
                CALL GETCARG(i, eqdsk_string, numargs)
            case ("-mgrid")
                i = i + 1
                lmgrid = .true.
                lcoil = .false.
                CALL GETCARG(i, mgrid_string, numargs)
            case ("-restart")
                i = i + 1
                lrestart_particles = .true.
                CALL GETCARG(i, restart_string, numargs)
            case ("-restart_grid")
                i = i + 1
                lrestart_grid = .true.
                CALL GETCARG(i, restart_string, numargs)
            case ("-coil")
                i = i + 1
                lcoil = .true.
                lmgrid = .false.
                CALL GETCARG(i, coil_string, numargs)
            case ("-vessel")
                i = i + 1
                lvessel = .true.
                CALL GETCARG(i, vessel_string, numargs)
            case ("-beamlet")
                i = i + 1
                lbbnbi = .true.
                CALL GETCARG(i, bbnbi_string, numargs)
            case ("-hitonly","-hit_only")
                lhitonly  = .true.
            case ("-depo")
                ldepo  = .true.
            case ("-raw")
                lraw  = .true.
            case ("-w7x")
                lw7x  = .true.
            case ("-beam_simple")
                lbeam_simple  = .true.
            case ("-collisions")
                lcollision = .true.
            case ("-plasma")
                lplasma_only = .true.
            case ("-rand")
                lrandomize = .true.
            case ("-suzuki")
                lsuzuki = .true.
            case ("-fusion")
                lfusion = .true.
            case ("-fusion_alpha")
                lfusion = .true.
                lfusion_alpha = .true.
            case ("-help", "-h") ! Output Help message
                write(6, *) ' Beam MC Code'
                write(6, *) ' Usage: xbeams3d <options>'
                write(6, *) '    <options>'
                write(6, *) '     -vmec ext:     VMEC input/wout extension'
                write(6, *) '     -hint ext:     HINT input/wout extension'
                write(6, *) '     -eqdsk in gf   BEAMS3D input file and gfile'
                !write(6,*)'     -pies ext:   PIES input extension (must have &INDATA namelist)'
                !write(6,*)'     -spec ext:     SPEC input extension (must have &INDATA namelist)'
                write(6, *) '     -vessel file:  Vessel File (for limiting)'
                write(6, *) '     -mgrid file:   MAKEGRID File (for vacuum)'
                write(6, *) '     -coil file:    Coils. File (for vacuum)'
                write(6, *) '     -restart ext:  BEAMS3D HDF5 extension for starting particles'
                write(6, *) '     -beamlet ext:  Beamlet file for beam geometry'
                write(6, *) '     -beam_simple:  Monoenergetic BEAMS'
                write(6, *) '     -ascot5:       Output data in ASCOT5 gyro-center format'
                write(6, *) '     -ascot5_fl:    Output data in ASCOT5 fieldline format'
                write(6, *) '     -ascot4:       Output data in ASCOT4 format'
                write(6, *) '     -raw:          Treat coil currents as raw (scale factors)'
                write(6, *) '     -vac:          Only vacuum field'
                write(6, *) '     -plasma:       Only plasma field'
                write(6, *) '     -depo:         Only Deposition'
                write(6, *) '     -collisions:   Force collision operator'
                write(6, *) '     -rand:         Randomize particle processor'
                write(6, *) '     -suzuki:       Force Suzuki NBI model'
                write(6, *) '     -fusion:       Fusion Reaction Rates for birth'
                write(6, *) '     -fusion_alpha: Fusion Reaction Rates for birth (alphas only)'
                write(6, *) '     -noverb:       Supress all screen output'
                write(6, *) '     -help:         Output help message'
            end select
            i = i + 1
        END DO
        DEALLOCATE(args)
      END IF
    ! Broadcast variables
#if defined(MPI_OPT)
      CALL MPI_BCAST(id_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(mgrid_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(coil_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(vessel_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(bbnbi_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(restart_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(eqdsk_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lvmec, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lpies, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lspec, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lhint, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(leqdsk, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lmgrid, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lcoil, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lvessel, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lvac, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lfidasim, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lfidasim2, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)  
      CALL MPI_BCAST(lascot, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lascotfl, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lascot4, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lbbnbi, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lrestart_grid, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lrestart_particles, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
      CALL MPI_BCAST(lhitonly,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lraw,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lplasma_only,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(ldepo,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lbeam_simple,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lcollision,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lw7x,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lrandomize,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lsuzuki,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lfusion,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(lfusion_alpha,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
      CALL MPI_BCAST(limas,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'beams3d_main',ierr_mpi)
#endif
      RETURN
      END SUBROUTINE beams3d_init_commandline

      SUBROUTINE beams3d_output_header
      IMPLICIT NONE
      INTEGER :: nshar
#if defined(MPI_OPT)
      CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ierr_mpi) ! MPI
#endif
      IF (lverb) THEN
        WRITE(6, '(a,f5.2)') 'BEAMS3D Version ', BEAMS3D_VERSION
#if defined(LHDF5)
        IF (h5par > 0) THEN
           WRITE(6,'(A)')      '-----  HDF5 (Parallel) Parameters  -----'
        ELSE
           WRITE(6,'(A)')      '-----  HDF5 Parameters  -----'
        ENDIF
        WRITE(6,'(A,I2,2(A,I2.2))')  '   HDF5_version:  ', h5major,'.',h5minor,' release: ',h5rel
#endif
        WRITE(6,'(A)')      '-----  MPI Parameters  -----'
        WRITE(6,'(A,I2,A,I2.2)')  '   MPI_version:  ', vmajor,'.',vminor
        WRITE(6,'(A,A)')  '   ', TRIM(mpi_lib_name(1:liblen))
        WRITE(6,'(A,I8)')  '   Nproc_total:  ', nprocs_beams
        WRITE(6,'(A,3X,I5)')  '   Nproc_shared: ', nshar
      END IF
      RETURN
      END SUBROUTINE beams3d_output_header

END MODULE BEAMS3D_INTERFACE_MOD