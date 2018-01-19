!-----------------------------------------------------------------------
!     Program:       BEAMS3D
!     Authors:       M. McMillan S. Lazerson
!     Date:          06/20/2012
!     Description:   The BEAMS3D code performs Monte-Carlo particle
!                    simulations on an R-phi-Z cylindrical grid.
!     References:
!-----------------------------------------------------------------------
PROGRAM BEAMS3D
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE beams3d_runtime
    USE wall_mod, ONLY: wall_free
    USE mpi_params ! MPI
    !-----------------------------------------------------------------------
    !     Local Variables
    !          numargs      Number of input arguments
    !          i            Index
    !          arg_len      Length of input strings
    !          arg1         Input file
    !          args         Input arguments
    !-----------------------------------------------------------------------
    IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
    INCLUDE 'mpif.h' ! MPI
!DEC$ ENDIF
    integer :: numargs, i, ier
    integer, parameter :: arg_len = 256
    character*(arg_len) :: arg1
    character*(arg_len), allocatable, dimension(:) :: args
    !-----------------------------------------------------------------------
    !     Begin Program
    !-----------------------------------------------------------------------
    myworkid = master
!DEC$ IF DEFINED (MPI_OPT)
    CALL MPI_INIT(ierr_mpi) ! MPI
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_INIT_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_COMM_RANK(MPI_COMM_BEAMS, myworkid, ierr_mpi) ! MPI
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_COMM_SIZE(MPI_COMM_BEAMS, nprocs_beams, ierr_mpi) ! MPI
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_SIZE_ERR, 'beams3d_main', ierr_mpi)
!DEC$ ENDIF
    pi = 4.0 * ATAN(1.0)
    pi2 = 8.0 * ATAN(1.0)
    mu0 = (16.0E-7) * ATAN(1.0)
    to3 = REAL(2)/REAL(3)
    lverb = .true.
    lread_input = .true.
    IF (myworkid == master) THEN
        !OPEN(6, RECL = 2**24)
        numargs = 0
        i = 0
        arg1 = ''
        lverb = .true.
        lvmec = .false.
        lpies = .false.
        lspec = .false.
        lcoil = .false.
        lmgrid = .false.
        lmu = .false.
        lvessel = .false.
        lvac = .false.
        lrestart = .false.
        lflux = .false.
        lhitonly  = .false.
        lplasma_only = .false.
        lraw = .false.
        ldepo = .false.
        lbeam_simple = .false.
        lcollision = .false.
        lw7x = .false.
        id_string = ''
        coil_string = ''
        mgrid_string = ''
        vessel_string = ''
        restart_string = ''

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
            case ("-plasma")
                lplasma_only = .true.
            case ("-vmec")
                i = i + 1
                lvmec = .true.
                lpies = .false.
                lspec = .false.
                CALL GETCARG(i, id_string, numargs)
            case ("-pies")
                i = i + 1
                lpies = .true.
                lvmec = .false.
                lspec = .false.
                CALL GETCARG(i, id_string, numargs)
            case ("-flux","-booz")
                lflux = .true.
            case ("-spec")
                i = i + 1
                lspec = .true.
                lpies = .false.
                lvmec = .false.
                CALL GETCARG(i, id_string, numargs)
            case ("-mgrid")
                i = i + 1
                lmgrid = .true.
                lcoil = .false.
                CALL GETCARG(i, mgrid_string, numargs)
            case ("-restart")
                i = i + 1
                lrestart = .true.
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
            case ("-hitonly","-hit_only")
                lhitonly  = .true.
            case ("-depo")
                ldepo  = .true.
                lplasma_only = .true.
            case ("-raw")
                lraw  = .true.
            case ("-w7x")
                lw7x  = .true.
            case ("-beam_simple")
                lbeam_simple  = .true.
            case ("-collisions")
                lcollision = .true.
            case ("-help", "-h") ! Output Help message
                write(6, *) ' Beam MC Code'
                write(6, *) ' Usage: xbeams3d <options>'
                write(6, *) '    <options>'
                write(6, *) '     -vmec ext:     VMEC input/wout extension'
                write(6, *) '     -booz file:    BOOZ_XFORM output file for mapping'
                !write(6,*)'     -pies ext:   PIES input extension (must have &INDATA namelist)'
                !write(6,*)'     -spec ext:     SPEC input extension (must have &INDATA namelist)'
                write(6, *) '     -vessel file:  Vessel File (for limiting)'
                write(6, *) '     -mgrid file:   MAKEGRID File (for vacuum)'
                write(6, *) '     -coil file:    Coils. File (for vacuum)'
                write(6, *) '     -beam_simple:  Monoenergetic BEAMS'
                write(6, *) '     -w7x:          W7-X beam model'
                !write(6,*)'     -restart ext:  FIELDLINES HDF5 extension.'
                write(6, *) '     -raw:          Treat coil currents as raw (scale factors)'
                write(6, *) '     -vac:          Only vacuum field'
                write(6, *) '     -plasma:       Only plasma field'
                write(6, *) '     -depo:         Only Deposition'
                write(6, *) '     -flux:         Calculate Flux Coordinates'
                write(6, *) '     -collisions:   Force collision operator'
                write(6, *) '     -noverb:       Supress all screen output'
                write(6, *) '     -help:         Output help message'
!DEC$ IF DEFINED (MPI_OPT)
                CALL MPI_FINALIZE(ierr_mpi)
                IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR, 'beams3d_main', ierr_mpi)
!DEC$ ENDIF
            end select
            i = i + 1
        END DO
        DEALLOCATE(args)
        WRITE(6, '(a,f5.2)') 'BEAMS3D Version ', BEAMS3D_VERSION
    ELSE
        lverb = .false. ! Shutup the workers
    END IF
    ! Broadcast variables
!DEC$ IF DEFINED (MPI_OPT)
    CALL MPI_BCAST(id_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(mgrid_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(coil_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(vessel_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(restart_string, 256, MPI_CHARACTER, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lvmec, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lpies, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lspec, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lmgrid, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lcoil, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lvessel, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lvac, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lmu, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lrestart, 1, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'beams3d_main', ierr_mpi)
    CALL MPI_BCAST(lhitonly,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(lraw,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(lflux,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(lplasma_only,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(ldepo,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(lbeam_simple,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(lcollision,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
    CALL MPI_BCAST(lw7x,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
!DEC$ ENDIF

    ! Initialize the Calculation
    CALL beams3d_init
    ! Follow Fieldlines
    CALL beams3d_follow
    CALL beams3d_write('TRAJECTORY_PARTIAL')
    !IF (myworkid == master) CALL wall_free(ier)
    ! Write some stuff
    CALL beams3d_diagnostics
    ! Clean up
    CALL beams3d_free
!DEC$ IF DEFINED (MPI_OPT)
    ierr_mpi=0
    CALL MPI_FINALIZE(ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR, 'beams3d_main', ierr_mpi)
!DEC$ ENDIF
    IF (lverb) WRITE(6, '(A)') '----- BEAMS3D DONE -----'

    !-----------------------------------------------------------------------
    !     End Program
    !-----------------------------------------------------------------------
END PROGRAM BEAMS3D
