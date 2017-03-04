      PROGRAM stellarator_optimizer
!
!     Driver routine (front end)
!
!     Invoke optimizer as follows:
!
!     xstellopt input_file_name  (input_file_name == name of vmec input_file)
!
      USE optim, ONLY: lone_step, lscale_only, lrestart
      USE mpi_params                                                            ! MPI

      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: arg_len = 200
      INTEGER :: numargs, iseq, ikey
      CHARACTER(len=arg_len) :: arg1_input
      CHARACTER(len=arg_len), ALLOCATABLE, DIMENSION(:) :: args
!-----------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
!
!     Call MPI Initialization routines:                                          ! MPI

      CALL MPI_INIT( ierr_mpi )                                                  ! MPI
      IF (ierr_mpi .ne. 0) STOP 'MPI_INIT error in stellopt'
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr_mpi )                       ! MPI
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in stellopt'
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr_mpi )                   ! MPI
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in stellopt'
!DEC$ ENDIF
      lscale_only = .false.
      lone_step = .false.
      lrestart = .false.
      
!
!     read input file name from command line
!      
!DEC$ IF DEFINED (MPI_OPT) .AND. DEFINED (LINUX)
      IF( myid .eq. master) THEN
!DEC$ ENDIF
         CALL getcarg(1, arg1_input, numargs)
         ALLOCATE (args(numargs))
         args(1) = arg1_input
         DO iseq = 2, numargs
            CALL getcarg(iseq, args(iseq), numargs)
         ENDDO

         IF (numargs .lt. 1 .or. args(1) .eq. '-h' .or. 
     1      args(1) .eq. '/h') THEN
            PRINT *,
     1      ' ENTER INPUT FILE NAME (OR LIST OF NAMES) OR INPUT-FILE',
     2      ' SUFFIX(ES) ON COMMAND LINE'
            PRINT *,' For example: '
            PRINT *,' xstellopt input.tftr (or tftr or ../input.tftr)'
            PRINT *,' Optional control parameters may be entered:'
            PRINT *,' xstellopt {(-o)ne_step} {(-s)cale_only}',
     1      '{(-r)estart} [list of files]'
            PRINT *,' where:'
            PRINT *,' (-o)ne_step:   runs input file(s) one step only'
            PRINT *,'          (effect is same as setting NITER_OPT=1'
            PRINT *,'               in the input file)'           
            PRINT *,
     1      ' (-s)cale_only: scales entries in input files according'
            PRINT *,'        to r00_scale, b00_scale values, writes'
            PRINT *,'        new input file, but does NOT execute'
            PRINT *,'        any optimization steps'
            PRINT *,' (-r)esume: resumes genetic algorithm or'
            PRINT *,' diff. evolution from previous population'
            STOP 'Invalid command line'
         ENDIF
         print *,'====================================================='
     1          ,'================='
         print *,'===================== S T E L L O P T  ( v2.0 )  ===='
     1          ,'================='
         print *,'====================================================='
     1          ,'================='

!DEC$ IF DEFINED (MPI_OPT) .AND. DEFINED (LINUX)
      ENDIF

C      PRINT *,"myid=",myid, " before first Bcast" ; CALL flush(6)
      CALL MPI_BCAST(numargs,1,MPI_INTEGER, master, MPI_COMM_WORLD, 
     1               ierr_mpi)

      IF( myid .ne. master) ALLOCATE (args(numargs))

C      PRINT *,myid, "before Bcast loop" ; CALL flush(6)
      DO iseq = 1, numargs
         CALL MPI_BCAST(args(iseq), arg_len, MPI_CHARACTER,master,
     1                  MPI_COMM_WORLD, ierr_mpi)
      ENDDO
!DEC$ ENDIF

C      DO iseq = 2, numargs
C         PRINT *, myid, iseq, args(iseq)
C      ENDDO
C     CALL flush(6)

      DO iseq = 1, numargs
         IF (args(iseq)(1:2).eq."-s" .or. args(iseq)(1:2).eq."-S"
     1   .or. args(iseq)(1:2).eq."/s" .or. args(iseq)(1:2).eq."/S")
     2   lscale_only = .true.          
         IF (args(iseq)(1:2).eq."-o" .or. args(iseq)(1:2).eq."-O"
     1   .or. args(iseq)(1:2).eq."/o" .or. args(iseq)(1:2).eq."/O")
     2   lone_step = .true.          
         IF (args(iseq)(1:2).eq."-r" .or. args(iseq)(1:2).eq."-R"
     1   .or. args(iseq)(1:2).eq."/r" .or. args(iseq)(1:2).eq."/R")
     2   lrestart = .true.
      END DO

!
!     read each input file name from command line & optimize it
!      
      DO iseq = 1, numargs
         IF (args(iseq)(1:1).ne.'-' .and. args(iseq)(1:1).ne.'/')
     1      CALL optimize(args(iseq))
      END DO

      DEALLOCATE(args)

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_FINALIZE(ierr_mpi)      !Close out MPI                             !MPI
      IF (ierr_mpi .ne. 0) STOP 'MPI_FINALIZE error in stellopt'
!DEC$ ENDIF
      END PROGRAM stellarator_optimizer
