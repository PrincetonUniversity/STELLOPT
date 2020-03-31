!-----------------------------------------------------------------------
!     Program:       TORLINES
!     Authors:       S. Lazerson
!     Date:          04/30/2013
!     Description:   
!     References:
!-----------------------------------------------------------------------
      PROGRAM TORLINES
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE torlines_runtime
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------  
      implicit none
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'
!DEC$ ENDIF  
      integer                                      :: numargs,i
      integer, parameter                           :: arg_len =256
      character*(arg_len)                          :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      myid = master
      numprocs = 1
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = MPI_SUCCESS
      CALL MPI_INIT( ierr_mpi )                                         ! MPI
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SPLIT( MPI_COMM_WORLD,0,myid,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_FIELDLINES, myid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, numprocs, ierr_mpi )          ! MPI
      CALL MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr_mpi)
!DEC$ ENDIF
      !OPEN(6,CARRIAGECONTROL='fortran')
      pi2 = 8 * ATAN(1.0)

      IF (myid == master) THEN
         numargs=0
         i=0
         arg1=''
         lverb=.true.
         lrestart=.false.
         lvmec    = .false.
         lpies    = .false.
         lspec    = .false.
         lcoil    = .false.
         lmgrid   = .false.
         lvac     = .false.
         lvessel  = .false.
         lraw     = .false.
         lemc3    = .false.
         lauto    = .false.
         ! First Handle the input arguments
         id_string = ''
         coil_string = ''
         mgrid_string = ''
         CALL GETCARG(1, arg1, numargs)
         ALLOCATE(args(numargs))
         args(1)=arg1
         ! Cycle through Arguments
         i = 1
         DO WHILE (i <= numargs)
            call GETCARG(i,args(i),numargs)
            select case (args(i))
               case ("-noverb")  ! No Verbose Output
                   lverb=.false.
               case ("-vac")  
                   lvac=.true.
               case ("-vmec")
                   i = i + 1
                   lvmec = .true.
                   lpies = .false.
                   lspec = .false.
                   CALL GETCARG(i,id_string,numargs)
               case ("-raw")
                   lraw = .true.
               case ("-emc3")
                   lemc3 = .true.
               case ("-auto")
                   lauto = .true.
               case ("-mgrid")
                   i = i + 1
                   lmgrid = .true.
                   lcoil  = .false.
                   CALL GETCARG(i,mgrid_string,numargs)
               case ("-coil","-coils")
                   i = i + 1
                   lcoil  = .true.
                   lmgrid = .false.
                   CALL GETCARG(i,coil_string,numargs)
               case ("-vessel")
                   i = i + 1
                   lvessel = .true.
                   CALL GETCARG(i,vessel_string,numargs)
               case ("-help","-h") ! Output Help message
                  write(6,*)' Usage: xtorlines input_suffix <options>'
                  write(6,*)'    input_file:  Output of VMEC or PIES'
                  write(6,*)'    <options>'
                  write(6,*)'     -emc3:   Generate EMC3 Grid'
                  write(6,*)'     -lauto:  Auto calculation of edge'
                  write(6,*)'     -raw:    Treat currents as raw'
                  write(6,*)'     -coil coil_file: Use coil for vacuum fields'
                  write(6,*)'     -mgrid mgrid_file:  Use mgrid for vacuum fields'
                  write(6,*)'     -vessel ves_file: Use vessel file'
                  write(6,*)'     -help:   Output help message'
                  stop
            end select
            i = i + 1
         enddo
         DEALLOCATE(args)
         WRITE(6,'(a,f4.2)') 'TORLINES Version ',TORLINES_VERSION
      ELSE IF (myid /= master) THEN
         lverb=.false.   ! Shutup the slaves
      END IF
      ! Broadcast variables

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(id_string,256,MPI_CHARACTER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(mgrid_string,256,MPI_CHARACTER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(coil_string,256,MPI_CHARACTER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(vessel_string,256,MPI_CHARACTER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lvmec,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lpies,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lspec,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lmgrid,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lcoil,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lvessel,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lraw,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lemc3,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lvac,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
      CALL MPI_BCAST(lauto,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'torlines_main',ierr_mpi)
!DEC$ ENDIF
      ! Initialize the Calculation
      CALL torlines_init
      IF (lverb) WRITE(6,*) ''
      IF (lverb) WRITE(6,'(a)')'---------- EXECUTION ----------'
      IF (lemc3) THEN
         CALL torlines_emc3
      ELSE IF (lauto) THEN
         CALL follow
      ELSE
         CALL torlines_follow
         CALL torlines_write('FIELDLINES')
      END IF
      
      IF (lverb) WRITE(6,*) ''
      IF (lverb) WRITE(6,'(a)')'---------- TORLINES DONE ----------'
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'torlines_main',ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM TORLINES
