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
      USE mpi_sharmem
      USE mpi_params
      USE mpi_inc
#if defined(LHDF5)
    USE hdf5
#endif
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------  
      implicit none
      integer                                      :: numargs,i,ier, vmajor, vminor, liblen, nshar
      integer                                      :: h5major, h5minor, h5rel, h5par
      integer, parameter                           :: arg_len =256
      character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_lib_name
      character*(arg_len)                          :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      myworkid = master
      ierr_mpi = MPI_SUCCESS
#if defined(MPI_OPT)
      CALL MPI_INIT(ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_INIT_ERR, 'fieldlines_main', ierr_mpi)
      CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_FIELDLINES, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'fieldlines_main', ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_FIELDLINES, myworkid, ierr_mpi )              ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'fieldlines_main', ierr_mpi)
      CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, nprocs_torlines, ierr_mpi )          ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_SIZE_ERR, 'fieldlines_main', ierr_mpi)
      CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_FIELDLINES, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ierr_mpi) ! MPI
      CALL MPI_GET_VERSION(vmajor,vminor,ier)
      CALL MPI_GET_LIBRARY_VERSION(mpi_lib_name,liblen,ier)
      CALL MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr_mpi)
#endif

#if defined(LHDF5)
    CALL H5GET_LIBVERSION_F(h5major, h5minor, h5rel, ier)
    h5par = 0
#endif

#if defined(HDF5_PAR)
    h5par = 1
#endif

      !OPEN(6,CARRIAGECONTROL='fortran')
      pi2 = 8 * ATAN(1.0)
      lverb = .false.

      IF (myworkid == master) THEN
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
         WRITE(6,'(A,I8)')  '   Nproc_total:  ', nprocs_torlines
         WRITE(6,'(A,3X,I5)')  '   Nproc_shared: ', nshar
      END IF
      ! Broadcast variables

#if defined(MPI_OPT)
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
#endif
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
#if defined(MPI_OPT)
      CALL MPI_COMM_FREE(MPI_COMM_SHARMEM,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'torlines_main',ierr_mpi)
#endif
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM TORLINES
