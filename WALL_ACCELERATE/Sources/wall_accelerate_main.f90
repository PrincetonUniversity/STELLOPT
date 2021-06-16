PROGRAM WALL_ACCELERATE
   !-----------------------------------------------------------------------
   !     Libraries
   !-----------------------------------------------------------------------
         USE safe_open_mod
         USE wall_mod

         USE mpi_sharmem
         USE mpi_params
         USE mpi_inc
   !-----------------------------------------------------------------------
   !     Local Variables
   !          numargs      Number of input arguments
   !          i            Index
   !          arg_len      Length of input strings
   !          arg1         Input file
   !          args         Input arguments
   !-----------------------------------------------------------------------
         IMPLICIT NONE
         integer                                      :: numargs,i,ier, vmajor, vminor, liblen, nshar, nprocs_wallacc
         integer, parameter                           :: arg_len =256
         character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_lib_name
         character*(arg_len)                          :: arg1
         character*(arg_len),allocatable,dimension(:) :: args
         CHARACTER(arg_len)                           :: vessel_string
         INTEGER                                      :: nface_des = 0

   !-----------------------------------------------------------------------
   !     Error codes (same as FIELDLINES)
   !-----------------------------------------------------------------------
         INTEGER, PARAMETER ::  MPI_ERR            = 80
         INTEGER, PARAMETER ::  MPI_INIT_ERR       = 801
         INTEGER, PARAMETER ::  MPI_RANK_ERR       = 802
         INTEGER, PARAMETER ::  MPI_SIZE_ERR       = 803
         INTEGER, PARAMETER ::  MPI_DUP_ERR        = 804
         INTEGER, PARAMETER ::  MPI_BARRIER_ERR    = 81
         INTEGER, PARAMETER ::  MPI_SEND_ERR       = 821
         INTEGER, PARAMETER ::  MPI_RECV_ERR       = 822
         INTEGER, PARAMETER ::  MPI_BCAST_ERR      = 83
         INTEGER, PARAMETER ::  MPI_FINE_ERR       = 89

   !-----------------------------------------------------------------------
   !     Parameters
   !-----------------------------------------------------------------------
         REAL, PARAMETER :: WALL_ACC_VERSION = 1.00
   !-----------------------------------------------------------------------
   !     Begin Program
   !-----------------------------------------------------------------------
         myworkid = master
         ierr_mpi = MPI_SUCCESS
#if defined(MPI_OPT)
         CALL MPI_INIT(ierr_mpi) ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_INIT_ERR, 'wall_acc_main', ierr_mpi)
         CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_WALLACC, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'wall_acc_main', ierr_mpi)
         CALL MPI_COMM_RANK( MPI_COMM_WALLACC, myworkid, ierr_mpi )              ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'wall_acc_main', ierr_mpi)
         CALL MPI_COMM_SIZE( MPI_COMM_WALLACC, nprocs_wallacc, ierr_mpi )          ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_SIZE_ERR, 'wall_acc_main', ierr_mpi)
         CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WALLACC, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
         CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
         CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ierr_mpi) ! MPI
         CALL MPI_GET_VERSION(vmajor,vminor,ier)
         CALL MPI_GET_LIBRARY_VERSION(mpi_lib_name,liblen,ier)
         CALL MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr_mpi)
#endif

         IF (myworkid == master) THEN
            numargs=0
            i=0
            arg1=''
            
            ! First Handle the input arguments
            CALL GETCARG(1, arg1, numargs)
            ALLOCATE(args(numargs))
            ! Cycle through Arguments
            i=1
            DO WHILE (i <= numargs)
               CALL GETCARG(i,args(i),numargs)
               SELECT CASE (args(i))
                  CASE ("-vessel")
                        i = i + 1
                        CALL GETCARG(i,vessel_string,numargs)
                  CASE ("-nface")
                        i = i + 1
                        CALL GETCARG(i,args(i),numargs)
                        READ(args(i),*,IOSTAT=ier) nface_des
                  CASE ("-help","-h") ! Output Help message
                     WRITE(6,*) 'Program to generate acceleration wall meshes. Input original mesh with -vessel'
                     WRITE(6,*)'    <options>'
                     WRITE(6,*)'     -vessel file:  Vessel File to speed up'
                     WRITE(6,*)'     -nface: Number of faces per block desired'
                     WRITE(6,*)'     -help          Output help message'
#if defined(MPI_OPT)
                     CALL MPI_FINALIZE(ierr_mpi)   
                     IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'wall_acc_main',ierr_mpi)
#endif
               end SELECT
               i = i + 1
            END DO
            DEALLOCATE(args)

            WRITE(6,'(a,f5.2)') 'WALL_ACC Version ', WALL_ACC_VERSION
            WRITE(6,'(A)')      '-----  MPI Parameters  -----'
            WRITE(6,'(A,I2,A,I2.2)')  '   MPI_version:  ', vmajor,'.',vminor
            WRITE(6,'(A,A)')  '   ', TRIM(mpi_lib_name(1:liblen))
            WRITE(6,'(A,I8)')  '   Nproc_total:  ', nprocs_wallacc
            WRITE(6,'(A,3X,I5)')  '   Nproc_shared: ', nshar


            vessel_string = TRIM(vessel_string)
            vessel_string = ADJUSTL(vessel_string)
         END IF

#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_WALLACC,ierr_mpi)
         CALL MPI_BCAST(vessel_string,256,MPI_CHARACTER, master, MPI_COMM_WALLACC,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'wall_acc_main',ierr_mpi)
         CALL MPI_BCAST(nface_des,1,MPI_INTEGER, master, MPI_COMM_WALLACC,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'wall_acc_main',ierr_mpi)
#endif
   
         IF (nface_des .NE. 0)  CALL SET_NFACE(nface_des)
   
#if defined(MPI_OPT)
         CALL wall_load_txt(vessel_string, ier, MPI_COMM_WALLACC)  
         CALL MPI_BARRIER(MPI_COMM_WALLACC,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'wall_acc_main',ierr_mpi)
#else
         CALL wall_load_txt(vessel_string, ier) 
#endif
            
         ! Clean up
#if defined(MPI_OPT)
         CALL wall_free(ier, MPI_COMM_WALLACC)
#else
         CALL wall_free(ier)
#endif
   
#if defined(MPI_OPT)
         CALL MPI_FINALIZE(ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'wall_acc_main',ierr_mpi)
#endif

         IF (myworkid == master) WRITE(6,'(A)')'----- WALL_ACCELERATE DONE -----'
         
   !-----------------------------------------------------------------------
   !     End Program
   !-----------------------------------------------------------------------
         END PROGRAM WALL_ACCELERATE
      