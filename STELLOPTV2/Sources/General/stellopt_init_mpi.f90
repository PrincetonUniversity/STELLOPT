!-----------------------------------------------------------------------
!     Subroutine:    stellopt_init_mpi
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/04/2019
!     Description:   This subroutine initializes the vmec MPI communicators.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_init_mpi
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod, ONLY: noptimizers
      USE mpi_params
      USE mpi_inc
!----------------------------------------------------------------------
!     Local Variables
!        nprocs_total :: total number of processors         
!----------------------------------------------------------------------
      INTEGER :: nprocs_total, vmajor, vminor, color, key, comm_shar, num_shar, numjoin
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
#if defined(MPI_OPT)
      ierr_mpi = 0
      ! Get total number of processors
      CALL MPI_COMM_SIZE( MPI_COMM_STEL, nprocs_total, ierr_mpi )
      CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi)
      CALL MPI_BCAST( noptimizers, 1, MPI_INTEGER, master, MPI_COMM_STEL, ierr_mpi)

      ! Create a local communicator
      IF (noptimizers <=0 ) THEN ! Every process an optimizer
         noptimizers = numprocs + 1
         color = MOD(myid,noptimizers)
         key = myid
         CALL MPI_COMM_SPLIT( MPI_COMM_STEL, color, key, MPI_COMM_MYWORLD, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi0',ierr_mpi)
      ELSE ! make is shared memory
         CALL MPI_COMM_SPLIT_TYPE( MPI_COMM_STEL, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_MYWORLD, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi1',ierr_mpi)
      END IF
      CALL MPI_COMM_RANK( MPI_COMM_MYWORLD, myworkid, ierr_mpi)
      CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, num_shar, ierr_mpi)

      ! Free Duplicate 
      CALL MPI_COMM_FREE( MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi2',ierr_mpi)

      ! Create MPI_COMM_STEL from masters
      color = MPI_UNDEFINED
      myid = -1
      IF (myworkid == master) color = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, color, myid, MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi3',ierr_mpi)
      IF (myworkid == master)THEN
         CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )              ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi4',ierr_mpi)
         CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )          ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi5',ierr_mpi)
      END IF

      ! Broadcast information to everyone in worker groups
      CALL MPI_BCAST(numprocs, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi6',ierr_mpi)

      ! Handle mismatch of numprocs to nopt
      IF (noptimizers > numprocs) THEN ! We Split
         ! We can't yet split the shared memory communicator
         noptimizers = numprocs 
      END IF

      IF (myid == master) THEN
         CALL MPI_GET_VERSION(vmajor,vminor,ier)
         WRITE(6,*) '-----  MPI Params.   -----'
         WRITE(6,*) '   MPI Version: ',vmajor,'.',vminor
         WRITE(6,*) '   Number of Processors: ',nprocs_total
         WRITE(6,*) '   Size of Shared Mem. Comm: ',num_shar
         WRITE(6,*) '   Number of Optimizer Threads:    ',numprocs
         CALL FLUSH(6)
      END IF
#endif
     
      RETURN
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_init_mpi

