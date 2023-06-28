c-----------------------------------------------------------------------
c     Module:        mpi_params
c     Authors:       Various
c     Date:          02/18/2013
c     Description:   This module constains the communicators and
c                    deffinitions for various parallel codes.
c                    Please use MPI_COMM_MYWORLD as a proxy for
c                    MPI_COMM_WORLD.  This way if a parallel code
c                    needs to call your code it can do so.
c                    STELLOPT will initialize MPI_COMM_WORLD then
c                    split off a subset of threads into MPI_COMM_MYWORLD
c                    for calling parallel codes, MPI_COMM_STEL are
c                    the processes running STELLOPT.  This is important
c                    because the parallel codes must be handled
c                    differently.
      MODULE mpi_params
      INTEGER, PARAMETER :: master=0
      INTEGER, PARAMETER :: WORKER_SPLIT_KEY = 3
      INTEGER :: myid=master, numprocs, ierr_mpi, numworkers
      INTEGER :: myworkid=master, my_master=master
      ! Communicators
      INTEGER :: MPI_COMM_WORKERS=-1, worker_id=-1       !communicator for worker processors only
      INTEGER :: MPI_COMM_WORKERS_OK=-1, worker_id_ok=-1 !communicator subgroup, vmec ran ok
      INTEGER :: MPI_COMM_SHARMEM = 718, myid_sharmem=-1 !communicator for shared memory
      INTEGER :: MPI_COMM_STEL = 327                     !communicator which is a copy of MPI_COMM_WORLD (user must set this up)
      INTEGER :: MPI_COMM_MYWORLD = 411                  !communicator 
      INTEGER :: MPI_COMM_FIELDLINES = 328               !communicator for FIELDLINES code
      INTEGER :: MPI_COMM_TORLINES = 329                 !communicator for TORLINES code
      INTEGER :: MPI_COMM_BEAMS = 330                    !communicator for BEAMS3D code
      INTEGER :: MPI_COMM_BOOZER = 331                   !communicator for BOOZ_XFORM code
      INTEGER :: MPI_COMM_DIAGNO = 332                   !communicator for DIAGNO code
      INTEGER :: MPI_COMM_WALLACC = 333                  !communicator for WALLACC code
      INTEGER :: MPI_COMM_THRIFT = 334                   !communicator for THRIFT code
      INTEGER :: MPI_COMM_PARVMEC = 101                  !communicator for PARVMEC code

      CONTAINS
      
      SUBROUTINE mpi_stel_abort(error)
#if defined(MPI_OPT)
      USE MPI
#endif
      IMPLICIT NONE        
      INTEGER, INTENT(in)                 :: error
      INTEGER                             :: length, temp
      CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
#if defined(MPI_OPT)
      CALL MPI_ERROR_STRING(error,message,length,temp)
      WRITE(6,*) '!!!!!!!!!!!!MPI_ERROR DETECTED!!!!!!!!!!!!!!'
      WRITE(6,*) '  MESSAGE: ',message(1:length)
      WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      CALL FLUSH(6) 
#else
      WRITE(6,*) '!!!!!!!!!!!!MPI_ERROR DETECTED!!!!!!!!!!!!!!'
      WRITE(6,*) '  MPI_STEL_ABORT CALLED BUT NO MPI'
      WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#endif
      !CALL MPI_ABORT(MPI_COMM_STEL,1,temp)
      END SUBROUTINE mpi_stel_abort

      SUBROUTINE MPI_CALC_MYRANGE(comm,n1,n2,mystart,myend)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in) :: n1, n2
      INTEGER, INTENT(out) :: mystart, myend
      INTEGER :: delta, local_size, local_rank, istat, maxend, k, i
      mystart = n1; myend = n2
#if defined(MPI_OPT)
      CALL MPI_COMM_SIZE( comm, local_size, istat)
      CALL MPI_COMM_RANK( comm, local_rank, istat )
      delta = CEILING(DBLE(n2-n1+1)/DBLE(local_size))
      mystart = n1 + local_rank*delta
      myend   = mystart + delta - 1
      maxend = local_size*delta
      IF (maxend>n2) THEN
         k = maxend-n2
         DO i = (local_size-k), local_size-1
            IF (local_rank > i) THEN
                  mystart = mystart - 1
                  myend   = myend - 1
            ELSEIF (local_rank==i) THEN
                  myend = myend - 1
            END IF
         END DO
      END IF
#endif
      RETURN
      END SUBROUTINE MPI_CALC_MYRANGE

      END MODULE mpi_params
