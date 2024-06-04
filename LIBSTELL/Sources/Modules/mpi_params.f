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

      !> Distribute the workload of operating on over n1:n2 (inclusive)
      !> over the compute ranks available in the given communicator
      !> and return the local ranges to be worked on in mystart and myend.
      !
      !> This routine must __always__ run,
      !> hence no `STOP` statements or similar are allowed here.
      !> If more ranks than work items are available in the communicator,
      !> this routine returns `myend` > `mystart`,
      !> which implies that loops like `DO i = mystart, myend` are simply skipped
      !> in ranks that do not get a share of the workload.
      SUBROUTINE MPI_CALC_MYRANGE(comm,n1,n2,mystart,myend)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: comm    !< communicator to distribute work over
      INTEGER, INTENT(in)    :: n1      !< lower bound of range to work on
      INTEGER, INTENT(in)    :: n2      !< upper bound of range to work on (inclusive)
      INTEGER, INTENT(out)   :: mystart !< lower bound of chunk this rank should work on
      INTEGER, INTENT(out)   :: myend   !< upper bound of chunk this rank should work on (inclusive)

      INTEGER :: local_size, local_rank, istat
      INTEGER :: total_work, work_per_rank, work_remainder

      ! Default if not using MPI: just work on full range
      mystart = n1
      myend   = n2

#if defined(MPI_OPT)

      ! `local_size` is the number of available ranks.
      ! We assume it is always > 0, i.e., 1, 2, 3, ...
      CALL MPI_COMM_SIZE(comm, local_size, istat)

      ! `local_rank` is the ID of the rank to compute `mystart` and `myend` for.
      ! We assume it is always >= 0, i.e., 0, 1, 2, ...
      CALL MPI_COMM_RANK(comm, local_rank, istat)

      ! Total number of items to work on.
      ! NOTE: n2 is the upper range bound, inclusive!
      total_work = n2 - n1 + 1

      ! size of chunks that are present in all ranks
      work_per_rank = total_work / local_size

      ! number of work items that remain after distributing
      ! equal chunks of work to all ranks
      work_remainder = MODULO(total_work, local_size)

      ! ranges corresponding to working on evenly distributed chunks
      ! `myend` is inclusive, i.e., the indices to work on are
      ! { mystart, mystart+1, ..., myend-1, myend }.
      ! Thus, one can use code like `DO i = mystart, myend`.
      mystart = n1 +  local_rank      * work_per_rank
      myend   = n1 + (local_rank + 1) * work_per_rank - 1

      IF (local_rank .lt. work_remainder) THEN
         ! The first `work_remainder` ranks get one additional item to work on.
         ! This takes care of the additional `work_remainder` items
         ! that need to be worked on, on top of the evenly distributed chunks.
         mystart = mystart + local_rank
         myend   = myend   + local_rank + 1
      ELSE
         ! All following ranks after the first `work_remainder` ones
         ! get their ranges just shifted by a constant offset,
         ! since they don't do any additional work.
         mystart = mystart + work_remainder
         myend   = myend   + work_remainder
      END IF

#endif
      RETURN
      END SUBROUTINE MPI_CALC_MYRANGE

      END MODULE mpi_params
