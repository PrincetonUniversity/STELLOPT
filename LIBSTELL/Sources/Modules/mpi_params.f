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
      INTEGER :: MPI_COMM_WORKERS=-1, worker_id=-1                 !communicator for worker processors only
      INTEGER :: MPI_COMM_WORKERS_OK=-1, worker_id_ok=-1           !communicator subgroup, vmec ran ok
      INTEGER :: MPI_COMM_STEL = 327                               !communicator which is a copy of MPI_COMM_WORLD (user must set this up)
      INTEGER :: MPI_COMM_MYWORLD = 411                            !communicator 
      INTEGER :: MPI_COMM_FIELDLINES = 328                         !communicator for FIELDLINES code
      INTEGER :: MPI_COMM_TORLINES = 329                         !communicator for FIELDLINES code
      INTEGER :: MPI_COMM_BEAMS = 330                            !communicator for BEAMS3D code
      INTEGER :: MPI_COMM_BOOZER = 331                           !communicator for BOOZ_XFORM code
      INTEGER :: MPI_COMM_DIAGNO = 332                           !communicator for BOOZ_XFORM code
      INTEGER :: MPI_COMM_PARVMEC = 101                           !communicator for PARVMEC code

!DEC$ IF DEFINED (MPI_OPT)
      CONTAINS
      
      SUBROUTINE mpi_stel_abort(error)
      USE MPI
      IMPLICIT NONE        
      INTEGER, INTENT(in)                 :: error
      INTEGER                             :: length, temp
      CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
      CALL MPI_ERROR_STRING(error,message,length,temp)
      WRITE(6,*) '!!!!!!!!!!!!MPI_ERROR DETECTED!!!!!!!!!!!!!!'
      WRITE(6,*) '  MESSAGE: ',message(1:length)
      WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      CALL FLUSH(6) 
      !CALL MPI_ABORT(MPI_COMM_STEL,1,temp)
      END SUBROUTINE
!DEC$ ENDIF  
      END MODULE mpi_params
