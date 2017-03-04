      SUBROUTINE ga_fitness_mpi (np, fvec, nopt, fcn, nfev, funcval)
      USE ga_mod
      USE mpi_params, ONLY: master, myid, MPI_COMM_STEL
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
      REAL(rprec), DIMENSION(nparam) :: x
!DEC$ ENDIF

      INTEGER :: np, nopt, nfev
      REAL(rprec), DIMENSION(nopt) :: fvec
      REAL(rprec) :: funcval(np)
      EXTERNAL fcn

!DEC$ IF DEFINED (MPI_OPT)
      INTEGER :: j, iflag, istat

      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: numprocs, i                           !mpi stuff
      INTEGER :: numsent, sender, ierr                       !mpi stuff
      INTEGER :: anstype, column                             !mpi stuff

!******************************************
!
!  mpi setup calls; set barrier so ALL processors get here before starting
!
      CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr )       !mpi stuff
      CALL MPI_COMM_size( MPI_COMM_STEL, numprocs, ierr )   !mpi stuff
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr)                 !mpi stuff

!******************************************
!
!     ****Master portion of the code****
!
      IF (myid .eq. master) THEN
         numsent = 0    !numsent is a counter used to track how many
                        !jobs have been sent to workers
c
c     SEND forward difference displacements from master to each
c           worker process. Tag with these with the column number.
c
         DO j = 1,MIN(numprocs-1,np)
            x(:) = parent(:,j)
            CALL MPI_SEND(x, nparam, MPI_REAL8, j,
     1                  j, MPI_COMM_STEL, ierr)
            IF (ierr .ne. 0) STOP 'MPI_SEND error(1) in ga_fitness_mpi'
            numsent = numsent+1
         END DO          !j = 1,MIN(numprocs-1,n)
c
c      Looping through the columns, collect answers from the workers.
c      As answers are received, new uncalculated columns are sent
c      out to these same workers.
c
         DO j = 1,np
            CALL MPI_RECV(fvec, nopt, MPI_REAL8,
     1           MPI_any_SOURCE, MPI_any_TAG,
     2           MPI_COMM_STEL, status, ierr)
            IF (ierr .ne. 0) STOP 'MPI_RECV error(1) in ga_fitness_mpi'
            sender     = status(MPI_SOURCE)
            anstype    = status(MPI_TAG)       ! column is tag value
            IF (anstype .gt. np) STOP 'ANSTYPE > NP IN ga_fitness_mpi'

            funcval(anstype) = -SUM(fvec(:nopt)**2)
            WRITE(6,'(a,1pe10.3,a,i3,a,i3)')' FUNCVAL = ',
     1       -funcval(anstype),
     2      ' for iteration ', anstype+nfev,' processor = ', sender

c
c           If more columns are left, THEN sEND another column to the worker(sender)
c           that just sent in an answer
c
            IF (numsent .lt. np) THEN
               numsent = numsent+1
               x(:) = parent(:,numsent)

               CALL MPI_SEND(x, nparam, MPI_REAL8,
     1                       sender, numsent, MPI_COMM_STEL, ierr)
               IF (ierr .ne. 0)
     1            STOP 'MPI_SEND error(2) in ga_fitness_mpi'

            ELSE                ! Tell worker that there is no more work to DO

               CALL MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8,
     1                       sender, 0, MPI_COMM_STEL, ierr)
               IF (ierr .ne. 0)STOP 'MPI_end error(3) in ga_fitness_mpi'
            ENDIF      ! IF( myid .eq. master ) THEN
         END DO     ! DO j = 1,n
c
c     ****Worker portion of the code****
c        Skip this when processor id exceeds work to be done
c
      ELSE IF (myid .le. np) THEN        ! i.e., IF( myid .ne. master )
c
c        Otherwise accept the next available column, check the tag,
c        and IF the tag is non-zero CALL SUBROUTINE fcn.
c        If the tag is zero, there are no more columns
c        and worker skips to the END.
c
 90      CALL MPI_RECV(x, nparam, MPI_REAL8, master,
     1                 MPI_any_TAG, MPI_COMM_STEL, status, ierr)
         IF (ierr .ne. 0) STOP 'MPI_RECV error(2) in ga_fitness_mpi'

         column = status(MPI_TAG)                !!ID of pseudo-processor issuing this message
         IF (column .eq. 0) THEN
            GOTO 200
         ELSE
            iflag = column
c           CALL the chisq fcn for the portion of displacement vector which
c           was just received. Note that WA stores the local fvec_min array

            CALL fcn(nopt, nparam, x, fvec, iflag, nfev)
            IF (iflag.ne.0) GOTO 300
c
c           Send this function evaluation back to the master process tagged
c           with the column number so the master knows where to put it
c
            CALL MPI_SEND(fvec, nopt, MPI_REAL8, master,
     1                    column, MPI_COMM_STEL, ierr)
            IF (ierr .ne. 0) STOP 'MPI_SEND error(4) in ga_fitness_mpi'
            GOTO 90    !Return to 90 and check IF master process has sent ANY more jobs
         END IF
 200     CONTINUE
      ENDIF       ! IF( myid .ne. master )

!
!     Broadcast the funcval array to ALL processors FROM master
!
      CALL MPI_BCAST(funcval, np, MPI_REAL8, master,
     1     MPI_COMM_STEL, ierr)
      IF (ierr .ne. 0) GOTO 100

      RETURN

 100  CONTINUE
      WRITE (6,*) ' MPI_BCAST error in ga_fitness_mpi: IERR=', ierr

      RETURN

 300  CONTINUE
      WRITE (6,*) ' IFLAG = ', iflag, ' in ga_fitness_mpi CALL to fcn'
      STOP

!DEC$ ENDIF
      END SUBROUTINE ga_fitness_mpi
