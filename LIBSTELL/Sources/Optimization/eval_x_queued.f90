!-----------------------------------------------------------------------
!     Subroutine:    EVAL_X_QUEUED
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/02/2017
!     Description:   This subroutine evaluates an array of xvectors
!                    and returns them in fvec.  It uses a queued
!                    execution to achieve this.
!-----------------------------------------------------------------------
      SUBROUTINE EVAL_X_QUEUED(fcn, m, n, NP, xvec, fvec, iter, &
                               HYPER_COMM)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Variables
!        fcn     Target function one wishes to mimizie
!        m       Dimensionality of the function
!        n       Dimensionality of the function space
!        NP      Population size
!        xvec    Vector specifying mimium location in function space
!        fvec    Vector of function values at minimum
!        HYPER_COMM  MPI Communicator
!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: m, n, NP
      INTEGER, INTENT(inout) :: iter
      REAL(rprec), INTENT(inout) :: xvec(n,NP)
      REAL(rprec), INTENT(inout) :: fvec(m,NP)
      INTEGER, INTENT(inout) :: HYPER_COMM
      EXTERNAL  fcn
      
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      LOGICAL :: lsent
      INTEGER :: myid, rank, ierr_mpi, numsent, numrecd, idex, i, &
                 sender, anstype, istat
      REAL(rprec)    :: fnorm
      REAL, ALLOCATABLE :: r_vec(:)
      REAL(rprec), ALLOCATABLE :: x_temp(:), fvec_temp(:)
      CHARACTER(256) ::  map_file
#if defined(MPI_OPT)
      INTEGER :: mpi_stat(MPI_STATUS_SIZE)
#endif
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
#if defined(MPI_OPT)
      ierr_mpi = 0
      CALL MPI_COMM_RANK(HYPER_COMM, myid, ierr_mpi)
      CALL MPI_COMM_SIZE(HYPER_COMM, rank, ierr_mpi)
#endif

      ! ALLOCATIONS
      ALLOCATE(x_temp(n), fvec_temp(m), r_vec(n),STAT=istat)

      numsent = 0
      numrecd = 0
      idex    = 0

#if defined(MPI_OPT)
      ! Work Loop
      IF (myid == 0) THEN ! Master Loop
         ! Send initial datasets
         sender = MIN(rank-1,NP)
         DO i = 1, sender
            idex = idex + 1
            x_temp=xvec(:,idex)
            CALL MPI_SSEND(x_temp,n,MPI_DOUBLE,i,idex,HYPER_COMM,ierr_mpi)
            numsent = numsent+1
         END DO

         ! Handle NP < rank-1
         IF (sender < rank-1) THEN
            DO i = sender+1, rank-1
               CALL MPI_SSEND(MPI_BOTTOM,0,MPI_DOUBLE,i,0, &
                                HYPER_COMM,ierr_mpi)
            END DO
         END IF
          
         ! Result/More work loop
         DO
            CALL MPI_RECV(fvec_temp,m,MPI_DOUBLE,MPI_ANY_SOURCE, &
                          MPI_ANY_TAG,HYPER_COMM,mpi_stat,ierr_mpi)
            numrecd         = numrecd + 1
            sender          = mpi_stat(MPI_SOURCE)
            anstype         = mpi_stat(MPI_TAG)
            fvec(:,anstype) = fvec_temp
            fnorm           = SQRT(SUM(fvec_temp*fvec_temp))
            WRITE(6,'(2X,I6,8X,I3,7X,1ES12.4E2)') anstype,sender,&
                                                  fnorm**2
            CALL FLUSH(6)

            ! Send Work
            IF (numsent .lt. NP) THEN  ! More work to do 
               idex = idex + 1
               x_temp = xvec(:,idex)
               CALL MPI_SSEND(x_temp,n,MPI_DOUBLE,sender,idex, &
                                HYPER_COMM,ierr_mpi)
               numsent = numsent + 1
            ELSE ! No more work to do
               CALL MPI_SSEND(MPI_BOTTOM,0,MPI_DOUBLE,sender,0, &
                                HYPER_COMM,ierr_mpi)
            END IF 
               
            ! Exit if all work received
            IF (numrecd .eq. NP) EXIT
         END DO
      ELSE ! Worker Loop 
         DO
            CALL MPI_RECV(x_temp,n,MPI_DOUBLE,0, &
                          MPI_ANY_TAG,HYPER_COMM,mpi_stat,ierr_mpi)
            anstype    = mpi_stat(MPI_TAG)
            IF (anstype .ne. 0) THEN
               istat = anstype
               fvec_temp = 0
               CALL fcn(m,n,x_temp,fvec_temp,istat,iter)
               iter = iter + 1
               CALL MPI_SSEND(fvec_temp,m,MPI_DOUBLE,0, &
                             anstype,HYPER_COMM,ierr_mpi)
            ELSE
               EXIT
            END IF   
         END DO            
      END IF
      CALL MPI_BARRIER(HYPER_COMM, ierr_mpi)
!      iter = iter + NP
      CALL MPI_BCAST(iter,1,MPI_INTEGER,0,HYPER_COMM,ierr_mpi)
#endif

      ! DEALLOCATE
      DEALLOCATE(x_temp,fvec_temp)
      
      
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE EVAL_X_QUEUED
