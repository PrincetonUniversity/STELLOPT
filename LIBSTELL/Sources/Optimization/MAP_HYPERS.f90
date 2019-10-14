!-----------------------------------------------------------------------
!     Subroutine:    MAP_HYPERS
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/02/2017
!     Description:   This subroutine maps a hyperspace from a point.
!                    This is achieved by computing a set of points
!                    a fixed distance from a point in hyperspace.
!                    This set of points is randomly sampled but
!                    normalized to be a given distance from the point.
!                    Each generation steps a fixed distance away from
!                    the initial point.
!                    Note we use SSEND so communication is blocking
!                    for short messages.
!-----------------------------------------------------------------------
      SUBROUTINE MAP_HYPERS(fcn, m, n, NP, XCmin, XCmax, x, fvec, &
                            drho,erho,maxfev,HYPER_COMM)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod, ONLY: safe_open
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Variables
!        fcn     Target function one wishes to mimizie
!        m       Dimensionality of the function
!        n       Dimensionality of the function space
!        NP      Population size
!        XCmin   Minimum values of function space
!        XCmax   Maximum values of function space
!        x       Vector specifying mimium location in function space
!        fvec    Vector of function values at minimum
!        drho    Distance for each delta step
!        erho    Scaling factor for delta
!        maxfev  Maximum number of function evaluations.
!        HYPER_COMM  MPI Communicator
!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: m, n, NP, maxfev
      REAL(rprec), INTENT(in)  :: drho, erho
      REAL(rprec), INTENT(in)  :: XCmin(n), XCmax(n)
      REAL(rprec), INTENT(inout) :: x(n)
      REAL(rprec), INTENT(out) :: fvec(m)
      INTEGER, INTENT(inout) :: HYPER_COMM
      EXTERNAL  fcn
      
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      INTEGER :: rank, myid, ierr_mpi, iter, istat, nfeval, iunit, &
                 idex, i, j, mpi_req, ntot
      REAL(rprec)    :: rho, fnorm
      REAL, ALLOCATABLE :: r_vec(:)
      REAL(rprec), ALLOCATABLE :: x_temp(:), fvec_temp(:)
      REAL(rprec), ALLOCATABLE :: x_global(:,:), fvec_global(:,:)
      CHARACTER(256) ::  map_file
#if defined(MPI_OPT)
      INTEGER :: mpi_stat(MPI_STATUS_SIZE)
#endif
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      fvec = 0
#if defined(MPI_OPT)
      ierr_mpi = 0
      CALL MPI_COMM_RANK(HYPER_COMM, myid, ierr_mpi)
      CALL MPI_COMM_SIZE(HYPER_COMM, rank, ierr_mpi)
#endif

      ! Error check
      IF (erho .le. 0) THEN
         IF (myid .eq. 0) WRITE(6,'(A)') '!!!!! ERROR: ERHO <= 0 !!!!!'
         RETURN
      END IF
      IF (drho .le. 0) THEN
         IF (myid .eq. 0) WRITE(6,'(A)') '!!!!! ERROR: DRHO <= 0 !!!!!'
         RETURN
      END IF

      ! ALLOCATIONS
      ntot = NP*maxfev
      ALLOCATE(x_global(1:n,0:ntot),fvec_global(1:m,0:ntot),STAT=istat)
      ALLOCATE(x_temp(n), fvec_temp(m), r_vec(n),STAT=istat)

      ! Do first iteration
      IF (myid .eq. 0) THEN
         x_temp        = x
         x_global(:,0) = x
         fvec_temp(:)  = 0
         istat         = -1
         nfeval        = 0
         CALL fcn(m, n, x_temp, fvec_temp, istat, nfeval)
         fvec_global(:,0) = fvec_temp
         fnorm            = SQRT(SUM(fvec_temp*fvec_temp))

         WRITE(6,'(70("="))')
         WRITE(6,'(A)') '      ITER       PROC         CHISQ'
         WRITE(6,'(70("="))')
         WRITE(6,'(2(6X,I5),7X,1ES12.4E3)') nfeval,myid,fnorm**2

      END IF
#if defined(MPI_OPT)
      CALL MPI_BCAST(istat,1,MPI_INTEGER,0,HYPER_COMM,ierr_mpi)
      CALL MPI_BARRIER(HYPER_COMM, ierr_mpi)
#endif
      IF (istat < 0) THEN
         IF (myid == 0) WRITE(6,'(A)') '!!!!! First Eq failed !!!!!'
         DEALLOCATE(x_temp,fvec_temp)
         IF (ALLOCATED(x_global)) DEALLOCATE(x_global)
         IF (ALLOCATED(fvec_global)) DEALLOCATE(fvec_global)
         RETURN
      END IF

      ! Cleanup after first iteration
      istat = -100
      IF (myid .eq. 0) istat = -101
      call fcn (m, n, x_temp, fvec_temp, istat, nfeval)
      nfeval = 1

      ! Setup Random numbers
      CALL RANDOM_SEED()  ! Processor reinitializes the seed

      ! Setup Global array
      idex = 1
      DO iter = 1, maxfev
         rho = drho*(REAL(iter)/REAL(maxfev))**erho
         DO i = 1, NP
            CALL RANDOM_NUMBER(r_vec)
            x_global(:,i) = x+rho*r_vec/SQRT(SUM(r_vec*r_vec))
            idex = idex + 1
         END DO
      END DO

      ! Evaluate
      CALL eval_x_queued(fcn,m,n,ntot,x_global(1:n,1:ntot), &
                         fvec_global(1:m,1:ntot),nfeval,HYPER_COMM)

      ! Write File   
      IF (myid .eq. 0) THEN
         map_file = 'map_hypers.dat'
         istat = 0
         CALL safe_open(iunit,istat,TRIM(map_file),'new','formatted')
         WRITE(iunit,'(4(2X,i6))') m,n,ntot+1
         DO i = 0, ntot
            WRITE(iunit,'(1p,4ES22.12E3)') (x_global(j,i), j=1,n)
         END DO
         DO i = 0, ntot
            WRITE(iunit,'(1p,4ES22.12E3)') (fvec_global(j,i), j=1,m)
         END DO
         CLOSE(iunit)
      END IF

      ! DEALLOCATE
      DEALLOCATE(x_temp,fvec_temp,r_vec)
      IF (ALLOCATED(x_global)) DEALLOCATE(x_global)
      IF (ALLOCATED(fvec_global)) DEALLOCATE(fvec_global)
      
      
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE MAP_HYPERS
