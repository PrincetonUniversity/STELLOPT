!-----------------------------------------------------------------------
!     Subroutine:    MAP
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/01/2014
!     Description:   This subroutine performs and n-dimensional mapping
!                    of a parameter hyperspace.
!-----------------------------------------------------------------------
      SUBROUTINE MAP(obj, n, m, XCmin, XCmax, NP, iwrite, ndiv, MAP_COMM)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE gade_mod, ONLY: pso_cleanup
      USE safe_open_mod, ONLY: safe_open
      USE mpi_inc
!-----------------------------------------------------------------------
!     Arguments
!        obj         Functional to evaluate
!        n           Length of X-Vector
!        m           Length of f-vector
!        XCmin/max   Min/max values of x-vector elements
!        NP          Population size
!        iwrite      Unit number to write to
!        ndiv        Number of subdivisions to evaluate
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: n,m, ndiv, iwrite, NP
      INTEGER, INTENT(inout) :: MAP_COMM
      DOUBLE PRECISION, DIMENSION(n), INTENT(in) :: XCmin, XCmax
      EXTERNAL  obj
      
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: myid, numprocs
!DEC$ ELSE
      INTEGER, PARAMETER :: myid = 0
      INTEGER, PARAMETER :: numprocs = 1
!DEC$ ENDIF
      INTEGER :: num_search, i, j, itemp, ierr_mpi, numsent, iflag, &
                 num_recd, dex1, dex2, ierr, iunit
      INTEGER , ALLOCATABLE :: bin_step(:)
      DOUBLE PRECISION :: temp_norm
      DOUBLE PRECISION, ALLOCATABLE :: x_temp(:), val_temp(:), dx(:)
      DOUBLE PRECISION, ALLOCATABLE :: grid(:,:),vals(:,:)
      CHARACTER(256) ::  map_file
      
      INTEGER, PARAMETER :: master = 0
      INTEGER, PARAMETER :: FLAG_SINGLETASK = -1
      INTEGER, PARAMETER :: FLAG_CLEANUP = -100
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BARRIER(MAP_COMM, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: BARRIER 1"
      CALL MPI_COMM_RANK( MAP_COMM, myid, ierr_mpi )              ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: COMM_RANK"
      CALL MPI_COMM_SIZE( MAP_COMM, numprocs, ierr_mpi )          ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: COMM_SIZE"
!DEC$ ENDIF
      
      ! The extent of the search space is
      !       ndiv^N
      !   N    :  Number of variables
      !   ndiv :  Number of steps
      num_search = ndiv**n
      ALlOCATE(x_temp(n),val_temp(m),dx(n))
      ALLOCATE(grid(num_search,n),vals(num_search,m), &
               bin_step(n))
      
      ! Set up Grid
      DO i = 1, n
         dx(i) = (XCmax(i)-XCmin(i))/(ndiv-1)
      END DO
      bin_step = 0
      DO i = 1, num_search
         bin_step = 0
         j = n
         itemp = i - 1
         DO WHILE(itemp > 0)
            bin_step(j) = MOD(itemp,ndiv)
            j=j-1
            itemp = itemp / ndiv
         END DO
         grid(i,1:n) = XCmin(1:n) + dx(1:n) * bin_step(1:n)
      END DO
      
      ! Have master write to screen.  Note that stellopt_fcn will attempt
      ! to abort the run if the first evaluation doesn't work so we'll need
      ! to check iflag and handle that case.
      IF (myid == master) THEN
         x_temp(1:n) = grid(1,1:n)
         iflag = FLAG_SINGLETASK
         i=0
         val_temp = 0.0
         CALL obj(m,n,x_temp,val_temp,iflag,i) ! Run 0,0,0 case
         IF (iflag /= 0) THEN
            val_temp(1:m) = 10*SQRT(1.0E10/m)
            vals(1,1:m) = val_temp(1:m)
            temp_norm = SQRT(SUM(vaL_temp(1:m)**2))
         ELSE
            vals(1,1:m) = val_temp(1:m)
            temp_norm = SQRT(SUM(vaL_temp(1:m)**2))
            iflag = pso_cleanup
            i=0
            CALL obj(m,n,x_temp,val_temp,iflag,i)    ! Clean-up 0,0,0 case
         END IF
         WRITE(6, '(2x,A,3x,A,7x,A)') 'Iteration','Processor','Chi-Sq'
         WRITE(6, '(2x,i6,8x,i3,7x,1ES22.12E3)') 1, myid, temp_norm*temp_norm
         CALL FLUSH(6)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MAP_COMM, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: BARRIER 2"
!DEC$ ENDIF
      
!DEC$ IF DEFINED (MPI_OPT)
      ! Queued Workload
      IF (myid .eq. master) THEN
         numsent = 1
         ! First send out work
         DO j = 1, MIN(numprocs-1,num_search)
            numsent = numsent+1
            x_temp(1:n) = grid(numsent,1:n)
            CALL MPI_SEND(x_temp,n,MPI_REAL8,j,numsent,&
                          MAP_COMM,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: MPI_SEND 1"
         END DO
         num_recd = 1
         ! Now Recieve and send work
         DO WHILE (num_recd < num_search)
            ! Get Data
            ierr_mpi = 0
            val_temp = 0.0
            CALL MPI_RECV(val_temp,m,MPI_REAL8, &
                          MPI_ANY_SOURCE, MPI_ANY_TAG, &
                          MAP_COMM, status, ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: MPI_RECV 1"
            num_recd = num_recd + 1
            ! Sort Data
            dex1 = status(MPI_SOURCE)
            dex2 = status(MPI_TAG)
            vals(dex2,1:m) = val_temp(1:m)
            temp_norm = SQRT(SUM(vaL_temp(1:m)**2))
            WRITE(6, '(2x,i6,8x,i3,7x,1ES22.12E3)') dex2, dex1, temp_norm*temp_norm
            CALL FLUSH(6)
            IF (numsent < num_search) THEN
                ! Send More Work
                numsent = numsent + 1
                x_temp(1:n) = grid(numsent,1:n)
                ierr_mpi = 0
                CALL MPI_SEND(x_temp, n, MPI_REAL8, dex1, numsent, MAP_COMM, ierr_mpi)
                IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: MPI_SEND 2"
            ELSE
                ! Send ALL DONE
                ierr_mpi = 0
                CALL MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8, dex1, 0, MAP_COMM, ierr_mpi) 
                IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: MPI_SEND 3"
            END IF
         END DO
      ELSE
         dex2 = -1
         DO WHILE (dex2 /= 0)
            ! Get work
            ierr_mpi = 0
            x_temp = 0.0
            CALL MPI_RECV(x_temp,n,MPI_REAL8,master,MPI_ANY_TAG,MAP_COMM, status, ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: MPI_RECV 2"
            dex2 = status(MPI_TAG)
            IF (dex2 /= 0) THEN
               iflag = dex2
               iunit = dex2
               val_temp = 0.0
               CALL obj(m,n,x_temp,val_temp,iflag,iunit)
               ierr_mpi = 0
               CALL MPI_SEND(val_temp,m,MPI_REAL8,master,dex2,MAP_COMM,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: MPI_SEND 4"
               dex2 = -1
            END IF
         END DO
      END IF
!DEC$ ELSE
      STOP 'LIBSTELL must be compile in parallel to use this feature!'
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MAP_COMM, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) STOP "MAP: BARRIER 3"
!DEC$ ENDIF
   
      IF (myid .eq. master) THEN
         WRITE(6,'(A)') '------- Outputing to map.dat --------'
         CALL FLUSH(6)
         map_file = 'map.dat'
         CALL safe_open(iunit,ierr,TRIM(map_file),'new','formatted')
         WRITE(iunit,'(4(2X,i6))') m,n,ndiv,num_search
         DO i = 1, num_search
            WRITE(iunit,'(1p,4ES22.12E3)') (grid(i,j), j=1,n)
         END DO
         DO i = 1, num_search
            WRITE(iunit,'(1p,4ES22.12E3)') (vals(i,j), j=1,m)
         END DO
         CLOSE(iunit)
      END IF
      
      DEALLOCATE(x_temp, val_temp, dx)
      DEALLOCATE(grid, vals, bin_step)
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MAP_COMM, ierr_mpi)                 !mpi stuff
!DEC$ ENDIF

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE MAP
