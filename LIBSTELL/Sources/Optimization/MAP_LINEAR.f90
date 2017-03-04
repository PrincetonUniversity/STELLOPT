      SUBROUTINE MAP_LINEAR(obj, n, m,  XC0, XCmin, XCmax, NP, iwrite, ndiv)
!.......................................................................
!
! PARAMETER SPACE MAPPING FUNCTION
!
!.......................................................................
!
!  This Fortran 90 PROGRAM is written by Dr. Samuel Lazerson
!.........................................................................
!                obj : The User provided file for evlauting the objective function.
!                      SUBROUTINE obj(xc,fitness)
!                      WHERE "xc" is the REAL decision PARAMETER vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the REAL decision parameters.
!      XCMIN(Dim_XC) : The lower bound of the REAL decision parameters.
!      XCMAX(Dim_XC) : The upper bound of the REAL decision parameters.
!                 NP : Population SIZE.
!             iWRITE : The unit specfier for writing to an EXTERNAL data file.

!    best_XC(Dim_XC) : The best REAL decision parameters.  If non-zero, this
!                      array also CONTAINS an initial guess that is included
!                      in the population
!

      USE stel_kinds
      USE mpi_params
      USE gade_mod, ONLY: pso_cleanup
      USE safe_open_mod, ONLY: safe_open
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF

      INTEGER, INTENT(in) :: n,m, ndiv, iwrite, NP
      REAL(rprec), DIMENSION(n), INTENT(in) :: XCmin, XCmax, XC0
      EXTERNAL  obj
      
      INTEGER :: i,j,k,nbin,dex_base,dex1,dex2,itemp,num_search,&
                 numsent, ierr,bin_val, iflag, iflag2, ierr_flag, &
                 iunit, num_recd
      REAL(rprec) :: temp_norm, enorm
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER , ALLOCATABLE :: bin_step(:)
      REAL(rprec), ALLOCATABLE :: x_temp(:), val_temp(:), dx(:)
      REAL(rprec), ALLOCATABLE :: grid(:,:),vals(:,:)
      CHARACTER(256) ::  map_file
      
      INTEGER, PARAMETER :: FLAG_SINGLETASK = -1
      INTEGER, PARAMETER :: FLAG_CLEANUP = -100

!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
!DEC$ ENDIF
      
      ! The extent of the search space is
      !       ndiv^N
      !   N    :  Number of variables
      !   ndiv :  Number of steps
      num_search = 2*n + 1
      !IF (myid == master) PRINT *,'Number of searches',num_search
      ALlOCATE(x_temp(n),val_temp(m),dx(n))
      ALLOCATE(grid(num_search,n),vals(num_search,m))
      
      ! Set up Grid
      DO i = 1, n
         dx(i) = (XCmax(i)-XCmin(i))/(ndiv-1)
      END DO
      grid(1,1:n) = XC0(1:n)
      itemp = 2
      DO i = 1, n
         grid(itemp,1:n) = XC0(1:n)
         grid(itemp,i)   = XC0(i) - dx(i)
         itemp = itemp + 1
      END DO
      DO i = 1, n
         grid(itemp,1:n) = XC0(1:n)
         grid(itemp,i)   = XC0(i) + dx(i)
         itemp = itemp + 1
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
            temp_norm = enorm(m,val_temp)
         ELSE
            vals(1,1:m) = val_temp(1:m)
            temp_norm = enorm(m,val_temp)
            iflag = pso_cleanup
            i=0
            CALL obj(m,n,x_temp,val_temp,iflag,i)    ! Clean-up 0,0,0 case
         END IF
         WRITE(6, '(2x,A,3x,A,7x,A)') 'Iteration','Processor','Chi-Sq'
         WRITE(6, '(2x,i6,8x,i3,7x,1ES22.12E3)') 1, myid, temp_norm*temp_norm
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      
      ! Queued Workload
      IF (myid .eq. master) THEN
         numsent = 1
         ! First send out work
         DO j = 1, MIN(numprocs-1,num_search)
            numsent = numsent+1
            x_temp(1:n) = grid(numsent,1:n)
            CALL MPI_SEND(x_temp,n,MPI_REAL8,j,numsent,&
                          MPI_COMM_STEL,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
         END DO
         num_recd = 1
         ! Now Recieve and send work
         DO WHILE (num_recd < num_search)
            ! Get Data
            ierr_mpi = 0
            val_temp = 0.0
            CALL MPI_RECV(val_temp,m,MPI_REAL8, &
                          MPI_ANY_SOURCE, MPI_ANY_TAG, &
                          MPI_COMM_STEL, status, ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
            num_recd = num_recd + 1
            ! Sort Data
            dex1 = status(MPI_SOURCE)
            dex2 = status(MPI_TAG)
            vals(dex2,1:m) = val_temp(1:m)
            temp_norm = enorm(m,val_temp)
            WRITE(6, '(2x,i6,8x,i3,7x,1ES22.12E3)') dex2, dex1, temp_norm*temp_norm
            IF (numsent < num_search) THEN
                ! Send More Work
                numsent = numsent + 1
                x_temp(1:n) = grid(numsent,1:n)
                ierr_mpi = 0
                CALL MPI_SEND(x_temp, n, MPI_REAL8, dex1, numsent, MPI_COMM_STEL, ierr_mpi)
                IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
                !PRINT *,'sent ', numsent,'/',num_search
            ELSE
                ! Send ALL DONE
                ierr_mpi = 0
                CALL MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8, dex1, 0, MPI_COMM_STEL, ierr_mpi) 
                IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
            END IF
         END DO
      ELSE
         dex2 = -1
         DO WHILE (dex2 /= 0)
            ! Get work
            ierr_mpi = 0
            x_temp = 0.0
            CALL MPI_RECV(x_temp,n,MPI_REAL8,master,MPI_ANY_TAG,MPI_COMM_STEL, status, ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
            dex2 = status(MPI_TAG)
            IF (dex2 /= 0) THEN
               iflag = dex2
               iunit = dex2
               val_temp = 0.0
               CALL obj(m,n,x_temp,val_temp,iflag,iunit)
               ierr_mpi = 0
               CALL MPI_SEND(val_temp,m,MPI_REAL8,master,dex2,MPI_COMM_STEL,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
               dex2 = -1
            END IF
         END DO
      END IF


!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
   
      IF (myid .eq. master) THEN
         map_file = 'map_linear.dat'
         CALL safe_open(iunit,ierr,TRIM(map_file),'new','formatted')
         WRITE(iunit,'(4(2X,i6))') m,n,num_search
         DO i = 1, num_search
            WRITE(iunit,'(1p,4ES22.12E3)') (grid(i,j), j=1,n)
         END DO
         DO i = 1, num_search
            WRITE(iunit,'(1p,4ES22.12E3)') (vals(i,j), j=1,m)
         END DO
         CLOSE(iunit)
      END IF
      
      DEALLOCATE(x_temp, val_temp, dx)
      DEALLOCATE(grid, vals)
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
!DEC$ ENDIF
      
      END SUBROUTINE MAP_LINEAR
