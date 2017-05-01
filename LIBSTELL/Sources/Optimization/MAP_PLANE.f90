      SUBROUTINE MAP_PLANE(obj, n, m, XC0, XC1, XC2, factor, NP, iwrite, ndiv)
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
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
!DEC$ ENDIF

      INTEGER, INTENT(in) :: n,m, ndiv, iwrite, NP
      REAL(rprec), INTENT(inout) :: factor
      REAL(rprec), DIMENSION(n), INTENT(in) :: XC0, XC1, XC2
      EXTERNAL  obj
      
      LOGICAL :: lexist
      INTEGER :: i,j,k,nbin,dex_base,dex1,dex2,itemp,num_search,&
                 numsent, ierr,bin_val, iflag, iflag2, ierr_flag, &
                 iunit, num_recd, n1, n2
      REAL(rprec) :: temp_norm, enorm
      REAL(rprec) :: XVMIN(2),XVMAX(2), factor2(2)
      REAL(rprec), ALLOCATABLE :: x_temp(:), val_temp(:), dx1(:), dx2(:)
      REAL(rprec), ALLOCATABLE :: grid(:,:),vals(:,:)
      CHARACTER(256) ::  map_file
      
      INTEGER, PARAMETER :: FLAG_SINGLETASK = -1
      INTEGER, PARAMETER :: FLAG_CLEANUP = -100

!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
!DEC$ ENDIF


      ! This is done to allow customized search domains

      n1=0; n2=0
      INQUIRE(FILE='MAP_VECTOR',EXIST=lexist)
      IF (lexist) THEN
         ierr = 0
         CALL safe_open(iunit,ierr,'MAP_VECTOR','old','formatted')
         READ(iunit,*) n1, XVMIN(1), XVMAX(1)
         READ(iunit,*) n2, XVMIN(2), XVMAX(2)
         CLOSE(IUNIT)
         factor = 0
      ELSE
         n1 = ndiv
         n2 = ndiv
      END IF
      

      ! The extent of the search space is
      !       ndiv*ndiv
      num_search = n1*n2
      ALLOCATE(x_temp(n),val_temp(m),dx1(n),dx2(n))
      ALLOCATE(grid(num_search,n),vals(num_search,m))

      IF (myid == master) THEN
        ! Debugging
        DO i=1, n
           WRITE(6,'(I4,1X,ES20.10,1X,ES20.10,1X,ES20.10)') i,XC0(i),XC1(i),XC2(i)
        END DO
        ! Initialize starting point
        IF (factor < 0) THEN
           x_temp = XC0 + factor*(XC1-XC0) &
                        + factor*(XC2-XC0)
           factor2(1) = 2*ABS(factor)
           factor2(2) = 2*ABS(factor)
        ELSEIF (lexist) THEN
           x_temp = XC0 + XVMIN(1)*(XC1-XC0) &
                        + XVMIN(2)*(XC2-XC0)
           factor2(1) = XVMAX(1)-XVMIN(1)
           factor2(2) = XVMAX(2)-XVMIN(2)
           WRITE(6,'(A)') '!!!!!! USING DATA FROM MAP_VECTOR !!!!!!!'
           WRITE(6,'(A,I6)') 'NUMSEARCH: ',num_search 
           WRITE(6,'(A,I4,1X,F5.2,1X,F5.2)') 'N1, XVMIN(1), XVMAX(1): ',n1, XVMIN(1), XVMAX(1)
           WRITE(6,'(A,I4,1X,F5.2,1X,F5.2)') 'N2, XVMIN(2), XVMAX(2): ',n2, XVMIN(2), XVMAX(2)
        ELSE
           x_temp = XC0
           factor2(1) = factor
           factor2(2) = factor
        ENDIF

        ! Set up Grid
        dx1=factor2(1)*(XC1-XC0)/(n1-1)
        dx2=factor2(2)*(XC2-XC0)/(n2-1)
        k=1
        DO i = 1, n1
           DO j = 1, n2
              grid(k,:) = x_temp + dx1*REAL(i-1) &
                              + dx2*REAL(j-1)
              k=k+1
           END DO
        END DO
      END IF

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
            iflag = -110
            i=0
            CALL obj(m,n,x_temp,val_temp,iflag,1)    ! Create the 00001 file
         END IF
         WRITE(6, '(2x,A,3x,A,7x,A)') 'Iteration','Processor','Chi-Sq'
         WRITE(6, '(2x,i6,8x,i3,7x,1ES22.12E3)') 1, myid, temp_norm*temp_norm
         CALL FLUSH(6)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      
!DEC$ IF DEFINED (MPI_OPT)
      ! Queued Workload
      IF (myid .eq. master) THEN
         numsent = 1
         ! First send out work
         DO j = 1, MIN(numprocs-1,num_search-1)
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
            CALL FLUSH(6)
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
      ELSEIF (myid < num_search) THEN
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
               iflag = -110 ! JUST Input
               CALL obj(m,n,x_temp,val_temp,iflag,iunit)    ! Output the input files and 
               ierr_mpi = 0
               CALL MPI_SEND(val_temp,m,MPI_REAL8,master,dex2,MPI_COMM_STEL,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
               dex2 = -1
            END IF
         END DO
      END IF
!DEC$ ELSE
      STOP 'LIBSTELL must be compile for parallel to use this feature!'
!DEC$ ENDIF


!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
   
      IF (myid .eq. master) THEN
         map_file = 'map_plane.dat'
         WRITE(6,'(A)') '------- Outputing to map_plane.dat --------'
         CALL FLUSH(6)
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
      
      DEALLOCATE(x_temp, val_temp, dx1, dx2)
      DEALLOCATE(grid, vals)
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
!DEC$ ENDIF
      
      END SUBROUTINE MAP_PLANE
