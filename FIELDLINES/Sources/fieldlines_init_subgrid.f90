!-----------------------------------------------------------------------
!     Module:        fieldlines_init_subgrid
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/29/2012
!     Description:   This subroutine initializes fieldline starting
!                    points using the good flux surfaces.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_subgrid
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl
      USE fieldlines_lines
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i,j,ier, n1,n2,mystart
      REAL(rprec) :: r1, r2, r0, z0, z1, z2, rinner, zinner, phiend_temp,&
                     phi0, ra, za, phia
      REAL(rprec), ALLOCATABLE :: delta_r(:), delta_z(:), delta_l(:)
!-----------------------------------------------------------------------
!     External Functions
!          A00ADF               NAG Detection
!-----------------------------------------------------------------------
!      EXTERNAL A00ADF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      phiend_temp = pi2
      ra = r_start(1); za = z_start(1); phia=phi_start(1)
      mystart = LBOUND(R_lines,1)
      ! Get inner point
      !STOP 'DISABLED Due to upgrades'
      IF (myid == master) THEN
         WRITE(6,'(A)') '===========EDGE SEARCH=========='
         CALL FLUSH(6)
      END IF
      ! Get outer point
      i   = myend
      DO WHILE (i >= mystart)
         ier = 0
         j   = 1
         DO WHILE (j < nsteps)
            IF ((R_lines(i,j) == R_lines(i,j+1)) .and.    &
               (Z_lines(i,j) == Z_lines(i,j+1))) ier = 1
            IF (ier == 1) EXIT
            j = j + 1
         END DO
         IF (ier == 0) EXIT
         i = i - 1
      END DO
      IF (i < mystart) i = 0 ! Do this so MPI_ALLREDUCE WORKS
      n2 = i + 1
      IF (n2 > nlines) n2 = nlines
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,n2,1,MPI_INTEGER,MPI_MAX,MPI_COMM_FIELDLINES,ierr_mpi)
      n1 = n2 - 2
      r1 = 0; z1 = 0; r2 = 0; z2 = 0
      IF (n1 >= mystart .and. n1<= myend) THEN
         r1 = R_lines(n1,0) ! Should be outter good line
         z1 = Z_lines(n1,0) ! Should be outter good line
      END IF
      IF (n2 >= mystart .and. n2<= myend) THEN
         r2 = R_lines(n2,0)   ! First bad line
         z2 = Z_lines(n2,0)   ! First bad line
      END IF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,z1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,r2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,z2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FIELDLINES,ierr_mpi)

         
      ! Refine grid from outter good to outter bad line
      !nlines = nr
      phiend_temp = MAXVAL(phi_end)
      r_start = -1
      z_start = -1
      phi_start = 0.0
      phi_end = 0.0
      IF (nlines > MAXLINES) nlines = MAXLINES-2
      DO i = 1, nlines
         r_start(i) = (i-1)*(r2-r1)/(nlines-1) + r1
         z_start(i) = (i-1)*(z2-z1)/(nlines-1) + z1
         phi_start(i) = phia
         phi_end(i)   = phiend_temp
      END DO
      IF (myid == master) THEN
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   AXIS_GUESS[R,Z]   = [',ra,',',za,'];'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_OUTER[R,Z]   = [',r1,',',z1,'];'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_LOST[R,Z]    = [',r2,',',z2,'];'
         WRITE(6,'(A,I4)')            '   # of Fieldlines: ',nlines
         CALL FLUSH(6)
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
      CALL MPI_BCAST(nlines,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, numprocs, ierr_mpi )          ! MPI
!DEC$ ENDIF
   
      ! We follow the edge lines
      CALL fieldlines_follow
      mystart = LBOUND(R_lines,1)
      
      ! Refine outer point
      i   = myend
      DO WHILE (i >= mystart)
         ier = 0
         j   = 1
         DO WHILE (j < nsteps)
            IF ((R_lines(i,j) == R_lines(i,j+1)) .and.    &
                (Z_lines(i,j) == Z_lines(i,j+1))) ier = 1
            IF (ier == 1) EXIT
            j = j + 1
         END DO
         IF (ier == 0) EXIT
         i = i - 1
      END DO
      IF (i < mystart) i = 0 ! Do this so MPI_ALLREDUCE WORKS
      n2 = i + 1
      IF (n2 > nlines) n2 = nlines
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,n2,1,MPI_INTEGER,MPI_MAX,MPI_COMM_FIELDLINES,ierr_mpi)
      r2 = 0; z2 = 0
      IF (n2 >= mystart .and. n2<= myend) THEN
         r2 = R_lines(n2,0)   ! First bad line
         z2 = Z_lines(n2,0)   ! First bad line
      END IF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,r2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,z2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FIELDLINES,ierr_mpi)
         
      IF (myid == master) THEN
         ! Get the axis
         r1 = r0
         r0 = ra
         z0 = za
         phi0 = phia
         CALL fieldlines_find_axis(r0,z0,phi0)
      END IF
      
      IF (COUNT(r_hc > 0) > 0) THEN
         IF (myid == master) THEN
            CALL fieldlines_find_hc
         END IF

!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
         CALL MPI_BCAST(nlines,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
         CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
         CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
         CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
         CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
         CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, numprocs, ierr_mpi )          ! MPI
!DEC$ ENDIF
      
         ! Follow homoclines
         CALL fieldlines_follow  
         mystart = LBOUND(R_lines,1)
        
         ! Store the Homoclines
         ALLOCATE(Rhc_lines(mystart:myend,0:nsteps),Zhc_lines(mystart:myend,0:nsteps))
         Rhc_lines = R_lines
         Zhc_lines = Z_lines
      END IF
         
      ! Final Grid
      !nlines = nr
      IF (nlines > MAXLINES) nlines = MAXLINES-2
      r_start = -1
      z_start = -1
      phi_start = 0.0
      phi_end = 0.0
      DO i = 1, nlines
         r_start(i) = (i-1)*(r2-r0)/(nlines-1) + r0
         z_start(i) = (i-1)*(z2-z0)/(nlines-1) + z0
         phi_start(i) = phia
         phi_end(i)   = phiend_temp
      END DO
      IF (myid == master) THEN ! Final grid from axis to refined grid
         WRITE(6,'(A)') '===========FINAL GRID=========='
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   AXIS[R,Z]   = [',r0,',',z0,']'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE[R,Z]   = [',r2,',',z2,']'
         WRITE(6,'(A,I4)')            '   # of Fieldlines: ',nlines
         CALL FLUSH(6)
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
       CALL MPI_BCAST(nlines,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, numprocs, ierr_mpi )          ! MPI
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
      
      
      CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
      CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_main',ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_subgrid
