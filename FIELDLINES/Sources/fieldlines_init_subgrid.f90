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
      USE mpi_params                                                    ! MPI
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF  
      INTEGER :: i,j,ier, n1,n2
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
      ! Get inner point
      STOP 'DISABLED Due to upgrades'
      IF (myid == master) THEN
         WRITE(6,'(A)') '===========EDGE SEARCH=========='
         CALL FLUSH(6)
         ! Get inner point (gives a better axis guess)
         !i   = 1
         !DO WHILE (i <= nlines)
         !   ier = 0
         !   j   = 1
         !   DO WHILE (j < nsteps)
         !     IF ((R_lines(i,j) == R_lines(i,j+1)) .and.    &
         !         (Z_lines(i,j) == Z_lines(i,j+1))) ier = 1
         !     IF (ier == 1) EXIT
         !     j = j + 1
         !  END DO
         !  IF (ier == 0) EXIT
         !  i = i + 1
         !END DO
         !PRINT *,i,j
         !n1 = i - 1
         !if (n1 < 1) n1 = 1
         !if (n1 > nlines-1) n1 = nlines-1
         !r1 = R_lines(n1,0)   ! First bad line (for a better axis guess)
         !z1 = Z_lines(n1,0)   ! First bad line (for a better axis guess)
         !rinner = r1
         !zinner = z1
         !WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_INNER[R,Z]   = [',rinner,',',zinner,'];'
         
         ! Get outer point
         i   = nlines
         DO WHILE (i >= 1)
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
         n2 = i + 1
         if (n2 > nlines) n2 = nlines
         r2 = R_lines(n2,0)   ! First bad line
         z2 = Z_lines(n2,0)   ! First bad line
         r1 = R_lines(n2-2,0) ! Should be outter good line
         z1 = Z_lines(n2-2,0) ! Should be outter good line
         
         ! Find closest point to axis
         !r_start(1),z_start(1),phi_start(1) is axis guess
         !ALLOCATE(delta_r(nlines),delta_z(nlines),delta_l(nlines))
         !delta_r(:) = 1.0E30; delta_z(:) = 1.0E30; delta_l(:) = 1.0E30
         !delta_r(n1+1:n2-1) = R_lines(n1+1:n2-1,npoinc+1) - R_lines(n1+1:n2-1,0)
         !delta_z(n1+1:n2-1) = Z_lines(n1+1:n2-1,npoinc+1) - Z_lines(n1+1:n2-1,0)
         !delta_l = SQRT(delta_r*delta_r+delta_z*delta_z)
         !n1 = MINLOC(delta_l,DIM=1)
         !r0 = R_lines(n1,0)
         !z0 = Z_lines(n1,0)
         !DEALLOCATE(delta_r,delta_z,delta_l)
         
         ! Refine grid from outter good to outter bad line
         nlines = nr
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
      
      IF (myid == master) THEN! Get outer point
         i   = nlines
         DO WHILE (i >= 1)
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
         n2 = i + 1
         IF (n2 > nlines) n2 = nlines
         r2 = R_lines(n2,0)   ! First bad line
         z2 = Z_lines(n2,0)   ! First bad line
         
         ! Get the axis
         r1 = r0
         r0 = ra
         z0 = za
         phi0 = pi2/10
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
        
         IF (myid == master) THEN ! Store the Homoclines
            ALLOCATE(Rhc_lines(nlines,0:nsteps),Zhc_lines(nlines,0:nsteps))
            Rhc_lines = R_lines
            Zhc_lines = Z_lines
         END IF
      END IF
         
      IF (myid == master) THEN ! Final grid from axis to refined grid
         nlines = nr
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
         WRITE(6,'(A)') '===========FINAL GRID=========='
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   AXIS[R,Z]   = [',r0,',',z0,']'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE[R,Z]   = [',r2,',',z2,']'
         WRITE(6,'(A,I4)')            '   # of Fieldlines: ',nlines
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
