!-----------------------------------------------------------------------
!     Module:        fieldlines_gridgen
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/07/2023
!     Description:   This subroutine attempts to generate a grid.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_gridgen
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
      USE mpi_sharmem
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, j, k, i1, i2, j2, k2, l, l1, l2, n1, ier, mystart
      REAL(rprec) :: r1, r2, r0, z0, z1, z2, rinner, zinner, phiend_temp,&
                     phi0, ra, za, phia, dr, dz
      LOGICAL, DIMENSION(:,:,:), POINTER :: lgoodline
      INTEGER :: win_lgoodline
      INTEGER :: MPI_COMM_LOCAL
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! First search for the magnetic axis
      ra = r_start(1); za = z_start(1); phia=phi_start(1)
      IF (myworkid == master) THEN
         ! Get the axis
         r1 = r0
         r0 = ra
         z0 = za
         phi0 = phia
         !WRITE(6,'(A,3(F8.5,A))')     '   AXIS_GUESS[R,PHI,Z] = [',ra,',',phia,',',za,'];'
         CALL fieldlines_find_axis(r0,z0,phi0)
         WRITE(6,'(A,3(F8.5,A))') '         AXIS[R,PHI,Z] = [',r0,',',phi0,',',z0,']'
         ra = r0; za = z0; phia = phi0
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_gridgen',ierr_mpi)
      ! Now allocate the helper
      CALL mpialloc(lgoodline, nr, nphi, nz, myid_sharmem, master, MPI_COMM_SHARMEM, win_lgoodline)
      IF (myid_sharmem == master) lgoodline = .FALSE.
#endif

      ! Now do a grid search
      IF (myworkid == master) THEN
         r_start = -1; z_start = -1; phi_start = 0; phi_end = 0;
         n1 = 1
      !  Do just Rmajor at Zaxis
         DO i = 2, nr-1
            r_start(n1) = raxis(i)
            z_start(n1) = za
            phi_start(n1)   = phiaxis(1)
            phi_end(n1)     = phiaxis(1)+10*pi2
            n1 = n1 + 1
         END DO
      !  Do the entire PHI=0 Plane
      !   DO i = 2, nr-1
      !      DO j = 2, nz-1
      !         r_start(n1)     = raxis(i)
      !         z_start(n1)     = zaxis(j)
      !         phi_start(n1)   = phiaxis(1)
      !         phi_end(n1)     = phiaxis(1)+10*pi2
      !         n1 = n1 + 1
      !      END DO
      !   END DO
         nlines = n1 - 1
         WRITE(6,'(A)') '===========PHI = 0 Grid=========='
      END IF
      dr = (raxis(nr)-raxis(1))/DBLE(nr-1)
      dz = (zaxis(nz)-zaxis(1))/DBLE(nz-1)
      npoinc = nphi - 1 ! so that steps are in phiaxis
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(nlines,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, nprocs_fieldlines, ierr_mpi )          ! MPI
#endif
      CALL fieldlines_follow  ! This call on subgrid grid

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_gridgen',ierr_mpi)
#endif

      ! Now find lines still in the grid
      mystart = LBOUND(R_lines,1)
      myend   = UBOUND(R_lines,1)
      l1      = LBOUND(R_lines,2)
      l2      = UBOUND(R_lines,2)
      DO i = mystart,myend
         IF ((R_lines(i,l2-1) < raxis(nr)) .AND. (R_lines(i,l2-1) > raxis(1)) .AND. &
            (Z_lines(i,l2-1) < zaxis(nz)) .AND. (Z_lines(i,l2-1) > zaxis(1))) THEN
               DO l = l1,l2-1
                  !i2 = COUNT(raxis<=R_lines(i,l))
                  j2 = MOD(l,npoinc)+1
                  !k2 = COUNT(zaxis<=Z_lines(i,l))
                  i2 = COUNT((raxis-0.5*dr)   <= R_lines(i,l))
                  !j = COUNT((phiaxis-0.5*dp) <= zeta_1D(s))
                  k2 = COUNT((zaxis-0.5*dz)   <= Z_lines(i,l))
                  i2 = MIN(MAX(i2,1),nr)
                  j2 = MIN(MAX(j2,1),nphi)
                  k2 = MIN(MAX(k2,1),nz)
                  lgoodline(i2,j2,k2) = .TRUE.
               END DO
         END IF
      END DO
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_gridgen',ierr_mpi)
      ! Have master nodes reduce the lgood variable
      i = MPI_UNDEFINED
      IF (myid_sharmem == master) i = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_FIELDLINES,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
      IF (myid_sharmem == master) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,lgoodline,nr*nphi*nz,MPI_LOGICAL,MPI_LOR,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_gridgen',ierr_mpi)
#endif
      IF (myworkid == master) THEN
         WRITE(6,'(A,I8,A,I8)') '  PHI=0 NGOOD/NTOTAL = ',COUNT(lgoodline(2:nr,1,2:nz)),'/',(nr-2)*(nz-2)
         WRITE(6,'(A,I8,A,I8)') '        NGOOD/NTOTAL = ',COUNT(lgoodline(2:nr,:,2:nz)),'/',(nr-2)*(nz-2)*nphi
      END IF

      ! Now do a new search for missing lines
      IF (myworkid == master) THEN
         WRITE(6,'(A)') '=========== Creating Grid =========='
         ! Axis
         r_start(1) = ra
         z_start(1) = za
         phi_start(1) = phiaxis(1)
         n1 = 2
         ! PHI = 0
         DO l = 1, nr*nz
            i = MOD(l-1,nr)+1
            j = MOD(l-1,nr*nz)
            j = FLOOR(REAL(j) / REAL(nr))+1
            IF (.not.lgoodline(i,1,j)) CYCLE
            r_start(n1)     = raxis(i)
            z_start(n1)     = zaxis(j)
            phi_start(n1)   = phiaxis(1)
            n1 = n1 + 1
         END DO
         ! Now find missing points
         DO l = 1, nr*nphi*nz
            i = MOD(l-1,nr)+1
            j = MOD(l-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(l) / REAL(nr*nphi))
            i = MAX(MIN(i,nr-1),2)
            k = MAX(MIN(k,nz-1),2)
            IF (lgoodline(i,j,k)) CYCLE
            i2 = 0
            IF (lgoodline(i-1, j, k-1)) i2 = i2+1
            IF (lgoodline(i-1, j, k  )) i2 = i2+1
            IF (lgoodline(i-1, j, k+1)) i2 = i2+1
            IF (lgoodline(i  , j, k-1)) i2 = i2+1
            IF (lgoodline(i  , j, k+1)) i2 = i2+1
            IF (lgoodline(i+1, j, k-1)) i2 = i2+1
            IF (lgoodline(i+1, j, k  )) i2 = i2+1
            IF (lgoodline(i+1, j, k+1)) i2 = i2+1
            IF (i2 > 5) THEN
               r_start(n1)     = raxis(i)
               phi_start(n1)   = phiaxis(j)
               z_start(n1)     = zaxis(k)
               n1 = n1 + 1
            END IF
         END DO
         nlines = n1 - 1
      END IF
      phi_end(1:nlines) = pi2*100.0
      npoinc = nphi - 1 ! so that steps are in phiaxis
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(nlines,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_gridgen',ierr_mpi)
      CALL MPI_COMM_SIZE( MPI_COMM_FIELDLINES, nprocs_fieldlines, ierr_mpi )          ! MPI
      IF (ASSOCIATED(lgoodline))   CALL mpidealloc(lgoodline,   win_lgoodline)
#endif
      CALL fieldlines_follow  ! This call on subgrid grid
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_gridgen
