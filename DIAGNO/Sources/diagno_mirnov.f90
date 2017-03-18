!-----------------------------------------------------------------------
!     Module:        diagno_mirnov
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/16/2012
!     Description:   This subroutine calculates the response of the
!                    Mirnov arrays.  Here the array is assumed to be
!                    calculated via:
!                    signal[T]=\frac{\int\vec{b}\cdot d\vec{l}}{\int |d\vec{l}|}
!                    which is just the average magnetic field integrated
!                    along each array.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_mirnov
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime, pi2_diag => pi2
      USE biotsavart
      USE safe_open_mod
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
#endif
      INTEGER(KIND=BYTE_8) :: icount, chunk
      INTEGER :: ier, iunit, nseg, npts,i,j,k,ig, ncg
      REAL(rprec) :: nx, ny, nz, int_fac, xp1, yp1, zp1, xp2, yp2, zp2
      REAL(rprec) :: xvec(3),bvec(3)
      REAL(rprec), ALLOCATABLE :: flux(:)
      REAL(rprec), ALLOCATABLE :: xp(:,:), yp(:,:), rp(:,:), phip(:,:),&
                                  zp(:,:), dx(:,:), dy(:,:), dz(:,:), &
                                  xmid(:,:), ymid(:,:), zmid(:,:), &
                                  eff_area(:,:), flux_mut(:,:)
      
!-----------------------------------------------------------------------
!     External Functions
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      INTERFACE
         REAL FUNCTION db_bode(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_bode
         
         REAL FUNCTION db_midpoint(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         real(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_midpoint
      
         REAL FUNCTION db_simpson(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_simpson
      END INTERFACE
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Basic copy of MPI_COMM_DIANGO
      mylocalid = myworkid
      mylocalmaster = master
      MPI_COMM_LOCAL = MPI_COMM_DIAGNO
      numprocs_local = nprocs_diagno
      ncg = SIZE(coil_group)
      IF (ncg == 0) ncg = 1 ! Just so we don't allocate a zero size array
 

      if(lverb) write(6,*)' --Calculating Mirnov Array Values'
      int_fac=1./int_step
      IF (lverb) THEN
         WRITE(6,'(5X,A,2(9X,A,9X))') 'Loop','nseg','flux'
      END IF
      
      ! Read in diagnostic file
      IF (myworkid == master) THEN
         iunit = 36
         CALL safe_open(iunit,ier,TRIM(mirnov_file),'old','formatted')
         READ(iunit,*) nseg,npts
         ALLOCATE(xp(nseg,1:npts),yp(nseg,1:npts),zp(nseg,1:npts),rp(nseg,1:npts),phip(nseg,1:npts))
         ALLOCATE(dx(nseg,1:npts),dy(nseg,1:npts),dz(nseg,1:npts))
         ALLOCATE(xmid(nseg,0:npts),ymid(nseg,0:npts),zmid(nseg,0:npts))
         ALLOCATE(flux(nseg),eff_area(nseg,npts),flux_mut(nseg,ncg))
         DO i = 1, nseg
            DO j = 1, npts
               READ(iunit,*) rp(i,j), zp(i,j), phip(i,j), eff_area(i,j)
            END DO
         END DO
         CLOSE(iunit)
         phip = phip * onerad ! Input in theta
         xp = rp * cos(phip)
         yp = rp * sin(phip)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_mirnov1',ierr_mpi)
      CALL MPI_BCAST(nseg,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov1',ierr_mpi)
      CALL MPI_BCAST(npts,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov1',ierr_mpi)
#endif
      IF (myworkid /= master) THEN
         ALLOCATE(xp(nseg,1:npts),yp(nseg,1:npts),zp(nseg,1:npts),rp(nseg,1:npts),phip(nseg,1:npts))
         ALLOCATE(dx(nseg,1:npts),dy(nseg,1:npts),dz(nseg,1:npts))
         ALLOCATE(xmid(nseg,0:npts),ymid(nseg,0:npts),zmid(nseg,0:npts))
         ALLOCATE(flux(nseg),eff_area(nseg,npts),flux_mut(nseg,ncg))
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_mirnov2',ierr_mpi)
      CALL MPI_BCAST(xp,nseg*npts,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov2',ierr_mpi)
      CALL MPI_BCAST(yp,nseg*npts,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov2',ierr_mpi)
      CALL MPI_BCAST(zp,nseg*npts,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov2',ierr_mpi)
      CALL MPI_BCAST(rp,nseg*npts,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov2',ierr_mpi)
      CALL MPI_BCAST(phip,nseg*npts,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov2',ierr_mpi)
      CALL MPI_BCAST(eff_area,nseg*npts,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_mirnov2',ierr_mpi)
#endif
      
      ! Now create calculate helper arrays
      dx(:,2:npts-1) = xp(:,3:npts) - xp(:,1:npts-2)
      dy(:,2:npts-1) = yp(:,3:npts) - yp(:,1:npts-2)
      dz(:,2:npts-1) = zp(:,3:npts) - zp(:,1:npts-2)
      xmid(:,2:npts-1) = 0.5 * ( xp(:,3:npts) + xp(:,1:npts-2) )
      ymid(:,2:npts-1) = 0.5 * ( yp(:,3:npts) + yp(:,1:npts-2) )
      zmid(:,2:npts-1) = 0.5 * ( zp(:,3:npts) + zp(:,1:npts-2) )

      ! Divide up the work
      chunk = FLOOR(REAL(nseg) / REAL(numprocs_local))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1
#if defined(MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs_local), moffsets(numprocs_local))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(numprocs_local)+mnum(numprocs_local)-1) == nseg) EXIT
         IF (i == numprocs_local) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs_local) = moffsets(i+1:numprocs_local) + 1
         i=i+1
      END DO
      mystart = moffsets(mylocalid+1)
      chunk  = mnum(mylocalid+1)
      myend   = mystart + chunk - 1
#endif
      
      ! Read the mutual induction file if present
      iunit = 36
      IF (luse_mut) THEN
         CALL safe_open(iunit,ier,TRIM(mir_mut_file),'old','formatted')
         READ(iunit,'(2I4)') nseg, nextcur
         IF (ALLOCATED(flux_mut)) DEALLOCATE(flux_mut)
         ALLOCATE(flux_mut(nseg,nextcur), STAT=ier)
         DO i = 1, nseg
            DO ig = 1, nextcur
               READ(iunit,'(1X,I4,1X,I4,E22.14)') j,k,flux_mut(i,ig)
            END DO
            IF (i<mystart .or. i>myend) flux_mut(i,:)=0
         END DO
         CLOSE(iunit)
         DO ig = 1, nextcur
            IF (luse_extcur(ig)) THEN
               flux_mut(:,ig) = flux_mut(:,ig) * extcur(ig)
            ELSE
               flux_mut(:,ig) = 0
            END IF
         END DO
         ncg = nextcur
      END IF

      
      ! Calculate Field
      flux = 0
      DO i = mystart, myend
         SELECT CASE (int_type)
            CASE('midpoint')
               DO j = 2, npts-1
                  DO k = 1, int_step
                     xp1 = xp(i,j) + (k-1)*int_fac*dx(i,j)
                     yp1 = yp(i,j) + (k-1)*int_fac*dy(i,j)
                     zp1 = zp(i,j) + (k-1)*int_fac*dz(i,j)
                     xp2 = xp(i,j) + k*int_fac*dx(i,j)
                     yp2 = yp(i,j) + k*int_fac*dy(i,j)
                     zp2 = zp(i,j) + k*int_fac*dz(i,j)
                     IF (lcoil) THEN
                        DO ig = 1, SIZE(coil_group)
                          IF (.not.luse_extcur(ig)) CYCLE
                          flux_mut(i,ig)=flux_mut(i,ig) + db_midpoint(xp1,yp1,zp1,xp2,yp2,zp2,ig,LDB = .TRUE.)
                        END DO
                     END IF
                     IF (lvac) CYCLE
                     ier = 0
                     flux(i)=flux(i) + db_midpoint(xp1,yp1,zp1,xp2,yp2,zp2,ier,LDB = .TRUE.)
                  END DO
               END DO
            CASE('simpson')
               DO j = 2, npts-1
                  DO k = 1, int_step
                     xp1 = xp(i,j) + (k-1)*int_fac*dx(i,j)
                     yp1 = yp(i,j) + (k-1)*int_fac*dy(i,j)
                     zp1 = zp(i,j) + (k-1)*int_fac*dz(i,j)
                     xp2 = xp(i,j) + k*int_fac*dx(i,j)
                     yp2 = yp(i,j) + k*int_fac*dy(i,j)
                     zp2 = zp(i,j) + k*int_fac*dz(i,j)
                     IF (lcoil) THEN
                        DO ig = 1, SIZE(coil_group)
                           IF (.not.luse_extcur(ig)) CYCLE
                           flux_mut(i,ig)=flux_mut(i,ig) + db_simpson(xp1,yp1,zp1,xp2,yp2,zp2,ig,LDB = .TRUE.)
                        END DO
                     END IF
                     IF (lvac) CYCLE
                     ier = 0
                     flux(i)=flux(i) + db_simpson(xp1,yp1,zp1,xp2,yp2,zp2,ier,LDB = .TRUE.)
                  END DO
               END DO
            CASE('bode')
               DO j = 2, npts-1
                  DO k = 1, int_step
                     xp1 = xp(i,j) + (k-1)*int_fac*dx(i,j)
                     yp1 = yp(i,j) + (k-1)*int_fac*dy(i,j)
                     zp1 = zp(i,j) + (k-1)*int_fac*dz(i,j)
                     xp2 = xp(i,j) + k*int_fac*dx(i,j)
                     yp2 = yp(i,j) + k*int_fac*dy(i,j)
                     zp2 = zp(i,j) + k*int_fac*dz(i,j)
                     IF (lcoil) THEN
                        DO ig = 1, SIZE(coil_group)
                          IF (.not.luse_extcur(ig)) CYCLE
                          flux_mut(i,ig)=flux_mut(i,ig) + db_bode(xp1,yp1,zp1,xp2,yp2,zp2,ig,LDB = .TRUE.)
                        END DO
                     END IF
                     IF (lvac) CYCLE
                     ier = 0
                     flux(i)=flux(i) + db_bode(xp1,yp1,zp1,xp2,yp2,zp2,ier,LDB = .TRUE.)
                  END DO
               END DO
            CASE('adaptive')
              
         END SELECT
            !IF (lverb) THEN
            !   WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, npts, flux(i)
            !END IF
      END DO

      ! Now handle the arrays
      IF (lmut) THEN
#if defined(MPI_OPT)
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,flux_mut,nseg*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(flux_mut,flux_mut,nseg*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         END IF
         flux = SUM(flux_mut,DIM=2)
#endif
      ELSE
         ! Plasma + vacuum
         flux = flux + SUM(flux_mut, DIM=2)
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_rogo3',ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        flux,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
#endif
      END IF
      
      ! Write files
      IF (myworkid == master) THEN
         iunit = 30
         IF (lmut) THEN
            CALL safe_open(iunit,ier,TRIM(rog_mut_file),'replace','formatted')
            WRITE(iunit,'(2(I4,1X))') nseg,ncg
            DO i = 1, nseg
               DO ig = 1, ncg
                  WRITE(iunit,'(2(I6,1X),1ES22.12E3)') i,ig, flux_mut(i,ig)
               END DO
            END DO
         ELSE
            CALL safe_open(iunit,ier,'diagno_mirnov.'//TRIM(id_string),'replace','formatted')
            WRITE(iunit,'(i6,1p,1ES22.12E3)')(i,flux(i),i=1,nseg)
            WRITE(iunit,'(A)')'   #   flux [Wb]'
         END IF
         DO i = 1, nseg
            IF (lverb) WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, npts, flux(i)
         END DO
         CLOSE(iunit)
      END IF
      
      ! Cleanup
      DEALLOCATE(xp,yp,zp,rp,phip)
      DEALLOCATE(dx,dy,dz,xmid,ymid,zmid)
      DEALLOCATE(flux,eff_area)
      
      
 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_mirnov
