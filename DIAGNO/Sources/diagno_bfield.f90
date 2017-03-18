!-----------------------------------------------------------------------
!     Module:        diagno_bfield
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the field at a point in
!                    space.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_bfield
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime
      USE virtual_casing_mod
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
      INTEGER :: ier, iunit, ncoils, i, ig, iunit_out, ncg
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: xp, yp, rp, phip, zp, bx, by, br, bphi, bz, modb,&
                     bxp, byp, bzp
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: bx_mut, by_mut, bz_mut
      REAL(rprec) :: xvec(3),bvec(3)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Basic copy of MPI_COMM_DIANGO
      mylocalid = myworkid
      mylocalmaster = master
      MPI_COMM_LOCAL = MPI_COMM_DIAGNO
      numprocs_local = nprocs_diagno
      ncg = SIZE(coil_group)

      ! Read File
      if(lverb) write(6,*)' ---Calculating Magnetic Field Test Points'
      IF(lverb .and. lrphiz) THEN
         WRITE(6,'(14X,A,7(6X,A,6X))') 'No','r','phi','z','Br','Bphi','Bz','|B|'
      ELSE IF (lverb) THEN
         WRITE(6,'(14X,A,7(6X,A,6X))') 'No','x','y','z','Bx','By','Bz','|B|'
      END IF
      iunit = 26; iunit_out = 27;
      IF (myworkid == master) THEN
         CALL safe_open(iunit,ier,TRIM(bfield_points_file),'old','formatted')
         READ(iunit,*) ncoils
         ALLOCATE(xp(ncoils), yp(ncoils), rp(ncoils), phip(ncoils), zp(ncoils), bx(ncoils), &
                  by(ncoils), br(ncoils), bphi(ncoils), bz(ncoils), modb(ncoils),&
                  bxp(ncoils), byp(ncoils), bzp(ncoils))
         DO i = 1, ncoils
            IF (lrphiz) THEN
               READ(iunit,*) rp(i),phip(i),zp(i)
               xp(i) = rp(i) * cos(phip(i))
               yp(i) = rp(i) * sin(phip(i))
            ELSE
               READ(iunit,*) xp(i),yp(i),zp(i)
            END IF
         END DO
         CLOSE(iunit)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_bfield1',ierr_mpi)
      CALL MPI_BCAST(ncoils,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_bfield1',ierr_mpi)
#endif
      IF (myworkid /= master) ALLOCATE(xp(ncoils), yp(ncoils), rp(ncoils), phip(ncoils), zp(ncoils), bx(ncoils), &
                                       by(ncoils), br(ncoils), bphi(ncoils), bz(ncoils), modb(ncoils),&
                                       bxp(ncoils), byp(ncoils), bzp(ncoils))
      ALLOCATE(bx_mut(ncoils,ncg), by_mut(ncoils,ncg), bz_mut(ncoils,ncg))

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_bfield2',ierr_mpi)
      CALL MPI_BCAST(xp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_bfield2',ierr_mpi)
      CALL MPI_BCAST(yp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_bfield2',ierr_mpi)
      CALL MPI_BCAST(zp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_bfield2',ierr_mpi)
      CALL MPI_BCAST(rp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_bfield2',ierr_mpi)
      CALL MPI_BCAST(phip,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_bfield2',ierr_mpi)
#endif

      ! Divide up the work
      chunk = FLOOR(REAL(ncoils) / REAL(numprocs_local))
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
         IF ((moffsets(numprocs_local)+mnum(numprocs_local)-1) == ncoils) EXIT
         IF (i == numprocs_local) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs_local) = moffsets(i+1:numprocs_local) + 1
         i=i+1
      END DO
      mystart = moffsets(mylocalid+1)
      chunk  = mnum(mylocalid+1)
      myend   = mystart + chunk - 1
#endif

      ! Do Work
      bx = 0; by = 0; bz = 0;
      bx_mut = 0; by_mut = 0; bz_mut = 0;
      DO i = mystart, myend
         xvec(1) = xp(i); xvec(2)=yp(i); xvec(3)=zp(i)
         IF (lcoil) THEN
            DO ig = 1, ncg
               bvec=0
               CALL bsc_b(coil_group(ig),xvec,bvec)
               bx_mut(i,ig) = bvec(1)
               by_mut(i,ig) = bvec(2)
               bz_mut(i,ig) = bvec(3)
            END DO 
         END IF
         IF (lmut .or. lvac) CYCLE
         CALL bfield_vc(xp(i),yp(i),zp(i),bx(i),by(i),bz(i),ier)
      END DO

      ! Now handle the arrays
      IF (lmut) THEN
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,bx_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,by_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,bz_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(bx_mut,bx_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(by_mut,by_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(bz_mut,bz_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         END IF
      ELSE
         ! Plasma + vacuum
         bx = bx + SUM(bx_mut, DIM=2)
         by = by + SUM(by_mut, DIM=2)
         bz = bz + SUM(bz_mut, DIM=2)
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_bfield3',ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        bx,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        by,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        bz,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
#endif
      END IF
      
      ! PRINT TO SCREEN
      IF (myworkid == master) THEN
         IF (lmut) THEN
         ELSE
            CALL safe_open(iunit_out,ier,'diagno_btest.'//TRIM(id_string),'replace','formatted')
            modb = sqrt(bx*bx+by*by+bz*bz)
            br   = bx * cos(phip) + by * sin(phip)
            bphi = by * cos(phip) - bx * sin(phip)
            IF (lrphiz) THEN
               DO i = 1, ncoils
                  IF (lverb) WRITE(6,'(13X,I8,1X,7E14.5)') i,rp(i),phip(i),zp(i),br(i),bphi(i),bz(i),modb(i)
                  WRITE(iunit_out,'(13X,I8,1X,7E14.5)') i,rp(i),phip(i),zp(i),br(i),bphi(i),bz(i),modb(i)
               END DO
               WRITE(iunit_out,'(A)')'   #   rp[m]    phip[deg]    zp[m]      B_R[T]    B_PHI[T]      B_Z[T]      |B|[T]'
            ELSE
               DO i = 1, ncoils
                  IF (lverb) WRITE(6,'(13X,I8,1X,7E14.5)') i,xp(i),yp(i),zp(i),bx(i),by(i),bz(i),modb(i)
                  WRITE(iunit_out,'(13X,I8,1X,7E14.5)') i,xp(i),yp(i),zp(i),bx(i),by(i),bz(i),modb(i)
               END DO
               WRITE(iunit_out,'(A)')'   #   xp[m]      yp[m]      zp[m]      B_X[T]      B_Y[T]      B_Z[T]      |B|[T]'
            END IF
            CLOSE(iunit_out)
         END IF
      END IF

      ! Clean up
      DEALLOCATE(xp, yp, rp, phip, zp, bx, by, br, bphi, bz, modb,bxp, byp, bzp)
      DEALLOCATE(bx_mut,by_mut,bz_mut) 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_bfield
