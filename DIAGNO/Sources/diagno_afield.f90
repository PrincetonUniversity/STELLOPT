!-----------------------------------------------------------------------
!     Module:        diagno_afield
!     Authors:       J. Schilling (jonathan.schilling@mail.de),
!                    S. Lazerson (lazerson@pppl.gov)
!     Date:          05/12/2024
!     Description:   This subroutine calculates the vector potential
!                    at a point in space.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_afield
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime
      USE virtual_casing_mod
      USE biotsavart
      USE safe_open_mod
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
#endif
      INTEGER(KIND=BYTE_8) :: icount, chunk
      INTEGER :: ier, iunit, ncoils, i, ig, iunit_out, ncg
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: xp, yp, rp, phip, zp, ax, ay, ar, aphi, az, moda,&
                     axp, ayp, azp
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: ax_mut, ay_mut, az_mut
      REAL(rprec) :: xvec(3),avec(3)

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
      if(lverb) write(6,*)' ---Calculating Vector Potential Test Points'
      IF(lverb .and. lrphiz) THEN
         WRITE(6,'(14X,A,7(6X,A,6X))') 'No','r','phi','z','Ar','Aphi','Az','|A|'
      ELSE IF (lverb) THEN
         WRITE(6,'(14X,A,7(6X,A,6X))') 'No','x','y','z','Ax','Ay','Az','|A|'
      END IF
      iunit = 26; iunit_out = 27;
      IF (myworkid == master) THEN
         CALL safe_open(iunit,ier,TRIM(afield_points_file),'old','formatted')
         READ(iunit,*) ncoils
         ALLOCATE(xp(ncoils), yp(ncoils), rp(ncoils), phip(ncoils), zp(ncoils), &
                  ax(ncoils), ay(ncoils), ar(ncoils), aphi(ncoils), az(ncoils), &
                  moda(ncoils), axp(ncoils), ayp(ncoils), azp(ncoils))
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
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_afield1',ierr_mpi)
      CALL MPI_BCAST(ncoils,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_afield1',ierr_mpi)
#endif
      IF (myworkid /= master) ALLOCATE(xp(ncoils), yp(ncoils), rp(ncoils), phip(ncoils), zp(ncoils), &
                                       ax(ncoils), ay(ncoils), ar(ncoils), aphi(ncoils), az(ncoils), &
                                       moda(ncoils), axp(ncoils), ayp(ncoils), azp(ncoils))
      ALLOCATE(ax_mut(ncoils,ncg), ay_mut(ncoils,ncg), az_mut(ncoils,ncg))

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_afield2',ierr_mpi)
      CALL MPI_BCAST(xp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_afield2',ierr_mpi)
      CALL MPI_BCAST(yp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_afield2',ierr_mpi)
      CALL MPI_BCAST(zp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_afield2',ierr_mpi)
      CALL MPI_BCAST(rp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_afield2',ierr_mpi)
      CALL MPI_BCAST(phip,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_afield2',ierr_mpi)
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
      ax = 0; ay = 0; az = 0;
      ax_mut = 0; ay_mut = 0; az_mut = 0;
      DO i = mystart, myend
         xvec(1) = xp(i); xvec(2)=yp(i); xvec(3)=zp(i)
         IF (lcoil) THEN
            DO ig = 1, ncg
               avec=0
               CALL bsc_a(coil_group(ig),xvec,avec)
               ax_mut(i,ig) = avec(1)
               ay_mut(i,ig) = avec(2)
               az_mut(i,ig) = avec(3)
            END DO
         END IF
         IF (lmut .or. lvac) CYCLE
         CALL afield_vc(xp(i),yp(i),zp(i),ax(i),ay(i),az(i),ier)
      END DO

      ! Now handle the arrays
      IF (lmut) THEN
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,ax_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,ay_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,az_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(ax_mut,ax_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(ay_mut,ay_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(az_mut,az_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         END IF
      ELSE
         ! Plasma + vacuum
         ax = ax + SUM(ax_mut, DIM=2)
         ay = ay + SUM(ay_mut, DIM=2)
         az = az + SUM(az_mut, DIM=2)
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_afield3',ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        ax,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        ay,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        az,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
#endif
      END IF

      ! PRINT TO SCREEN
      IF (myworkid == master) THEN
         IF (lmut) THEN
         ELSE
            CALL safe_open(iunit_out,ier,'diagno_atest.'//TRIM(id_string),'replace','formatted')
            moda = sqrt(ax*ax+ay*ay+az*az)
            ar   = ax * cos(phip) + ay * sin(phip)
            aphi = ay * cos(phip) - ax * sin(phip)
            IF (lrphiz) THEN
               DO i = 1, ncoils
                  IF (lverb) WRITE(6,'(13X,I8,1X,7E14.5)') i,rp(i),phip(i),zp(i),ar(i),aphi(i),az(i),moda(i)
                  WRITE(iunit_out,'(13X,I8,1X,7E14.5)') i,rp(i),phip(i),zp(i),ar(i),aphi(i),az(i),moda(i)
               END DO
               WRITE(iunit_out,'(A)')'   #   rp[m]    phip[deg]    zp[m]      A_R[Tm]    A_PHI[Tm]      A_Z[Tm]      |A|[Tm]'
            ELSE
               DO i = 1, ncoils
                  IF (lverb) WRITE(6,'(13X,I8,1X,7E14.5)') i,xp(i),yp(i),zp(i),ax(i),ay(i),az(i),moda(i)
                  WRITE(iunit_out,'(13X,I8,1X,7E14.5)') i,xp(i),yp(i),zp(i),ax(i),ay(i),az(i),moda(i)
               END DO
               WRITE(iunit_out,'(A)')'   #   xp[m]      yp[m]      zp[m]      A_X[Tm]      A_Y[Tm]      A_Z[Tm]      |A|[Tm]'
            END IF
            CLOSE(iunit_out)
         END IF
      END IF

      ! Clean up
      DEALLOCATE(xp, yp, rp, phip, zp, ax, ay, ar, aphi, az, moda, axp, ayp, azp)
      DEALLOCATE(ax_mut,ay_mut,az_mut)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE diagno_afield
