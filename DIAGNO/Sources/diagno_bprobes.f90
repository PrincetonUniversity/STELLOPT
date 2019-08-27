!-----------------------------------------------------------------------
!     Module:        diagno_bprobes
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/10/2012
!     Description:   This subroutine calculates the response of a b-dot
!                    probe.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_bprobes
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime, pi2_diag => pi2
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
      INTEGER :: ier, iunit, ncoils,i,ig, ncg, j, k
      REAL(rprec) :: nx, ny, nz, bxp, byp, bzp
      REAL(rprec) :: xvec(3),bvec(3)
      REAL(rprec), ALLOCATABLE :: xp(:), yp(:), rp(:), phip(:), zp(:), bx(:), by(:),&
                     br(:), bphi(:), bz(:), modb(:), th_inc(:),&
                     phi_inc(:), eff_area(:), bn(:), flux(:)
      REAL(rprec), ALLOCATABLE :: bx_mut(:,:), by_mut(:,:), bz_mut(:,:), br_mut(:,:),&
                     bphi_mut(:,:), modb_mut(:,:), bn_mut(:,:)
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
      
      ! Read in diagnostic file
      if(lverb) write(6,*)' --Calculating Magnetic Probe Values'
      IF(lverb .and. lrphiz) THEN
         WRITE(6,'(5X,A,7(9X,A,8X))') 'No','x','y','z','Br','Bphi','Bz','|B|'
      ELSE IF (lverb) THEN
         WRITE(6,'(5X,A,7(9X,A,8X))') 'No','x','y','z','Bx','By','Bz','|B|'
      END IF
      IF (myworkid == master) THEN
         iunit = 29
         CALL safe_open(iunit,ier,TRIM(bprobes_file),'old','formatted')
         READ(iunit,*) ncoils
         ALLOCATE(xp(ncoils),yp(ncoils),zp(ncoils),rp(ncoils),phip(ncoils))
         ALLOCATE(bx(ncoils),by(ncoils),bz(ncoils),br(ncoils),bphi(ncoils),modb(ncoils))
         ALLOCATE(th_inc(ncoils), phi_inc(ncoils), eff_area(ncoils))
         ALLOCATE(flux(ncoils))
         IF (lmut) THEN
            ALLOCATE(bx_mut(ncoils,ncg), by_mut(ncoils,ncg), bz_mut(ncoils,ncg))
            !ALLOCATE(bx_ind(ncoils,ncg),by_ind(ncoils,ncg),&
            !         bz_ind(ncoils,ncg),br_ind(ncoils,ncg),&
            !         bphi_ind(ncoils,ncg),modb_ind(ncoils,ncg),&
            !         bn_ind(ncoils,ncg))
         END IF
         DO i = 1, ncoils
            IF (lrphiz) THEN
               READ(iunit,*) rp(i), phip(i), zp(i), th_inc(i), phi_inc(i), eff_area(i)
               xp(i) = rp(i) * cos(phip(i))
               yp(i) = rp(i) * sin(phip(i))
            ELSE
               READ(iunit,*) xp(i), yp(i) ,zp(i), th_inc(i), phi_inc(i), eff_area(i)
               rp(i) = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
               IF (xp(i) == 0.0 .and. yp(i) == 0.0) THEN
                  phip(i) = 0.0
               ELSE 
                  phip(i) = ATAN2(yp(i),xp(i))
                  IF (phip(i) < 0) phip(i) = phip(i) + pi2
               END IF
            END IF
         END DO
         CLOSE(iunit)
         th_inc = th_inc * onerad
         phi_inc = phi_inc * onerad
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_probes1',ierr_mpi)
      CALL MPI_BCAST(ncoils,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes1',ierr_mpi)
#endif
      IF (myworkid /= master) ALLOCATE(xp(ncoils), yp(ncoils), rp(ncoils), phip(ncoils), zp(ncoils), bx(ncoils), &
                                       by(ncoils), br(ncoils), bphi(ncoils), bz(ncoils), modb(ncoils),&
                                       th_inc(ncoils), phi_inc(ncoils), eff_area(ncoils),flux(ncoils))
      ALLOCATE(bx_mut(ncoils,ncg), by_mut(ncoils,ncg), bz_mut(ncoils,ncg),br_mut(ncoils,ncg),bphi_mut(ncoils,ncg),&
               bn_mut(ncoils,ncg), modb_mut(ncoils,ncg))

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(xp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(yp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(zp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(rp,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(phip,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(th_inc,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(phi_inc,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
      CALL MPI_BCAST(eff_area,ncoils,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_probes2',ierr_mpi)
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
      
      ! Read the mutual induction file if present
      iunit = 36
      IF (luse_mut) THEN
         CALL safe_open(iunit,ier,TRIM(bprobes_mut_file),'old','formatted')
         READ(iunit,'(2I4)') ncoils, nextcur
         IF (ALLOCATED(bx_mut)) DEALLOCATE(bx_mut)
         IF (ALLOCATED(by_mut)) DEALLOCATE(by_mut)
         IF (ALLOCATED(bz_mut)) DEALLOCATE(bz_mut)
         ALLOCATE(bx_mut(ncoils,nextcur),by_mut(ncoils,nextcur),bz_mut(ncoils,nextcur), STAT=ier)
         bx_mut = 0; by_mut = 0; bz_mut = 0;
         DO i = 1, ncoils
            DO ig = 1, nextcur
               READ(iunit,'(2(I6,1X),3ES22.12E3)') j,k, bx_mut(i,ig), by_mut(i,ig), bz_mut(i,ig)
            END DO
            IF (i<mystart .or. i>myend) THEN
               bx_mut(i,:) = 0
               by_mut(i,:) = 0
               bz_mut(i,:) = 0
            END IF
         END DO
         CLOSE(iunit)
         DO ig = 1, nextcur
            IF (luse_extcur(ig)) THEN
               bx_mut(:,ig) = bx_mut(:,ig) * extcur(ig)
               by_mut(:,ig) = by_mut(:,ig) * extcur(ig)
               bz_mut(:,ig) = bz_mut(:,ig) * extcur(ig)
            ELSE
               bx_mut(:,ig) = 0
               by_mut(:,ig) = 0
               bz_mut(:,ig) = 0
            END IF
         END DO
         ncg = nextcur
      END IF

      ! Do Work
      bx = 0; by = 0; bz = 0;
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
         IF (lvac) CYCLE
         ier = 0
         CALL bfield_vc(xp(i),yp(i),zp(i),bx(i),by(i),bz(i),ier)
      END DO


      ! Now handle the arrays
      IF (lmut) THEN
#if defined(MPI_OPT)
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,bx_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,by_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,bz_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(bx_mut,bx_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(by_mut,by_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
            CALL MPI_REDUCE(bz_mut,bz_mut,ncoils*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         END IF
#endif
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

      
      ! Write files
      IF (myworkid == master) THEN
         iunit = 30
         IF (lmut) THEN
            CALL safe_open(iunit,ier,TRIM(bprobes_mut_file),'replace','formatted')
            WRITE(iunit,'(2(I4,1X))') ncoils,ncg
            DO i = 1, ncoils
               DO ig = 1, SIZE(coil_group)
                  WRITE(iunit,'(2(I6,1X),3ES22.12E3)') i,ig, bx_mut(i,ig), by_mut(i,ig), bz_mut(i,ig)
               END DO
            END DO
            WRITE(iunit,'(A)')'#  Coilgroup  BX  BY  BZ'
         ELSE
            modb = sqrt(bx*bx+by*by+bz*bz)
            IF(lverb .and. lrphiz) THEN
               br   = bx * cos(phip) + by * sin(phip)
               bphi = by * cos(phip) - bx * sin(phip)
               DO i = 1, ncoils
                  WRITE(6,'(i6,1p,7e18.8)') i, xp(i), yp(i), zp(i), br(i), bphi(i), bz(i), modb(i)
               END DO
            ELSE IF (lverb) THEN
               DO i = 1, ncoils
                  WRITE(6,'(i6,1p,7e18.8)') i, xp(i), yp(i), zp(i), bx(i), by(i), bz(i), modb(i)
               END DO
            END IF
            CALL FLUSH(6)
            ! Adjust for diagnostic
            bx = bx * sin(phi_inc)*cos(th_inc)
            by = by * sin(phi_inc)*sin(th_inc)
            bz = bz * cos(phi_inc)
            flux = bx + by + bz
            flux = eff_area * flux * bprobe_turns(1:ncoils)
            CALL safe_open(iunit,ier,'diagno_bth.'//TRIM(id_string),'replace','formatted')
            WRITE(iunit,'(i6,1p,5ES22.12E3)')(i,xp(i),yp(i),zp(i),modb(i),flux(i),i=1,ncoils)
            WRITE(iunit,'(A)')'   #   xp[m]      yp      zp      |B| [T]      flux [Wb]'
         END IF
         CLOSE(iunit)
      END IF
      
      ! Cleanup
      DEALLOCATE(xp,yp,zp,rp,phip)
      DEALLOCATE(bx,by,bz,br,bphi,modb)
      DEALLOCATE(th_inc,phi_inc,eff_area)
      DEALLOCATE(flux)
      DEALLOCATE(bx_mut,by_mut,bz_mut,br_mut,&
                  bphi_mut,modb_mut,bn_mut)
      
 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_bprobes
