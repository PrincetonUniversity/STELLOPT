!-----------------------------------------------------------------------
!     Module:        beams3d_divB
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          07/14/2025
!     Description:   This subroutine cleans divergence B from the
!                    grid.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_divb
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 hr, hp, hz, hri, hpi, hzi, &
                                 phimax, B_R, B_Z, B_PHI, S_ARR
      USE mpi_params
      USE mpi_inc
      USE mpi_sharmem

!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, j, k, s, l, mystart, myend, jm1, jp1
      REAL(rprec) :: br, dBRdR, dBPdP, dBZdZ
      REAL(rprec) :: dPOTdR2, dPOTdP2, dPOTdZ2, dPOTdR
      REAL(rprec) :: dr,dp,dz,fact, scale_divb

      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL

      INTEGER :: win_drinv, win_dpinv, win_dzinv, win_divb, win_pot,&
                 win_pot2, win_divb2, win_rinv
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: rinv
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: drinv
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: dpinv
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: dzinv
      DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: DIVB, POT, POT2, &
                                                    DIVB2

      REAL(rprec), PARAMETER :: omega = 0.50_rprec

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      ! Allocate the helpers
      CALL mpialloc(rinv, nr, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rinv)
      CALL mpialloc(drinv, nr, myid_sharmem, 0, MPI_COMM_SHARMEM, win_drinv)
      CALL mpialloc(dpinv, nphi, myid_sharmem, 0, MPI_COMM_SHARMEM, win_dpinv)
      CALL mpialloc(dzinv, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_dzinv)
      CALL mpialloc(DIVB, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_divb)
      CALL mpialloc(POT, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_pot)
      CALL mpialloc(POT2, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_pot2)
      CALL mpialloc(DIVB2, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_divb2)


      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nphi*nz, mystart, myend)

      ! Setup ICT for values and derivatives
      IF (lverb) THEN
         WRITE(6,'(A)')   '----- Fixing DIV(B) -----'
         CALL FLUSH(6)
      ENDIF

      !mystart = 1
      !myend = nr*nphi*nz

      ! Initialize helpers
      IF (mylocalid == mylocalmaster) THEN
         DIVB(:,:,:) = 0.0_rprec
         DIVB2(:,:,:) = 0.0_rprec
         POT(:,:,:) = 0.0_rprec
         POT2(:,:,:) = 0.0_rprec
         FORALL(i=1:nr) rinv(i) = 1.0_rprec/raxis(i)
         FORALL(i=2:nr-1) drinv(i) = 2.0_rprec/(raxis(i+1)-raxis(i-1))
         FORALL(j=2:nphi-1) dpinv(j) = 2.0_rprec/(phiaxis(j+1)-phiaxis(j-1))
         FORALL(k=2:nz-1) dzinv(k) = 2.0_rprec/(zaxis(k+1)-zaxis(k-1))
         dpinv(1) = 2.0_rprec/(phiaxis(2)+phiaxis(2))
         dpinv(nphi) = dpinv(1)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif

      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Computing DIVB [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Compute divergence of B
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         ! CYCLE if boundary is hit
         IF ((i==1) .or. (i==nr)) CYCLE
         IF ((k==1) .or. (k==nz)) CYCLE
         IF ((j==nphi)) CYCLE
         IF (S_ARR(i,j,k) > 2.0) CYCLE ! This is done to avoid coil regions
         ! This is done to deal with peridic boundary conditions
         jm1 = j-1
         jp1 = j+1
         IF (j==1) jm1 = nphi-1
         ! Manually calcualte DIVB (finite difference)
         br    = B_R(i,j,k)
         dBRdR = 0.5*(  B_R(i+1, j  , k  ) -   B_R(i-1, j  , k  ))*drinv(i)
         dBPdP = 0.5*(B_PHI(i  , jp1, k  ) - B_PHI(i  , jm1, k  ))*dpinv(j)
         dBZdZ = 0.5*(  B_Z(i  , j  , k+1) -   B_Z(i  , j  , k-1))*dzinv(k)
         ! Divergence in Cyl div(B) = (1/R)*d(RB_R)/dR + (1/R)*d(B_PHI)/dPHI + d(B_Z)/dZ
         DIVB(i,j,k)  = dBRdR + (br+dBPdP)*rinv(i) + dBZdZ
         ! Screen output
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            END IF
         END IF
         CALL FLUSH(6)
      END DO

      !scale_divb = 1E3*MAXVAL(DIVB)

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif

      IF (mylocalid == mylocalmaster) THEN
         DIVB(:,nphi,:) = DIVB(:,1,:)
         !DIVB = DIVB/scale_divb
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif

      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Solving Poissons Equation [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Compute the potential
      DO l = 1, 1000
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            ! CYCLE if boundary is hit
            IF ((i==1) .or. (i==nr)) CYCLE
            IF ((k==1) .or. (k==nz)) CYCLE
            IF ((j==nphi)) CYCLE
            jm1 = j-1
            jp1 = j+1
            IF (j==1) jm1 = nphi-1
            ! This only works for equidistant grids
            dPOTdR2 = (POT(i+1, j  , k  ) + POT(i-1, j  , k  ))*drinv(i)*drinv(i)
            dPOTdP2 = (POT(i  , jp1, k  ) + POT(i  , jm1, k  ))*dpinv(j)*dpinv(j)*rinv(i)*rinv(i)
            dPOTdZ2 = (POT(i  , j  , k+1) + POT(i  , j  , k-1))*dzinv(k)*dzinv(k)
            dPOTdR  = 0.5*(POT(i+1, j  , k  ) - POT(i-1, j  , k  ))*drinv(i)*rinv(i)
            dr      = 1.0_rprec/drinv(i)
            dp      = 1.0_rprec/dpinv(j)
            dz      = 1.0_rprec/dzinv(k)
            fact    = raxis(i)*raxis(i)*dr*dr*dp*dp*dz*dz/(2.0*raxis(i)*raxis(i)*dp*dp*dz*dz &
               + 2.0*dr*dr*dz*dz + 2.0*raxis(i)*raxis(i)*dr*dr*dp*dp)
            POT2(i,j,k) = (dPOTdR2 + dPOTdP2 + dPOTdZ2 + dPOTdR - DIVB(i,j,k))*fact
         END DO
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif
         ! Check for convergence
         IF (ALL((ABS(POT-POT2)/ABS(POT2)) < 1.0E-3)) EXIT

         ! Apply BC to POT2
         IF (mylocalid == mylocalmaster) THEN
            POT2(1,:,:)  = POT2(2,:,:)
            POT2(nr,:,:) = POT2(nr-1,:,:)
            POT2(:,:,1)  = POT2(:,:,2)
            POT2(:,:,nz) = POT2(:,:,nz-1)
            POT2(:,nphi,:) = POT2(:,1,:)
            POT = (1.0_rprec - omega)*POT + omega*POT2
            POT2 = POT
            PRINT *,l,POT(128,45,128)
         END IF
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif
         ! Screen output
         IF (lverb) THEN
            CALL backspace_out(6,6)
            WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*l)/(1000-1+1)),']%'
         END IF
         CALL FLUSH(6)
      END DO

      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Correcting DIVB [',0,']%'
         CALL FLUSH(6)
      END IF


      !IF (mylocalid == mylocalmaster) POT = POT * scale_divb

      ! Remove the div(B) component -nabla(POT)
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         ! CYCLE if boundary is hit
         IF ((i==1) .or. (i==nr)) CYCLE
         IF ((k==1) .or. (k==nz)) CYCLE
            IF ((j==nphi)) CYCLE
         ! This is done to deal with peridic boundary conditions
         jm1 = j-1
         jp1 = j+1
         IF (j==1) jm1 = nphi-1
         !IF (j==nphi) jp1 = 2
         ! Manually calcualte DIVB (finite difference)
         B_R(i,j,k)   = B_R(i,j,k)   - 0.5*(POT(i+1, j  , k  ) - POT(i-1, j  , k  ))*drinv(i)
         B_PHI(i,j,k) = B_PHI(i,j,k) - 0.5*(POT(i  , jp1, k  ) - POT(i  , jm1, k  ))*dpinv(j)*rinv(i)
         B_Z(i,j,k)   = B_Z(i,j,k)   - 0.5*(POT(i  , j  , k+1) - POT(i  , j  , k-1))*dzinv(k)
         ! Screen output
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            END IF
         END IF
         CALL FLUSH(6)
      END DO
      
      IF (mylocalid == mylocalmaster) THEN 
         B_R(:,nphi,:)   = B_R(:,1,:)
         B_PHI(:,nphi,:) = B_PHI(:,1,:)
         B_Z(:,nphi,:)   = B_Z(:,1,:)
      END IF

#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif

      ! Compute DIVB2
      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Recomputing DIVB [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Compute divergence of B-final
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         ! CYCLE if boundary is hit
         IF ((i==1) .or. (i==nr)) CYCLE
         IF ((k==1) .or. (k==nz)) CYCLE
         IF ((j==nphi)) CYCLE
         ! This is done to deal with peridic boundary conditions
         jm1 = j-1
         jp1 = j+1
         IF (j==1) jm1 = nphi-1
         !IF (j==nphi) jp1 = 2
         ! Manually calcualte DIVB (finite difference)
         br    = B_R(i,j,k)
         dBRdR = 0.5*(  B_R(i+1, j  , k  ) -   B_R(i-1, j  , k  ))*drinv(i)
         dBPdP = 0.5*(B_PHI(i  , jp1, k  ) - B_PHI(i  , jm1, k  ))*dpinv(j)
         dBZdZ = 0.5*(  B_Z(i  , j  , k+1) -   B_Z(i  , j  , k-1))*dzinv(k)
         ! Divergence in Cyl div(B) = (1/R)*d(RB_R)/dR + (1/R)*d(B_PHI)/dPHI + d(B_Z)/dZ
         DIVB2(i,j,k)  = dBRdR + (br+dBPdP)*rinv(i) + dBZdZ
         ! Screen output
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            END IF
         END IF
         CALL FLUSH(6)
      END DO


      IF (mylocalid == mylocalmaster) DIVB2(:,nphi,:) = DIVB2(:,1,:)

#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif
      ! Final
      IF (mylocalid == mylocalmaster) THEN
         DO s = 1, nr*nphi*nz
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            WRITE(327,*) s, i, j, k, DIVB(i,j,k), DIVB2(i,j,k)
         END DO
      END IF

      IF (ASSOCIATED(rinv)) CALL mpidealloc(rinv,win_rinv)
      IF (ASSOCIATED(drinv)) CALL mpidealloc(drinv,win_drinv)
      IF (ASSOCIATED(dpinv)) CALL mpidealloc(dpinv,win_dpinv)
      IF (ASSOCIATED(dzinv)) CALL mpidealloc(dzinv,win_dzinv)
      IF (ASSOCIATED(POT)) CALL mpidealloc(POT,win_POT)
      IF (ASSOCIATED(POT2)) CALL mpidealloc(POT,win_POT2)
      IF (ASSOCIATED(DIVB)) CALL mpidealloc(DIVB,win_DIVB)
      IF (ASSOCIATED(DIVB2)) CALL mpidealloc(DIVB2,win_DIVB2)

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_divb