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

!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! For splines
      INTEGER :: i, j, k, s, l, mystart, myend, jm1, jp1
      REAL(rprec) :: br, dBRdR, dBPdP, dBZdZ
      REAL(rprec) :: dPOTdR2, dPOTdP2, dPOTdZ2, dPOTdR, rinv, omega
      REAL(rprec) :: rt,pt,zt,fact
      REAL(rprec), DIMENSION(nr) :: drinv
      REAL(rprec), DIMENSION(nphi) :: dpinv
      REAL(rprec), DIMENSION(nz) :: dzinv
      REAL(rprec), DIMENSION(nr,nphi,nz) :: DIVB, POT, POT2, DIVB2

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
#if defined(MPI_OPT)
      !CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      !CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      !CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      !mylocalmaster = master

      ! Break up the Work
      !CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nphi*nz, mystart, myend)

      ! Setup ICT for values and derivatives
      !ict = (/1,1,1,1,0,0,0,0/)
      IF (lverb) THEN
         WRITE(6,'(A)')   '----- Fixing DIV(B) -----'
         !WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'DIV(B) Correction [',0,']%'
         CALL FLUSH(6)
      ENDIF

      mystart = 1
      myend = nr*nphi*nz

      DIVB(:,:,:) = 0
      FORALL(i=2:nr-1) drinv(i) = 1.0_rprec/(raxis(i+1)-raxis(i-1))
      FORALL(j=2:nphi-1) dpinv(j) = 1.0_rprec/(phiaxis(j+1)-phiaxis(j-1))
      FORALL(k=2:nz-1) dzinv(k) = 1.0_rprec/(zaxis(k+1)-zaxis(k-1))
      dpinv(1) = 1.0_rprec/(phiaxis(2)+phiaxis(2))
      dpinv(nphi) = dpinv(1)

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
         dBRdR = (  B_R(i+1, j  , k  ) -   B_R(i-1, j  , k  ))*drinv(i)
         dBPdP = (B_PHI(i  , jp1, k  ) - B_PHI(i  , jm1, k  ))*dpinv(j)
         dBZdZ = (  B_Z(i  , j  , k+1) -   B_Z(i  , j  , k-1))*dzinv(k)
         ! Divergence in Cyl div(B) = (1/R)*d(RB_R)/dR + (1/R)*d(B_PHI)/dPHI + d(B_Z)/dZ
         DIVB(i,j,k)  = dBRdR + (br+dBPdP)/raxis(i) + dBZdZ
         ! Screen output
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            END IF
         END IF
         CALL FLUSH(6)
      END DO

      ! IF (lverb) THEN
      !    WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Thresholding DIVB [',0,']%'
      !    CALL FLUSH(6)
      ! END IF

      ! ! Deal with coil singularities
      ! DO s = mystart, myend
      !    i = MOD(s-1,nr)+1
      !    j = MOD(s-1,nr*nphi)
      !    j = FLOOR(REAL(j) / REAL(nr))+1
      !    k = CEILING(REAL(s) / REAL(nr*nphi))
      !    ! CYCLE if boundary is hit
      !    IF ((i==1) .or. (i==nr)) CYCLE
      !    IF ((k==1) .or. (k==nz)) CYCLE
      !    IF ((j==nphi)) CYCLE
      !    ! This is done to deal with peridic boundary conditions
      !    jm1 = j-1
      !    jp1 = j+1
      !    IF (j==1) jm1 = nphi-1
      !    br    = SQRT(B_R(i,j,k)*B_R(i,j,k) + B_PHI(i,j,k)*B_PHI(i,j,k) + B_Z(i,j,k)*B_Z(i,j,k))
      !    IF (br > 10.0) DIVB(i,j,k) = 0.0_rprec
      !    ! Screen output
      !    IF (MOD(s,nr) == 0) THEN
      !       IF (lverb) THEN
      !          CALL backspace_out(6,6)
      !          WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
      !       END IF
      !    END IF
      !    CALL FLUSH(6)
      ! END DO

      DIVB(:,nphi,:) = DIVB(:,1,:)

      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Solving Poissons Equation [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Compute the potential
      POT  = 1.0_rprec
      POT2 = 1.0_rprec
      omega = 0.2_rprec
      DO l = 1, 100
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
            !IF (j==nphi) jp1 = 2
            ! This only works for equidistant grids
            rinv   = 1.0_rprec / raxis(i)
            dPOTdR2 = 4.0_rprec*(POT(i+1, j  , k  ) + POT(i-1, j  , k  ))*drinv(i)*drinv(i)
            dPOTdP2 = 4.0_rprec*(POT(i  , jp1, k  ) + POT(i  , jm1, k  ))*dpinv(j)*dpinv(j)*rinv*rinv
            dPOTdZ2 = 4.0_rprec*(POT(i  , j  , k+1) + POT(i  , j  , k-1))*dzinv(k)*dzinv(k)
            dPOTdR  = (POT(i+1, j  , k  ) - POT(i-1, j  , k  ))*drinv(i)
            rt      = 2.0_rprec/drinv(i)
            pt      = 2.0_rprec/dpinv(j)
            zt      = 2.0_rprec/dzinv(k)
            fact    = raxis(i)*raxis(i)*rt*rt*pt*pt*zt*zt/(2.0*raxis(i)*raxis(i)*pt*pt*zt*zt &
               + 2.0*rt*rt*zt*zt + 2.0*rt*rt*pt*pt)
            POT2(i,j,k) = (dPOTdR2 + dPOTdP2 + dPOTdZ2 + dPOTdR - DIVB(i,j,k))*fact
            !POT2(i,j,k) = dPOTdR2 + dPOTdZ2 + (dPOTdP2*rinv + dPOTdR)*rinv - DIVB(i,j,k)
            !WRITE(328,*) l,i,j,k,dPOTdR2,dPOTdP2,dPOTdZ2,dPOTdR,POT2(i,j,k)
         END DO
         ! Apply BC to POT2
         POT2(1,:,:) = POT2(2,:,:)
         POT2(nr,:,:) = POT2(nr-1,:,:)
         POT2(:,:,1) = POT2(:,:,2)
         POT2(:,:,nz) = POT2(:,:,nz-1)
         POT2(:,nphi,:) = POT2(:,1,:)
         POT = (1.0_rprec - omega)*POT + omega*POT2
         PRINT *,l,POT(:,1,128)
         ! Screen output
         IF (lverb) THEN
            CALL backspace_out(6,6)
            WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*l)/(100-1+1)),']%'
         END IF
         CALL FLUSH(6)
      END DO

      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Correcting DIVB [',0,']%'
         CALL FLUSH(6)
      END IF

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
         B_R(i,j,k)   = B_R(i,j,k)   - (POT(i+1, j  , k  ) - POT(i-1, j  , k  ))*drinv(i)
         B_PHI(i,j,k) = B_PHI(i,j,k) - (POT(i  , jp1, k  ) - POT(i  , jm1, k  ))*dpinv(j)/raxis(i)
         B_Z(i,j,k)   = B_Z(i,j,k)   - (POT(i  , j  , k+1) - POT(i  , j  , k-1))*dzinv(k)
         ! Screen output
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            END IF
         END IF
         CALL FLUSH(6)
      END DO

      B_R(:,nphi,:)   = B_R(:,1,:)
      B_PHI(:,nphi,:) = B_PHI(:,1,:)
      B_Z(:,nphi,:)   = B_Z(:,1,:)

      ! Compute DIVB2
      DIVB2(:,:,:) = 0

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
         dBRdR = (  B_R(i+1, j  , k  ) -   B_R(i-1, j  , k  ))*drinv(i)
         dBPdP = (B_PHI(i  , jp1, k  ) - B_PHI(i  , jm1, k  ))*dpinv(j)
         dBZdZ = (  B_Z(i  , j  , k+1) -   B_Z(i  , j  , k-1))*dzinv(k)
         ! Divergence in Cyl div(B) = (1/R)*d(RB_R)/dR + (1/R)*d(B_PHI)/dPHI + d(B_Z)/dZ
         DIVB2(i,j,k)  = dBRdR + (br+dBPdP)/raxis(i) + dBZdZ
         ! Screen output
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            END IF
         END IF
         CALL FLUSH(6)
      END DO
      DIVB2(:,nphi,:) = DIVB2(:,1,:)

      ! Final
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         WRITE(327,*) s, i, j, k, DIVB(i,j,k), DIVB2(i,j,k)
      END DO



      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_divb