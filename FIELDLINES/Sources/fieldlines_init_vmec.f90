!-----------------------------------------------------------------------
!     Subroutine:    fieldlines_init_vmec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/25/2012
!     Description:   This subroutine reads the VMEC wout file and
!                    initializes the plasma fields.  This is achieved
!                    through utilization of a virtual casing principle.
!                    This is stil under development but working.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_vmec
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, extcur_in => extcur, nextcur_vmec => nextcur
      USE vmec_input, ONLY: nzeta_vmec => nzeta
      USE vmec_utils
      USE virtual_casing_mod, pi2_vc => pi2
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, vc_adapt_tol, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl
      USE fieldlines_lines, ONLY: nlines
      USE wall_mod, ONLY: wall_load_mn, wall_info,vertex,face
      USE mpi_params                                                    ! MPI
!      USE mpi
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
!DEC$ IF DEFINED (MPI_OPT)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
!DEC$ ENDIF  
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, s, i, j, k, nu, nv, mystart,myend
      INTEGER, ALLOCATABLE :: xn_temp(:), xm_temp(:)
      REAL :: br_vc, bphi_vc, bz_vc, xaxis_vc, yaxis_vc, zaxis_vc,&
              bx_vc, by_vc
      REAL(rprec) :: br, bphi, bz, sflx
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_temp(:,:),zmns_temp(:,:),&
                           bumnc_temp(:,:),bvmnc_temp(:,:),&
                           rmns_temp(:,:),zmnc_temp(:,:),&
                           bumns_temp(:,:),bvmns_temp(:,:)
      LOGICAL :: lvolint = .false.
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      IF ((nprocs_fieldlines) > nlocal) THEN
         i = myworkid/nlocal
         CALL MPI_COMM_SPLIT( MPI_COMM_FIELDLINES,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
         CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
         mylocalmaster = master
      ELSE
         ! Basic copy of MPI_COMM_FIELDLINES
         CALL MPI_COMM_DUP( MPI_COMM_FIELDLINES, MPI_COMM_LOCAL, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'fieldlines_init_vmec: MPI_COMM_DUP',ierr_mpi)
         mylocalid = myworkid
         mylocalmaster = master
         numprocs_local = nprocs_fieldlines
      END IF

      ! Open VMEC file
      ! Initialize Virtual Casing
      nu = 8 * mpol
      nu = 2 ** CEILING(log(DBLE(nu))/log(2.0_rprec))
      IF (nu < 128) nu = 128
      nv = 8 * ntor + 1
      nv = 2 ** CEILING(log(DBLE(nv))/log(2.0_rprec))
      IF (nv < 128) nv = 128
      ALLOCATE(xm_temp(mnmax),xn_temp(mnmax), STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM_TEMP XN_TEMP',ier)
      
      ! Initialize VC
      IF (.not.lvolint .and. .not.lplasma_only) THEN
         ! Load Variables
         ALLOCATE(rmnc_temp(mnmax,2),zmns_temp(mnmax,2))
         ALLOCATE(bumnc_temp(mnmax,1),bvmnc_temp(mnmax,1))
         xm_temp=INT(xm)
         xn_temp=-INT(xn)/nfp
         rmnc_temp(:,1)=rmnc(:,ns-1)
         rmnc_temp(:,2)=rmnc(:,ns)
         zmns_temp(:,1)=zmns(:,ns-1)
         zmns_temp(:,2)=zmns(:,ns)
         bumnc_temp(:,1) = (1.5*bsupumnc(:,ns) - 0.5*bsupumnc(:,ns-1))
         bvmnc_temp(:,1) = (1.5*bsupvmnc(:,ns) - 0.5*bsupvmnc(:,ns-1))
         IF (lasym) THEN
            ALLOCATE(rmns_temp(mnmax,2),zmnc_temp(mnmax,2))
            ALLOCATE(bumns_temp(mnmax,1),bvmns_temp(mnmax,1))
            rmns_temp(:,1)=rmns(:,ns-1)
            rmns_temp(:,2)=rmns(:,ns)
            zmnc_temp(:,1)=zmnc(:,ns-1)
            zmnc_temp(:,2)=zmnc(:,ns)
            bumns_temp(:,1) = 1.5*bsupumns(:,ns) - 0.5*bsupumns(:,ns-1)
            bvmns_temp(:,1) = 1.5*bsupvmns(:,ns) - 0.5*bsupvmns(:,ns-1)
            CALL init_virtual_casing(mnmax,nu,nv,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         RMNS=rmns_temp, ZMNC=zmnc_temp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         BUMNS=bumns_temp,BVMNS=bvmns_temp)
            DEALLOCATE(rmns_temp,zmnc_temp)
            DEALLOCATE(bumns_temp,bvmns_temp)
         ELSE
            CALL init_virtual_casing(mnmax,nu,nv,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp)
         END IF
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)
      ELSE IF (.not. lplasma_only) THEN
         xm_temp = INT(xm)
         xn_temp = -INT(xn)
         IF (lasym) THEN
             CALL init_volint(mnmax,nu,nv,ns,xm_temp,xn_temp,rmnc,zmns,nfp,&
                              JUMNC=isigng*currumnc, JVMNC=isigng*currvmnc,&
                              RMNS=rmns,ZMNC=zmnc,&
                              JUMNS=isigng*currumns, JVMNS=isigng*currvmns)
         ELSE
             CALL init_volint(mnmax,nu,nv,ns,xm_temp,xn_temp,rmnc,zmns,nfp,&
                              JUMNC=isigng*currumnc, JVMNC=isigng*currvmnc)
         END IF
      END IF
      
      DEALLOCATE(xm_temp,xn_temp)
      adapt_tol = 1.0E-8
      adapt_rel = vc_adapt_tol
      IF (vc_adapt_tol < 0) adapt_tol = adapt_rel
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_vmec',ierr_mpi)
!DEC$ ENDIF
      
      IF (lverb) THEN
         IF (.not.lplasma_only) CALL virtual_casing_info(6)
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Plasma Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      ! Break up the Work
      chunk = FLOOR(REAL(nr*nphi*nz) / REAL(numprocs_local))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1

      ! This section sets up the work so we can use ALLGATHERV
!DEC$ IF DEFINED (MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs_local), moffsets(numprocs_local))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(numprocs_local)+mnum(numprocs_local)-1) == nr*nphi*nz) EXIT
         IF (i == numprocs_local) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs_local) = moffsets(i+1:numprocs_local) + 1
         i=i+1
      END DO
      mystart = moffsets(mylocalid+1)
      chunk  = mnum(mylocalid+1)
      myend   = mystart + chunk - 1
!DEC$ ENDIF
         IF (lafield_only) THEN
            DO s = mystart, myend
               i = MOD(s-1,nr)+1
               j = MOD(s-1,nr*nphi)
               j = FLOOR(REAL(j) / REAL(nr))+1
               k = CEILING(REAL(s) / REAL(nr*nphi))
               sflx = 0.0
               CALL GetAcyl(raxis_g(i),phiaxis(j),zaxis_g(k),&
                            br, bphi, bz, SFLX=sflx,info=ier)
               IF (ier == 0 .and. bphi /= 0 .and. sflx<=1) THEN
                  B_R(i,j,k)   = br
                  B_PHI(i,j,k) = bphi
                  B_Z(i,j,k)   = bz
               ELSE IF (lplasma_only) THEN
                  B_R(i,j,k)   = 0.0
                  B_PHI(i,j,k) = 1.0
                  B_Z(i,j,k)   = 0.0
               ELSE IF (ier == -3 .or. bphi == 0 .or. s>1) THEN
                  xaxis_vc = raxis_g(i)*cos(phiaxis(j))
                  yaxis_vc = raxis_g(i)*sin(phiaxis(j))
                  zaxis_vc = zaxis_g(k)
                  ier = 1
                  CALL vecpot_vc(xaxis_vc,yaxis_vc,zaxis_vc,bx_vc,by_vc,bz_vc,ier)
                  IF (ier == 0) THEN
                     br_vc   = bx_vc * cos(phiaxis(j)) + by_vc * sin(phiaxis(j))
                     bphi_vc = by_vc * cos(phiaxis(j)) - bx_vc * sin(phiaxis(j))
                     IF (ABS(br_vc) > 0)   B_R(i,j,k)   = B_R(i,j,k) + br_vc
                     IF (ABS(bphi_vc) > 0) B_PHI(i,j,k) = B_PHI(i,j,k) + bphi_vc
                     IF (ABS(bz_vc) > 0)   B_Z(i,j,k)   = B_Z(i,j,k) + bz_vc
                  END IF
               ELSE
                  ! This is an error code check
                  PRINT *,'ERROR in GetAcyl Detected'
                  PRINT *,'R,PHI,Z',raxis_g(i),phiaxis(j),zaxis_g(k)
                  print *,'br,bphi,bz,myworkid',br,bphi,bz,myworkid
                  stop 'ERROR in GetAcyl'
               END IF
               IF (lverb .and. (MOD(s,nr) == 0)) THEN
                  CALL backspace_out(6,6)
                  WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
                  CALL FLUSH(6)
               END IF
            END DO
         ELSE
            DO s = mystart, myend
               i = MOD(s-1,nr)+1
               j = MOD(s-1,nr*nphi)
               j = FLOOR(REAL(j) / REAL(nr))+1
               k = CEILING(REAL(s) / REAL(nr*nphi))
               sflx = 0.0
               ! The GetBcyl Routine returns -3 if cyl2flx thinks s>1
               ! however, if cyl2flx fails to converge then s may be
               ! greater than 1 but cyl2flux won't throw the -3 code.
               ! In this case GetBcyl returns br,bphi,bz = 0.  So
               ! bphi == 0 or ier ==-3 indicate that a point is
               ! outside the VMEC domain.
               CALL GetBcyl(raxis_g(i),phiaxis(j),zaxis_g(k),&
                                  br, bphi, bz, SFLX=sflx,info=ier)
               IF (ier == 0 .and. bphi /= 0 .and. sflx<=1) THEN
                  B_R(i,j,k)   = br
                  B_PHI(i,j,k) = bphi
                  B_Z(i,j,k)   = bz
               ELSE IF (lplasma_only) THEN
                  B_R(i,j,k)   = 0.0
                  B_PHI(i,j,k) = 1.0
                  B_Z(i,j,k)   = 0.0
               ELSE IF (ier == -3 .or. bphi == 0 .or. s>1) THEN
                  xaxis_vc = raxis_g(i)*cos(phiaxis(j))
                  yaxis_vc = raxis_g(i)*sin(phiaxis(j))
                  zaxis_vc = zaxis_g(k)
                  ier = 1
                  CALL bfield_vc(xaxis_vc,yaxis_vc,zaxis_vc,bx_vc,by_vc,bz_vc,ier)
                  IF (ier == 0) THEN
                     br_vc   = bx_vc * cos(phiaxis(j)) + by_vc * sin(phiaxis(j))
                     bphi_vc = by_vc * cos(phiaxis(j)) - bx_vc * sin(phiaxis(j))
                     IF (ABS(br_vc) > 0)   B_R(i,j,k)   = B_R(i,j,k) + br_vc
                     IF (ABS(bphi_vc) > 0) B_PHI(i,j,k) = B_PHI(i,j,k) + bphi_vc
                     IF (ABS(bz_vc) > 0)   B_Z(i,j,k)   = B_Z(i,j,k) + bz_vc
                  END IF
               ELSE
                  ! This is an error code check
                  PRINT *,'ERROR in GetBcyl Detected'
                  PRINT *,'R,PHI,Z',raxis_g(i),phiaxis(j),zaxis_g(k)
                  print *,'br,bphi,bz,myworkid',br,bphi,bz,myworkid
                  stop 'ERROR in GetBcyl'
               END IF
               IF (lverb .and. (MOD(s,nr) == 0)) THEN
                  CALL backspace_out(6,6)
                  WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
                  CALL FLUSH(6)
               END IF
            END DO
         END IF
      
      ! Free variables
      IF (.not. lplasma_only) CALL free_virtual_casing
      IF (.not. (lemc3 .or. ledge_start)) CALL read_wout_deallocate
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         CALL FLUSH(6)
      END IF    
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_vmec',ierr_mpi)
!       ! Adjust indexing to send 2D arrays
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        B_R,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        B_PHI,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        B_Z,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_LOCAL,ierr_mpi)
       DEALLOCATE(mnum)
       DEALLOCATE(moffsets)
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
      !IF (nprocs_fieldlines > nlocal) THEN
         ierr_mpi=0
      !   For John Schmitt
      !   CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      !   IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'fieldlines_init_vmec: MPI_COMM_FREE',ierr_mpi)
      !END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_vmec',ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_vmec
