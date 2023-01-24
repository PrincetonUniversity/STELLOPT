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
      USE read_wout_mod, phi_vmec => phi
      USE vmec_utils
      USE virtual_casing_mod, pi2_vc => pi2
      USE fieldlines_runtime, extcur_fieldlines => extcur
      USE fieldlines_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, vc_adapt_tol, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl, PRES_G
      USE wall_mod, ONLY: wall_load_mn, wall_info,vertex,face
      USE mpi_params                                                    ! MPI
      USE EZspline_obj
      USE EZspline
!      USE mpi
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
      LOGICAL :: lnyquist, luse_vc, lcreate_wall
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, s, i, j, k, nu, nv, mystart, myend, mnmax_temp, u, v
      INTEGER, ALLOCATABLE :: xn_temp(:), xm_temp(:)
      REAL :: br_vc, bphi_vc, bz_vc, xaxis_vc, yaxis_vc, zaxis_vc,&
              bx_vc, by_vc
      REAL(rprec) :: br, bphi, bz, sflx, uflx
      DOUBLE PRECISION, ALLOCATABLE :: mfact(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_temp(:,:),zmns_temp(:,:),&
                           bumnc_temp(:,:),bvmnc_temp(:,:),&
                           rmns_temp(:,:),zmnc_temp(:,:),&
                           bumns_temp(:,:),bvmns_temp(:,:)
      LOGICAL :: lvolint = .false.
      TYPE(EZspline1_r8) :: p_spl
      INTEGER :: bcs1(2)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/0,0/)
      ! Handle what we do
      luse_vc = ((lcoil .or. lmgrid) .and. .not.lplasma_only)

      ! Divide up Work
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      ! Open VMEC file
      IF (myworkid == master) THEN
         IF (ALLOCATED(extcur)) DEALLOCATE(extcur) ! From reading MGRID
         CALL read_wout_file(TRIM(id_string),ier)
         IF (ier /= 0) CALL handle_err(VMEC_WOUT_ERR,'beams3d_init_vmec',ier)
      END IF
      
#if defined(MPI_OPT)
      ! We do this to avoid multiple opens of wout file
      CALL MPI_BCAST(ns,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(mpol,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(ntor,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(nfp,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(mnyq,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(nnyq,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(mnmax,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(mnmax_nyq,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(lasym,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(lthreed,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(lwout_opened,1,MPI_LOGICAL, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (myworkid /= master) THEN
         ALLOCATE(phi_vmec(ns),presf(ns),iotaf(ns),phipf(ns))
         ALLOCATE(xm(mnmax),xn(mnmax),xm_nyq(mnmax_nyq),xn_nyq(mnmax_nyq))
         ALLOCATE(rmnc(mnmax,ns),zmns(mnmax,ns),lmns(mnmax,ns),bsupumnc(mnmax_nyq,ns),bsupvmnc(mnmax_nyq,ns))
         IF (lasym) ALLOCATE(rmns(mnmax,ns),zmnc(mnmax,ns),lmnc(mnmax,ns),bsupumns(mnmax_nyq,ns),bsupvmns(mnmax_nyq,ns))
      END IF
      CALL MPI_BCAST(phi_vmec,ns,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(presf,ns,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(phipf,ns,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(iotaf,ns,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(xm,mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(xn,mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(xm_nyq,mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(xn_nyq,mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(rmnc,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(zmns,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(lmns,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(bsupumnc,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(bsupvmnc,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (lasym) THEN
         CALL MPI_BCAST(rmns,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         CALL MPI_BCAST(zmnc,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         CALL MPI_BCAST(lmnc,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         CALL MPI_BCAST(bsupumns,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
         CALL MPI_BCAST(bsupvmns,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      END IF
#endif

      ! Write info to screen
      IF (lverb) THEN
         WRITE(6,'(A)')               '----- VMEC Information -----'
         WRITE(6,'(A,A)')             '   FILE: ',TRIM(id_string)
         WRITE(6,'(A,F9.5,A,F9.5,A)') '   R       = [',rmin_surf,',',rmax_surf,']'
         WRITE(6,'(A,F8.5,A,F8.5,A)') '   Z       = [',-zmax_surf,',',zmax_surf,']'
         IF (ABS(Itor) > 1E8) THEN
            WRITE(6,'(A,F7.3,A,F7.3,A)') '   BETA    = ',betatot,';  I  = ',Itor*1E-9,' [GA]'
         ELSEIF (ABS(Itor) > 1E5) THEN
            WRITE(6,'(A,F7.3,A,F7.3,A)') '   BETA    = ',betatot,';  I  = ',Itor*1E-6,' [MA]'
         ELSEIF (ABS(Itor) > 1E2) THEN
            WRITE(6,'(A,F7.3,A,F7.3,A)') '   BETA    = ',betatot,';  I  = ',Itor*1E-3,' [kA]'
         ELSE
            WRITE(6,'(A,F7.3,A,F7.3,A)') '   BETA    = ',betatot,';  I  = ',Itor,' [A]'
         END IF
         WRITE(6,'(A,F7.3,A)')        '   AMINOR  = ',Aminor,' [m]'
         WRITE(6,'(A,F7.3,A)')        '   PHIEDGE = ',phi_vmec(ns),' [Wb]'
         WRITE(6,'(A,F7.3,A)')        '   VOLUME  = ',Volume,' [m^3]'
      END IF

      IF (luse_vc) THEN
         nu = 8 * mpol + 1 
         nu = 2 ** CEILING(log(DBLE(nu))/log(2.0_rprec))
         nv = 8 * ntor + 1
         nv = 2 ** CEILING(log(DBLE(nv))/log(2.0_rprec))
         IF (nv < 128) nv = 128

         ! Handle Nyquist issues
         IF (SIZE(xm_nyq) > SIZE(xm)) THEN
            mnmax_temp = SIZE(xm_nyq)
            lnyquist = .true.
         ELSE
            mnmax_temp = mnmax
            lnyquist = .false.
         END IF
         ALLOCATE(xm_temp(mnmax_temp),xn_temp(mnmax_temp), STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM_TEMP XN_TEMP',ier)
         ALLOCATE(rmnc_temp(mnmax_temp,2),zmns_temp(mnmax_temp,2),&
                  bumnc_temp(mnmax_temp,1),bvmnc_temp(mnmax_temp,1), STAT = ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RMNC_TEMP ZMNS_TEMP BUMNC_TEMP BVMNC_TEMP',ier)
         IF (lasym) ALLOCATE(rmns_temp(mnmax_temp,2),zmnc_temp(mnmax_temp,2),&
                     bumns_temp(mnmax_temp,1),bvmns_temp(mnmax_temp,1), STAT = ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RMNS_TEMP ZMNC_TEMP BUMNS_TEMP BVMNS_TEMP',ier)
         IF (lnyquist) THEN
            xm_temp = xm_nyq
            xn_temp = -xn_nyq/nfp  ! Because init_virtual_casing uses (mu+nv) not (mu-nv*nfp)
            IF(lverb) WRITE(6,'(A)')        '   NYQUIST DETECTED IN WOUT FILE!'
            DO u = 1,mnmax_temp
               DO v = 1, mnmax
                  IF ((xm(v) .eq. xm_nyq(u)) .and. (xn(v) .eq. xn_nyq(u))) THEN
                     rmnc_temp(u,1) = rmnc(v,ns-1)
                     zmns_temp(u,1) = zmns(v,ns-1)
                     rmnc_temp(u,2) = rmnc(v,ns)
                     zmns_temp(u,2) = zmns(v,ns)
                     IF (lasym) THEN
                        rmns_temp(u,1) = rmns(v,ns-1)
                        zmnc_temp(u,1) = zmnc(v,ns-1)
                        rmns_temp(u,2) = rmns(v,ns)
                        zmnc_temp(u,2) = zmnc(v,ns)
                     END IF
                  END IF
               END DO
            END DO
         ELSE
            xm_temp = xm
            xn_temp = -xn/nfp  ! Because init_virtual_casing uses (mu+nv) not (mu-nv*nfp)
            rmnc_temp(:,1) = rmnc(:,ns-1)
            zmns_temp(:,1) = zmns(:,ns-1)
            rmnc_temp(:,2) = rmnc(:,ns)
            zmns_temp(:,2) = zmns(:,ns)
            IF (lasym) THEN
               rmns_temp(:,1) = rmns(:,ns-1)
               zmnc_temp(:,1) = zmnc(:,ns-1)
               rmns_temp(:,2) = rmns(:,ns)
               zmnc_temp(:,2) = zmnc(:,ns)
            END IF
         ENDIF
         ! Get B onto full grid
         ALLOCATE(mfact(mnmax_temp,2))
         WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
            mfact(:,1)= 1.5
            mfact(:,2)=-0.5
         ELSEWHERE
            mfact(:,1)= 1.5*SQRT((ns-1.0)/(ns-1.5))
            mfact(:,2)=-0.5*SQRT((ns-1.0)/(ns-2.5))
         ENDWHERE
         bumnc_temp(:,1) = mfact(:,1)*bsupumnc(:,ns) + mfact(:,2)*bsupumnc(:,ns-1)
         bvmnc_temp(:,1) = mfact(:,1)*bsupvmnc(:,ns) + mfact(:,2)*bsupvmnc(:,ns-1)
         IF (lasym) THEN
            bumns_temp(:,1) = mfact(:,1)*bsupumns(:,ns) + mfact(:,2)*bsupumns(:,ns-1)
            bvmns_temp(:,1) = mfact(:,1)*bsupvmns(:,ns) + mfact(:,2)*bsupvmns(:,ns-1)
            CALL init_virtual_casing(mnmax_temp,nu,nv,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         RMNS=rmns_temp, ZMNC=zmnc_temp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         BUMNS=bumns_temp,BVMNS=bvmns_temp,&
                                         COMM=MPI_COMM_FIELDLINES)
            DEALLOCATE(rmns_temp,zmnc_temp)
            DEALLOCATE(bumns_temp,bvmns_temp)
         ELSE
            CALL init_virtual_casing(mnmax_temp,nu,nv,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         COMM=MPI_COMM_FIELDLINES)
         END IF
         DEALLOCATE(mfact)
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)

         IF (lpres) THEN
            CALL EZspline_init(p_spl,ns,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,&
                               'EZspline_init/fieldlines_init_vmec',ier)
            p_spl%isHermite = 1
            p_spl%x1 = phi_vmec/phi_vmec(ns) !Spline over normalized toroidal flux
            CALL EZspline_setup(p_spl,presf,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,&
                               'EZspline_setup/fieldlines_init_vmec',ier)
         END IF 

         adapt_tol = 0.0
         adapt_rel = vc_adapt_tol
         DEALLOCATE(xm_temp,xn_temp)
      END IF
      
      IF (lverb) THEN
         IF (.not.lplasma_only) CALL virtual_casing_info(6)
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Plasma Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)


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
               ELSE IF (ier == -3 .or. bphi == 0 .or. sflx>1) THEN
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
                  B_R(i,j,k)   = B_R(i,j,k) + br
                  B_PHI(i,j,k) = B_PHI(i,j,k)+ bphi
                  B_Z(i,j,k)   = B_Z(i,j,k) + bz
                  IF (lpres) THEN 
                     CALL EZspline_interp(p_spl,sflx,PRES_G(i,j,k),ier)
                  END IF
               ELSE IF (lplasma_only) THEN
                  B_R(i,j,k)   = 0.0
                  B_PHI(i,j,k) = 1.0
                  B_Z(i,j,k)   = 0.0
               ELSE IF (ier == -3 .or. bphi == 0 .or. sflx>1) THEN
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
                  IF (lpres) THEN ! constant presssure outsise LCFS 
                     sflx = 1.0
                     CALL EZspline_interp(p_spl,sflx,PRES_G(i,j,k),ier)
                  END IF
               ! ELSE
               !    ! This is an error code check
               !    PRINT *,'ERROR in GetBcyl Detected'
               !    PRINT *,'R,PHI,Z',raxis_g(i),phiaxis(j),zaxis_g(k)
               !    print *,'br,bphi,bz,myworkid',br,bphi,bz,myworkid
               !    stop 'ERROR in GetBcyl'
               END IF
               IF (lverb .and. (MOD(s,nr) == 0)) THEN
                  CALL backspace_out(6,6)
                  WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
                  CALL FLUSH(6)
               END IF
            END DO
         END IF
      
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      
      ! Free variables
      IF (luse_vc) CALL free_virtual_casing(MPI_COMM_FIELDLINES)
      IF (myworkid == master) THEN
         CALL read_wout_deallocate
      ELSE
         lwout_opened = .FALSE.
         DEALLOCATE(phi_vmec,presf,phipf,iotaf)
         DEALLOCATE(xm,xn,xm_nyq,xn_nyq)
         DEALLOCATE(rmnc,zmns,lmns,bsupumnc,bsupvmnc)
         IF (lasym) DEALLOCATE(rmns,zmnc,lmnc,bsupumns,bsupvmns)
      END IF

      IF (EZspline_allocated(p_spl)) CALL EZspline_free(p_spl,ier)
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         CALL FLUSH(6)
      END IF    
      
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_vmec',ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_vmec: MPI_COMM_LOCAL',ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_vmec',ierr_mpi)
#endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_vmec
