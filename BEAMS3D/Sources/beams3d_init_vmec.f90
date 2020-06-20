!-----------------------------------------------------------------------
!     Subroutine:    beams3d_init_vmec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/25/2012
!     Description:   This subroutine reads the VMEC wout file and
!                    initializes the plasma fields.  This is achieved
!                    through utilization of a virtual casing principle.
!                    This is stil under development but working.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_vmec
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod
      USE vmec_utils
      USE virtual_casing_mod, pi2_vc => pi2
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, vc_adapt_tol, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl, TE_spl_s, NE_spl_s, TI_spl_s, &
                                 nte, nne, nti, TE, NE, TI, Vp_spl_s, S_ARR,&
                                 U_ARR, POT_ARR, POT_spl_s, nne, nte, nti, npot, &
                                 ZEFF_spl_s, nzeff, ZEFF_ARR, req_axis, zeq_axis, &
                                 phiedge_eq, reff_eq
      USE wall_mod, ONLY: wall_load_mn, wall_info,vertex,face
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
      LOGICAL :: lnyquist, luse_vc, lcreate_wall
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, s, i, j, k, nu, nv, mystart, myend, mnmax_temp, u, v
      INTEGER :: bcs1_s(2)
      INTEGER, ALLOCATABLE :: xn_temp(:), xm_temp(:)
      REAL :: br_vc, bphi_vc, bz_vc, xaxis_vc, yaxis_vc, zaxis_vc,&
              bx_vc, by_vc, dr_temp
      REAL(rprec) :: br, bphi, bz, sflx, uflx, xaxis, yaxis
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_temp(:,:),zmns_temp(:,:),&
                           bumnc_temp(:,:),bvmnc_temp(:,:),&
                           rmns_temp(:,:),zmnc_temp(:,:),&
                           bumns_temp(:,:),bvmns_temp(:,:)
      LOGICAL :: lold_field = .true.  ! If false we attempt to subtract the vacuum field from the VMEC total field
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle what we do
      luse_vc = (lcoil .and. .not.lplasma_only)
      lcreate_wall = (lplasma_only .and. .not. lvessel)

      ! Divide up Work
      mylocalid = myworkid
      numprocs_local = 1
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master
      bx_vc = 0.0; by_vc = 0.0; bz_vc = 0.0

      ! Open VMEC file
      IF (myworkid == master) THEN
         CALL read_wout_file(TRIM(id_string),ier)
         IF (ier /= 0) CALL handle_err(VMEC_WOUT_ERR,'beams3d_init_vmec',ier)
      END IF
      
#if defined(MPI_OPT)
      ! We do this to avoid multiple opens of wout file
      CALL MPI_BCAST(ns,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mpol,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(ntor,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nfp,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mnyq,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nnyq,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mnmax,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mnmax_nyq,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(lasym,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(lthreed,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(lwout_opened,1,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (myworkid /= master) THEN
         ALLOCATE(vp(ns),phi(ns))
         ALLOCATE(xm(mnmax),xn(mnmax),xm_nyq(mnmax_nyq),xn_nyq(mnmax_nyq))
         ALLOCATE(rmnc(mnmax,ns),zmns(mnmax,ns),bsupumnc(mnmax_nyq,ns),bsupvmnc(mnmax_nyq,ns))
         IF (lasym) ALLOCATE(rmns(mnmax,ns),zmnc(mnmax,ns),bsupumns(mnmax_nyq,ns),bsupvmns(mnmax_nyq,ns))
      END IF
      CALL MPI_BCAST(vp,ns,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(phi,ns,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(xm,mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(xn,mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(xm_nyq,mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(xn_nyq,mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(rmnc,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(zmns,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(bsupumnc,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(bsupvmnc,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (lasym) THEN
         CALL MPI_BCAST(rmns,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
         CALL MPI_BCAST(zmnc,ns*mnmax,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
         CALL MPI_BCAST(bsupumns,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
         CALL MPI_BCAST(bsupvmns,ns*mnmax_nyq,MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS,ierr_mpi)
      END IF
      phiedge_eq = phi(ns)
      reff_eq = Aminor
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
         WRITE(6,'(A,F7.3,A)')        '   AMINOR  = ',reff_eq,' [m]'
         WRITE(6,'(A,F7.3,A)')        '   PHIEDGE = ',phiedge_eq,' [Wb]'
         WRITE(6,'(A,F7.3,A)')        '   VOLUME  = ',Volume,' [m^3]'
      END IF

      ! Load the Vp Spline if using the beams
      ! Note in VMEC dV/ds = \int\int sqrt(g) dtheta dphi
      ! However VMEC integrates over u and v
      ! So Vp is missing a factor of 4*pi*pi
      ! Also is on the half grid
      vp(1) = 1.5*vp(2) - 0.5*vp(3) ! Both on half grid
      vp(2:ns-1) = 0.5*(vp(2:ns-1)+vp(3:ns)) ! Average to full grid
      vp(ns) = 2*vp(ns) - vp(ns-1) ! ns on half grid, ns-1 on full grid
      vp = vp*pi2*pi2
      bcs1_s=(/ 0, 0 /)
      CALL EZspline_init(Vp_spl_s,ns,bcs1_s,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_vmec',ier)
      Vp_spl_s%isHermite   = 1
      CALL EZspline_setup(Vp_spl_s,vp(1:ns),ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_vmec',ier)

      ! Calculate the axis values (need negative since we haven't flipped Kernel yet)
      IF (mylocalid == master) THEN
         DO i = 1, nphi
            req_axis(i) = 0
            zeq_axis(i) = 0
            DO u = 1, mnmax
               req_axis(i) = req_axis(i)+rmnc(u,1)*COS(-xn(u)*phiaxis(i))
               zeq_axis(i) = zeq_axis(i)+zmns(u,1)*SIN(-xn(u)*phiaxis(i))
            END DO
         END DO
         IF (lasym) THEN
            DO i = 1, nphi
               DO u = 1, mnmax
                  req_axis(i) = req_axis(i)+rmns(u,1)*SIN(-xn(u)*phiaxis(i))
                  zeq_axis(i) = zeq_axis(i)+zmnc(u,1)*COS(-xn(u)*phiaxis(i))
               END DO
            END DO
         END IF
      END IF


      ! If we ask for a plasma-only run and do not provide a vessel then make one.
      IF (lcreate_wall) THEN
         lvessel = .TRUE.  ! Do this so the other parts of the code know there is a vessel
         k = ns
         CALL wall_load_mn(DBLE(rmnc(1:mnmax,k)),DBLE(zmns(1:mnmax,k)),DBLE(xm),-DBLE(xn),mnmax,120,120,COMM=MPI_COMM_BEAMS)
         IF (lverb) CALL wall_info(6)
         !IF (mylocalid /= master) DEALLOCATE(vertex,face)
      END IF

      ! Initialize Virtual Casing
      IF (luse_vc) THEN
         nu = 8 * mpol + 1 
         nu = 2 ** CEILING(log(DBLE(nu))/log(2.0_rprec))
         nv = 8 * ntor + 1
         nv = 2 ** CEILING(log(DBLE(nv))/log(2.0_rprec))
         IF (nv < 128) nv = 128
         dr_temp    = 1./REAL(ns)

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
         rmnc_temp = 0; zmns_temp =0
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
         bumnc_temp(:,1) = (1.5*bsupumnc(:,ns) - 0.5*bsupumnc(:,ns-1))
         bvmnc_temp(:,1) = (1.5*bsupvmnc(:,ns) - 0.5*bsupvmnc(:,ns-1))
         IF (lasym) THEN
            bumns_temp(:,1) = 1.5*bsupumns(:,ns) - 0.5*bsupumns(:,ns-1)
            bvmns_temp(:,1) = 1.5*bsupvmns(:,ns) - 0.5*bsupvmns(:,ns-1)
            CALL init_virtual_casing(mnmax_temp,nu,nv,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         RMNS=rmns_temp, ZMNC=zmnc_temp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         BUMNS=bumns_temp,BVMNS=bvmns_temp,&
                                         COMM=MPI_COMM_BEAMS)
            DEALLOCATE(rmns_temp,zmnc_temp)
            DEALLOCATE(bumns_temp,bvmns_temp)
         ELSE
            CALL init_virtual_casing(mnmax_temp,nu,nv,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         COMM=MPI_COMM_BEAMS)
         END IF
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)
         
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
      chunk = FLOOR(REAL(nr*nphi*nz) / REAL(numprocs_local))
      mystart = mylocalid*chunk + 1
      myend = mystart + chunk - 1

      ! This section sets up the work so we can use ALLGATHERV
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
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
      DEALLOCATE(mnum)
      DEALLOCATE(moffsets)
#endif

      IF (mylocalid == mylocalmaster) THEN
         TE = 0; NE = 0; TI=0; S_ARR=1.5; U_ARR=0; POT_ARR=0; ZEFF_ARR = 1;
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         sflx = 0.0
         ! The GetBcyl Routine returns -3 if cyl2flx thinks s>1
         ! however, if cyl2flx fails to converge then s may be
         ! greater than 1 but cyl2flux will not throw the -3 code.
         ! In this case GetBcyl returns br,bphi,bz = 0.  So
         ! bphi == 0 or ier ==-3 indicate that a point is
         ! outside the VMEC domain.
         CALL GetBcyl(raxis_g(i),phiaxis(j),zaxis_g(k),&
                      br, bphi, bz, SFLX=sflx,UFLX=uflx,info=ier)
         IF (ier == 0 .and. bphi /= 0) THEN ! We have field data
            ! Save Grid data
            S_ARR(i,j,k) = MAX(sflx,0.0)
            IF (uflx<0)  uflx = uflx+pi2
            U_ARR(i,j,k) = uflx
            ! Handle equilibrium data
            IF (sflx <=1.0) THEN ! Inside equilibrium
               B_R(i,j,k)   = br
               B_PHI(i,j,k) = bphi
               B_Z(i,j,k)   = bz
            ELSE IF (.not. luse_vc) THEN  ! Overwrite data outside
               B_R(i,j,k)   = br
               B_PHI(i,j,k) = bphi
               B_Z(i,j,k)   = bz
               sflx = 1.0 ! Assume s=1 for lplasma_only
            END IF
            IF (sflx <= 1) THEN
               IF (nte > 0) CALL EZspline_interp(TE_spl_s,sflx,TE(i,j,k),ier)
               IF (nne > 0) CALL EZspline_interp(NE_spl_s,sflx,NE(i,j,k),ier)
               IF (nti > 0) CALL EZspline_interp(TI_spl_s,sflx,TI(i,j,k),ier)
               IF (npot > 0) CALL EZspline_interp(POT_spl_s,sflx,POT_ARR(i,j,k),ier)
               IF (nzeff > 0) CALL EZspline_interp(ZEFF_spl_s,sflx,ZEFF_ARR(i,j,k),ier)
            END IF
         ELSE IF (.not. luse_vc) THEN
            B_R(i,j,k)   = 0
            B_PHI(i,j,k) = 1
            B_Z(i,j,k)   = 0
         END IF
         ! Virtual Casing
         IF (luse_vc .and. sflx > 1) THEN
            xaxis_vc = raxis_g(i)*cos(phiaxis(j))
            yaxis_vc = raxis_g(i)*sin(phiaxis(j))
            zaxis_vc = zaxis_g(k)
            adapt_rel = 0.0
            adapt_tol = vc_adapt_tol*sqrt(B_R(i,j,k)**2+B_PHI(i,j,k)**2+B_Z(i,j,k)**2)
            ier = 1
            CALL bfield_vc(xaxis_vc,yaxis_vc,zaxis_vc,bx_vc,by_vc,bz_vc,ier)
            IF (ier == 0) THEN
               br_vc   = bx_vc * cos(phiaxis(j)) + by_vc * sin(phiaxis(j))
               bphi_vc = by_vc * cos(phiaxis(j)) - bx_vc * sin(phiaxis(j))
               IF (ABS(br_vc) > 0)   B_R(i,j,k)   = B_R(i,j,k) + br_vc
               IF (ABS(bphi_vc) > 0) B_PHI(i,j,k) = B_PHI(i,j,k) + bphi_vc
               IF (ABS(bz_vc) > 0)   B_Z(i,j,k)   = B_Z(i,j,k) + bz_vc
            ELSE IF (ier > 1) THEN ! Only report real errors not failure to find tolerance errors
               WRITE(6,*) myworkid,mylocalid,i,j,k,ier
            END IF
         !ELSE
            ! This is an error code check
            !PRINT *,'ERROR in GetBcyl Detected'
            !PRINT *,'R,PHI,Z',raxis_g(i),phiaxis(j),zaxis_g(k)
            !print *,'br,bphi,bz,myworkid',br,bphi,bz,mylocalid
            !CALL FLUSH(6)
            !stop 'ERROR in GetBcyl'
         END IF
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END IF
      END DO
      
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      
      ! Free variables
      IF (luse_vc) CALL free_virtual_casing(MPI_COMM_BEAMS)
      IF (myworkid == master) THEN
         CALL read_wout_deallocate
      ELSE
         lwout_opened = .FALSE.
         DEALLOCATE(vp,phi)
         DEALLOCATE(xm,xn,xm_nyq,xn_nyq)
         DEALLOCATE(rmnc,zmns,bsupumnc,bsupvmnc)
         IF (lasym) DEALLOCATE(rmns,zmnc,bsupumns,bsupvmns)
      END IF
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         WRITE(6,*)
         CALL FLUSH(6)
      END IF    

      ! Fix ZEFF
      IF (mylocalid == mylocalmaster) WHERE(ZEFF_ARR < 1) ZEFF_ARR = 1

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_vmec',ierr_mpi)
#endif

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_vmec
