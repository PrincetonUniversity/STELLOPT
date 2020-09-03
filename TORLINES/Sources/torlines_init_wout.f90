!-----------------------------------------------------------------------
!     Subroutine:    pies_init_wout
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine initizlizes the PIES background
!                    coordinates from a VMEC wout file.  This code
!                    assumes the user wishes to extend the VMEC domain
!                    by extra surfaces.  First, the VMEC surfaces are
!                    transformed to real space (along with the
!                    magnetic field).  They are then interpolated in
!                    real space from their flux space representation
!                    to a coordinate grid equidistant in rho.  Then
!                    the extra surfaces are added extrapolated in real
!                    space.  Virtual casing is employed to calculate the
!                    plasma field outside the VMEC domain along with
!                    interpolation from the mgrid file used to produce
!                    the VMEC equilibrium.  Then all quantities are
!                    Fourier transformed back to Fourier space on the
!                    new background coordinates.
!-----------------------------------------------------------------------
      SUBROUTINE torlines_init_wout
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, extcur_in => extcur, nextcur_in => nextcur,&
                         lasym_in => lasym, mgrid_file_in => mgrid_file
      USE torlines_realspace
      USE torlines_runtime
      USE torlines_background, ONLY: rho, nfp_m => nfp, bound_separation,&
                               rmnc_sav, rmns_sav, zmns_sav, zmnc_sav, &
                               xm_sav, xn_sav, vc_adapt_tol, isgn, &
                               rb_ws, zb_ws, dr_m, ixm_m, ixn_m, &
                               rmnc_m, rmns_m, zmnc_m, zmns_m, u0, v0,&
                               mnmax_m, isgn
      USE virtual_casing_mod, pi2_vc => pi2, nuvp_vc => nuvp
      USE wall_mod, ONLY: wall_free, collide
      USE mpi_params
      USE EZspline_obj
      USE EZspline
      use mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          uv          Dummy index over real-space
!          ik          Dummy index over radial surfaces
!          vsurf       VMEC surface in PIES background coordinates
!          rho_in      VMEC Rho
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, ier, mn, im, in, ik , i, j, ns1, k1,mn0,&
                 u,v
      INTEGER :: bcs1(2)
      INTEGER, ALLOCATABLE :: xn_temp(:), xm_temp(:)
      REAL(rprec) :: val1, val2, dval, scale, rhomax, dr, f0_temp, dz,&
                     alvb, b1, c1, re, ae, r1, alub
      REAL(rprec), ALLOCATABLE :: rreal_in(:,:,:), zreal_in(:,:,:)
      REAL(rprec), ALLOCATABLE :: bureal_in(:,:,:), bvreal_in(:,:,:)
      REAL(rprec), ALLOCATABLE :: rho_vmec(:), rho_ext(:), ftemp(:),fmn(:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumnc_temp(:,:), bsupvmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumns_temp(:,:), bsupvmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bumnc_temp(:,:), bvmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bmnc_temp(:,:), bmns_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_tempr(:,:), zmns_tempr(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rmns_tempr(:,:), zmnc_tempr(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: bsupumnc_tempr(:,:), bsupvmnc_tempr(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: bsupumns_tempr(:,:), bsupvmns_tempr(:,:)
      REAL(rprec), ALLOCATABLE :: rmnc1(:),zmns1(:),rmnc2(:),zmns2(:)
      REAL(rprec), ALLOCATABLE :: rmns1(:),zmnc1(:),rmns2(:),zmnc2(:)
      TYPE(EZspline1_r8) :: f_spl
      LOGICAL :: lhit
      DOUBLE PRECISION :: x0, y0, z0, x2, y2, z2, xw, yw, zw

      EXTERNAL :: map_v
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/0,0/)
      CALL read_wout_file(TRIM(id_string),ierr=ier)
      IF (ier /= 0) CALL handle_err(VMEC_WOUT_ERR,id_string,ier)
      ! Readjust EXTCUR
      IF (ALLOCATED(extcur)) DEALLOCATE(extcur)
      nextcur = nextcur_in
      ALLOCATE(extcur(1:nextcur))
      extcur(1:nextcur) = extcur_in(1:nextcur)
      lasym = lasym_in
      nfp_m = nfp
      ! Output some stuff to the screen
      if (lverb) THEN
         write(6,*) '-----VMEC File Parameters-----'
         write(6,'(A,A)') '    file: ',TRIM(id_string)
         write(6,'(A,I3,A,I3)') '       m: ',mpol,'   nu: ',nu
         write(6,'(A,I3,A,I3)') '       n: ',ntor,'   nv: ',nv
         write(6,'(A,I4,A)') '   mnmax: ',mnmax
         write(6,'(A,I3)') '     nfp: ',nfp
         write(6,'(A,I3)') '      ns: ',ns
         IF (ABS(Itor) .ge. 1.0E6) THEN
            write(6,'(A,F8.3,A)')             '   Total Current: ',Itor/1.0E6,' [MA]'
         ELSE IF (ABS(Itor) .ge. 1.0E3) THEN
            write(6,'(A,F8.3,A)')             '   Total Current: ',Itor/1.0E3,' [kA]'
         ELSE
            write(6,'(A,F8.3,A)')             '   Total Current: ',Itor,' [A]'
         END IF
         CALL FLUSH(6)
      END IF
      ! Get fields on full mesh
      bsupumnc(:,1) = 1.5*bsupumnc(:,2) - 0.5*bsupumnc(:,3)
      bsupvmnc(:,1) = 1.5*bsupvmnc(:,2) - 0.5*bsupvmnc(:,3)
      bsupumnc(:,2:ns-1) = 0.5*(bsupumnc(:,2:ns-1)+bsupumnc(:,3:ns))
      bsupvmnc(:,2:ns-1) = 0.5*(bsupvmnc(:,2:ns-1)+bsupvmnc(:,3:ns))
      bsupumnc(:,ns) = 2.0*bsupumnc(:,ns) - bsupumnc(:,ns-1) ! Note at this point ns on half grid ns-1 on full grid
      bsupvmnc(:,ns) = 2.0*bsupvmnc(:,ns) - bsupvmnc(:,ns-1)
      IF (lasym) THEN
         bsupumns(:,1) = 1.5*bsupumns(:,2) - 0.5*bsupumns(:,3)
         bsupvmns(:,1) = 1.5*bsupvmns(:,2) - 0.5*bsupvmns(:,3)
         bsupumns(:,2:ns-1) = 0.5*(bsupumns(:,2:ns-1)+bsupumns(:,3:ns))
         bsupvmns(:,2:ns-1) = 0.5*(bsupvmns(:,2:ns-1)+bsupvmns(:,3:ns))
         bsupumns(:,ns) = 2.0*bsupumns(:,ns) - bsupumns(:,ns-1) ! Note at this point ns on half grid ns-1 on full grid
         bsupvmns(:,ns) = 2.0*bsupvmns(:,ns) - bsupvmns(:,ns-1)
      END IF
      ! ----- Interpolate VMEC arrays to the temp array
      CALL EZspline_init(f_spl,ns,bcs1,ier)
      ALLOCATE(rho_vmec(ns),fmn(nrho),ftemp(1:ns))
      ALLOCATE(rmnc_sav(mnmax,nrho),zmns_sav(mnmax,nrho))
      ALLOCATE(bsupumnc_temp(mnmax_nyq,nrho),bsupvmnc_temp(mnmax_nyq,nrho))
      ALLOCATE(bmnc_temp(mnmax_nyq,nrho))
      IF (lasym) THEN
         ALLOCATE(rmns_sav(mnmax,nrho),zmnc_sav(mnmax,nrho))
         ALLOCATE(bsupumns_temp(mnmax_nyq,nrho),bsupvmns_temp(mnmax_nyq,nrho))
         ALLOCATE(bmns_temp(mnmax_nyq,nrho))
      END IF
      ! Setup Rho_VMEC
      FORALL(ik = 1:ns) rho_vmec(ik) =  SQRT(REAL(ik-1) / REAL(ns-1))
      f_spl%isHermite = 1
      f_spl%x1 = rho_vmec
      ! Setup RHO (we do this so rho=1 is always a surface and equidistant grid)
      ae = SUM(rmnc(:,1)) ! RAXIS
      re = SUM(rmnc(:,ns)) ! REDGE
      !IF (lverb) PRINT *,re,ae
      r1 = (re-ae)+(bound_separation-1) ! This is the real scale
      !IF (lverb) PRINT *,r1,re-ae
      IF (myid_sharmem==master) THEN
         FORALL(ik = 1:nrho) rho(ik) =  REAL(ik-1) / REAL(nrho-1)
         rho = rho * r1  ! Scale rho
         j=nrho
         DO ik = 1, nrho
            IF (rho(ik) >= (re-ae)) THEN
               rho(ik) = 1.0
               j = ik
               nexternal = nrho-j
               EXIT
            END IF
         END DO
         FORALL(ik = 1:nrho) rho(ik) =  REAL(ik-1) / REAL(j-1)
      END IF
#if defined(MPI_OPT)
      CALL MPI_BCAST(nexternal,1,MPI_INTEGER, master, MPI_COMM_SHARMEM,ierr_mpi)
      CALL MPI_BCAST(j,1,MPI_INTEGER, master, MPI_COMM_SHARMEM,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
#endif
      !IF (lverb) PRINT *,rho
      ! Spline to rho axis (all quantities now on full grid)
      DO mn = 1, mnmax
         IF (xm(mn) == 0 .and. xn(mn) == 0) mn0 = mn
         ! RMNC
         f0_temp = rmnc(mn,1)
         ftemp(1:ns) = (rmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/2.+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_setup_wout1',ier)
         CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
         rmnc_sav(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) rmnc_sav(mn,1) = f0_temp
         ! ZMNS
         f0_temp = zmns(mn,1)
         ftemp(1:ns) = (zmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/2.+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_setup_wout2',ier)
         CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
         zmns_sav(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) zmns_sav(mn,1) = f0_temp
         IF (lasym) THEN
            ! RMNS
            f0_temp = rmns(mn,1)
            ftemp(1:ns) = (rmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init_wout',ier)
            CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
            rmns_sav(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
            IF (xm(mn) < 2) rmns_sav(mn,1) = f0_temp
            ! ZMNC
            f0_temp = zmnc(mn,1)
            ftemp(1:ns) = (zmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init_wout',ier)
            CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
            zmnc_sav(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
            IF (xm(mn) < 2) zmnc_sav(mn,1) = f0_temp
         END IF
      END DO
      ! Nyquist sized arrays
      DO mn = 1, mnmax_nyq
         ! BSUPUMNC
         f0_temp = bsupumnc(mn,1)
         ftemp(1:ns) = (bsupumnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_setup_wout3',ier)
         CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
         bsupumnc_temp(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
         IF (xm_nyq(mn) < 2) bsupumnc_temp(mn,1) = f0_temp
         ! BSUPVMNC
         f0_temp = bsupvmnc(mn,1)
         ftemp(1:ns) = (bsupvmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_setup_wout4',ier)
         CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
         bsupvmnc_temp(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
         IF (xm_nyq(mn) < 2) bsupvmnc_temp(mn,1) = f0_temp
         ! BMNC
         f0_temp = bmnc(mn,1)
         ftemp(1:ns) = (bmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_setup_wout4',ier)
         CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
         bmnc_temp(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
         IF (xm_nyq(mn) < 2) bmnc_temp(mn,1) = f0_temp
         IF (lasym) THEN
            ! BSUPUMNS
            f0_temp = bsupumns(mn,1)
            ftemp(1:ns) = (bsupumns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init_wout',ier)
            CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
            bsupumns_temp(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
            IF (xm_nyq(mn) < 2) bsupumns_temp(mn,1) = f0_temp
            ! BSUPVMNS
            f0_temp = bsupvmns(mn,1)
            ftemp(1:ns) = (bsupvmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init_wout',ier)
            CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
            bsupvmns_temp(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
            IF (xm_nyq(mn) < 2) bsupvmns_temp(mn,1) = f0_temp
            ! BMNS
            f0_temp = bmns(mn,1)
            ftemp(1:ns) = (bmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm_nyq(mn))/2.+1)
            ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init_wout',ier)
            CALL EZspline_interp(f_spl,nrho,rho,fmn,ier)
            bmns_temp(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm_nyq(mn))-2)
            IF (xm_nyq(mn) < 2) bmns_temp(mn,1) = f0_temp
         END IF
      END DO
      DEALLOCATE(ftemp,fmn)
      CALL EZspline_free(f_spl,ier)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'torlines_main',ierr_mpi)
#endif
      IF ((nrho > ns) .and. (bound_separation <= 1)) THEN
         nrho = ns
         vsurf = nrho
      ELSE
         dr = 1./(ns-1)
         DO ik = 1, nrho
            IF (rho(ik) <= 1.0) vsurf = ik
         END DO
      END IF
      isgn = isigng
      !IF (lverb) PRINT *,'test'
      ! Setup VC if needed
      IF (.not.lvac .and. lvc_field .and. (bound_separation > 1)) THEN
         ! Initialize the Virtual Casing
         ALLOCATE(xm_temp(mnmax_nyq),xn_temp(mnmax_nyq), STAT=ier)
         ALLOCATE(rmnc_tempr(mnmax_nyq,2),zmns_tempr(mnmax_nyq,2))
         ALLOCATE(bsupumnc_tempr(mnmax_nyq,1),bsupvmnc_tempr(mnmax_nyq,1))
         IF (lasym) THEN
            ALLOCATE(rmns_tempr(mnmax_nyq,2),zmnc_tempr(mnmax_nyq,2))
            ALLOCATE(bsupumns_tempr(mnmax_nyq,1),bsupvmns_tempr(mnmax_nyq,1))
         END IF
         xm_temp=xm_nyq
         xn_temp=-xn_nyq/nfp
         DO u = 1, mnmax_nyq
            DO v = 1, mnmax
               IF ((xm(v) .eq. xm_nyq(u)) .and. (xn(v) .eq. xn_nyq(u))) THEN
                  rmnc_tempr(u,1)=rmnc(v,ns-1)
                  rmnc_tempr(u,2)=rmnc(v,ns)
                  zmns_tempr(u,1)=zmns(v,ns-1)
                  zmns_tempr(u,2)=zmns(v,ns)
                  IF (lasym) THEN
                     rmns_tempr(u,1)=rmns(v,ns-1)
                     rmns_tempr(u,2)=rmns(v,ns)
                     zmnc_tempr(u,1)=zmnc(v,ns-1)
                     zmnc_tempr(u,2)=zmnc(v,ns)
                  END IF
               END IF
            END DO
         END DO
         bsupumnc_tempr(:,1) = bsupumnc(:,ns)
         bsupvmnc_tempr(:,1) = bsupvmnc(:,ns)
         IF (lasym) THEN
            bsupumns_tempr(:,1) = bsupumns(:,ns)
            bsupvmns_tempr(:,1) = bsupvmns(:,ns)
            CALL init_virtual_casing(mnmax_nyq,nu_vc,nv_vc,xm_temp,xn_temp,&
                                     rmnc_tempr,zmns_tempr,nfp,&
                                     RMNS=rmns_tempr, ZMNC=zmnc_tempr,&
                                     BUMNC=bsupumnc_tempr,BVMNC=bsupvmnc_tempr,&
                                     BUMNS=bsupumns_tempr,BVMNS=bsupvmns_tempr,&
                                     COMM=MPI_COMM_FIELDLINES)
            DEALLOCATE(rmns_tempr,zmnc_tempr)
            DEALLOCATE(bsupumns_tempr,bsupvmns_tempr)
         ELSE
            CALL init_virtual_casing(mnmax_nyq,nu_vc,nv_vc,xm_temp,xn_temp,&
                                     rmnc_tempr,zmns_tempr,nfp,&
                                     BUMNC=bsupumnc_tempr,BVMNC=bsupvmnc_tempr,&
                                     COMM=MPI_COMM_FIELDLINES)
         END IF
         DEALLOCATE(rmnc_tempr,zmns_tempr)
         DEALLOCATE(bsupumnc_tempr,bsupvmnc_tempr)
         DEALLOCATE(xm_temp,xn_temp)
         adapt_tol = 0.0
         adapt_rel = vc_adapt_tol
         IF (adapt_rel < 0.0) adapt_tol = -1.0
         MIN_CLS = 0
      ELSE IF (.not.lvac .and. .not.lvc_field .and. (bound_separation > 1)) THEN
         ALLOCATE(xm_temp(mnmax),xn_temp(mnmax), STAT=ier)
         xm_temp = xm
         xn_temp = -xn
         MIN_CLS = 0
         IF (lasym) THEN
             CALL init_volint(mnmax,nu_vc,nv_vc,ns,xm_temp,xn_temp,rmnc,zmns,nfp,&
                              JUMNC=isigng*currumnc, JVMNC=isigng*currvmnc,&
                              RMNS=rmns,ZMNC=zmnc,&
                              JUMNS=isigng*currumns, JVMNS=isigng*currvmns)
         ELSE
             CALL init_volint(mnmax,nu_vc,nv_vc,ns,xm_temp,xn_temp,rmnc,zmns,nfp,&
                              JUMNC=isigng*currumnc, JVMNC=isigng*currvmnc)
         END IF
         adapt_tol = 0.0
         adapt_rel = vc_adapt_tol
      END IF
      ! Now transform everthing into realspace
      !ALLOCATE(xu(1:nu),xv(1:nv))
      IF (myid_sharmem==master) THEN
         rreal  = 0
         zreal  = 0
         bsreal = 0
         bureal = 0
         bvreal = 0
         ! Do Nyquist first so we have harmonic for jacobian later
         ALLOCATE(xm_sav(mnmax_nyq),xn_sav(mnmax_nyq), STAT=ier)
         xm_sav = INT(xm_nyq)
         xn_sav = -INT(xn_nyq)/nfp
         CALL mntouv(1,nrho,mnmax_nyq,nu,nv,xu,xv,bsupumnc_temp,xm_sav,xn_sav,bureal,0,1)
         CALL mntouv(1,nrho,mnmax_nyq,nu,nv,xu,xv,bsupvmnc_temp,xm_sav,xn_sav,bvreal,0,0)
         CALL mntouv(1,nrho,mnmax_nyq,nu,nv,xu,xv,bmnc_temp,xm_sav,xn_sav,breal,0,0)
         IF (lasym) THEN
            CALL mntouv(1,nrho,mnmax_nyq,nu,nv,xu,xv,bsupumns_temp,xm_sav,xn_sav,bureal,1,0)
            CALL mntouv(1,nrho,mnmax_nyq,nu,nv,xu,xv,bsupvmns_temp,xm_sav,xn_sav,bvreal,1,0)
            CALL mntouv(1,nrho,mnmax_nyq,nu,nv,xu,xv,bmns_temp,xm_sav,xn_sav,breal,1,0)
         END IF
         DEALLOCATE(xm_sav,xn_sav)
         ALLOCATE(xm_sav(mnmax),xn_sav(mnmax), STAT=ier)
         xm_sav = INT(xm)
         xn_sav = -INT(xn)/nfp
         FORALL(ik = 1:nu) xu(ik) = REAL(ik-1)/REAL(nu-1)
         FORALL(ik = 1:nv) xv(ik) = REAL(ik-1)/REAL(nv-1)
         CALL mntouv(1,nrho,mnmax,nu,nv,xu,xv,rmnc_sav,xm_sav,xn_sav,rreal,0,1)
         CALL mntouv(1,nrho,mnmax,nu,nv,xu,xv,zmns_sav,xm_sav,xn_sav,zreal,1,0)
         IF (lasym) THEN
            CALL mntouv(1,nrho,mnmax,nu,nv,xu,xv,rmns_sav,xm_sav,xn_sav,rreal,1,0)
            CALL mntouv(1,nrho,mnmax,nu,nv,xu,xv,zmnc_sav,xm_sav,xn_sav,zreal,0,0)
         END IF
      ELSE
         ALLOCATE(xm_sav(mnmax),xn_sav(mnmax), STAT=ier)
         xm_sav = INT(xm)
         xn_sav = -INT(xn)/nfp
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
#endif
      DEALLOCATE(bsupumnc_temp,bsupvmnc_temp, bmnc_temp)
      IF (ALLOCATED(rmns_temp)) DEALLOCATE(bsupumns_temp,bsupvmns_temp, bmns_temp)
      ! Calculated the external surfaces
      IF (bound_separation > 1) THEN
         ! Try recalcing the exterior surfaces
         alub = pi2/(nu-1)
         alvb = pi2/(nv-1)
         re = 0
         re = 10*epsilon(re)
         ae = re
         dr_m = bound_separation-1
         mnmax_m = mnmax
         IF (ALLOCATED(rmnc_m)) DEALLOCATE(rmnc_m)
         IF (ALLOCATED(zmns_m)) DEALLOCATE(zmns_m)
         IF (ALLOCATED(ixm_m)) DEALLOCATE(ixm_m)
         IF (ALLOCATED(ixn_m)) DEALLOCATE(ixn_m)
         IF (ALLOCATED(rmns_m)) DEALLOCATE(rmns_m)
         IF (ALLOCATED(zmnc_m)) DEALLOCATE(zmnc_m)
         ALLOCATE(rmnc_m(mnmax))
         ALLOCATE(zmns_m(mnmax))
         ALLOCATE(ixm_m(mnmax))
         ALLOCATE(ixn_m(mnmax))
         ALLOCATE(rmns_m(mnmax))
         ALLOCATE(zmnc_m(mnmax))
         rmnc_m = 0.0; zmns_m = 0.0; rmns_m = 0.0; zmnc_m = 0.0
         rmnc_m = rmnc_sav(:,vsurf)
         zmns_m = zmns_sav(:,vsurf)
         IF (lasym) THEN
            rmns_m = rmns_sav(:,vsurf)
            zmnc_m = zmnc_sav(:,vsurf)
         END IF
         ixm_m = xm
         ixn_m = -xn/nfp
         ! Check sign of Jacobian
         u0 = pi2/360
         v0 = 0.0; rb_ws = 0.0; zb_ws = 0.0
         isgn = 1
         DO mn = 1, mnmax_m
            IF (ixm_m(mn) == 0) CYCLE
            rb_ws = rb_ws + rmnc_m(mn)*cos(ixm_m(mn)*u0)&
                       + rmns_m(mn)*sin(ixm_m(mn)*u0)
            zb_ws = zb_ws + zmns_m(mn)*sin(ixm_m(mn)*u0)&
                       + zmnc_m(mn)*cos(ixm_m(mn)*u0)
         END DO
         IF (rb_ws*zb_ws < 0) isgn = -1
         IF (myid_sharmem==master) THEN
            DO v = 1, nv ! So we do nv=1 always
               DO u = 1, nu
                  u0 = alub*(u - 1)
                  v0 = alvb*(v - 1)  !!Np*(Real toroidal angle)
                  b1 = v0 - pi2/6
                  c1 = v0 + pi2/6
                  r1 = b1
                  ier = 0
                  DO
                    call fzero(Map_v, b1, c1, r1, re, ae, ier)
                    IF (ier <= 1) THEN
                       EXIT
                    ELSE IF (ier <5) THEN !F(B) = 0 (but we need to recalc B)
                       CALL MAp_v(b1)
                       EXIT
                    ELSE IF (ier ==5) THEN
                       r1=b1
                       b1=r1 - pi2/6
                       c1=r1 + pi2/6
                    END IF
                  END DO
                  IF (lvessel) THEN
                     b1 = pi2*xv(v)/nfp
                     v0 = v0/nfp
                     x0 = rreal(vsurf,u,v)*cos(v0)
                     y0 = rreal(vsurf,u,v)*sin(v0)
                     z0 = zreal(vsurf,u,v)
                     x2 = rb_ws*cos(v0)
                     y2 = rb_ws*sin(v0)
                     z2 = zb_ws
                     xw = 0; yw = 0; zw = 0
                     CALL collide(x0,y0,z0,x2,y2,z2,xw,yw,zw,lhit)
                     IF (lhit) THEN
                        rb_ws = sqrt(xw*xw+yw*yw)
                        zb_ws = zw
                     END IF
                  END IF
                  rreal(nrho,u,v) = rb_ws
                  zreal(nrho,u,v) = zb_ws
                  dr = (rreal(nrho,u,v) - rreal(vsurf,u,v))/(nrho-vsurf)
                  dz = (zreal(nrho,u,v) - zreal(vsurf,u,v))/(nrho-vsurf)
                  DO ik = 1, nrho-vsurf-1
                     rreal(vsurf+ik,u,v) = ik*dr+rreal(vsurf,u,v)
                     zreal(vsurf+ik,u,v) = ik*dz+zreal(vsurf,u,v)
                  END DO
               END DO
            END DO
         END IF
         ier = 0
         IF (ALLOCATED(rmnc_m)) DEALLOCATE(rmnc_m)
         IF (ALLOCATED(zmns_m)) DEALLOCATE(zmns_m)
         IF (ALLOCATED(ixm_m)) DEALLOCATE(ixm_m)
         IF (ALLOCATED(ixn_m)) DEALLOCATE(ixn_m)
         IF (ALLOCATED(rmns_m)) DEALLOCATE(rmns_m)
         IF (ALLOCATED(zmnc_m)) DEALLOCATE(zmnc_m)
      END IF
      IF (lvessel) CALL wall_free(ier)
      IF (myid_sharmem==master) THEN
         rreal(:,:,nv) = rreal(:,:,1)
         zreal(:,:,nv) = zreal(:,:,1)
         bureal(:,:,nv) = bureal(:,:,1)
         bvreal(:,:,nv) = bvreal(:,:,1)
         ! Redefine rho from [0, 1]
         !FORALL(ik = 1:nrho) rho(ik) = REAL(ik-1)/REAL(k-1)       ! Equidistant
         !FORALL(ik = 1:nrho) rho(ik) = (rreal(ik,1,1)-rreal(1,1,1))/(rreal(k,1,1)-rreal(1,1,1)) ! Rout
         FORALL(ik = 1:nrho) rho(ik) = sqrt((rreal(ik,1,1)-rreal(1,1,1))**2+(zreal(ik,1,1)-zreal(1,1,1))**2) ! rho_out
         rho = rho / rho(nrho)
      END IF
      ! Now do some cleanup
      ALLOCATE(rs(1:nrho,1:nu,1:nv),ru(1:nrho,1:nu,1:nv),rv(1:nrho,1:nu,1:nv),&
               zs(1:nrho,1:nu,1:nv),zu(1:nrho,1:nu,1:nv),zv(1:nrho,1:nu,1:nv),&
               g11(1:nrho,1:nu,1:nv),g12(1:nrho,1:nu,1:nv),g13(1:nrho,1:nu,1:nv),&
               g22(1:nrho,1:nu,1:nv),g23(1:nrho,1:nu,1:nv),g33(1:nrho,1:nu,1:nv),&
               detg(1:nrho,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'METRICS',ier)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
#endif
      CALL calc_metrics
      ! Deallocations
      IF (ALLOCATED(rho_vmec)) DEALLOCATE(rho_vmec)
      IF (ALLOCATED(rmnc_sav)) DEALLOCATE(rmnc_sav,zmns_sav)
      IF (ALLOCATED(rmns_sav)) DEALLOCATE(rmns_sav,zmnc_sav)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_init_wout
