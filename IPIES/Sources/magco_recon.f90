!-----------------------------------------------------------------------
!     Subroutine:    magco_recon
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/6/2011
!     Description:   This subroutine reconstructs the magnetic
!                    coordinates from the field line data.
!-----------------------------------------------------------------------
      SUBROUTINE magco_recon
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_magco
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          ik          Radial dummy index
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          ierr        Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,u,v,ik,j,mn,ngood,nbad,ikg,ikb,npoints,ik2,ik3
      INTEGER :: bad_surf(0:k)
      INTEGER :: bcs1(2)
      INTEGER, ALLOCATABLE :: bad_index(:), good_index(:)
      REAL(rprec) :: pi2, x1, y1, rho1, x2, y2, rho2, pi2_inv,dxu
      REAL(rprec) :: rho_local(0:k), theta_local(0:k), phi_local(0:k)
      REAL(rprec) :: theta_old(0:k)
      REAL(rprec) :: xsi_local(0:k), eta_local(0:k)
      REAL(rprec) :: R_local(0:k), Z_local(0:k), BS_local(0:k), BU_local(0:k), BV_local(0:k)
      REAL(rprec), ALLOCATABLE :: rho_good(:), rho_bad(:)
      REAL(rprec), ALLOCATABLE :: r_temp(:), z_temp(:), bs_temp(:),bu_temp(:),bv_temp(:)
      REAL(rprec), ALLOCATABLE :: fnuv(:)
      TYPE(EZspline1_r8) :: rm_spl,zm_spl,bsm_spl,bum_spl,bvm_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      pi2 = 8 * ATAN(1._rprec)
      ALLOCATE(rholn(0:k,0:nintw),thetaln(0:k,0:nintw),philn(0:k,0:nintw), STAT=ier)
      ! Extract theta,phi,R,and Z of the Field lines
      theta_old = 0.0
      xsi_local   = xsiln(0:k,0)
      eta_local   = etaln(0:k,0)
      rho_local   = sqrt(xsi_local*xsi_local+eta_local*eta_local)
      theta_local = 0.0
      thetaln(0:k,0) = theta_local
      philn(0:k,0)   = 0.0
      ! Now do rest of points
      DO j = 1, nintw
         xsi_local   = xsiln(0:k,j)-xsiln(0,j)
         eta_local   = etaln(0:k,j)-etaln(0,j)
         rho_local   = sqrt(xsi_local*xsi_local+eta_local*eta_local)
         theta_local = 0.0
         WHERE(rho_local > 0.0) theta_local=atan2(eta_local,xsi_local)
         WHERE (theta_local < 0.0) theta_local = theta_local + pi2
         rholn(0:k,j)   = rho_local
         thetaln(0:k,j) = theta_local
         philn(0:k,j) = j*dphi
      END DO
      !DEALLOCATE(xsiln,etaln)  ! For now we don't deallocate them so we can write them in write_pies_hdf5
      ! Fourier transform rho(theta)
      mnmax_polar = m + 1
      nu_polar = nintw+1
      nv_polar = 1
      xv_polar(1) = 0.0
      ALLOCATE(xu_polar(1:nintw+1),xm_polar(1:mnmax_polar),xn_polar(1:mnmax_polar),STAT=ier)
      ALLOCATE(rholn_polar(0:k,1:nintw+1,1:1),STAT=ier)
      FORALL(u=0:m) xm_polar(u+1) = u
      xn_polar(:) = 0
      rholn_polar(:,:,1) = rholn(:,:)
      pi2_inv=1/pi2
      ! Construct magnetic coordinates using rho(theta)
      ALLOCATE(rho_polar(0:k,1:nu,1:nv),STAT=ier)
      dxu = 0.5*pi2/nu
      DO ik = 0, k-1
         DO u = 1, nu
            xu_polar = thetaln(ik,:)
            rho_polar(ik,u,:)=SUM(rholn(ik,:),MASK = ((xu_polar > pi2*xu(u)-dxu) .and. (xu_polar < pi2*xu(u)+dxu)))
            npoints=COUNT((xu_polar > pi2*xu(u)-dxu) .and. (xu_polar < pi2*xu(u)+dxu))
            IF (npoints > 0) THEN
               rho_polar(ik,u,:)=rho_polar(ik,u,:)/npoints
            ELSE
               rho_polar(ik,u,:)=1.0
            END IF
         END DO
      END DO
      rho_polar(0,:,:)=0.0
      rho_polar(k,:,:)=1.0
      WHERE(rho_polar > 1.0) rho_polar=1.0
      WHERE(rho_polar < 0.0) rho_polar=0.0
      WRITE(63,*) rho_polar(:,:,1)
      DEALLOCATE(rholn)
      DEALLOCATE(xu_polar,xm_polar,xn_polar)
      DEALLOCATE(rholn_polar)
      ! Now extract spline quantities
      ALLOCATE(rreal_magco(0:k,1:nu,1:nv),zreal_magco(0:k,1:nu,1:nv),STAT=ier)
      ALLOCATE(bsreal_magco(0:k,1:nu,1:nv),bureal_magco(0:k,1:nu,1:nv),bvreal_magco(0:k,1:nu,1:nv),STAT=ier)
      DO v = 1, nv
         DO u = 1, nu
            rho_local(:) = rho_polar(:,u,v)
            theta_local(:) = xu(u)*pi2
            phi_local(:)   = xv(v)*pi2
            CALL EZspline_interp(R_spl,k+1,rho_local,theta_local,phi_local,R_local,ier)
            CALL EZspline_interp(Z_spl,k+1,rho_local,theta_local,phi_local,Z_local,ier)
            CALL EZspline_interp(BS_spl,k+1,rho_local,theta_local,phi_local,BS_local,ier)
            CALL EZspline_interp(BU_spl,k+1,rho_local,theta_local,phi_local,BU_local,ier)
            CALL EZspline_interp(BV_spl,k+1,rho_local,theta_local,phi_local,BV_local,ier)
            rreal_magco(0:k,u,v) = R_local
            zreal_magco(0:k,u,v) = Z_local
            bsreal_magco(0:k,u,v) = BS_local
            bureal_magco(0:k,u,v) = BU_local
            bvreal_magco(0:k,u,v) = BV_local
         END DO
      END DO
      DEALLOCATE(rho_polar)
      ! Now Fourier Transform surfaces
      ALLOCATE(rmnc_magco(1:mnmax,0:k),zmns_magco(1:mnmax,0:k),STAT=ier)
      ALLOCATE(bsmns_magco(1:mnmax,0:k),bumnc_magco(1:mnmax,0:k),bvmnc_magco(1:mnmax,0:k),STAT=ier)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmnc_magco,xm,xn,rreal_magco,0,1)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmns_magco,xm,xn,zreal_magco,1,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bsmns_magco,xm,xn,bsreal_magco,1,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bumnc_magco,xm,xn,bureal_magco,0,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bvmnc_magco,xm,xn,bvreal_magco,0,0)
      rmnc_magco(:,k) = rmnc(:,k)
      zmns_magco(:,k) = zmns(:,k)
      bsmns_magco(:,k) = bsmns(:,k)
      bumnc_magco(:,k) = bumnc(:,k)
      bvmnc_magco(:,k) = bvmnc(:,k)
      IF (lasym) THEN
         ALLOCATE(rmns_magco(1:mnmax,0:k),zmnc_magco(1:mnmax,0:k),STAT=ier)
         ALLOCATE(bsmnc_magco(1:mnmax,0:k),bumns_magco(1:mnmax,0:k),bvmns_magco(1:mnmax,0:k),STAT=ier)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmns_magco,xm,xn,rreal_magco,1,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmnc_magco,xm,xn,zreal_magco,0,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bsmnc_magco,xm,xn,bsreal_magco,0,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bumns_magco,xm,xn,bureal_magco,1,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bvmns_magco,xm,xn,bvreal_magco,1,0)
         rmns_magco(:,k) = rmns(:,k)
         zmnc_magco(:,k) = zmnc(:,k)
         bsmnc_magco(:,k) = bsmnc(:,k)
         bumns_magco(:,k) = bumns(:,k)
         bvmns_magco(:,k) = bvmns(:,k)
      END IF
      ! Now get Fourier Transformed Magnetic Coordinates
      rreal_magco = 0.0
      zreal_magco = 0.0
      bsreal_magco = 0.0
      bureal_magco = 0.0
      bvreal_magco = 0.0
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,rmnc_magco,xm,xn,rreal_magco,0,1)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,zmns_magco,xm,xn,zreal_magco,1,0)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bsmns_magco,xm,xn,bsreal_magco,1,0)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bumnc_magco,xm,xn,bureal_magco,0,0)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bvmnc_magco,xm,xn,bvreal_magco,0,0)
      IF (lasym) THEN
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,rmns_magco,xm,xn,rreal_magco,1,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,zmnc_magco,xm,xn,zreal_magco,0,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bsmnc_magco,xm,xn,bsreal_magco,0,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bumns_magco,xm,xn,bureal_magco,1,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bvmns_magco,xm,xn,bvreal_magco,1,0)
      END IF
      ! Search for overlapping surfaces and mark as bad
      bad_surf = 0
      ik = 1
      DO WHILE (ik .lt. k)
         DO u = 1, nu
            DO v = 1, nv
               x1 = rreal_magco(ik,u,v)-rreal_magco(0,u,v)
               y1 = zreal_magco(ik,u,v)-zreal_magco(0,u,v)
               rho1 = sqrt(x1*x1+y1*y1)
               x2 = rreal_magco(ik+1,u,v)-rreal_magco(0,u,v)
               y2 = zreal_magco(ik+1,u,v)-zreal_magco(0,u,v)
               rho2 = sqrt(x2*x2+y2*y2)
               IF (rho1 .ge. rho2) THEN
                  bad_surf(ik-1) = 1
                  bad_surf(  ik) = 1
                  bad_surf(ik+1) = 1
               END IF
            END DO
         END DO
         ik = ik + 1
      END DO
      ! Throw out single surfaces
      DO ik = 1, k - 1
         IF ((bad_surf(ik-1) == 1) .and. (bad_surf(ik+1) == 1)) bad_surf(ik) = 1
      END DO
      WHERE(hitwal == 1) bad_surf = 1
      bad_surf(0) = 0
      bad_surf(k) = 0
      nbad  = COUNT(bad_surf==1)
      ngood = k-nbad+1
      WRITE(6,'(3X,i5,3X)',ADVANCE='yes') nbad
      ! Create good/bad rho surface arrays
      ALLOCATE(rho_good(1:ngood),rho_bad(1:nbad),good_index(ngood),bad_index(nbad),STAT=ier)
      ALLOCATE(r_temp(1:ngood),z_temp(1:ngood),STAT=ier)
      ALLOCATE(bs_temp(1:ngood),bu_temp(1:ngood),bv_temp(1:ngood),STAT=ier)
      ikg=1
      ikb=1
      bcs1=(/0,0/)
      DO ik = 0, k
         IF (bad_surf(ik)==1) THEN
            rho_bad(ikb)=rho(ik)
            bad_index(ikb) = ik
            ikb=ikb+1
         ELSE
            rho_good(ikg)=rho(ik)
            good_index(ikg)=ik
            ikg=ikg+1
         END IF
      END DO
      ! Now create splines over radial indicies
      DO u = 1, nu
         DO v = 1, nv
            ! Init Splines
            CALL EZspline_init(rm_spl,ngood,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/magco_recon',ier)
            CALL EZspline_init(zm_spl,ngood,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/magco_recon',ier)
            CALL EZspline_init(bsm_spl,ngood,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/magco_recon',ier)
            CALL EZspline_init(bum_spl,ngood,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/magco_recon',ier)
            CALL EZspline_init(bvm_spl,ngood,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/magco_recon',ier)   
            rm_spl%isHermite = 1
            zm_spl%isHermite = 1
            bsm_spl%isHermite = 1
            bum_spl%isHermite = 1
            bvm_spl%isHermite = 1
            rm_spl%x1        = rho_good
            zm_spl%x1        = rho_good
            bsm_spl%x1        = rho_good
            bum_spl%x1        = rho_good
            bvm_spl%x1        = rho_good
            ! Get good surfaces locations as knots
            ikg = 1
            DO ik =0, k
               IF (bad_surf(ik) == 0) THEN
                  r_temp(ikg) = rreal_magco(ik,u,v)
                  z_temp(ikg) = zreal_magco(ik,u,v)
                  bs_temp(ikg) = bsreal_magco(ik,u,v)
                  bu_temp(ikg) = bureal_magco(ik,u,v)
                  bv_temp(ikg) = bvreal_magco(ik,u,v)
                  ikg=ikg+1
               END IF
            END DO
            CALL EZspline_setup(rm_spl,r_temp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/magco_recon',ier)
            CALL EZspline_setup(zm_spl,z_temp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/magco_recon',ier)
            CALL EZspline_setup(bsm_spl,bs_temp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/magco_recon',ier)
            CALL EZspline_setup(bum_spl,bu_temp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/magco_recon',ier)
            CALL EZspline_setup(bvm_spl,bv_temp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/magco_recon',ier)
            DO ik = 1, nbad
               CALL EZspline_interp(rm_spl,rho_bad(ik),rreal_magco(bad_index(ik),u,v),ier)
               CALL EZspline_interp(zm_spl,rho_bad(ik),zreal_magco(bad_index(ik),u,v),ier)
               CALL EZspline_interp(bsm_spl,rho_bad(ik),bsreal_magco(bad_index(ik),u,v),ier)
               CALL EZspline_interp(bum_spl,rho_bad(ik),bureal_magco(bad_index(ik),u,v),ier)
               CALL EZspline_interp(bvm_spl,rho_bad(ik),bvreal_magco(bad_index(ik),u,v),ier)
            END DO
            CALL EZspline_free(rm_spl,ier)
            CALL EZspline_free(zm_spl,ier)
            CALL EZspline_free(bsm_spl,ier)
            CALL EZspline_free(bum_spl,ier)
            CALL EZspline_free(bvm_spl,ier)
         END DO
      END DO
      ! Now Fourier Transform surfaces
      rmnc_magco = 0.0
      zmns_magco = 0.0
      bsmns_magco = 0.0
      bumnc_magco = 0.0
      bvmnc_magco = 0.0
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmnc_magco,xm,xn,rreal_magco,0,1)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmns_magco,xm,xn,zreal_magco,1,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bsmns_magco,xm,xn,bsreal_magco,1,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bumnc_magco,xm,xn,bureal_magco,0,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bvmnc_magco,xm,xn,bvreal_magco,0,0)
      IF (lasym) THEN
         rmns_magco = 0.0
         zmnc_magco = 0.0
         bsmnc_magco = 0.0
         bumns_magco = 0.0
         bvmns_magco = 0.0
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmns_magco,xm,xn,rreal_magco,1,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmnc_magco,xm,xn,zreal_magco,0,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bsmnc_magco,xm,xn,bsreal_magco,0,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bumns_magco,xm,xn,bureal_magco,1,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bvmns_magco,xm,xn,bvreal_magco,1,0)
      END IF
      ! DEALLOCATE ARRAYS
      DEALLOCATE(rho_good,rho_bad,good_index,bad_index)
      DEALLOCATE(r_temp,z_temp)
      DEALLOCATE(bs_temp,bu_temp,bv_temp)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE magco_recon
