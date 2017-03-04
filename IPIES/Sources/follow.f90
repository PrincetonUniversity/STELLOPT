!-----------------------------------------------------------------------
!     Subroutine:    follow
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/8/2011
!     Description:   This subroutine preforms field line following in
!                    the PIES background coordinates system.  After this
!                    subroutine the background coordinates are the 
!                    magnetic coordinates.
!-----------------------------------------------------------------------
      SUBROUTINE follow
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
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          ierr        Error flag
!          bcs1/2/3    Boundary Condition Arrays for EZspline
!          factor      Factor used to convert B^v B^phi
!          brho_real   bsreal/bvreal
!          btheta_real bureal/bvreal
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,u,v,ik,form,mn
      INTEGER :: bcs1(2),bcs2(2),bcs3(2)
      REAL(rprec) :: factor
      REAL(rprec) :: pi2,bxsi_sum,beta_sum
      REAL(rprec) :: bxsi_temp(1:nu_spline,1:nv_spline)
      REAL(rprec) :: beta_temp(1:nu_spline,1:nv_spline)
      REAL(rprec) :: brho_real(0:k,1:nu_spline,1:nv_spline)
      REAL(rprec) :: btheta_real(0:k,1:nu_spline,1:nv_spline)
      REAL(rprec) :: rhobth_real(0:k,1:nu_spline,1:nv_spline)
      REAL(rprec) :: bxsi_real(0:k,1:nu_spline,1:nv_spline)
      REAL(rprec) :: beta_real(0:k,1:nu_spline,1:nv_spline)
      REAL(rprec) :: bumnc_temp(1:mnmax,0:k)
      REAL(rprec) :: bumns_temp(1:mnmax,0:k)
      REAL(rprec) :: rbumnc_temp(1:mnmax,0:k)
      REAL(rprec) :: rbumns_temp(1:mnmax,0:k)
      REAL(rprec) :: rbumnc(1:mnmax,0:k)
      REAL(rprec) :: rbumns(1:mnmax,0:k)
      REAL(rprec) :: xu_spline(1:nu_spline)
      REAL(rprec) :: xv_spline(1:nv_spline)
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      form = 1
      pi2 = 8 * ATAN(1._rprec)
      factor = pi2/nfp
      hitsrf = k
      ! Remap the real values to u=-pi to pi
      FORALL(u=1:nu_spline) xu_spline(u)=REAL(u-1)/REAL(nu_spline-1)-0.5
      FORALL(v=1:nv_spline) xv_spline(v)=REAL(v-1)/REAL(nv_spline-1)
      ! Now recompute arrays on new grid
      IF (ALLOCATED(bsreal)) DEALLOCATE(bsreal,bureal,bvreal)
      IF (ALLOCATED(rreal)) DEALLOCATE(rreal,zreal)
      ALLOCATE(rreal(0:k,1:nu_spline,1:nv_spline),&
               zreal(0:k,1:nu_spline,1:nv_spline),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z (real) - FOLLOW',ier)
      ALLOCATE(bsreal(0:k,1:nu_spline,1:nv_spline),&
               bureal(0:k,1:nu_spline,1:nv_spline),&
               bvreal(0:k,1:nu_spline,1:nv_spline),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BS BU BV (real) - FOLLOW',ier)
      rreal = 0.0; zreal = 0.0; bsreal = 0.0; bureal = 0.0; bvreal = 0.0
      CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,rmnc,xm,xn,rreal,0,1)
      CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,zmns,xm,xn,zreal,1,0)
      CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bsmns,xm,xn,bsreal,1,0)
      CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bumnc,xm,xn,bureal,0,0)
      CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bvmnc,xm,xn,bvreal,0,0)
      IF (lasym) THEN
         CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,rmns,xm,xn,rreal,1,0)
         CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,zmnc,xm,xn,zreal,0,0)
         CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bsmnc,xm,xn,bsreal,0,0)
         CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bumns,xm,xn,bureal,1,0)
         CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bvmns,xm,xn,bvreal,1,0)   
      END IF
      ! First calculate some quantities for field-line following
      !   1) Create B^s/B^v and B^u/B^v in real space
      !   2) Transform B^u/B^v to Fourier Space
      !   3) Multiply by rho
      !   4) Transform back to real space
      brho_real   = bsreal/(bvreal*factor)                                       !B^s/B^v
      btheta_real = bureal/(bvreal*factor)                                       !B^u/B^v
      CALL uvtomn(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bumnc_temp,xm,xn,btheta_real,0,1)
      FORALL (ik=0:k) rbumnc_temp(:,ik) = rho(ik) * bumnc_temp(:,ik)
      rhobth_real = 0.0
      CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,rbumnc_temp,xm,xn,rhobth_real,0,1)
      IF (lasym) THEN
         CALL uvtomn(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,bumns_temp,xm,xn,btheta_real,1,0)
         FORALL (ik=0:k) rbumns_temp(:,ik) = rho(ik) * bumns_temp(:,ik)
         CALL mntouv(0,k,mnmax,nu_spline,nv_spline,xu_spline,xv_spline,rbumns_temp,xm,xn,rhobth_real,1,0)
      END IF
      ! Allocate helper arrays
      ALLOCATE(bxsi_array(0:k),beta_array(0:k),hitwal(0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BXSI_ARRAY',ier)
      ! Initialize Splines
      bcs1 = (/0,0/)   ! Not-a-knot
      bcs2 = (/-1,-1/) ! Periodic
      bcs3 = (/-1,-1/) ! Periodic
      CALL EZspline_init(R_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      CALL EZspline_init(Z_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      CALL EZspline_init(bxsi_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      CALL EZspline_init(beta_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      CALL EZspline_init(BS_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      CALL EZspline_init(BU_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      CALL EZspline_init(BV_spl,k+1,nu_spline,nv_spline,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/follow',ier)
      ! Setup Spline Axes
      R_spl%x1    = rho
      Z_spl%x1    = rho
      bxsi_spl%x1 = rho
      beta_spl%x2 = rho
      BS_spl%x1    = rho
      BU_spl%x1    = rho
      BV_spl%x1    = rho
      R_spl%x2    = xu_spline*pi2
      Z_spl%x2    = xu_spline*pi2
      bxsi_spl%x2 = xu_spline*pi2
      beta_spl%x2 = xu_spline*pi2
      BS_spl%x2    = xu_spline*pi2
      BU_spl%x2    = xu_spline*pi2
      BV_spl%x2    = xu_spline*pi2
      R_spl%x3    = xv_spline*pi2
      Z_spl%x3    = xv_spline*pi2
      bxsi_spl%x3 = xv_spline*pi2
      beta_spl%x3 = xv_spline*pi2
      BS_spl%x3    = xv_spline*pi2
      BU_spl%x3    = xv_spline*pi2
      BV_spl%x3    = xv_spline*pi2
      R_spl%isHermite    = 1
      Z_spl%isHermite    = 1
      bxsi_spl%isHermite = 1
      beta_spl%isHermite = 1
      BS_spl%isHermite    = 1
      BU_spl%isHermite    = 1
      BV_spl%isHermite    = 1
      ! Calculate Bxsi and Beta components
      ! bxsi_real = cos(theta)*brho_real-sin(theta)*rbth_real ~ BR/Bphi
      ! beta_real = sin(theta)*brho_real+cos(theta)*rbth_real ~ BZ/Bphi
      DO u = 1, nu_spline
         bxsi_real(:,u,:) = brho_real(:,u,:)*cos(pi2*xu_spline(u))-rhobth_real(:,u,:)*sin(pi2*xu_spline(u))
         beta_real(:,u,:) = brho_real(:,u,:)*sin(pi2*xu_spline(u))+rhobth_real(:,u,:)*cos(pi2*xu_spline(u))
      END DO
      bxsi_temp=bxsi_real(1,:,:)
      beta_temp=beta_real(1,:,:)
      DO v = 1, nv_spline
         bxsi_sum = SUM(bxsi_temp(:,v),DIM=1)
         beta_sum = SUM(beta_temp(:,v),DIM=1)
         bxsi_real(0,:,v) = bxsi_sum / nu_spline
         beta_real(0,:,v) = beta_sum / nu_spline
      END DO
      ! Calculate Splines
      CALL EZspline_setup(R_spl,rreal,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      CALL EZspline_setup(Z_spl,zreal,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      CALL EZspline_setup(bxsi_spl,bxsi_real,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      CALL EZspline_setup(beta_spl,beta_real,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      CALL EZspline_setup(BS_spl,bsreal,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      CALL EZspline_setup(BU_spl,bureal,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      CALL EZspline_setup(BV_spl,bvreal,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/follow',ier)
      IF (form == 1) THEN
         ! Calculate Magnetic Axis
         hitwal(:) = 0
         CALL magaxis
         WRITE(6,'(6X,f5.3,5X,f5.3)',advance='no') r_axis,z_axis
         CALL FLUSH(6)
         ! Follow Field lines
         CALL follow_fieldlines
         ! Calculate magnetic coordiantes
         CALL magco_recon
      ELSE IF (form == 2) THEN
         ! Calculate Magnetic Axis
         hitwal(:) = 0
         CALL magaxis
         WRITE(6,'(6X,f5.3,5X,f5.3)',advance='no') r_axis,z_axis
         ! Calculate invariant curve
         CALL invariant_surf
         PRINT *,'INVARIANT_SURF DONE'
         stop
      END IF
      ! Now calculate iota
      CALL pies_iotaslv
      DEALLOCATE(thetaln,philn)
      ! Now replace the background coordinates with the magnetic
      DEALLOCATE(rreal,zreal,bsreal,bureal,bvreal)
      ALLOCATE(rreal(0:k,1:nu,1:nv),zreal(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RREAL ZREAL (FOLLOW)',ier)
      ALLOCATE(bsreal(0:k,1:nu,1:nv),bureal(0:k,1:nu,1:nv),bvreal(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSREAL BUREAL BVREAL (FOLLOW)',ier)
      rreal  = rreal_magco
      zreal  = zreal_magco
      bsreal = bsreal_magco
      bureal = bureal_magco
      bvreal = bvreal_magco
      rmnc   = rmnc_magco
      zmns   = zmns_magco
      bsmns  = bsmns_magco
      bumnc  = bumnc_magco
      bvmnc  = bvmnc_magco
      IF (lasym) THEN
         rmns   = rmns_magco
         zmnc   = zmnc_magco
         bsmnc  = bsmnc_magco
         bumns  = bumns_magco
         bvmns  = bvmns_magco
         DEALLOCATE(rmns_magco,zmnc_magco,bsmnc_magco,bumns_magco,bvmns_magco)
      END IF
      ! Now calculated the metric elements
      CALL calc_metrics
      ! Deallocae Arrays
      CALL EZspline_free(R_spl,ier)
      CALL EZspline_free(Z_spl,ier)
      CALL EZspline_free(bxsi_spl,ier)
      CALL EZspline_free(beta_spl,ier)
      CALL EZspline_free(BS_spl,ier)
      CALL EZspline_free(BU_spl,ier)
      CALL EZspline_free(BV_spl,ier)
      DEALLOCATE(bxsi_array,beta_array, hitwal)
      DEALLOCATE(rreal_magco,zreal_magco,bsreal_magco,bureal_magco,bvreal_magco)
      DEALLOCATE(rmnc_magco,zmns_magco,bsmns_magco,bumnc_magco,bvmnc_magco)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE follow
