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
      SUBROUTINE pies_init_wout
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, nfp_in => nfp, lasym_in => lasym, &
                      lfreeb_in => lfreeb, mnmax_in => mnmax, &
                      rmnc_in => rmnc, zmns_in => zmns, &
                      rmns_in => rmns, zmnc_in => zmnc, &
                      bsupumnc_in => bsupumnc, bsupumns_in => bsupumns,&
                      bsupvmnc_in => bsupvmnc, bsupvmns_in => bsupvmns,&
                      xn_in => xn, xm_in => xm, mgrid_file_in => mgrid_file,&
                      nextcur_in => nextcur, extcur_in => extcur, &
                      iotaf_in => iotaf
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile, ONLY: press, iprime, p_spl, ip_spl,torflux_edge,&
                              curtor, p_cubspl, ip_cubspl, iotaf
      USE pies_fieldlines, ONLY: dkmin
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          uv          Dummy index over real-space
!          ik          Dummy index over radial surfaces
!          i           Dummy index
!          j           Dummy index
!          bcs1        Boundary Condition Array for EZspline
!          val1        Dummy value used in boundary extension
!          val2        Dummy value used in boundary extension
!          dval        Dummy value used in boundary extension
!          rhomax      Twice the maximum VMEC rho in each v section
!          rreal_in    VMEC R in (s,u,v) space
!          zreal_in    VMEC Z in (s,u,v) space
!          bureal_in   VMEC Bsupu in (s,u,v) space
!          bvreal_in   VMEC Bsupv in (s,u,v) space
!          rho_vmec    VMEC radial grid
!          f_spl       Dummy spline object for spline to PIES grid
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, ier, mn, im, in, ik , i, j
      INTEGER :: bcs1(2)
      REAL(rprec) :: val1, val2, dval, rhomax,pi2, scale, theta,r1, z1,&
                     r2, z2, theta2, dtheta
      REAL(rprec), ALLOCATABLE :: rreal_in(:,:,:), zreal_in(:,:,:)
      REAL(rprec), ALLOCATABLE :: bureal_in(:,:,:), bvreal_in(:,:,:)
      REAL(rprec), ALLOCATABLE :: rho_vmec(:), rho_vmec_half(:), rho_ext(:), ftemp(:), ftemp_half(:)
      REAL(rprec), ALLOCATABLE :: extrap_surfs(:,:,:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumnc_temp(:,:), bsupvmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bumnc_temp(:,:), bvmnc_temp(:,:)
      TYPE(EZspline1_r8) :: f_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/0,0/)
      pi2 = 8 * ATAN(1._rprec)
      CALL read_wout_file(TRIM(id_string),ierr=ier)
      IF (ier /= 0) CALL handle_err(VMEC_WOUT_ERR,id_string,ier)
      ! Now set some background values
      lasym = lasym_in
      lfreeb = lfreeb_in
      IF (lfreeb) THEN
         mgrid_file = TRIM(mgrid_file_in)
         nv_vmec = 32 ! Will need to read the INDATA namelist to get this
         nextcur = nextcur_in
         ALLOCATE(extcur(1:nextcur),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'EXTCUR',ier)
         extcur(1:nextcur) = extcur_in(1:nextcur_in)
      ENDIF
      IF (isigng >0) THEN
         signgs = -1
      ELSE
         signgs = 1
      END IF
      nfp = nfp_in
      m = mpol-1
      n = ntor
      if (k == -1) k = ns
      nu = 4*m+1
      nv = 4*n+1
      mnmax = (m+1)*(2*n+1)
      nuv = nu*nv
      nuvp = nuv*nfp
      vsurf = k - extsurfs
      iotamn = MINVAL(iotaf_in(:))
      iotamn = iotamn - iotamn*0.4
      iotamx = MAXVAL(iotaf_in(:))
      iotamx = iotamx + iotamx*0.4
      dkmin  = MIN(ABS(dkmin),ABS(iotamx))
      dkmin  = MIN(ABS(dkmin),ABS(iotamn))
      torflux_edge = phi(ns)
      curtor = Itor
      rbtor_pies=RBtor
      id_string = TRIM(input_extension)
      if (lverb) THEN
         write(6,*) '-----VMEC File Parameters-----'
         write(6,'(A,A)') '    file: ',TRIM(id_string)
         write(6,'(A,I3,A,I3)') '       m: ',m,'   nu: ',nu
         write(6,'(A,I3,A,I3)') '       n: ',n,'   nv: ',nv
         write(6,'(A,I4,A,I4,A,I5)') '   mnmax: ',mnmax,'  nuv: ',nuv,'   nuvp: ',nuvp
         write(6,'(A,I3)') '     nfp: ',nfp
         write(6,'(A,I3)') '      ns: ',ns
         write(6,'(A,L3)') '  lfreeb: ',lfreeb
         write(6,'(A,F6.3,A,F6.3,A)')    '    iota: [',iotamn,',',iotamx,']'  
         write(6,'(A,F6.3)')             'torflux_edge: ',torflux_edge 
         IF (ABS(curtor) .ge. 1.0E6) THEN
            write(6,'(A,F8.3,A)')             'Total Current: ',curtor/1.0E6,' [MA]'
         ELSE IF (ABS(curtor) .ge. 1.0E3) THEN
            write(6,'(A,F8.3,A)')             'Total Current: ',curtor/1.0E3,' [kA]'
         ELSE
            write(6,'(A,F8.3,A)')             'Total Current: ',curtor,' [A]'
         END IF
      END IF
      ! Now allocate 1D variables
      ALLOCATE(xm(1:mnmax),xn(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM XN',ier)
      ALLOCATE(rho(0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO',ier)
      ALLOCATE(rho_vmec(1:ns),rho_vmec_half(1:ns),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO VMEC',ier)
      ALLOCATE(xu(1:nu),xv(1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XU XV',ier)
      ! Now allocate 2D variables
      ALLOCATE(rmnc(1:mnmax,0:k),zmns(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      ALLOCATE(bsmns(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMNS',ier)
      ALLOCATE(bumnc(1:mnmax,0:k),bvmnc(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNC BVMNC',ier)
      IF (lasym) THEN
         ALLOCATE(rmns(1:mnmax,0:k),zmnc(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Rs Zc',ier)
         ALLOCATE(bsmnc(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMNC',ier)
         ALLOCATE(bumns(1:mnmax,0:k),bvmns(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNS BVMNS',ier)
      END IF
      ALLOCATE(rreal(0:k,1:nu,1:nv),zreal(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z (real)',ier)
      ALLOCATE(bsreal(0:k,1:nu,1:nv),bureal(0:k,1:nu,1:nv),bvreal(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'B (real)',ier)
      ! Now initizlie modes
      mn = 1
      DO i=0,m
         DO j=-n,n
           xm(mn) = i
           xn(mn) = j
           mn = mn + 1
         END DO
      END DO
      FORALL(ik = 0:k) rho(ik) = REAL(ik) / REAL(vsurf)
      rho = rho * rho
      DO ik = 1,ns
         rho_vmec(ik) = REAL(ik-1)/REAL(ns-1)
         rho_vmec_half(ik) = (REAL(ik)-0.5)/REAL(ns)
      END DO
      rho_vmec_half(1) = 0.0
      ! Now initialize xu and xv
      FORALL(i=1:nu) xu(i) = REAL(i-1)/REAL(nu)
      FORALL(i=1:nv) xv(i) = REAL(i-1)/REAL(nv)
      ! Extrapolate Fields to Full mesh
      !bsupumnc_in(:,1) = 0.0
      !bsupumnc_in(:,1) = 1.5*bsupumnc_in(:,2) - 0.5*bsupumnc_in(:,3)
      !bsupvmnc_in(:,1) = 1.5*bsupvmnc_in(:,2) - 0.5*bsupvmnc_in(:,3)
      !bsupumnc_in(:,ns) = 1.5*bsupumnc_in(:,ns-1) - 0.5 * bsupumnc_in(:,ns-2)
      !bsupvmnc_in(:,ns) = 1.5*bsupvmnc_in(:,ns-1) - 0.5 * bsupvmnc_in(:,ns-2)
      !DO ik = 2,ns-1
      !   bsupumnc_in(:,ik)=0.5*(bsupumnc_in(:,ik)+bsupumnc_in(:,ik+1))
      !   bsupvmnc_in(:,ik)=0.5*(bsupvmnc_in(:,ik)+bsupvmnc_in(:,ik+1))
      !END DO
      ! Now interpolate to PIES coordinate in flux space
      ALLOCATE(ftemp(1:ns),ftemp_half(1:ns),STAT=ier)
      ALLOCATE(rmnc_temp(1:mnmax_in,0:k),zmns_temp(1:mnmax_in,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      ALLOCATE(bsupumnc_temp(1:mnmax_in,0:k),bsupvmnc_temp(1:mnmax_in,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNC BVMNC',ier)
      DO mn = 1, mnmax_in
         !RMNC
         ! We want the fourier coefficients to go to zero linearly in toroidal flux
         ! so we multiply by flux**(-m/2 + 1)
         ! we then interpolate in flux space rho**2
         ! Finally we remove the factor in rho (PIES) space.
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec
         ftemp(:) = rmnc_in(mn,:)*rho_vmec**(-xm_in(mn)/2+1)
         IF (xm(mn) /= 0.0) ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
         DO ik = 1, vsurf
            CALL EZspline_interp(f_spl,rho(ik),rmnc_temp(mn,ik),ier)
         END DO
         rmnc_temp(mn,:) = rmnc_temp(mn,:)*SQRT(rho)**(xm_in(mn)-2)
         rmnc_temp(mn,0) = rmnc_in(mn,1)
         CALL EZspline_free(f_spl,ier)
         ! ZMNS
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec
         ftemp(:) = zmns_in(mn,:)*rho_vmec**(-xm_in(mn)/2+1)
         IF (xm(mn) /= 0.0) ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
         DO ik = 1, vsurf
            CALL EZspline_interp(f_spl,rho(ik),zmns_temp(mn,ik),ier)
         END DO
         zmns_temp(mn,:) = zmns_temp(mn,:)*SQRT(rho)**(xm_in(mn)-2)
         zmns_temp(mn,0) = zmns_in(mn,1)
         CALL EZspline_free(f_spl,ier)
         
         !BSUPUMNC
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec_half
         IF (xm(mn) /= 0) THEN
            ftemp_half(1:ns) = bsupumnc_in(mn,1:ns)*rho_vmec_half(1:ns)**(-xm_in(mn)/2+1)
            ftemp_half(1) = 0.0
         ELSE
            ftemp_half(1:ns) = bsupumnc_in(mn,1:ns)
            ftemp_half(1) = (4.0 * bsupumnc_in(mn,2) - bsupumnc_in(mn,3))/3.0
         END IF
         CALL EZspline_setup(f_spl,ftemp_half,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
         DO ik = 0, vsurf-1
            CALL EZspline_interp(f_spl,rho(ik),bsupumnc_temp(mn,ik),ier)
         END DO
         IF (xm(mn) /= 0) THEN
            bsupumnc_temp(mn,:) = bsupumnc_temp(mn,:)*SQRT(rho)**(xm_in(mn)-2)
         END IF
         bsupumnc_temp(mn,vsurf) = 1.5 * bsupumnc_in(mn,ns) - 0.5 * bsupumnc_in(mn,ns-1)
         CALL EZspline_free(f_spl,ier)
         
         
         !BSUPVMNC
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec_half
         IF (xm(mn) /= 0) THEN
            ftemp_half(1:ns) = bsupvmnc_in(mn,1:ns)*rho_vmec_half(1:ns)**(-xm_in(mn)/2+1)
            ftemp_half(1) = 0.0
         ELSE
            ftemp_half(1:ns) = bsupvmnc_in(mn,1:ns)
            ftemp_half(1) = (4.0 * bsupvmnc_in(mn,2) - bsupvmnc_in(mn,3))/3.0
         END IF
         CALL EZspline_setup(f_spl,ftemp_half,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
         DO ik = 0, vsurf-1
            CALL EZspline_interp(f_spl,rho(ik),bsupvmnc_temp(mn,ik),ier)
         END DO
         IF (xm(mn) /= 0) THEN
            bsupvmnc_temp(mn,:) = bsupvmnc_temp(mn,:)*SQRT(rho)**(xm_in(mn)-2)
         END IF
         bsupvmnc_temp(mn,vsurf) = 1.5 * bsupvmnc_in(mn,ns) - 0.5 * bsupvmnc_in(mn,ns-1)
         CALL EZspline_free(f_spl,ier)
      END DO
      DEALLOCATE(ftemp, ftemp_half)
      ! Now copy from rmnc_temp to rmnc  (mu-nv) -> (mu+nv)
      rmnc = 0.0
      zmns = 0.0
      bumnc = 0.0
      bvmnc = 0.0
      xn_in = - xn_in/nfp
      DO mn = 1, mnmax
         DO j = 1, mnmax_in
            IF ((xm(mn) .eq. xm_in(j)) .and. (xn(mn) .eq. xn_in(j))) THEN 
               rmnc(mn,0:vsurf) =  rmnc_temp(j,0:vsurf)
               zmns(mn,0:vsurf) =  zmns_temp(j,0:vsurf)
               bumnc(mn,0:vsurf) =  bsupumnc_temp(j,0:vsurf)
               bvmnc(mn,0:vsurf) =  bsupvmnc_temp(j,0:vsurf)
               !bumnc(mn,0:k) =  bsupumnc_temp(j,0:k)
               !bvmnc(mn,0:k) =  bsupvmnc_temp(j,0:k)
            END IF
         END DO
      END DO
      DEALLOCATE(rmnc_temp,zmns_temp)
      DEALLOCATE(bsupumnc_temp, bsupvmnc_temp)
      ! Now initialize profile splines
      ALLOCATE(press(0:k),iprime(0:k),iotaf(1:ns),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'press iprime iotaf',ier)
      CALL EZspline_init(p_spl,ns,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
      CALL EZspline_init(ip_spl,ns,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
      p_spl%isHermite = 1
      ip_spl%isHermite = 1                                               
      p_spl%x1 = phi/phi(ns)                                            !Spline over normalized toroidal flux
      ip_spl%x1 = phi/phi(ns)                                           !Spline over normalized toroidal flux
      CALL EZspline_setup(p_spl,presf,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
      CALL EZspline_setup(ip_spl,ABS(jcurv)/MAXVAL(ABS(jcurv)),ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
      DO ik = 0, vsurf
         CALL EZspline_interp(p_spl,rho(ik),press(ik),ier)
         CALL EZspline_interp(ip_spl,rho(ik),iprime(ik),ier)
      END DO
      press(vsurf+1:k)  = press(vsurf)
      iprime(vsurf+1:k) = 0.0
      iotaf(:) = iotaf_in(:)
      ! Now consturct PIES cubic splines using cubspl
      ALLOCATE(p_cubspl(1:4,1:ns),ip_cubspl(1:4,1:ns),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'CUBIC SPLINES',ier)
      p_cubspl  = 0.0
      ip_cubspl = 0.0
      FORALL(ik=1:ns) p_cubspl(1,ik) = presf(ik)
      FORALL(ik=1:ns) ip_cubspl(1,ik) = ABS(jcurv(ik))/MAXVAL(ABS(jcurv))
      IF (ANY(presf .gt. 0)) THEN
         CALL cubspl(rho_vmec,p_cubspl,ns,0,0)
      END IF
      IF (ABS(curtor) .gt. 0) THEN
         CALL cubspl(rho_vmec,ip_cubspl,ns,0,0)
      ElSE
         ip_cubspl = 0.0
      END IF
      ! Now transform to realspace
      rreal=0.0
      zreal=0.0
      bsreal=0.0
      bureal=0.0
      bvreal=0.0
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,rmnc,xm,xn,rreal,0,1)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,zmns,xm,xn,zreal,1,0)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bumnc,xm,xn,bureal,0,0)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bvmnc,xm,xn,bvreal,0,0)
      IF (lasym) THEN
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,rmns,xm,xn,rreal,1,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,zmnc,xm,xn,zreal,0,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bumns,xm,xn,bureal,1,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,bvmns,xm,xn,bvreal,1,0)
      END IF
      ! Free VMEC wout memory
      CALL read_wout_deallocate
      ! ****************************************************************
      ! ** Beyond this point the VMEC wout variables are deallocated  **
      ! ****************************************************************
      ! Calculate the external surfaces
      FORALL(ik = 0:k) rho(ik) = REAL(ik) / REAL(vsurf)
      IF (extsurfs > 0) THEN
         !scale = 2.0
         scale = 3.0*(MAXVAL(rho)-1.0)+1.0
         ALLOCATE(ftemp(0:vsurf+1))
         ALLOCATE(rmnc_temp(1:mnmax,0:k),zmns_temp(1:mnmax,0:k))
         ALLOCATE(rho_ext(0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO_EXT',ier)
         rho_ext = rho
         rho_ext(vsurf+1) = scale
         !rho_ext(vsurf+1) = 1+2*(rho(k)-1)
         DO j = 1, nv
            rhomax = scale * MAXVAL( SQRT( (rreal(vsurf,:,j) - rreal(0,1,j))**2 &
                                        +(zreal(vsurf,:,j) - zreal(0,1,j))**2))
            DO i = 1, nu
               r1 = rreal(vsurf,i,j)-rreal(0,1,j)
               z1 = zreal(vsurf,i,j)-zreal(0,1,j)
               r2 = rreal(vsurf-1,i,j)-rreal(0,1,j)
               z2 = zreal(vsurf-1,i,j)-zreal(0,1,j)
               IF ((r1 /=0 ) .and. (z1 /= 0)) THEN
                  theta = ATAN2(z1,r1)
                  IF (theta < 0) theta = theta + pi2
                  theta2 = ATAN2(z2,r2)
                  IF (theta2 < 0) theta2 = theta2 + pi2
                  theta = theta + (theta-pi2*xu(i))
               ELSE
                  theta = 0.0
               END IF
               rreal(vsurf+1,i,j) = rhomax * COS(pi2*xu(i))+rreal(0,1,j)
               zreal(vsurf+1,i,j) = rhomax * SIN(pi2*xu(i))+zreal(0,1,j)
               ! Spline R
               CALL EZspline_init(f_spl,vsurf+2,bcs1,ier)
               IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
               f_spl%isHermite = 1
               f_spl%x1 = rho_ext(0:vsurf+1)
               CALL EZspline_setup(f_spl,rreal(0:vsurf+1,i,j),ier)
               IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
               DO ik = vsurf+1, k
                  CALL EZspline_interp(f_spl,rho(ik),rreal(ik,i,j),ier)
               END DO
               CALL EZspline_free(f_spl,ier)
               ! Spline Z
               CALL EZspline_init(f_spl,vsurf+2,bcs1,ier)
               IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
               f_spl%isHermite = 1
               f_spl%x1 = rho_ext(0:vsurf+1)
               CALL EZspline_setup(f_spl,zreal(0:vsurf+1,i,j),ier)
               IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
               DO ik = vsurf+1, k
                  CALL EZspline_interp(f_spl,rho(ik),zreal(ik,i,j),ier)
               END DO
               CALL EZspline_free(f_spl,ier)
            END DO
         END DO
         rmnc_temp = 0; zmns_temp = 0
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmnc_temp,xm,xn,rreal,0,1)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmns_temp,xm,xn,zreal,1,0)
         rmnc(:,vsurf+1:k) = rmnc_temp(:,vsurf+1:k)
         zmns(:,vsurf+1:k) = zmns_temp(:,vsurf+1:k)
         IF (lasym) THEN
!            rmns = 0; zmnc = 0
            CALL uvtomn(vsurf+1,k,mnmax,nu,nv,xu,xv,rmns,xm,xn,rreal,0,0)
            CALL uvtomn(vsurf+1,k,mnmax,nu,nv,xu,xv,zmnc,xm,xn,zreal,1,0)
         END IF
         DEALLOCATE(ftemp)
         DEALLOCATE(rho_ext)
         DEALLOCATE(rmnc_temp,zmns_temp)
      END IF
      !DO ik = 0, k
      !   WRITE(27,*) rmnc(:,ik)
      !   WRITE(28,*) zmns(:,ik)
      !END DO
      ! Now we define our metric derivaitives
      ALLOCATE(rs(0:k,1:nu,1:nv),zs(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RS ZS',ier)
      ALLOCATE(ru(0:k,1:nu,1:nv),zu(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RU ZU',ier)
      ALLOCATE(rv(0:k,1:nu,1:nv),zv(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RV ZV',ier)
      ALLOCATE(g11(0:k,1:nu,1:nv),g12(0:k,1:nu,1:nv),g13(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'G11 G12 G13',ier)
      ALLOCATE(g22(0:k,1:nu,1:nv),g23(0:k,1:nu,1:nv),g33(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'G22 G23 G33',ier)
      ALLOCATE(detg(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'DETG',ier)
      CALL calc_metrics
      ! Now use virtual casing to calculate plasma response
      IF (lfreeb .and. extsurfs > 0) THEN
         ALLOCATE(bsmns_temp(1:mnmax,0:k),bumnc_temp(1:mnmax,0:k),bvmnc_temp(1:mnmax,0:k))
         !CALL virtual_casing(vsurf)
         CALL mgrid_field(vsurf+1,1)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bsmns_temp,xm,xn,bsreal,1,1)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bumnc_temp,xm,xn,bureal,0,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,bvmnc_temp,xm,xn,bvreal,0,0)
         bsmns(:,vsurf+1:k) = bsmns_temp(:,vsurf+1:k)
         bumnc(:,vsurf+1:k) = bumnc_temp(:,vsurf+1:k)
         bvmnc(:,vsurf+1:k) = bvmnc_temp(:,vsurf+1:k)
         DEALLOCATE(bsmns_temp,bumnc_temp,bvmnc_temp)
      END IF
      FORALL(ik = 0:k) rho(ik) = REAL(ik) / REAL(k)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE pies_init_wout
