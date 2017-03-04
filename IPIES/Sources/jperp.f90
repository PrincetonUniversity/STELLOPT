!-----------------------------------------------------------------------
!     Subroutine:    jperp
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/14/2011
!     Description:   This subroutine calculates the perpendicular
!                    current density from:
!                               B x grad(p)
!                     j_perp = ------------
!                                 |B|^2
!                    Note that we enforce j^s = 0.0 explicitly
!                    This is also the routine which calculates the
!                    pressure profile and derivative.
!-----------------------------------------------------------------------
      SUBROUTINE jperp
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error flag
!          ik           Radial dummy index
!          u            Poloidal dummy index
!          v            Toroidal dummy index
!          torflux_norm Normalized toroidal flux
!          dpdpsi       Pressure gradient dp/dpsi
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, ik, u, v, mn
      REAL(rprec) :: pi2, mu0, I_norm
      REAL(rprec) :: dpds(0:k), dsdr(0:k)
      REAL(rprec), ALLOCATABLE :: gpr(:,:,:),gpphi(:,:,:),gpz(:,:,:)
      REAL(rprec), ALLOCATABLE :: jr(:,:,:),jphi(:,:,:),jz(:,:,:)
      REAL(rprec), ALLOCATABLE :: br(:,:,:),bphi(:,:,:),bz(:,:,:)
      REAL(rprec), ALLOCATABLE :: bsq(:,:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      mu0 = 2*pi2*1e-7
      CALL toroidal_flux
      PRINT *,' '
      PRINT *,'torflux',torflux
      PRINT *,' '
      PRINT *,'torflux_norm',torflux_norm
      DO ik = 0, k
         IF (torflux_norm(ik) <= 1.0) THEN
            CALL EZspline_interp(p_spl,torflux_norm(ik),press(ik),ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_interp/jperp',ier)
            CALL EZspline_interp(ip_spl,torflux_norm(ik),iprime(ik),ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_interp/jperp',ier)
            CALL EZspline_derivative(p_spl,1,torflux_norm(ik),dpds(ik),ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_derivative/jperp',ier)
         ELSE
            press(ik)  = 0.0
            iprime(ik) = 0.0
            dpds(ik)   = 0.0
         END IF
      END DO
      WHERE(press == 0.0) dpds = 0.0
      ! Now we calculate jperp
      ALLOCATE(gpr(0:k,1:nu,1:nv),gpphi(0:k,1:nu,1:nv),gpz(0:k,1:nu,1:nv), STAT=ier) 
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'GPR GPPHI GPZ (jperp)',ier)
      ALLOCATE(jr(0:k,1:nu,1:nv),jphi(0:k,1:nu,1:nv),jz(0:k,1:nu,1:nv), STAT=ier) 
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'JR JPHI JZ (jperp)',ier)
      ALLOCATE(br(0:k,1:nu,1:nv),bphi(0:k,1:nu,1:nv),bz(0:k,1:nu,1:nv), STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BR BPHI BZ (jperp)',ier)
      ALLOCATE(bsq(0:k,1:nu,1:nv), STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSQ (jperp)',ier)
      PRINT *,' '
      PRINT *,'dpds',dpds
      dsdr(0)     = (torflux_norm(1) - torflux_norm(0))/(rho(1))
      dsdr(1:k-1) = (torflux_norm(2:k) - torflux_norm(0:k-2))/(rho(2:k)-rho(0:k-2))
      dsdr(k)     = (torflux_norm(k) - torflux_norm(k-1))/(rho(k)-rho(k-1))
      PRINT *,' '
      PRINT *,'dsdr',dsdr
      ! new dpds <- s is PIES radial coordiante now
      dpds(0)     = (press(1)   - press(0)    )/(rho(1))
      dpds(1:k-1) = (press(2:k) - press(0:k-2))/(rho(2:k)-rho(0:k-2))
      dpds(k)     = (press(k)   - press(k-1)  )/(rho(k)-rho(k-1))
      PRINT *,' '
      PRINT *,'dpds(new)',dpds
      
      br   = bsreal*rs+bureal*ru+bvreal*rv
      bphi = rreal*bvreal
      bz   = bsreal*zs+bureal*zu+bvreal*zv
      bsq  = br*br+bphi*bphi+bz*bz
      PRINT *,' '
      PRINT *,'br',br(0:k,1,1)
      PRINT *,' '
      PRINT *,'bphi',bphi(0:k,1,1)
      PRINT *,' '
      PRINT *,'bz',bz(0:k,1,1)
      DO u = 1, nu
         DO v = 1, nv
            gpr(:,u,v)   = dpds(:) * rs(:,u,v)
            gpphi(:,u,v) = 0.0
            gpz(:,u,v)   = dpds(:) * zs(:,u,v)
         END DO
      END DO
      jr = 0.0; jphi = 0.0; jz =0.0
      jr   = (bphi * gpz) / bsq
      jphi = (bz   * gpr   - br   * gpz) / bsq
      jz   =-(bphi * gpr) / bsq
      PRINT *,' '
      PRINT *,'jr',jr(0:k,1,1)
      PRINT *,' '
      PRINT *,'jphi',jphi(0:k,1,1)
      PRINT *,' '
      PRINT *,'jz',jz(0:k,1,1)
      stop
      ! Now add in Net current
      I_norm = curtor/SUM(iprime)
      IF (curtor == 0.0) I_norm = 0
      DO ik = 0, k
         jphi(ik,:,:) = jphi(ik,:,:) + iprime(ik)*jphi(ik,:,:)*I_norm
      END DO
      ! Now we need to transform to s,u,v space
      jsreal = 0.0; jureal=0.0; jvreal =0.0
      CALL cyl2suv(0,k,nu,nv,rreal,jr,jphi,jz,jsreal,jureal,jvreal,1._rprec)
      ! Now Fourier Transform the quantities
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,jsmns,xm,xn,jsreal,1,1)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,jumnc,xm,xn,jureal,0,0)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,jvmnc,xm,xn,jvreal,0,0)
      jsreal = 0.0; jsmns = 0.0
      IF (lasym) THEN
         jsmnc = 0.0
         !CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,jsmnc,xm,xn,jsreal,0,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,jumns,xm,xn,jureal,1,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,jvmns,xm,xn,jvreal,1,0)
      END IF
      PRINT *,'got here'
      ! Deallocate Arrays
      DEALLOCATE(jr, jphi, jz)
      DEALLOCATE(br, bphi, bz)
      DEALLOCATE(bsq)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jperp
