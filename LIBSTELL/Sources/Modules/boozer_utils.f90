!-----------------------------------------------------------------------
!     Module:        boozer_utils
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          07/02/2012
!     Description:   This module contains various routines for working
!                    with Boozer Coordinate data as generated by
!                    the BOOX_XFORM code.
!-----------------------------------------------------------------------
      MODULE boozer_utils
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_boozer_mod ! BOOZER Function in LIBSTELL
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      INTEGER     ::  domain_flag, nfp
      INTEGER     ::  bcs0(2) = (/ 0, 0/)
      INTEGER     ::  bcs1(2) = (/-1,-1/)
      REAL(rprec) ::  pi2, R_target, Z_target, PHI_target
      TYPE(EZspline3_r8) :: R_spl, Z_spl, P_spl, MODB_spl
      
!-----------------------------------------------------------------------
!     Subroutines
!         load_boozer:     Loads boozer variables
!         get_boozer_s:    Return boozer s,theta,phi coordiante
!-----------------------------------------------------------------------
      CONTAINS
        
      SUBROUTINE load_boozer(filename,iflag)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: filename
      INTEGER, INTENT(inout) :: iflag
      INTEGER ::  ier, iunit,nvar_in
      INTEGER ::  nu, nv, u, v, mn, dex, ns
      INTEGER, ALLOCATABLE :: im(:), in(:)
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:), rho(:)
      REAL(rprec), ALLOCATABLE :: r_temp(:,:,:), z_temp(:,:,:)
      REAL(rprec), ALLOCATABLE :: p_temp(:,:,:), modb_temp(:,:,:)
      REAL(rprec), ALLOCATABLE :: fmn_temp(:,:)
      pi2 = 8.0 * ATAN(1.0)
      IF (iflag < 0) RETURN
      ! Read the BOOZ_XFORM output
      CALL read_boozer_file(TRIM(filename),iflag)
      IF (iflag .ne. 0) RETURN
      ! Create the radial array
      dex = COUNT(idx_b == 1)
      ns = dex
      ALLOCATE(rho(ns))
      v=1
      DO u = 1, ns_b
         IF (idx_b(u) == 1) THEN
            rho(v) = phi_b(u)
            v = v + 1
         END IF
      END DO
      rho = rho / MAXVAL(rho)
      ! Get the realspace R and Z, p and |B|
      nu = 4* mboz_b + 1
      nv = 4* nboz_b + 1
      ALLOCATE(im(mnboz_b),in(mnboz_b))
      ALLOCATE(xu(nu),xv(nv))
      ALLOCATE(r_temp(nu,nv,ns),z_temp(nu,nv,ns))
      ALLOCATE(p_temp(nu,nv,ns),modb_temp(nu,nv,ns))
      ALLOCATE(fmn_temp(mnboz_b,ns))
      FORALL(u=1:nu) xu(u) = REAL(u-1)/REAL(nu-1)
      FORALL(v=1:nv) xv(v) = REAL(v-1)/REAL(nv-1)
      im =  ixm_b
      in = -ixn_b/nfp_b
      nfp = nfp_b
      r_temp =0; z_temp =0; p_temp=0; modb_temp = 0.0;
      CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,rmnc_b,im,in,r_temp,0,1)
      CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,zmns_b,im,in,z_temp,1,0)
      CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,pmns_b,im,in,p_temp,1,0)
      CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,bmnc_b,im,in,modb_temp,0,0)
      IF (lasym_b) THEN
         CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,rmns_b,im,in,r_temp,1,0)
         CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,zmnc_b,im,in,z_temp,0,0)
         CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,pmnc_b,im,in,p_temp,0,0)
         CALL mntouv(1,ns,mnboz_b,nu,nv,xu,xv,bmns_b,im,in,modb_temp,1,0)
      END IF
      ! Constructed 3D Splines
      IF (EZspline_allocated(R_spl)) CALL EZspline_free(R_spl,iflag)
      IF (EZspline_allocated(Z_spl)) CALL EZspline_free(Z_spl,iflag)
      IF (EZspline_allocated(P_spl)) CALL EZspline_free(P_spl,iflag)
      IF (EZspline_allocated(MODB_spl)) CALL EZspline_free(MODB_spl,iflag)
      CALL EZspline_init(R_spl,nu,nv,ns,bcs1,bcs1,bcs0,iflag)
      CALL EZspline_init(Z_spl,nu,nv,ns,bcs1,bcs1,bcs0,iflag)
      CALL EZspline_init(P_spl,nu,nv,ns,bcs1,bcs1,bcs0,iflag)
      CALL EZspline_init(MODB_spl,nu,nv,ns,bcs1,bcs1,bcs0,iflag)
      R_spl%x1 = xu*pi2
      R_spl%x2 = xv*pi2
      R_spl%x3 = rho
      R_spl%isHermite = 1
      Z_spl%x1 = xu*pi2
      Z_spl%x2 = xv*pi2
      Z_spl%x3 = rho
      Z_spl%isHermite = 1
      P_spl%x1 = xu*pi2
      P_spl%x2 = xv*pi2
      P_spl%x3 = rho
      P_spl%isHermite = 1
      MODB_spl%x1 = xu*pi2
      MODB_spl%x2 = xv*pi2
      MODB_spl%x3 = rho
      MODB_spl%isHermite = 1
      CALL EZspline_setup(R_spl,r_temp,iflag)
      CALL EZspline_setup(Z_spl,z_temp,iflag)
      CALL EZspline_setup(P_spl,p_temp,iflag)
      CALL EZspline_setup(MODB_spl,modb_temp,iflag)
      RETURN
      END SUBROUTINE load_boozer
      
      SUBROUTINE rzfunct_booz(m,n,x,fvec,fjac,ldfjac,iflag)
      IMPLICIT NONE
      INTEGER m,n,ldfjac,iflag, ier
      DOUBLE PRECISION x(n),fvec(m),fjac(ldfjac,n)
      REAL(rprec) :: R_temp, Z_temp, P_temp
      REAL(rprec) :: R_grad(3), Z_grad(3), P_grad(3)
      IF (x(2) < 0.0) x(2) = x(2) + pi2
      x(2) = MOD(x(2),pi2)
      IF (x(3) < 0.0) x(3) = x(3) + pi2
      x(3) = MOD(x(3),pi2)
      IF (x(1) < 0) THEN
         x(1) = -x(1)/2
         x(2) = x(2) +0.5*pi2
         x(2) = MOD(x(2),pi2)
      END IF
      ier = 0
      CALL EZspline_isInDomain(R_spl,x(2),x(3),x(1),ier)
      IF (ier .ne. 0) THEN
         iflag = -1
         domain_flag = -1
      END IF
      IF (iflag == 1) THEN
         CALL EZspline_interp(R_spl,x(2),x(3),x(1),R_temp,iflag)
         CALL EZspline_interp(Z_spl,x(2),x(3),x(1),Z_temp,iflag)
         CALL EZspline_interp(P_spl,x(2),x(3),x(1),P_temp,iflag)
         fvec(1) = (R_temp - R_target)
         fvec(2) = (Z_temp - Z_target)
         fvec(3) = (x(3) + P_temp - PHI_target) ! PHI(cyl) = PHI(booz)+P(s,u,v) 
         !PRINT *,R_temp,Z_temp,P_temp
      ELSE IF (iflag == 2) THEN
         CALL EZspline_gradient(R_spl,x(2),x(3),x(1),R_grad,iflag)
         CALL EZspline_gradient(Z_spl,x(2),x(3),x(1),Z_grad,iflag)
         CALL EZspline_gradient(P_spl,x(2),x(3),x(1),P_grad,iflag)
         fjac(1,1) = R_grad(3) !dR/ds
         fjac(1,2) = R_grad(1) !dR/du
         fjac(1,3) = R_grad(2) !dR/dv
         fjac(2,1) = Z_grad(3) !dZ/ds
         fjac(2,2) = Z_grad(1) !dZ/du
         fjac(2,3) = Z_grad(2) !dZ/dv
         fjac(3,1) = P_grad(3) !dP/ds
         fjac(3,2) = P_grad(1) !dP/du
         fjac(3,3) = 1+P_grad(2) !dP/dv
      END IF
      RETURN
      END SUBROUTINE rzfunct_booz
      
      SUBROUTINE get_booz_s(r_val,phi_val,z_val,s_val,ier,u_val,v_val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  r_val
      REAL(rprec), INTENT(in)    ::  phi_val
      REAL(rprec), INTENT(in)    ::  z_val
      REAL(rprec), INTENT(out)   ::  s_val
      REAL(rprec), INTENT(out), OPTIONAL   ::  u_val
      REAL(rprec), INTENT(out), OPTIONAL   ::  v_val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER, PARAMETER :: mfunct=3
      INTEGER, PARAMETER :: nvars=3
      INTEGER, PARAMETER :: ldfjac=3
      INTEGER :: ik, maxfev_local, nfev, info, njev, maxfev, nprint
      INTEGER, DIMENSION(nvars) :: ipvt
      REAL(rprec) :: ftol,xtol,gtol,factor, mode
      REAL(rprec), DIMENSION(nvars)  :: xc_opt, diag, qtf, wa1, wa2, wa3
      REAL(rprec), DIMENSION(mfunct) :: fval,wa4
      REAL(rprec), DIMENSION(ldfjac,nvars) :: fjac
      IF (ier < 0) RETURN
      IF (EZspline_allocated(R_spl) .and. EZspline_allocated(Z_spl) .and. &
          EZspline_allocated(P_spl)) THEN
         R_target = r_val
         Z_target = z_val
         PHI_target = phi_val
         IF (PHI_target < 0) PHI_target = PHI_target + pi2
         PHI_target = MOD(PHI_target,pi2/nfp)*nfp
         xc_opt(1) = 0.5
         xc_opt(2) = 0.0
         xc_opt(3) = PHI_target
         ftol = 1.0E-12
         xtol = 1.0E-12
         gtol = 1.0E-12
         maxfev_local = 5000
         diag(:) = 1.0
         factor = 0.1
         nprint = 0
         domain_flag = 0
         DO ik = 1, 4
            CALL lmder_serial(rzfunct_booz,mfunct,nvars,xc_opt,fval,fjac,ldfjac,ftol,xtol,gtol,&
                    maxfev_local,diag,mode,factor,nprint,info,nfev,njev,ipvt,qtf,&
                    wa1,wa2,wa3,wa4)
            IF (info < 4) EXIT
            xc_opt(2) = xc_opt(2) + 0.05*pi2
            info = 0
            nfev = 0
            njev = 0
         END DO
         ier = info
         IF (info < 4) ier = 0
         IF (domain_flag .ne. 0) THEN
            ier = 9
            s_val = 1.5
            IF (PRESENT(u_val)) u_val = 2*pi2
            IF (PRESENT(v_val)) v_val = 2*pi2
            RETURN
         END IF
         s_val = xc_opt(1)
         IF (PRESENT(u_val)) THEN
            u_val = xc_opt(2)
            u_val = MOD(u_val,pi2)
            DO WHILE (u_val < 0)
               u_val = u_val + pi2
            END DO
         END IF
         IF (PRESENT(v_val)) THEN
            v_val = xc_opt(3)
            v_val = MOD(v_val,pi2)
            DO WHILE (v_val < 0)
               v_val = v_val + pi2
            END DO
         END IF
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_booz_s
         
         
      SUBROUTINE mntouv(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: k1
      INTEGER, INTENT(in) :: k
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nu
      INTEGER, INTENT(in) :: nv
      REAL(rprec), INTENT(in) :: xu(1:nu)
      REAL(rprec), INTENT(in) :: xv(1:nv)           
      REAL(rprec), INTENT(in) :: fmn(1:mnmax,k1:k)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      REAL(rprec), INTENT(inout) :: f(1:nu,1:nv,k1:k)
      INTEGER, INTENT(in) :: signs
      INTEGER, INTENT(in) :: calc_trig
      INTEGER     :: mn, i, ier, ik
      REAL(rprec) :: xm_temp(1:mnmax,1)
      REAL(rprec) :: xn_temp(1:mnmax,1)
      REAL(rprec) :: pi2_l
      REAL(rprec) :: mt(1:mnmax,1:nu)
      REAL(rprec) :: nz(1:mnmax,1:nv)
      REAL(rprec) :: fmn_temp(1:mnmax,1:nu)
      REAL(rprec) :: xu_temp(1,1:nu)
      REAL(rprec) :: xv_temp(1,1:nv)
      REAL(rprec) :: fmn_help(1:mnmax)
      REAL(rprec), ALLOCATABLE, SAVE :: cosmt(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: sinmt(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: cosnz(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: sinnz(:,:)
      pi2_l = 8 * ATAN(1.)
      IF (calc_trig == 1) THEN
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         ALLOCATE(cosmt(1:mnmax,1:nu),sinmt(1:mnmax,1:nu),&
                  cosnz(1:mnmax,1:nv),sinnz(1:mnmax,1:nv),STAT=ier)
         FORALL(i=1:mnmax) xm_temp(i,1)=REAL(xm(i))
         FORALL(i=1:mnmax) xn_temp(i,1)=REAL(xn(i))
         FORALL(i=1:nu) xu_temp(1,i)=xu(i)
         FORALL(i=1:nv) xv_temp(1,i)=xv(i)
         mt = MATMUL(xm_temp,xu_temp)
         nz = MATMUL(xn_temp,xv_temp)
         FORALL(mn=1:mnmax,i=1:nu) cosmt(mn,i) = dcos(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nu) sinmt(mn,i) = dsin(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) cosnz(mn,i) = dcos(pi2_l*nz(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) sinnz(mn,i) = dsin(pi2_l*nz(mn,i))
      END IF
      IF (SIGNS == 0) THEN
         DO ik = k1,k
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik)  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),cosnz) &
                                   - MATMUL(TRANSPOSE((fmn_temp*sinmt)),sinnz)
         END DO
      ELSE IF (SIGNS == 1) THEN
         DO ik = k1,k
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik) + MATMUL(TRANSPOSE((fmn_temp*sinmt)),cosnz) &
                                  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),sinnz)
         END DO
      END IF
      END SUBROUTINE mntouv
      
      SUBROUTINE free_booz(ier)
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: ier
      IF (ier < 0) RETURN
      CALL read_boozer_deallocate
      CALL EZspline_free(R_spl,ier)
      CALL EZspline_free(Z_spl,ier)
      CALL EZspline_free(P_spl,ier)
      CALL EZspline_free(MODB_spl,ier)
      RETURN
      END SUBROUTINE free_booz
      
      END MODULE boozer_utils