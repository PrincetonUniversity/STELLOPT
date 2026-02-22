!-----------------------------------------------------------------------
!     Module:        spline_coils_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/10/2025
!     Description:   This module contains routines for representing
!                    coils as splines in 3D space.
!-----------------------------------------------------------------------
      MODULE spline_coils_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE bsc_T
      USE biotsavart, ONLY: coil_group, nfp => nfp_bs
      USE safe_open_mod
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER, PRIVATE :: pi2 = 6.283185482025146D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: zero = 0.0D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: one = 1.0D+00
      INTEGER, PRIVATE :: mnmax, ns, ncoilgroups, n_kts, k_kts, nw_coil, nh_coil
      DOUBLE PRECISION, PRIVATE :: factor
      DOUBLE PRECISION, PRIVATE :: curvature_mean, curvature_max, curvature_min
      DOUBLE PRECISION, PRIVATE :: torsion_mean, torsion_max, torsion_min
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE :: xm, xn, &
         rmnc, zmns, rmnc0, zmns0, t_kts
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: rho_kts, &
         theta_kts, zeta_kts
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: &
         curvature, torsion, dLength
      TYPE(EZspline1_r8), DIMENSION(20) :: RHO_spl, THETA_spl, ZETA_spl
      
!-----------------------------------------------------------------------
!     Module SUBROUTINES/FUNCTIONS
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE init_spline_coils(ns_in, ncoilgroups_in, n_in, nk_in, &
            rho_in, theta_in, zeta_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ns_in
      INTEGER, INTENT(in) :: ncoilgroups_in
      INTEGER, INTENT(in) :: n_in
      INTEGER, INTENT(in) :: nk_in
      DOUBLE PRECISION, INTENT(in) :: rho_in(ncoilgroups_in,n_in)
      DOUBLE PRECISION, INTENT(in) :: theta_in(ncoilgroups_in,n_in)
      DOUBLE PRECISION, INTENT(in) :: zeta_in(ncoilgroups_in,n_in)
      INTEGER :: i
      INTEGER :: bcs0(2), ier
      bcs0=(/-1,-1/)
      ns = ns_in
      ncoilgroups = ncoilgroups_in
      n_kts = n_in+1
      k_kts = nk_in-n_in
      IF (ALLOCATED(rho_kts)) DEALLOCATE(rho_kts)
      IF (ALLOCATED(theta_kts)) DEALLOCATE(theta_kts)
      IF (ALLOCATED(zeta_kts)) DEALLOCATE(zeta_kts)
      IF (ALLOCATED(t_kts)) DEALLOCATE(t_kts)
      ALLOCATE(rho_kts(ncoilgroups,n_kts), theta_kts(ncoilgroups,n_kts), &
            zeta_kts(ncoilgroups,n_kts), t_kts(n_kts))
      rho_kts(:,1:n_in) = rho_in
      theta_kts(:,1:n_in) = theta_in
      zeta_kts(:,1:n_in) = zeta_in
      rho_kts(:,n_kts) = rho_kts(:,1)
      theta_kts(:,n_kts) = theta_kts(:,1)+pi2
      zeta_kts(:,n_kts) = zeta_kts(:,1)
      FORALL(i=1:n_kts) t_kts(i) = DBLE(i-1)/DBLE(n_kts-1)
      ! Spline stuff
      DO i = 1, ncoilgroups
         ier = 0
         IF (EZspline_allocated(RHO_spl(i))) CALL EZspline_free(RHO_spl(i),ier)
         IF (EZspline_allocated(THETA_spl(i))) CALL EZspline_free(THETA_spl(i),ier)
         IF (EZspline_allocated(ZETA_spl(i))) CALL EZspline_free(ZETA_spl(i),ier)
         CALL EZspline_init(RHO_spl(i),n_kts,bcs0,ier)
         CALL EZspline_init(THETA_spl(i),n_kts,bcs0,ier)
         CALL EZspline_init(ZETA_spl(i),n_kts,bcs0,ier)
         RHO_spl(i)%x1          = t_kts
         THETA_spl(i)%x1        = t_kts
         ZETA_spl(i)%x1         = t_kts
         RHO_spl(i)%isHermite   = 1
         THETA_spl(i)%isHermite = 1
         ZETA_spl(i)%isHermite  = 1
         CALL EZspline_setup(RHO_spl(i),rho_kts(i,:),ier,EXACT_DIM=.true.)
         CALL EZspline_setup(THETA_spl(i),theta_kts(i,:),ier,EXACT_DIM=.true.)
         CALL EZspline_setup(ZETA_spl(i),zeta_kts(i,:),ier,EXACT_DIM=.true.)
      END DO
      RETURN
      END SUBROUTINE init_spline_coils

      SUBROUTINE init_boundary_spline_coils(mnmax_in,xm_in,xn_in,rmnc_in,zmns_in,rmnc_ax,zmns_ax)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mnmax_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: xm_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: xn_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: rmnc_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: zmns_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: rmnc_ax
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: zmns_ax
      mnmax = mnmax_in
      IF (ALLOCATED(xm)) DEALLOCATE(xm)
      IF (ALLOCATED(xn)) DEALLOCATE(xn)
      IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc)
      IF (ALLOCATED(zmns)) DEALLOCATE(zmns)
      IF (ALLOCATED(rmnc0)) DEALLOCATE(rmnc0)
      IF (ALLOCATED(zmns0)) DEALLOCATE(zmns0)
      ALLOCATE(xm(mnmax), xn(mnmax), rmnc(mnmax), zmns(mnmax), rmnc0(mnmax), zmns0(mnmax))
      xm = xm_in
      xn = xn_in
      rmnc = rmnc_in
      zmns = zmns_in
      rmnc0 = rmnc_ax
      zmns0 = zmns_ax
      nfp = MINVAL(xn, MASK = xn > 0)
      xn = xn / nfp
      factor = pi2/nfp
      RETURN
      END SUBROUTINE init_boundary_spline_coils

      SUBROUTINE rhothetazeta2xyz(rho_in,theta_in,zeta_in,x_out,y_out,z_out)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: rho_in, theta_in, zeta_in
      DOUBLE PRECISION, INTENT(out) :: x_out, y_out, z_out
      INTEGER :: mn
      DOUBLE PRECISION :: R, Z, RN, ZN, REDGE, ZEDGE, PHI, N
      DOUBLE PRECISION :: rho_ext, whi, wlo, wloo, whio, cop, sip
      R = zero; Z = zero
      RN = zero; ZN = zero
      REDGE = zero; ZEDGE = zero
      ! Extrapolation stuff (like VMEC)
      PHI = zeta_in/nfp
      rho_ext = rho_in + 1.0
      whi   = (rho_ext*rho_ext-1.0)*2.0
      wlo   = (2.0 - whi)/2.0
      wloo  = wlo*rho_ext
      whio  = whi*rho_ext/SQRT(2.0)
      DO mn = 1, mnmax
         cop = cos(xm(mn)*theta_in+xn(mn)*zeta_in)
         sip = sin(xm(mn)*theta_in+xn(mn)*zeta_in)
         REDGE = REDGE + rmnc(mn)*cop
         ZEDGE = ZEDGE + zmns(mn)*sip
         IF ((xm(mn) == 0) .and. (xn(mn) == 0)) THEN
            R =  R  + rmnc(mn)*cop
         ELSEIF (MOD(int(xm(mn)),2)==0) THEN
            R = R + rmnc(mn)*wlo*cop
            Z = Z + zmns(mn)*wlo*sip
         ELSE
            R = R + rmnc(mn)*wloo*cop
            Z = Z + zmns(mn)*wloo*sip
         END IF
         IF ((xm(mn)==1) .and. (xn(mn)==0)) THEN
            ! Note we use odd here since xm==1
            R    =  R + 4.0*whio*cop
            Z    =  Z + 4.0*whio*sip
         END IF
      END DO
      RN    = R - REDGE
      ZN    = Z - ZEDGE
      N     = SQRT(RN*RN+ZN*ZN)
      RN    = RN/N; ZN = ZN/N
      R     = (REDGE + rho_in*RN)
      x_out = R*COS(PHI)
      y_out = R*SIN(PHI)
      z_out = ZEDGE + rho_in*ZN
      RETURN
      END SUBROUTINE rhothetazeta2xyz

      SUBROUTINE spline_to_coils(normal_sign)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: normal_sign
      INTEGER :: i, j, mn, ns1, ier
      DOUBLE PRECISION :: rho, theta, zeta, X, Y, Z, cop, sip, L
      DOUBLE PRECISION, DIMENSION(ns) :: Rc,Zc,Pc
      DOUBLE PRECISION, DIMENSION(3,ns) :: xnod_in, xnod_ss, xnod_bb
      CHARACTER(len=100) :: s_name
      CHARACTER(len=100) :: l_name
      CHARACTER(len=100) :: c_name
      TYPE(bsc_coil)     :: coil_temp
      TYPE(bsc_rs)       :: rot_mat
      nw_coil = 1; nh_coil = 1
      ns1 = ns - 1
      ! Deallocated the coils if allocated
      IF (ALLOCATED(coil_group)) THEN
         DO i = 1, SIZE(coil_group)
            CALL bsc_destroy(coil_group(i))
         END DO
         DEALLOCATE(coil_group)
      END IF
      ! Allocate the coilgroups
      ALLOCATE(coil_group(ncoilgroups))
      ! Create each coil object
      DO i = 1,ncoilgroups
         WRITE(l_name,*) 'i = ',i
         WRITE(c_name,'(A,I2.2)') 'MODULAR_COIL_',i
         CALL bsc_construct_coilcoll(coil_group(i),TRIM(c_name),TRIM(l_name))
         DO j = 1, ns
            l = DBLE(j-1)/DBLE(ns-1)
            ier = 0
            CALL EZspline_interp(RHO_spl(i),l,rho,ier)
            CALL EZspline_interp(THETA_spl(i),l,theta,ier)
            CALL EZspline_interp(ZETA_spl(i),l,zeta,ier)
            CALL rhothetazeta2xyz(rho,theta,zeta,X,Y,Z)
            RC(j) = SQRT(X*X+Y*Y)
            PC(j) = ATAN2(Y,X)
            ZC(j) = Z
            xnod_in(1,j) = X
            xnod_in(2,j) = Y
            xnod_in(3,j) = Z
         END DO
         xnod_in(:,ns) = xnod_in(:,1)
         ! Now create the first coil
         WRITE(s_name, '(a4,i5.5)') 'ID ', 1
         CALL bsc_construct_coil(coil_temp,'fil_loop',TRIM(s_name),'',one,xnod_in(1:3,1:ns))
         CALL bsc_append(coil_group(i),coil_temp)
         ! Now create the stellarator symmetric coil
         Zc = -Zc
         Pc = factor - Pc
         xnod_ss(1,2:ns) = Rc(ns1:1:-1)*cos(Pc(ns1:1:-1))
         xnod_ss(2,2:ns) = Rc(ns1:1:-1)*sin(Pc(ns1:1:-1))
         xnod_ss(3,2:ns) = Zc(ns1:1:-1)
         xnod_ss(:,1) = xnod_ss(:,ns)
         CALL bsc_construct_coil(coil_temp,'fil_loop',TRIM(s_name),'',one,xnod_ss(1:3,1:ns))
         CALL bsc_append(coil_group(i),coil_temp)
         ! Now create rest of the coils
         DO j = 2, nfp
            cop  = cos((j-1)*factor)
            sip  = sin((j-1)*factor)
            xnod_bb(1,:) = xnod_in(1,:)*cop - xnod_in(2,:)*sip
            xnod_bb(2,:) = xnod_in(2,:)*cop + xnod_in(1,:)*sip
            xnod_bb(3,:) = xnod_in(3,:)
            CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_bb(1:3,1:ns))
            CALL bsc_append(coil_group(i),coil_temp)
            xnod_bb(1,:) = xnod_ss(1,:)*cop - xnod_ss(2,:)*sip
            xnod_bb(2,:) = xnod_ss(2,:)*cop + xnod_ss(1,:)*sip
            xnod_bb(3,:) = xnod_ss(3,:)
            CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_bb(1:3,1:ns))
            CALL bsc_append(coil_group(i),coil_temp)
         END DO
      END DO
      RETURN
      END SUBROUTINE spline_to_coils

      SUBROUTINE xyz2rhothetazeta(x_in,y_in,z_in,rho_out,theta_out,zeta_out)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x_in, y_in, z_in
      DOUBLE PRECISION, INTENT(inout) :: rho_out, theta_out, zeta_out
      INTEGER :: nfe
      DOUBLE PRECISION :: R_in, fnorm, fmin0, fmin, fact_local, &
         rho_min, theta_min, X, Y, Z, R, X1, Y1, Z1, R1, delta, &
         dRdrho, dZdrho, dRdtheta, dZdtheta, dR, dZ, tau, &
         delrho, deltheta
      zeta_out = ATAN2(y_in,x_in)*nfp
      rho_out = MAX(rho_out,0.0)
      R_in = SQRT(x_in*x_in + y_in*y_in)
      fnorm = one / SQRT(R_in*R_in+Z_in*Z_in)
      fmin0 = 1.0D+10
      fmin  = 1.0D+10
      delta = 1.0D-03
      fact_local = one
      rho_min = rho_out; theta_min = theta_out
      nfe = 0
      DO WHILE ((nfe .lt. 1000) .and. (fmin .gt. 1.0E-6))
         nfe = nfe + 1
         ! Compute R,Z
         CALL rhothetazeta2xyz(rho_out,theta_out,zeta_out,X,Y,Z)
         R = SQRT(X*X + Y*Y)
         ! Compute dR/drho and dZ/drho
         CALL rhothetazeta2xyz(rho_out+delta,theta_out,zeta_out,X1,Y1,Z1)
         R1 = SQRT(X1*X1 + Y1*Y1)
         dRdrho = (R1-R)/delta
         dZdrho = (Z1-Z)/delta
         ! Compute dR/dtheta and dZ/dtheta
         CALL rhothetazeta2xyz(rho_out,theta_out+delta,zeta_out,X1,Y1,Z1)
         R1 = SQRT(X1*X1 + Y1*Y1)
         dRdtheta = (R1-R)/delta
         dZdtheta = (Z1-Z)/delta
         ! Compute Function minimization (R0,Z0)
         dR = R - R_in
         dZ = Z - Z_in
         fmin = (dR*dR+dZ*dZ)*fnorm
         !PRINT *,nfe,rho_out,theta_out,R,R1,Z,Z1
         !PRINT *,'===',dRdrho,dZdrho,dRdtheta,dZdtheta
         !PRINT *,'===',dR,dZ,fmin
         ! Compute Descent Direction
         IF (fmin .gt. fmin0) THEN
            fact_local = (2*fact_local)/3
            rho_out = rho_min; theta_out = theta_min
            ! REDIRECT ALONG STEEPEST-DESCENT PATH
            IF (6*fact_local .lt. one) THEN
               !xu(1) = ru1; xu(3) = zu1
               !xs(1) = rs1; xs(3) = zs1
               !dels =-(s*rs1 + u*zs1)/(rs1**2 + zs1**2)
               !delu =-(x0(1)*xu(1) + x0(3)*xu(3))/(xu(1)**2 + xu(3)**2)
               delrho =-(dR*dRdrho + dZ*dZdrho)/(dRdrho*dRdrho + dZdrho*dZdrho)
               deltheta =-(dR*dRdtheta + dZ*dZdtheta)/(dRdtheta*dRdtheta + dZdtheta*dZdtheta)
            END IF
         ELSE
            fmin0 = fmin
            fact_local = one
            rho_min = rho_out
            theta_min = theta_out
            !NEWTON STEP
            tau = dRdtheta*dZdrho - dZdtheta * dRdrho
            !dels = ( x0(1)*xu(3) - x0(3)*xu(1))/tau
            !delu = (-x0(1)*xs(3) + x0(3)*xs(1))/tau
            delrho = ( dR*dZdtheta - dZ*dRdtheta)/tau
            deltheta = (-dR*dZdrho + dZ*dRdrho)/tau
            IF (fmin .gt. 1.0D-03) THEN
               delrho = delrho*0.5; deltheta = deltheta*0.5
            END IF
         END IF
         !PRINT *,'===',delrho,deltheta
         rho_out = MIN(MAX(rho_out + delrho*fact_local,1.0D-3),10.0)
         theta_out = MOD(theta_out + deltheta*fact_local,pi2)
      END DO
      RETURN
      END SUBROUTINE xyz2rhothetazeta

      SUBROUTINE set_currents(nextcur,extcur)
      INTEGER, INTENT(IN) :: nextcur
      DOUBLE PRECISION, DIMENSION(nextcur), INTENT(IN) :: extcur
      INTEGER i,j
      DO i = 1, nextcur
         DO j = 1, coil_group(i)%ncoil
            coil_group(i)%coils(j)%current = extcur(i)
         END DO
      END DO
      RETURN
      END SUBROUTINE set_currents

      SUBROUTINE coils_to_multifilament(nw,nh,width,height)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nw
      INTEGER, INTENT(IN) :: nh
      DOUBLE PRECISION, INTENT(IN) :: width
      DOUBLE PRECISION, INTENT(IN) :: height
      INTEGER :: i,j,ns1,l,k
      DOUBLE PRECISION :: xc,yc,zc,ntotal
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xn,yn,zn,nt,xb,yb,zb
      DOUBLE PRECISION, DIMENSION(3,ns) :: xnod
      TYPE(bsc_coil)     :: coil_temp
      TYPE (bsc_coilcoll), DIMENSION(:), ALLOCATABLE, TARGET :: coil_single
      nw_coil = nw; nh_coil = nh
      ! Save the original coil
      ALLOCATE(coil_single(ncoilgroups))
      DO i = 1, ncoilgroups
         CALL bsc_construct_coilcoll(coil_single(i),TRIM(coil_group(i)%s_name),TRIM(coil_group(i)%l_name))
         coil_single(i)%ncoil = coil_group(i)%ncoil
         coil_single(i)%coils = coil_group(i)%coils
      END DO
      ! Deallocated the coils if allocated
      IF (ALLOCATED(coil_group)) THEN
         DO i = 1, SIZE(coil_group)
            CALL bsc_destroy(coil_group(i))
         END DO
         DEALLOCATE(coil_group)
      END IF
      ! Allocate the coilgroups
      ALLOCATE(coil_group(ncoilgroups))
      ! Allocate the helpers
      ns1 = ns - 1
      ALLOCATE(xn(ns1),yn(ns1),zn(ns1),nt(ns1))
      ALLOCATE(xb(ns1),yb(ns1),zb(ns1))
      ! Now loop over each coil
      DO i = 1, ncoilgroups
         CALL bsc_construct_coilcoll(coil_group(i),TRIM(coil_single(i)%s_name),TRIM(coil_single(i)%l_name))
         DO j = 1, coil_single(i)%ncoil
            ! Compute geometric center
            ntotal = SIZE(coil_single(i)%coils(j)%xnod(1,:))-1
            xc = SUM(coil_single(i)%coils(j)%xnod(1,2:))/(ns-1)
            yc = SUM(coil_single(i)%coils(j)%xnod(2,2:))/(ns-1)
            zc = SUM(coil_single(i)%coils(j)%xnod(3,2:))/(ns-1)
            ! Compute Normal vector
            xn = coil_single(i)%coils(j)%xnod(1,:ns1)-xc
            yn = coil_single(i)%coils(j)%xnod(2,:ns1)-yc
            zn = coil_single(i)%coils(j)%xnod(3,:ns1)-zc
            nt = xn * coil_single(i)%coils(j)%ehnod(1,:) &
               + yn * coil_single(i)%coils(j)%ehnod(2,:) &
               + zn * coil_single(i)%coils(j)%ehnod(3,:)
            xn = xn - nt * coil_single(i)%coils(j)%ehnod(1,:)
            yn = yn - nt * coil_single(i)%coils(j)%ehnod(2,:)
            zn = zn - nt * coil_single(i)%coils(j)%ehnod(3,:)
            nt = SQRT(xn*xn + yn*yn + zn*zn)
            xn = xn / nt
            yn = yn / nt
            zn = zn / nt
            ! Compute bi-normal vector
            xb = coil_single(i)%coils(j)%ehnod(2,:) * zn &
               - coil_single(i)%coils(j)%ehnod(3,:) * yn
            yb = coil_single(i)%coils(j)%ehnod(3,:) * xn &
               - coil_single(i)%coils(j)%ehnod(1,:) * zn
            zb = coil_single(i)%coils(j)%ehnod(1,:) * yn &
               - coil_single(i)%coils(j)%ehnod(2,:) * xn
            nt = SQRT(xb*xb + yb*yb + zb*zb)
            xb = xb / nt
            yb = yb / nt
            zb = zb / nt
            DO l = 1, nw
               DO k = 1, nh
                  xnod = coil_single(i)%coils(j)%xnod
                  xnod(1,1:ns1) = xnod(1,1:ns1) - xb*width/2 - xn*height/2 &
                            + xb*width*(l-1)/(nw-1) + xn*height*(k-1)/(nh-1)
                  xnod(2,1:ns1) = xnod(2,1:ns1) - yb*width/2 - yn*height/2 &
                            + yb*width*(l-1)/(nw-1) + yn*height*(k-1)/(nh-1)
                  xnod(3,1:ns1) = xnod(3,1:ns1) - zb*width/2 - zn*height/2 &
                            + zb*width*(l-1)/(nw-1) + zn*height*(k-1)/(nh-1)
                  xnod(:,ns) = xnod(:,1)
                  CALL bsc_construct_coil(coil_temp,'fil_loop',TRIM(coil_single(i)%coils(j)%s_name),'',coil_single(i)%coils(j)%current/(nh*nw),xnod(1:3,1:ns))
                  CALL bsc_append(coil_group(i),coil_temp)
               END DO
            END DO
         END DO
      END DO
      ! Deallocated the coils if allocated
      IF (ALLOCATED(coil_single)) THEN
         DO i = 1, SIZE(coil_single)
            CALL bsc_destroy(coil_single(i))
         END DO
         DEALLOCATE(coil_single)
      END IF
      ! Deallocate helpers
      DEALLOCATE(xn,yn,zn,nt,xb,yb,zb)
      RETURN
      END SUBROUTINE coils_to_multifilament

      SUBROUTINE dump_coil_info(iunit)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: iunit
      INTEGER :: k,j
      WRITE(iunit,'(A,2X,I6)') ' Number of coil groups',SIZE(coil_group)
      WRITE(iunit,'(A)')       ' ======================================='
      DO k = 1, SIZE(coil_group)
         WRITE(iunit,'(I6)') k
         WRITE(iunit,'(A,2X,A)') ' S_NAME ',TRIM(coil_group(k)%s_name)
         WRITE(iunit,'(A,2X,A)') ' L_NAME ',TRIM(coil_group(k)%l_name)
         WRITE(iunit,'(A,2X,I6)') ' Number of coils ',coil_group(k)%ncoil
         WRITE(iunit,'(A)')       ' ---------------------------------------'
         DO j = 1, coil_group(k)%ncoil
            WRITE(iunit,'(I6)') j
            WRITE(iunit,'(A,2X,A)') ' S_NAME ',TRIM(coil_group(k)%coils(j)%c_type)
            WRITE(iunit,'(A,2X,A)') ' S_NAME ',TRIM(coil_group(k)%coils(j)%s_name)
            WRITE(iunit,'(A,2X,A)') ' L_NAME ',TRIM(coil_group(k)%coils(j)%l_name)
            WRITE(iunit,'(A,2X,ES20.10)') ' Current ',coil_group(k)%coils(j)%current
            WRITE(iunit,'(A,2X,I6)') ' Number of points ',SIZE(coil_group(k)%coils(j)%xnod)
         ENDDO
      WRITE(iunit,'(A)')       ' ======================================='

      END DO
      END SUBROUTINE dump_coil_info

      ! SUBROUTINE compute_coil_coil_distance(outext)
      ! IMPLICIT NONE
      ! CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: outext
      ! INTEGER :: nc1, nc2
      ! DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dist
      ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: x2d, y2d, z2d, d2d
      ! !--------------------------------------------------------------
      ! !    Output the coil_coil_distance
      ! !--------------------------------------------------------------
      ! loutput = .FALSE.
      ! IF (PRESENT(outext)) THEN
      !    loutput = .TRUE.
      !    CALL safe_open(iunit_out,ier,TRIM('coil_dist.'//TRIM(outext)),'unknown','formatted')
      !    WRITE(iunit_out,'(I6,2X,I6,2X,I6,2X,I6)') ncoilgroups,nw_coil,nh_coil,ns
      ! END IF
      ! !--------------------------------------------------------------
      ! !    First compare every coil to every coil inside a given group
      ! ! This is a mess
      ! !--------------------------------------------------------------
      ! ALLOCATE(dist(nw_coil*nh_coil*ncoilgroups,1))
      ! DO i1 = 1, ncoilgroups
      !    DO j1 = 1, nw_coil*nh_coil
      !       n1 = nw_coil*nh_coil + 1
      !       nc1 = SIZE(coil_group(i1)%coils(j1)%xnod,2)
      !       ALLOCATE(dist(nc1))
      !       DO j2 = n1, coil_group(i1)%ncoil
      !          nc2 = SIZE(coil_group(i1)%coils(j2)%xnod,2)
      !          ALLOCATE(x2d(nc1,nc2),y2d(nc1,nc2),z2d(nc1,nc2),d2d(nc1,nc2))
      !          FORALL(k=1:nc1) x2d(k,:) = coil_group(i1)%coils(j1)%xnod(1,k)
      !          FORALL(k=1:nc1) y2d(k,:) = coil_group(i1)%coils(j1)%xnod(2,k)
      !          FORALL(k=1:nc1) z2d(k,:) = coil_group(i1)%coils(j1)%xnod(3,k)
      !          FORALL(k=1:nc2) x2d(:,k) = x2d(:,k) - coil_group(i1)%coils(j2)%xnod(1,k)
      !          FORALL(k=1:nc2) y2d(:,k) = y2d(:,k) - coil_group(i1)%coils(j2)%xnod(2,k)
      !          FORALL(k=1:nc2) z2d(:,k) = z2d(:,k) - coil_group(i1)%coils(j2)%xnod(3,k)
      !          d2d = x2d*x2d+y2d*y2d+z2d*z2d
      !          WHERE(d2d < 1.0E-6) d2d = 1.0E6
      !          dist = SQRT(MINVAL(d2d,DIM=2))
      !          dist_min = MIN(SQRT(MINVAL(d2d)),dist_min)
      !          DEALLOCATE(x2d,y2d,z2d,d2d)
      !       END DO
      !       IF (loutput) THEN
      !          DO k = 1, nc1
      !             WRITE(iunit_out,'(2(2X,I3),14(2X,ES22.12))') &
      !                i1,j1,coil_group(i1)%coils(j1)%xnod(1,k), &
      !                coil_group(i1)%coils(j1)%xnod(2,k), &
      !                coil_group(i1)%coils(j1)%xnod(3,k), &
      !                dist(k)
      !          END DO
      !       END IF
      !       DEALLOCATE(dist)
      !    END DO
      ! END DO
      ! !--------------------------------------------------------------
      ! !    Now we compare differnt coil groups 
      ! !--------------------------------------------------------------
      ! DO i1 = 1, ncoilgroups
      !    n1 = i1+1
      !    DO i2 = n1,ncoilgroups
      !       DO j1 = 1, coil_group(i1)%ncoil
      !          nc1 = SIZE(coil_group(i1)%coils(j1)%xnod,2)
      !          ALLOCATE(dist(nc1))
      !          DO j2 = 1, coil_group(i2)%ncoil
      !             nc2 = SIZE(coil_group(i2)%coils(j2)%xnod,2)
      !             ALLOCATE(x2d(nc1,nc2),y2d(nc1,nc2),z2d(nc1,nc2),d2d(nc1,nc2))
      !             FORALL(k=1:nc1) x2d(k,:) = coil_group(i1)%coils(j1)%xnod(1,k)
      !             FORALL(k=1:nc1) y2d(k,:) = coil_group(i1)%coils(j1)%xnod(2,k)
      !             FORALL(k=1:nc1) z2d(k,:) = coil_group(i1)%coils(j1)%xnod(3,k)
      !             FORALL(k=1:nc2) x2d(:,k) = x2d(:,k) - coil_group(i2)%coils(j2)%xnod(1,k)
      !             FORALL(k=1:nc2) y2d(:,k) = y2d(:,k) - coil_group(i2)%coils(j2)%xnod(2,k)
      !             FORALL(k=1:nc2) z2d(:,k) = z2d(:,k) - coil_group(i2)%coils(j2)%xnod(3,k)
      !             d2d = x2d*x2d+y2d*y2d+z2d*z2d
      !             WHERE(d2d < 1.0E-6) d2d = 1.0E6
      !             dist_min = MIN(SQRT(MINVAL(d2d)),dist_min)
      !             DEALLOCATE(x2d,y2d,z2d,d2d)
      !          END DO
      !       END DO
      !    END DO
      ! END DO
      ! IF (loutput) CLOSE(iunit_out)
      ! RETURN
      ! END SUBROUTINE compute_coil_coil_distance

      SUBROUTINE compute_coil_curvature(outext)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: outext
      LOGICAL :: loutput
      INTEGER :: i, j, k, l, nc, nc1, iunit_out, ier, ntotal_coils
      DOUBLE PRECISION :: hs
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
         xc, yc, zc, xcp, ycp, zcp, xcpp, ycpp, zcpp, &
         xcppp, ycppp, zcppp
      loutput = .FALSE.
      IF (PRESENT(outext)) THEN
         loutput = .TRUE.
         CALL safe_open(iunit_out,ier,TRIM('coil_curvature.'//TRIM(outext)),'unknown','formatted')
         WRITE(iunit_out,'(I6,2X,I6,2X,I6,2X,I6)') ncoilgroups,nw_coil,nh_coil,ns
      END IF
      IF (ALLOCATED(curvature)) DEALLOCATE(curvature)
      IF (ALLOCATED(torsion)) DEALLOCATE(torsion)
      IF (ALLOCATED(dLength)) DEALLOCATE(dLength)
      i = ncoilgroups*nw_coil*nh_coil
      ALLOCATE(curvature(i,ns))
      ALLOCATE(torsion(i,ns))
      ALLOCATE(dLength(i,ns))
      curvature = zero
      torsion   = zero
      curvature_min = 0.0; curvature_max = 0.0; curvature_mean = 0.0
      torsion_min = 0.0; torsion_max = 0.0; torsion_mean = 0.0
      ntotal_coils = 0
      DO i = 1, ncoilgroups
         DO j = 1, nw_coil*nh_coil
            ntotal_coils = ntotal_coils + 1
            !nc = SIZE(coil_group(i)%coils(j)%xnod,2)
            nc = ns
            nc1 = nc - 1
            hs  = 1.0D+00/nc1
            ALLOCATE(xc(nc),yc(nc),zc(nc))
            ALLOCATE(xcp(nc),ycp(nc),zcp(nc))
            ALLOCATE(xcpp(nc),ycpp(nc),zcpp(nc))
            ALLOCATE(xcppp(nc),ycppp(nc),zcppp(nc))
            xc = coil_group(i)%coils(j)%xnod(1,:)
            yc = coil_group(i)%coils(j)%xnod(2,:)
            zc = coil_group(i)%coils(j)%xnod(3,:)
            xcp(1:nc1) = xc(2:nc) - xc(1:nc1)
            ycp(1:nc1) = yc(2:nc) - yc(1:nc1)
            zcp(1:nc1) = zc(2:nc) - zc(1:nc1)
            ! Note that x(1)=x(nc) so we need
            xcp(nc) = xcp(1)
            ycp(nc) = ycp(1)
            zcp(nc) = zcp(1)
            xcp = xcp * hs; ycp = ycp * hs; zcp = zcp * hs;
            xcpp(1:nc1) = xcp(2:nc) - xcp(1:nc1)
            ycpp(1:nc1) = ycp(2:nc) - ycp(1:nc1)
            zcpp(1:nc1) = zcp(2:nc) - zcp(1:nc1)
            xcpp(nc) = xcpp(1)
            ycpp(nc) = ycpp(1)
            zcpp(nc) = zcpp(1)
            xcpp = xcpp * hs; ycpp = ycpp * hs; zcpp = zcpp * hs;
            xcppp(1:nc1) = xcpp(2:nc) - xcpp(1:nc1)
            ycppp(1:nc1) = ycpp(2:nc) - ycpp(1:nc1)
            zcppp(1:nc1) = zcpp(2:nc) - zcpp(1:nc1)
            xcppp(nc) = xcppp(1)
            ycppp(nc) = ycppp(1)
            zcppp(nc) = zcppp(1)
            xcppp = xcppp * hs; ycppp = ycppp * hs; zcppp = zcppp * hs;
            dLength(ntotal_coils,:) = SQRT(xcp*xcp+ycp*ycp+zcp*zcp)/hs
            curvature(ntotal_coils,:)   = &
                 SQRT((zcpp*ycp-ycpp*zcp)**2 &
                    + (xcpp*zcp-zcpp*xcp)**2 &
                    + (ycpp*xcp-xcpp*ycp)**2) &
                    / (xcp*xcp+ycp*ycp+zcp*zcp)**(3.0/2.0)
            torsion(ntotal_coils,:) = &
                      ((zcpp*ycp-ycpp*zcp)*xcppp &
                    +  (xcpp*zcp-zcpp*xcp)*ycppp &
                    +  (ycpp*xcp-xcpp*ycp)*zcppp) &
                    / ((zcpp*ycp-ycpp*zcp)**2 &
                    +  (xcpp*zcp-zcpp*xcp)**2 &
                    +  (ycpp*xcp-xcpp*ycp)**2)
            curvature_min = curvature_min + MINVAL(curvature(ntotal_coils,:))
            curvature_max = curvature_max + MAXVAL(curvature(ntotal_coils,:))
            curvature_mean = curvature_mean + SUM(curvature(ntotal_coils,:))/nc
            torsion_min = torsion_min + MINVAL(torsion(ntotal_coils,:))
            torsion_max = torsion_max + MAXVAL(torsion(ntotal_coils,:))
            torsion_mean = torsion_mean + SUM(torsion(ntotal_coils,:))/nc
            IF (loutput) THEN
               DO k = 1, nc
                  WRITE(iunit_out,'(2(2X,I3),14(2X,ES22.12))') &
                     i,j,xc(k),yc(k),zc(k),xcp(k),ycp(k),zcp(k),&
                     xcpp(k),ycpp(k),zcpp(k),xcppp(k),ycppp(k),zcppp(k), &
                     curvature(ntotal_coils,k),torsion(ntotal_coils,k)
               END DO
            END IF
            DEALLOCATE(xc,yc,zc)
            DEALLOCATE(xcp,ycp,zcp)
            DEALLOCATE(xcpp,ycpp,zcpp)
            DEALLOCATE(xcppp,ycppp,zcppp)
         END DO
      END DO
      curvature_min  = curvature_min/ntotal_coils
      curvature_max  = curvature_max/ntotal_coils
      curvature_mean = curvature_mean/ntotal_coils
      torsion_min    = torsion_min/ntotal_coils
      torsion_max    = torsion_max/ntotal_coils
      torsion_mean   = torsion_mean/ntotal_coils
      IF (loutput) CLOSE(iunit_out)
      RETURN
      END SUBROUTINE compute_coil_curvature

      SUBROUTINE compute_coil_energy(nw,nh,width,height,icoil,E)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nw
      INTEGER, INTENT(IN) :: nh
      DOUBLE PRECISION, INTENT(IN) :: width
      DOUBLE PRECISION, INTENT(IN) :: height
      INTEGER, INTENT(IN) :: icoil
      DOUBLE PRECISION, INTENT(out) :: E
      INTEGER :: i, ic, ix, j, jc, jx, nc, nc1
      DOUBLE PRECISION :: hs, Lij, &
         xi, yi, zi, xip, yip, zip, &
         xj, yj, zj, xjp, yjp, zjp, &
         dx, dy, dz, delta, a, b, aob, boa, kcoef
      E=0.0D+00
      i = icoil
      nc = ns
      nc1 = nc - 1
      hs = 1.0E+00/nc1
      ! Compute the normalizaton factor
      a = width / nw_coil
      b = height / nh_coil
      aob = a/b
      boa = b/a
      !kcoef = (4.0/3.0)*boa*ATAN(aob) + (4.0/3.0)*aob*ATAN(boa) + (1.0/6.0)*boa*boa*LOG(boa) + (1.0/6.0)*aob*aob*LOG(aob)
      kcoef = (8.0*boa*ATAN(aob) + 8.0*aob*ATAN(boa) + boa*boa*LOG(boa) + aob*aob*LOG(aob))/6.0
      kcoef = kcoef - (a**4 - 6*a*a*b*b + b**4)*LOG(aob+boa)/(6*a*a*b*b)
      delta = EXP(kcoef - 25.0/6.0)
      ! Compute the Mutual inductance
      DO ic = 1, nw_coil*nh_coil
         DO ix = 1, nc1
            xi  = 0.5*(coil_group(icoil)%coils(ic)%xnod(1,ix) + coil_group(icoil)%coils(ic)%xnod(1,ix+1))
            yi  = 0.5*(coil_group(icoil)%coils(ic)%xnod(2,ix) + coil_group(icoil)%coils(ic)%xnod(2,ix+1))
            zi  = 0.5*(coil_group(icoil)%coils(ic)%xnod(3,ix) + coil_group(icoil)%coils(ic)%xnod(3,ix+1))
            xip = (coil_group(icoil)%coils(ic)%xnod(1,ix+1) - coil_group(icoil)%coils(ic)%xnod(1,ix))
            yip = (coil_group(icoil)%coils(ic)%xnod(2,ix+1) - coil_group(icoil)%coils(ic)%xnod(2,ix))
            zip = (coil_group(icoil)%coils(ic)%xnod(3,ix+1) - coil_group(icoil)%coils(ic)%xnod(3,ix))
            DO j = 1, ncoilgroups
               Lij = 0.0D+00
               DO jc = 1, nw_coil*nh_coil
                  DO jx = 1, nc1
                     xj  = 0.5*(coil_group(j)%coils(jc)%xnod(1,jx) + coil_group(j)%coils(jc)%xnod(1,jx+1))
                     yj  = 0.5*(coil_group(j)%coils(jc)%xnod(2,jx) + coil_group(j)%coils(jc)%xnod(2,jx+1))
                     zj  = 0.5*(coil_group(j)%coils(jc)%xnod(3,jx) + coil_group(j)%coils(jc)%xnod(3,jx+1))
                     xjp = (coil_group(j)%coils(jc)%xnod(1,jx+1) - coil_group(j)%coils(jc)%xnod(1,jx))
                     yjp = (coil_group(j)%coils(jc)%xnod(2,jx+1) - coil_group(j)%coils(jc)%xnod(2,jx))
                     zjp = (coil_group(j)%coils(jc)%xnod(3,jx+1) - coil_group(j)%coils(jc)%xnod(3,jx))
                     dx = xi - xj
                     dy = yi - yj
                     dz = zi - zj
                     Lij   = Lij + (xip * xjp + yip * yjp + zip * zjp) / SQRT(dx*dx+dy*dy+dz*dz + delta*a*b)
                  END DO
               END DO
               E = E + Lij * coil_group(j)%coils(1)%current
            END DO
         END DO
      END DO
      E = 0.5 * E * coil_group(icoil)%coils(1)%current * 1.0E-7 * pi2 * pi2 * hs *hs
      RETURN
      END SUBROUTINE compute_coil_energy

      SUBROUTINE get_coil_ns(ns_out)
      IMPLICIT NONE
      INTEGER, INTENT(out) :: ns_out
      ns_out = ns
      RETURN
      END SUBROUTINE get_coil_ns

      SUBROUTINE get_coil_curvature(coil_filament,coil_seg,curvature_out)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: coil_filament
      INTEGER, INTENT(in) :: coil_seg
      DOUBLE PRECISION, INTENT(out) :: curvature_out
      curvature_out = curvature(coil_filament,coil_seg)
      RETURN
      END SUBROUTINE get_coil_curvature

      SUBROUTINE get_coil_curvature_avg(curve_mean,curve_max,curve_min)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(out) :: curve_mean
      DOUBLE PRECISION, INTENT(out) :: curve_max
      DOUBLE PRECISION, INTENT(out) :: curve_min
      curve_mean = curvature_mean
      curve_max  = curvature_max
      curve_min  = curvature_min
      RETURN
      END SUBROUTINE get_coil_curvature_avg

      SUBROUTINE get_coil_torsion(coil_filament,coil_seg,torsion_out)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: coil_filament
      INTEGER, INTENT(in) :: coil_seg
      DOUBLE PRECISION, INTENT(out) :: torsion_out
      torsion_out = torsion(coil_filament,coil_seg)
      RETURN
      END SUBROUTINE get_coil_torsion

      SUBROUTINE get_coil_torsion_avg(tor_mean,tor_max,tor_min)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(out) :: tor_mean
      DOUBLE PRECISION, INTENT(out) :: tor_max
      DOUBLE PRECISION, INTENT(out) :: tor_min
      tor_mean = torsion_mean
      tor_max  = torsion_max
      tor_min  = torsion_min
      RETURN
      END SUBROUTINE get_coil_torsion_avg

      SUBROUTINE get_coil_dl(coil_filament,coil_seg,dl_out)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: coil_filament
      INTEGER, INTENT(in) :: coil_seg
      DOUBLE PRECISION, INTENT(out) :: dl_out
      dl_out = dLength(coil_filament,coil_seg)
      RETURN
      END SUBROUTINE get_coil_dl

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE spline_coils_mod
      