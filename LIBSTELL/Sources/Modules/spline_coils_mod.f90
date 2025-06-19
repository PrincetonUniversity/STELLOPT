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
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE :: xm, xn, rmnc, zmns, t_kts
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: rho_kts, theta_kts, zeta_kts
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: curvature, torsion, dLength
      
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
      RETURN
      END SUBROUTINE init_spline_coils

      SUBROUTINE init_boundary_spline_coils(mnmax_in,xm_in,xn_in,rmnc_in,zmns_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mnmax_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: xm_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: xn_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: rmnc_in
      DOUBLE PRECISION, DIMENSION(mnmax_in), INTENT(in) :: zmns_in
      mnmax = mnmax_in
      IF (ALLOCATED(xm)) DEALLOCATE(xm)
      IF (ALLOCATED(xn)) DEALLOCATE(xn)
      IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc)
      IF (ALLOCATED(zmns)) DEALLOCATE(zmns)
      ALLOCATE(xm(mnmax), xn(mnmax), rmnc(mnmax), zmns(mnmax))
      xm = xm_in
      xn = xn_in
      rmnc = rmnc_in
      zmns = zmns_in
      nfp = MINVAL(xn, MASK = xn > 0)
      xn = xn / nfp
      factor = pi2/nfp
      RETURN
      END SUBROUTINE init_boundary_spline_coils

      SUBROUTINE spline_to_coils(normal_sign)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: normal_sign
      INTEGER :: i, j, mn, ns1
      DOUBLE PRECISION :: AX, AY, AZ, BX, BY, BZ, NX, NY, NZ, N, &
            R, Z, RU, ZU, RV, ZV, rho, theta, zeta, cop, sip, l, &
            X, Y, phi
      DOUBLE PRECISION, DIMENSION(ns) :: Rc,Zc,Pc
      DOUBLE PRECISION, DIMENSION(3,ns) :: xnod_in, xnod_ss, xnod_bb
      CHARACTER(len=100) :: s_name
      CHARACTER(len=100) :: l_name
      CHARACTER(len=100) :: c_name
      TYPE(bsc_coil)     :: coil_temp
      TYPE(bsc_rs)       :: rot_mat
      !INTERFACE
      !   REAL FUNCTION bvalue( t, bcoef, n, k, x, jderiv )
      !      integer jderiv,k,n
      !      double precision bcoef(n),t(n+k),x
      !   END FUNCTION bvalue
      !END INTERFACE
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
         CALL bsc_construct_coilcoll(coil_group(i),TRIM(c_name),l_name)
         DO j = 1, ns
            l = DBLE(j-1)/DBLE(ns-1)
            CALL spline_it(n_kts,t_kts,rho_kts(i,:),1,l,rho,0)
            CALL spline_it(n_kts,t_kts,theta_kts(i,:),1,l,theta,0)
            CALL spline_it(n_kts,t_kts,zeta_kts(i,:),1,l,zeta,0)
            !PRINT *,rho,theta,zeta
            !rho = bvalue(t_kts,rho_kts(i,:),n_kts,k_kts,l,0)
            !theta = bvalue(t_kts,theta_kts(i,:),n_kts,k_kts,l,0)
            !zeta = bvalue(t_kts,zeta_kts(i,:),n_kts,k_kts,l,0)
            R = zero; Z = zero; RU = zero; ZU = zero; RV = zero; ZV=zero
            phi = zeta/nfp
            DO mn = 1, mnmax
               cop = cos(xm(mn)*theta+xn(mn)*zeta)
               sip = sin(xm(mn)*theta+xn(mn)*zeta)
               R  =  R + rmnc(mn)*cop
               Z  =  Z + zmns(mn)*sip
               RU = RU - rmnc(mn)*sip*xm(mn)
               ZU = ZU + zmns(mn)*cop*xm(mn)
               RV = RV - rmnc(mn)*sip*xn(mn) ! dR/dzeta
               ZV = ZV + zmns(mn)*cop*xn(mn) ! dZ/dzeta
            END DO
            cop = cos(phi)
            sip = sin(phi)
            Ax = RU * cop; Ay = RU * sip; Az = ZU
            ! dR/dzeta
            Bx = RV * cop - R * sip/nfp; By = RV * sip + R * cop/nfp; Bz = ZV
            Nx = Ay*Bz - Az*By
            Ny = Az*Bx - Ax*Bz
            Nz = Ax*By - Ay*Bx
            N  = SQRT(Nx*Nx+Ny*Ny+Nz*Nz)*normal_sign
            Nx = Nx/N; Ny = Ny/N; Nz = Nz/N
            X  = R*cop + rho*Nx
            Y  = R*sip + rho*Ny
            Z  = Z    + rho*Nz
            Rc(j) = SQRT(X*X + Y*Y)
            Zc(j) = Z
            Pc(j) = ATAN2(Y,X)
            xnod_in(1,j) = X
            xnod_in(2,j) = Y 
            xnod_in(3,j) = Z 
         END DO
         xnod_in(:,ns) = xnod_in(:,1)
         ! Now create the first coil
         WRITE(s_name, '(a4,i5.5)') 'ID #', 1
         CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_in(1:3,1:ns))
         CALL bsc_append(coil_group(i),coil_temp)
         ! Now create the stellarator symmetric coil
         Zc = -Zc
         Pc = factor - Pc
         xnod_ss(1,2:ns) = Rc(ns1:1:-1)*cos(Pc(ns1:1:-1))
         xnod_ss(2,2:ns) = Rc(ns1:1:-1)*sin(Pc(ns1:1:-1))
         xnod_ss(3,2:ns) = Zc(ns1:1:-1)
         xnod_ss(:,1) = xnod_ss(:,ns)
         CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_ss(1:3,1:ns))
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
         CALL bsc_construct_coilcoll(coil_single(i),coil_group(i)%s_name,coil_group(i)%l_name)
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
         CALL bsc_construct_coilcoll(coil_group(i),coil_group(i)%s_name,coil_group(i)%l_name)
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
                  xnod(1,:) = xnod(1,:) - xb*width/2 - xn*height/2 &
                            + xb*width*(l-1)/(nw-1) + xn*height*(k-1)/(nh-1)
                  xnod(2,:) = xnod(2,:) - yb*width/2 - yn*height/2 &
                            + yb*width*(l-1)/(nw-1) + yn*height*(k-1)/(nh-1)
                  xnod(3,:) = xnod(3,:) - zb*width/2 - zn*height/2 &
                            + zb*width*(l-1)/(nw-1) + zn*height*(k-1)/(nh-1)
                  xnod(:,ns) = xnod(:,1)
                  CALL bsc_construct_coil(coil_temp,'fil_loop',coil_single(i)%coils(j)%s_name,'',coil_single(i)%coils(j)%current/(nh*nw),xnod(1:3,1:ns))
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
            dLength(ntotal_coils,:) = SQRT(xcp*xcp+ycp*ycp+zcp*zcp)
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
      