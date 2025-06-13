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
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER, PRIVATE :: pi2 = 6.283185482025146D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: zero = 0.0D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: one = 1.0D+00
      INTEGER, PRIVATE :: mnmax, ns, ncoilgroups, n_kts, k_kts
      DOUBLE PRECISION, PRIVATE :: factor
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE :: xm, xn, rmnc, zmns, t_kts
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: rho_kts, theta_kts, zeta_kts
      
!-----------------------------------------------------------------------
!     Module SUBROUTINES/FUNCTIONS
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE init_spline_coils(ns_in, ncoilgroups_in, n_in, nk_in, &
            rho_in, theta_in, zeta_in, t_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ns_in
      INTEGER, INTENT(in) :: ncoilgroups_in
      INTEGER, INTENT(in) :: n_in
      INTEGER, INTENT(in) :: nk_in
      DOUBLE PRECISION, INTENT(in) :: rho_in(ncoilgroups_in,n_in)
      DOUBLE PRECISION, INTENT(in) :: theta_in(ncoilgroups_in,n_in)
      DOUBLE PRECISION, INTENT(in) :: zeta_in(ncoilgroups_in,n_in)
      DOUBLE PRECISION, INTENT(in) :: t_in(n_in)
      ns = ns_in
      ncoilgroups = ncoilgroups_in
      n_kts = n_in
      k_kts = nk_in-n_in
      IF (ALLOCATED(rho_kts)) DEALLOCATE(rho_kts)
      IF (ALLOCATED(theta_kts)) DEALLOCATE(theta_kts)
      IF (ALLOCATED(zeta_kts)) DEALLOCATE(zeta_kts)
      IF (ALLOCATED(t_kts)) DEALLOCATE(t_kts)
      ALLOCATE(rho_kts(ncoilgroups,n_kts), theta_kts(ncoilgroups,n_kts), &
            zeta_kts(ncoilgroups,n_kts), t_kts(n_kts))
      rho_kts = rho_in
      theta_kts = theta_in
      zeta_kts = zeta_in
      t_kts = t_in
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
      TYPE(bsc_coil)     :: coil_temp
      TYPE(bsc_rs)       :: rot_mat
      !INTERFACE
      !   REAL FUNCTION bvalue( t, bcoef, n, k, x, jderiv )
      !      integer jderiv,k,n
      !      double precision bcoef(n),t(n+k),x
      !   END FUNCTION bvalue
      !END INTERFACE
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
         CALL bsc_construct_coilcoll(coil_group(i),'MODULAR_COIL',l_name)
         DO j = 1, ns1
            l = DBLE(j-1)/DBLE(ns)
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
         ! Here is where we'd put code to make it into a multi-filament coil maybe
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
                  CALL bsc_construct_coil(coil_temp,'fil_loop',coil_single(i)%coils(j)%s_name,'',one,xnod(1:3,1:ns))
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

      END SUBROUTINE coils_to_multifilament

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE spline_coils_mod