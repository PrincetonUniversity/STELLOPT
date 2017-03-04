!-----------------------------------------------------------------------
!     Subroutine:    toroidal_flux
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/14/2011
!     Description:   This subroutine calculates the toroidal flux.
!-----------------------------------------------------------------------
      SUBROUTINE toroidal_flux
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile
!-----------------------------------------------------------------------
!     Local Variables
!          ier         Error flag
!          ik          Radial dummy index
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          nu_calc     Number of Poloidal points to use
!          nv_calc     Number of Toroidal points to use
!          pi2         2*pi
!          r1          Inner semi-minor R
!          r2          Outer semi-minor R
!          z1          Inner semi-minor Z
!          z2          Outer semi-minor Z
!          rho1        Inner semi-minor radius
!          rho2        Outer semi-minor radius
!          drho        rho2-rho1
!          dtheta      2*pi/nu_calc
!          xu_calc     Poloidal index array
!          xv_calc     Toroidal index array
!          r_temp      Realspace R
!          z_temp      Realspace Z
!          bphi_temp   Realspace Toroidal Field
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, ik, u, v, nu_calc, nv_calc, mn
      REAL(rprec) :: pi2,r1,z1,rho1,r2,z2,rho2,drho,dtheta,jac
      REAL(rprec) :: fmn_temp(1:mnmax,0:k)
      REAL(rprec), ALLOCATABLE :: xu_calc(:), xv_calc(:)
      REAL(rprec), ALLOCATABLE :: r_temp(:,:,:),z_temp(:,:,:), bphi_temp(:,:,:)
      REAL(rprec), ALLOCATABLE :: rs_temp(:,:,:),ru_temp(:,:,:), zs_temp(:,:,:), zu_temp(:,:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      ! Default some values
      pi2 = 8 * ATAN(1._rprec)
      nu_calc = nu*10
      nv_calc = 1
      ALLOCATE(xu_calc(1:nu_calc),xv_calc(1:nv_calc),STAT=ier)
      FORALL(u=1:nu_calc) xu_calc(u) = REAL(u-1)/REAL(nu_calc)
      xv_calc(1) = 0.0
      ! Now fourier transform R, Z, and Bphi on this
      ALLOCATE(r_temp(0:k,1:nu_calc,1:nv_calc),&
               z_temp(0:k,1:nu_calc,1:nv_calc),&
               bphi_temp(0:k,1:nu_calc,1:nv_calc),&
               rs_temp(0:k,1:nu_calc,1:nv_calc),&
               ru_temp(0:k,1:nu_calc,1:nv_calc),&
               zs_temp(0:k,1:nu_calc,1:nv_calc),&
               zu_temp(0:k,1:nu_calc,1:nv_calc),STAT=ier)
      r_temp = 0.0
      z_temp = 0.0
      bphi_temp = 0.0
      CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu_calc,xv_calc,rmnc,xm,xn,r_temp,0,1)
      CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu_calc,xv_calc,zmns,xm,xn,z_temp,1,0)
      CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu_calc,xv_calc,bvmnc,xm,xn,bphi_temp,0,0)
      FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = -rmnc(mn,0:k)*xm(mn)
      CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu,xv,fmn_temp,xm,xn,ru_temp,1,0) 
      FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = zmns(mn,0:k)*xm(mn)
      CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu,xv,fmn_temp,xm,xn,zu_temp,0,0) 
      IF (lasym) THEN
         CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu_calc,xv_calc,rmns,xm,xn,r_temp,0,1)
         CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu_calc,xv_calc,zmnc,xm,xn,z_temp,1,0)
         CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu_calc,xv_calc,bvmns,xm,xn,bphi_temp,1,0)
         FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = rmns(mn,0:k)*xm(mn)
         CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu,xv,fmn_temp,xm,xn,ru_temp,0,0) 
         FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = -zmnc(mn,0:k)*xm(mn)
         CALL mntouv(0,k,mnmax,nu_calc,nv_calc,xu,xv,fmn_temp,xm,xn,zu_temp,1,0) 
      END IF
      DO u=1,nu_calc
         DO v=1,nv_calc
            rs_temp(1:k-1,u,v) = 0.5*(r_temp(2:k,u,v)-r_temp(0:k-2,u,v))
            zs_temp(1:k-1,u,v) = 0.5*(z_temp(2:k,u,v)-z_temp(0:k-2,u,v))
            rs_temp(0,u,v)     = (r_temp(1,u,v) - r_temp(0,u,v))
            zs_temp(0,u,v)     = (z_temp(1,u,v) - z_temp(0,u,v))
            rs_temp(k,u,v)     = (r_temp(k,u,v) - r_temp(k-1,u,v))
            zs_temp(k,u,v)     = (z_temp(k,u,v) - z_temp(k-1,u,v))
         END DO
      END DO
      ! Now Bphi = R*Bv*dv/dphi
      bphi_temp = bphi_temp * r_temp
      ! Now we integrate
      IF (ALLOCATED(torflux_pies)) DEALLOCATE(torflux_pies)
      IF (ALLOCATED(torflux_norm)) DEALLOCATE(torflux_norm)
      ALLOCATE(torflux_pies(0:k),torflux_norm(0:k),STAT=ier)
      torflux_pies = 0.0
      dtheta = pi2 / nu_calc
      torflux_pies(0) = 0.0
      DO v = 1, nv_calc
         DO u = 1, nu_calc
            r1 = r_temp(1,u,v)-r_temp(0,u,v)
            z1 = z_temp(1,u,v)-z_temp(0,u,v)
            rho1 = sqrt(r1*r1+z1*z1)
            jac = rs_temp(1,u,v)*zu_temp(1,u,v)-ru_temp(1,u,v)*zs_temp(1,u,v)
            torflux_pies(1) = torflux_pies(1) + bphi_temp(0,u,v)*rho1*rho1*dtheta*jac
         END DO
      END DO
      DO ik = 2, k
         torflux_pies(ik) = torflux_pies(ik-1)
         DO v = 1, nv_calc
            DO u = 1, nu_calc
               r1 = r_temp(ik,u,v)-r_temp(0,u,v)
               z1 = z_temp(ik,u,v)-z_temp(0,u,v)
               rho1 = sqrt(r1*r1+z1*z1)
               r2 = r_temp(ik-1,u,v)-r_temp(0,u,v)
               z2 = z_temp(ik-1,u,v)-z_temp(0,u,v)
               rho2 = sqrt(r2*r2+z2*z2)
               drho = rho1 - rho2
               jac = rs_temp(ik-1,u,v)*zu_temp(ik-1,u,v)-ru_temp(ik-1,u,v)*zs_temp(ik-1,u,v)
               torflux_pies(ik) = torflux_pies(ik) + bphi_temp(ik,u,v)*rho2*drho*dtheta*jac
            END DO
         END DO
      END DO
      torflux_norm = torflux_pies/torflux_edge
      ! DEALLOCATE ARRAYS
      DEALLOCATE(xu_calc,xv_calc)
      DEALLOCATE(r_temp,z_temp,bphi_temp)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE toroidal_flux
