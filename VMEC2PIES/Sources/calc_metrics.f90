!-----------------------------------------------------------------------
!     Subroutine:    calc_metrics
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/15/2011
!     Description:   This subroutine calculates the metric elements from
!                    the values currently stored in the background
!                    and realspace grids.
!
!                    We define our fuctions to be of the form
!                    f(rho,theta,phi) = sum{n=-ntor,ntor} sum{m=0,mpol) fmn(s) * cos(m*theta+n*N*phi)
!                    So we must chain rule to u and v where
!                    u = theta/2*pi  and  v = phi / 2*pi*N
!                    Thus
!                    df/du = df/dtheta * dtheta/du
!                    df/dv = df/dphi   * dphi/dv
!-----------------------------------------------------------------------
      SUBROUTINE calc_metrics
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY : rprec
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile, ONLY : torflux_norm
!-----------------------------------------------------------------------
!     Local Variables
!          ier         Error flag
!          u           Poloidal dummy index
!          v           Toroidal dummy index
!          mn          Fourier dummy index
!          factor      ratio dphi/dv
!          fmn_temp    Dummy array to hold derivative coefficients (ex. rmnc*m)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, u, v, mn
      REAL(rprec) :: factor,pi2
      REAL(rprec), ALLOCATABLE :: fmn_temp(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      ! Get Toroidal Flux
      CALL toroidal_flux
      ! Calculate deriviatives
      ALLOCATE(fmn_temp(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'FMN_TEMP (calc_metrics)',ier)
      rs = 0; zs = 0; ru=0; zu=0; rv=0; zv=0; 
      FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = -rmnc(mn,0:k)*xm(mn)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,ru,1,1)
      FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = -rmnc(mn,0:k)*xn(mn)*nfp
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,rv,1,0)  
      FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = zmns(mn,0:k)*xm(mn)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zu,0,0)  
      FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = zmns(mn,0:k)*xn(mn)*nfp
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zv,0,0)  
      IF (lasym) THEN
         FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = rmns(mn,0:k)*xm(mn)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,ru,0,0)
         FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = rmns(mn,0:k)*xn(mn)*nfp
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,rv,0,0)  
         FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = -zmnc(mn,0:k)*xm(mn)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zu,1,0)  
         FORALL(mn = 1:mnmax) fmn_temp(mn,0:k) = -zmnc(mn,0:k)*xn(mn)*nfp
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zv,1,0)  
      END IF
      DO u=1,nu
         DO v=1,nv
            rs(1:k-1,u,v) = (rreal(2:k,u,v)-rreal(0:k-2,u,v))/(rho(2:k)-rho(0:k-2))
            zs(1:k-1,u,v) = (zreal(2:k,u,v)-zreal(0:k-2,u,v))/(rho(2:k)-rho(0:k-2))
            rs(0,u,v)     = (rreal(1,u,v) - rreal(0,u,v))/(rho(1)-rho(0))
            zs(0,u,v)     = (zreal(1,u,v) - zreal(0,u,v))/(rho(1)-rho(0))
            rs(k,u,v)     = (rreal(k,u,v) - rreal(k-1,u,v))/(rho(k)-rho(k-1))
            zs(k,u,v)     = (zreal(k,u,v) - zreal(k-1,u,v))/(rho(k)-rho(k-1))
            !rs(1:k-1,u,v) = (rreal(2:k,u,v)-rreal(0:k-2,u,v))/(torflux_norm(2:k)-torflux_norm(0:k-2))
            !zs(1:k-1,u,v) = (zreal(2:k,u,v)-zreal(0:k-2,u,v))/(torflux_norm(2:k)-torflux_norm(0:k-2))
            !rs(0,u,v)     = (rreal(1,u,v) - rreal(0,u,v))/(torflux_norm(1)-torflux_norm(0))
            !zs(0,u,v)     = (zreal(1,u,v) - zreal(0,u,v))/(torflux_norm(1)-torflux_norm(0))
            !rs(k,u,v)     = (rreal(k,u,v) - rreal(k-1,u,v))/(torflux_norm(k)-torflux_norm(k-1))
            !zs(k,u,v)     = (zreal(k,u,v) - zreal(k-1,u,v))/(torflux_norm(k)-torflux_norm(k-1))
         END DO
      END DO
      ! Calculate metric elements
      g11 = rs*rs+zs*zs
      g12 = rs*ru+zs*zu
      g13 = rs*rv+zs*zv
      g22 = ru*ru+zu*zu
      g23 = ru*rv+zu*zv
      g33 = rv*rv+zv*zv + rreal*rreal
      detg= g11*(g22*g33-g23*g23) + g12*(g23*g13-g12*g33) + g13*(g12*g23-g22*g13)
      DEALLOCATE(fmn_temp)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE calc_metrics
