!-----------------------------------------------------------------------
!     Function:      fblin
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following
!-----------------------------------------------------------------------
      SUBROUTINE fblin(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_fieldlines
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (xsi,eta)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(2), qdot(2)
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier, ik
      REAL(rprec) :: beta_temp, bxsi_temp
      REAL(rprec) :: rho_local, theta_local, phi_local
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, PARAMETER :: ict(8)=(/1,0,0,0,0,0,0,0/)
      !INTEGER, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      rho_local = q(1)
      theta_local = q(2)
      phi_local = phi
      IF (rho_local > 1) rho_local = 1
      ! Handle Angles
      theta_local = MOD(theta_local,thmx)
      phi_local = MOD(phi,phmx)
      IF (theta_local < 0) theta_local = theta_local + thmx
      IF (phi_local < 0) phi_local = phi_local + phmx
      qdot = 0
      IF ((rho_local >= 0-eps1) .and. (rho_local <= 1+eps1) .and. &
          (theta_local >= 0-eps2) .and. (theta_local <= pi2+eps2) .and. &
          (phi_local >= 0-eps3) .and. (phi_local <= phmx+eps3)) THEN
         ! Get the gridpoint info
         i = MIN(MAX(COUNT(rho < rho_local),1),nrho-1)
         j = MIN(MAX(COUNT(xu < theta_local),1),nu-1)
         k = MIN(MAX(COUNT(xv < phi_local),1),nv-1)
         hx     = rho(i+1) - rho(i)
         hy     = xu(j+1) - xu(j)
         hz     = xv(k+1) - xv(k)
         hxi    = one / hx
         hyi    = one / hy
         hzi    = one / hz
         xparam = (rho_local - rho(i)) * hxi
         yparam = (theta_local - xu(j)) * hyi
         zparam = (phi_local - xv(k)) * hzi
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BXSI4D(1,1,1,1),nrho,nu,nv)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                BXSI4D(1,1,1,1),nrho,nu,nv)
         bxsi_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BETA4D(1,1,1,1),nrho,nu,nv)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                BETA4D(1,1,1,1),nrho,nu,nv)
         beta_temp = fval(1)
      END IF

      ! Spline Bxsi
      !CALL EZspline_interp(bxsi_spl,rho_local,theta_local,phi_local,bxsi_temp,ier)
      ! Spline Beta
      !CALL EZspline_interp(beta_spl,rho_local,theta_local,phi_local,beta_temp,ier)
      !IF (rho_local == 1) THEN
      !   qdot(1) = 0
      !   qdot(2) = 1
      !ELSE
         qdot(1) = bxsi_temp
         qdot(2) = beta_temp
      !END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin
