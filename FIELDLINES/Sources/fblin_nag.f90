!-----------------------------------------------------------------------
!     Function:      fblin_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    field line following (NAG)
!-----------------------------------------------------------------------
      SUBROUTINE fblin_nag(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
!      USE EZspline_obj
!      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(2), qdot(2)
!-----------------------------------------------------------------------
!     Local Variables
!          ier        Error flag
!          r_temp     R
!          z_temp     Z
!          phi_temp   phi (radians)
!          br_temp    br/bphi evaluated from splines
!          bz_temp    bz/bphi evaluated from splines
!          br_dot     Random force for diffusion
!-----------------------------------------------------------------------
      INTEGER :: ier
      DOUBLE PRECISION :: r_temp, phi_temp, z_temp, br_temp, bz_temp
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      r_temp = q(1)
      z_temp = q(2)
      phi_temp = MOD(phi,delta_phi)
      IF (phi_temp < 0) phi_temp = delta_phi + phi_temp
      br_temp = 0.0; bz_temp = 0.0
      !CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
      !IF (ier == 0) THEN
      IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
          (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
          (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
         ! Get the gridpoint info
         i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
         j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
         k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
         hx     = raxis(i+1) - raxis(i)
         hy     = phiaxis(j+1) - phiaxis(j)
         hz     = zaxis(k+1) - zaxis(k)
         hxi    = one / hx
         hyi    = one / hy
         hzi    = one / hz
         xparam = (r_temp - raxis(i)) * hxi
         yparam = (phi_temp - phiaxis(j)) * hyi
         zparam = (z_temp - zaxis(k)) * hzi
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1)
      END IF
      qdot(1) = br_temp
      qdot(2) = bz_temp
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fblin_nag
