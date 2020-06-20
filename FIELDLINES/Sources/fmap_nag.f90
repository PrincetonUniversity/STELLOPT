!-----------------------------------------------------------------------
!     Function:      fmap_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          09/01/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    return map calculation.
!-----------------------------------------------------------------------
      SUBROUTINE fmap_nag(phi,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE fieldlines_runtime, ONLY: lmu, mu, ladvanced
!      USE EZspline_obj
!      USE EZspline
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(6), qdot(6)
!-----------------------------------------------------------------------
!     Local Variables
!          ier        Error flag
!          r_temp     R
!          z_temp     Z
!          phi_temp   phi (radians)
!          br_temp    br/bphi evaluated from splines
!          bz_temp    bz/bphi evaluated from splines
!          br_r       dR_dot/dR for return map
!          br_z       dR_dot/dZ for return map
!          bz_r       dZ_dot/dR for return map
!          bz_z       dZ_dot/dZ for return map
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, br_temp, bz_temp
      REAL(rprec) :: br_r,br_z,bz_r,bz_z
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1,3) ! So weird behavior but this must match the sum of ict
      INTEGER, parameter :: ict(8)=(/1,1,0,1,0,0,0,0/)
      !INTEGER, parameter :: ict(10)=(/1,1,0,1,0,0,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      r_temp = q(1)
      z_temp = q(2)
      phi_temp = MOD(phi,delta_phi)
      IF (phi_temp < 0) phi_temp = delta_phi + phi_temp
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
         ! Evaluate the Splines
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BR4D(1,1,1,1),nr,nphi,nz)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1,1)
         br_r = fval(1,2)  !dBR/dR
         br_z = fval(1,3)  !dBR/dZ
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1,1)
         bz_r = fval(1,2)  !dBR/dR
         bz_z = fval(1,3)  !dBR/dZ
      ELSE
         br_temp = 0.0
         bz_temp = 0.0
      END IF
      qdot(1) = br_temp
      qdot(2) = bz_temp
      qdot(3) = br_r
      qdot(4) = br_z
      qdot(5) = bz_r
      qdot(6) = bz_z
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fmap_nag
