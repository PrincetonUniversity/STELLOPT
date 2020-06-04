!-----------------------------------------------------------------------
!     Function:      jacobian_lsode
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/16/2012
!     Description:   This subroutine calculates the Jacobian of the ODE for 
!                    field line following (LSODE).  Currently just a
!                    dummy routine.
!-----------------------------------------------------------------------
      SUBROUTINE jacobian_lsode(neq,phi,q,ml,mp,pd,nrpd)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      !USE EZspline_obj
      !USE EZspline                                            ! MPI
!-----------------------------------------------------------------------
!     Input Variables
!          phi        phi angle
!          q          (q(1),q(2)) = (R,Z)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: neq, nrpd, ml, mp
      DOUBLE PRECISION :: phi, q(2), pd(2,2)
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, br_temp, bz_temp, br_dot, bz_dot
      REAL(rprec) :: x,y,z,vx, vy, vz
      !DOUBLE PRECISION :: R_grad(3), Z_grad(3)
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1,2) ! So weird behavior but this must match the sum of ict
      INTEGER, PARAMETER :: ict(8)=(/0,1,0,1,0,0,0,0/)
      !INTEGER, parameter :: ict(10)=(/0,1,0,1,0,0,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier    = 0
      r_temp = q(1)
      z_temp = q(2)
      phi_temp = MOD(phi,delta_phi)
      IF (phi_temp < 0) phi_temp = delta_phi + phi_temp
      pd = 0
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
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                BR4D(1,1,1,1),nr,nphi,nz)
         pd(1,1) = fval(1,1) !dBR/dR
         pd(1,2) = fval(1,2) ! dBR/dZ
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                BZ4D(1,1,1,1),nr,nphi,nz)
         pd(2,1) = fval(1,1) ! dBZ/dR
         pd(2,2) = fval(1,2) ! dBZ/dZ
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jacobian_lsode
