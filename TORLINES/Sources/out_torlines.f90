!-----------------------------------------------------------------------
!     Function:      out_torlines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/2/2011
!     Description:   Save output from field line following while running.
!-----------------------------------------------------------------------
      SUBROUTINE out_torlines(phi,q)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_realspace
      USE torlines_fieldlines
      USE torlines_runtime
      USE torlines_background
!-----------------------------------------------------------------------
!     Input Parameters
!          phi          Location along fieldline in phi
!          q            (q(1),q(2)) = (R,Z)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(inout) :: phi
      DOUBLE PRECISION, INTENT(inout) :: q(2)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER             :: jint, ik, ier
      REAL(rprec) :: r_temp, z_temp, b_temp
      REAL(rprec) :: rho_local, theta_local, phi_local
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, PARAMETER :: ict(8)=(/1,0,0,0,0,0,0,0/)
      !INTEGER, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
!-----------------------------------------------------------------------
!     Begin Function
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
      r_temp = -1; z_temp=0; b_temp=0;
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
                         R4D(1,1,1,1),nrho,nu,nv)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                R4D(1,1,1,1),nrho,nu,nv)
         r_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         Z4D(1,1,1,1),nrho,nu,nv)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                Z4D(1,1,1,1),nrho,nu,nv)
         z_temp = fval(1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         B4D(1,1,1,1),nrho,nu,nv)
         !CALL R8FVTRICUB(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,&
         !                B4D(1,1,1,1),nrho,nu,nv)
         b_temp = fval(1)
      END IF
      R_lines(myline,myldex) = r_temp
      Z_lines(myline,myldex) = z_temp
      PHI_lines(myline,myldex) = phi
      B_lines(myline,myldex) = b_temp
      U_lines(myline,myldex) = theta_local
      IF (lverb .and. (MOD(myldex,1000)==0)) THEN
         CALL backspace_out(6,6)
         WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*((myline-1)*nsteps+myldex)/(nsteps*myend))),']%'
         CALL FLUSH(6)
      END IF

      ! Update position
      phi = phi + dphi
      myldex = myldex + 1
      IF (myldex == nsteps) phi = phi_end(myline)
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------
      END SUBROUTINE out_torlines
