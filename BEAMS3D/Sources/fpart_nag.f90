!-----------------------------------------------------------------------
!     Function:      fpart_nag
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          10/19/2021
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    particle orbit following.  
!
!-----------------------------------------------------------------------
      SUBROUTINE fpart_nag(t,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_grid
      USE beams3d_lines, ONLY: mycharge, mymass
      USE mpi_params, ONLY: myworkid
!-----------------------------------------------------------------------
!     Input Variables
!          t          time coordinates
!          q          (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: t, q(6), qdot(6)
!-----------------------------------------------------------------------
!     Local Variables
!          ier        Error flag
!          r_temp     R
!          z_temp     Z
!          phi_temp   phi (radians)
!          br_temp    br evaluated from splines
!          bz_temp    bz evaluated from splines
!          bphi_temp  bphi evaluated from splines
!          modb_temp  |B| evaluated from splines
!!         br_dot     Random force for diffusion
!          moment     mu, or moment of the particle; .5mVperp^2/|B|
!          A, B       constants to speed up the calculations; see equations below
!          gradb      vector with gradient of |B| evaluatedfrom splines
!          gradbr     vector with gradient of Br evaluated from splines
!          gradbz     vector with gradient of Bz evaluated from splines
!          gradbphi   vector with gradient of Bphi evaluated from splines
!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, br_temp, bz_temp, bphi_temp, rinv
      REAL(rprec) :: Efield(3),VxB(3)
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam
      REAL*8 :: fval(1,1), fvalE(1,3)
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ictE(8)=(/0,1,1,1,0,0,0,0/)
      REAL*8, PARAMETER :: one = 1
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

!     Do we want to convert the above vars to doubles before doing the calculations below? q and qdot are
!     doubles so it seems like that would be an option.
      ier      = 0
      r_temp   = q(1)
      phi_temp = MOD(q(2), phimax)
      IF (phi_temp < 0) phi_temp = phi_temp + phimax
      z_temp   = q(3)
      rinv = one/r_temp
!      CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
      IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
          (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
          (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
!      IF (ier == 0) THEN
         ! Get the gridpoint info
         i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
         j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
         k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
         xparam = (r_temp - raxis(i)) * hri(i)
         yparam = (phi_temp - phiaxis(j)) * hpi(j)
         zparam = (z_temp - zaxis(k)) * hzi(k)
         !CALL R8HERM3xyz(r_temp,phi_temp,z_temp,&
         !                BR_spl%x1(1),BR_spl%n1,&
         !                BR_spl%x2(1),BR_spl%n2,&
         !                BR_spl%x3(1),BR_spl%n3,&
         !                BR_spl%ilin1,BR_spl%ilin2,BR_spl%ilin3,&
         !                i,j,k,xparam,yparam,zparam,&
         !                hx,hxi,hy,hyi,hz,hzi,ier)
         ! Evaluate the Splines
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BR4D(1,1,1,1),nr,nphi,nz)
         br_temp = fval(1,1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
         bphi_temp = fval(1,1)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1,1)
         CALL R8HERM3FCN(ictE,1,1,fvalE,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         POT4D(1,1,1,1),nr,nphi,nz)
         Efield(1:3) =-fvalE(1,1:3)
         ! Fix gradients
         Efield(2)   = Efield(2)*rinv
         ! VxB
         VxB(1) = q(5)*bz_temp   - q(6)*bphi_temp
         VxB(2) = q(6)*br_temp   - q(4)*bz_temp
         VxB(3) = q(4)*bphi_temp - q(5)*br_temp
         ! Equations
         qdot(1:3) = q(4:6)
         qdot(4:6) = mycharge*(Efield(1:3)+VxB)/mymass

         ! dA_phi/dt = F_phi/R
         qdot(2) = qdot(2)*rinv
         qdot(5) = qdot(5)*rinv
      ELSE
         qdot(1:6) = 0
      END IF
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fpart_nag
