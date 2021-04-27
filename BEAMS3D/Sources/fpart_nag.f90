!-----------------------------------------------------------------------
!     Function:      fpart_nag
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    particle orbit following.  The ODE's are similar
!                    to those used by the DESORBS code (eqs. 16-20)
!                    R.H. Fowler, R.N. Morris, J.A. Rome,
!                    and K. Hantani, "Neutral Beam Injection Benchmark
!                    Studies for Stellarators/Heliotrons" Nuclear Fusion
!                    30, p.997 (1990)
!
!-----------------------------------------------------------------------
      SUBROUTINE fpart_nag(t,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_grid
      USE beams3d_runtime, ONLY: lneut
      USE beams3d_lines, ONLY: moment, mycharge, mymass, myv_neut
      USE mpi_params, ONLY: myworkid
!-----------------------------------------------------------------------
!     Input Variables
!          t          time coordinates
!          q          (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
!          qdot       dq/dt
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: t, q(4), qdot(4)
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
      REAL(rprec) :: r_temp, phi_temp, z_temp, modb_temp, br_temp, bz_temp, bphi_temp,&
                      vll, A, B, rinv, binv, cinv
      REAL(rprec) :: gradb(3),gradbr(3),gradbz(3),gradbphi(3),Efield(3),normb(3),bdgB(3),bxbdgB(3),ExB(3)
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam
      REAL*8 :: fval(1,4), fvalE(1,3)
      INTEGER, parameter :: ict(8)=(/1,1,1,1,0,0,0,0/)
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
      vll      = q(4)
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
         br_temp = fval(1,1); gradbr(1:3) = fval(1,2:4)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
         bphi_temp = fval(1,1); gradbphi(1:3) = fval(1,2:4)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         BZ4D(1,1,1,1),nr,nphi,nz)
         bz_temp = fval(1,1); gradbz(1:3) = fval(1,2:4)
         CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         MODB4D(1,1,1,1),nr,nphi,nz)
         modb_temp = fval(1,1); gradb(1:3) = fval(1,2:4)
         CALL R8HERM3FCN(ictE,1,1,fvalE,i,j,k,xparam,yparam,zparam,&
                         hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                         POT4D(1,1,1,1),nr,nphi,nz)
         Efield(1:3) =-fvalE(1,1:3)
         ! Fix gradients
         gradb(2)    = gradb(2)*rinv
         gradbr(2)   = gradbr(2)*rinv
         gradbphi(2) = gradbphi(2)*rinv
         gradbz(2)   = gradbz(2)*rinv
         Efield(2)   = Efield(2)*rinv
         ! Calculated some helpers
         binv    = one/modb_temp
         cinv    = one/mycharge
         ! Normalization
         normb(1) = br_temp*binv
         normb(2) = bphi_temp*binv
         normb(3) = bz_temp*binv
         ! (b.grad)B (https://mathworld.wolfram.com/ConvectiveOperator.html)
         bdgB(1) = normb(1)*gradbr(1)   + normb(2)*gradbr(2)   + normb(3)*gradbr(3)   - normb(2)*bphi_temp*rinv
         bdgB(2) = normb(1)*gradbphi(1) + normb(2)*gradbphi(2) + normb(3)*gradbphi(3) - normb(2)*br_temp*rinv
         bdgB(3) = normb(1)*gradbz(1)   + normb(2)*gradbz(2)   + normb(3)*gradbz(3)
         bxbdgB(1) = normb(2)*bdgB(3) - normb(3)*bdgB(2)
         bxbdgB(2) = normb(3)*bdgB(1) - normb(1)*bdgB(3)
         bxbdgB(3) = normb(1)*bdgB(2) - normb(2)*bdgB(1)
         ! ExB/(B*B)
         ExB(1) = Efield(2)*normb(3) - Efield(3)*normb(2)
         ExB(2) = Efield(3)*normb(1) - Efield(1)*normb(3)
         ExB(3) = Efield(1)*normb(2) - Efield(2)*normb(1)
         ExB    = ExB*binv
         ! Equations
         A       = moment*binv*cinv
         B       = mymass*vll*vll*binv*binv*cinv

         qdot(1) = normb(2)*gradb(3)-normb(3)*gradb(2)
         qdot(2) = normb(3)*gradb(1)-normb(1)*gradb(3)
         qdot(3) = normb(1)*gradb(2)-normb(2)*gradb(1)
         qdot(1:3) = A*qdot(1:3) + B*bxbdgB(1:3) + vll*normb(1:3) + ExB(1:3)
          
         qdot(4) = -moment*( normb(1)*gradb(1) + normb(2)*gradb(2) + normb(3)*gradb(3) ) &
                  +mycharge*(normb(1)*Efield(1) + normb(2)*Efield(2) + normb(3)*Efield(3) )

         qdot(2) = qdot(2)*rinv
         qdot(4) = qdot(4)/mymass
      ELSE
         qdot(1:4) = 0
      END IF
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fpart_nag
