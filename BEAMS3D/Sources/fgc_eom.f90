!-----------------------------------------------------------------------
!     Function:      fgc_eom
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          03/19/2024
!     Description:   This subroutine calculates the RHS of the ODE for 
!                    gyrocenter orbit following.  The ODE's are similar
!                    to those used by the ASCOT5 and LOCUST codes. For
!                    the derivation refer to
!                    Littlejohn R.G. 1983 Variational principles of 
!                    guiding centre motion J. Plasma Phys. 29 111–25
!              https://mathworld.wolfram.com/CylindricalCoordinates.html
!-----------------------------------------------------------------------
      SUBROUTINE fgc_eom(t,q,qdot)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_grid
      USE beams3d_runtime, ONLY: lneut
      USE beams3d_lines, ONLY: moment, mycharge, mymass, myv_neut,mytdex
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
!          gradb      vector with gradient of |B| evaluatedfrom splines
!          gradbr     vector with gradient of Br evaluated from splines
!          gradbz     vector with gradient of Bz evaluated from splines
!          gradbphi   vector with gradient of Bphi evaluated from splines
!          Bstar     vector with B-related quantities from lagrange formulation
!          Estar     vector with E-related quantities from lagrange formulation      

!-----------------------------------------------------------------------
      INTEGER :: ier
      REAL(rprec) :: r_temp, phi_temp, z_temp, modb_temp, br_temp, bz_temp, bphi_temp,&
                      vll, rinv, binv, cinv, BhatDotBstar, A, B
      REAL(rprec) :: Bhat(3),Bstar(3),&
                    gradB(3),gradbr(3),gradbz(3),gradbphi(3),gradBcrossB(3),curlB(3),&
                    Efield(3),Estar(3),EstarcrossBhat(3),ExB(3)
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam
      REAL*8 :: fval(1,4), fvalE(1,3)
      INTEGER, parameter :: ict(8)  = (/1,1,1,1,0,0,0,0/)
      INTEGER, parameter :: ictE(8) = (/0,1,1,1,0,0,0,0/)
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
      IF ((r_temp >= rmin-eps1) .and. (r_temp <= rmax+eps1) .and. &
          (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
          (z_temp >= zmin-eps3) .and. (z_temp <= zmax+eps3)) THEN
         ! Get the gridpoint info
         i = MIN(MAX(COUNT(raxis < r_temp),1),nr-1)
         j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
         k = MIN(MAX(COUNT(zaxis < z_temp),1),nz-1)
         xparam = (r_temp - raxis(i)) * hri(i)
         yparam = (phi_temp - phiaxis(j)) * hpi(j)
         zparam = (z_temp - zaxis(k)) * hzi(k)
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
         ! Fix grad(B) & grad(E)
         gradB(2)    = gradB(2)*rinv
         Efield(2)   = Efield(2)*rinv
         ! Simplification for Conv. Operator below
         gradbr(2)   = gradbr(2)*rinv
         gradbphi(2) = gradbphi(2)*rinv
         gradbz(2)   = gradbz(2)*rinv
         ! Calculated some helpers
         binv    = one/modb_temp
         cinv    = one/mycharge
         ! Normalization \vec{B}/|B|
         Bhat(1) = br_temp*binv
         Bhat(2) = bphi_temp*binv
         Bhat(3) = bz_temp*binv
         
         ! Curl(B) in cylindical coordiantes
         curlB(1)= gradbz(2) - gradbphi(3)
         curlB(2)= gradbr(3) - gradbz(1)
         curlB(3)=bphi_temp*rinv-gradbr(2)+gradbphi(1)
         
         ! grad(B)xB
         gradBcrossB(1) = gradB(2) * bz_temp   - gradB(3) * bphi_temp
         gradBcrossB(2) = gradB(3) * br_temp   - gradB(1) * bz_temp
         gradBcrossB(3) = gradB(1) * bphi_temp - gradB(2) * br_temp

         ! B* = B + vll*m*(curl(B) - grad(B)xb)/(B*q)
         A = vll * mymass * cinv *binv
         Bstar(1) = br_temp   + A * (curlB(1) - gradBcrossB(1) * binv )
         Bstar(2) = bphi_temp + A * (curlB(2) - gradBcrossB(2) * binv )
         Bstar(3) = br_temp   + A * (curlB(3) - gradBcrossB(3) * binv )

         ! E* = -grad(Phi) - mu*grad(B)/q
         B = moment * cinv
         Estar = Efield - B * gradB

         ! Note here we define this as 1/b.B*
         BhatDotBstar = Bhat(1)*Bstar(1) + Bhat(2)*BStar(2) + Bhat(3)*BStar(3)
         BhatDotBstar = one / BhatDotBstar

         ! dX/dt = (vll * B* + cross(E*,b))/(b.B*)
         ! dvll/dt = (q/m * B*.E*)/(b.B*)
         qdot(1) = vll * Bstar(1) + Estar(2) * Bhat(3) - Estar(3)*Bhat(2)
         qdot(2) = vll * Bstar(2) + Estar(3) * Bhat(1) - Estar(1)*Bhat(3)
         qdot(3) = vll * Bstar(3) + Estar(1) * Bhat(2) - Estar(2)*Bhat(1)
         qdot(4) = mycharge * ( Bstar(1) * Estar(1) + Bstar(2) * Estar(2) + Bstar(3) * Estar(3) ) / mymass

         qdot = qdot * BhatDotBstar

        ! Because dphi/dt = vphi/R !rad/s
         qdot(2) = qdot(2)*rinv
      ELSE
         qdot(1:4) = 0
      END IF
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fgc_eom

