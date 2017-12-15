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
      USE beams3d_grid, ONLY: BR_spl, BZ_spl, delta_t, BPHI_spl, MODB_spl, &
                              POT_spl, phimax, rmin, rmax
      USE beams3d_runtime, ONLY: lmu, mu, lneut
      USE beams3d_lines, ONLY: moment, mycharge, mymass, myv_neut, B_temp
      USE EZspline_obj
      USE EZspline
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
                      vll, A, B, rinv, binv, cinv, pot_temp
      REAL(rprec) :: gradb(3),gradbr(3),gradbz(3),gradbphi(3),Efield(3)
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(8)=(/1,1,1,1,0,0,0,0/)
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
      !IF (.not.lneut) THEN
         CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
         IF (ier == 0) THEN
            ! This is the old way, however, since every grid is the same
            ! it usues may redundant calls.
            !CALL EZspline_interp(BR_spl,r_temp,phi_temp,z_temp,br_temp,ier)
            !CALL EZspline_interp(BZ_spl,r_temp,phi_temp,z_temp,bz_temp,ier)
            !CALL EZspline_interp(BPHI_spl,r_temp,phi_temp,z_temp,bphi_temp,ier)
            !CALL EZspline_interp(MODB_spl,r_temp,phi_temp,z_temp,modb_temp,ier)
            !CALL EZspline_gradient(BR_spl,r_temp,phi_temp,z_temp,gradbr,ier)
            !CALL EZspline_gradient(BPHI_spl,r_temp,phi_temp,z_temp,gradbphi,ier)
            !CALL EZspline_gradient(BZ_spl,r_temp,phi_temp,z_temp,gradbz,ier)
            !CALL EZspline_gradient(MODB_spl,r_temp,phi_temp,z_temp,gradb,ier)
            ! We could calculate modb_temp but because the SQRT is so slow
            ! the many multiplications in spline win out in the end.
            !modb_temp = SQRT(br_temp*br_temp+bphi_temp*bphi_temp+bz_temp*bz_temp)
            !binv     = 1/modb_temp
            !gradb(1) = (br_temp*gradbr(1)+bphi_temp*gradbphi(1)+bz_temp*gradbz(1))*binv
            !gradb(2) = (br_temp*gradbr(2)+bphi_temp*gradbphi(2)+bz_temp*gradbz(2))*binv
            !gradb(3) = (br_temp*gradbr(3)+bphi_temp*gradbphi(3)+bz_temp*gradbz(3))*binv
            ! Get the gridpoint info
            CALL R8HERM3xyz(r_temp,phi_temp,z_temp,&
                            BR_spl%x1(1),BR_spl%n1,&
                            BR_spl%x2(1),BR_spl%n2,&
                            BR_spl%x3(1),BR_spl%n3,&
                            BR_spl%ilin1,BR_spl%ilin2,BR_spl%ilin3,&
                            i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,ier)
            ! Evaluate the Splines
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            BR_spl%fspl(1,1,1,1),BR_spl%n1,BR_spl%n2,BR_spl%n3)
            br_temp = fval(1); gradbr(1:3) = fval(2:4)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            BPHI_spl%fspl(1,1,1,1),BPHI_spl%n1,BPHI_spl%n2,BPHI_spl%n3)
            bphi_temp = fval(1); gradbphi(1:3) = fval(2:4)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            BZ_spl%fspl(1,1,1,1),BZ_spl%n1,BZ_spl%n2,BZ_spl%n3)
            bz_temp = fval(1); gradbz(1:3) = fval(2:4)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            MODB_spl%fspl(1,1,1,1),MODB_spl%n1,MODB_spl%n2,MODB_spl%n3)
            modb_temp = fval(1); gradb(1:3) = fval(2:4)
            CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,&
                            POT_spl%fspl(1,1,1,1),POT_spl%n1,POT_spl%n2,POT_spl%n3)
            pot_temp = fval(1); Efield(1:3) = fval(2:4)
            ! Calculated some helpers
            binv    = one/modb_temp
            cinv    = one/mycharge
            A       = moment*binv*binv*cinv
            B       = mymass*vll*vll*binv*binv*binv*binv*cinv
            qdot(1) = A*( bphi_temp*gradb(3)-bz_temp*rinv*gradb(2) )  -  (bphi_temp*Efield(3)-bz_temp*rinv*Efield(2))*binv*binv + &
                      B*( bphi_temp*br_temp*gradbz(1) + bphi_temp*bphi_temp*rinv*gradbz(2) + bphi_temp*bz_temp*gradbz(3) &
                          - bz_temp*br_temp*gradbphi(1) - bz_temp*bphi_temp*rinv*gradbphi(2) - bz_temp*bz_temp*gradbphi(3) &
                          - bz_temp*bphi_temp*br_temp*rinv) + vll*br_temp*binv
             
            qdot(2) = rinv*A*( bz_temp*gradb(1)-br_temp*gradb(3) )    - rinv*(bz_temp*Efield(1)-br_temp*Efield(3))*binv*binv +  &
                      rinv*B*( bz_temp*br_temp*gradbr(1) + bz_temp*bphi_temp*rinv*gradbr(2) + bz_temp*bz_temp*gradbr(3) &
                           - bz_temp*bphi_temp*bphi_temp*rinv - br_temp*br_temp*gradbz(1) - br_temp*bphi_temp*rinv*gradbz(2) &
                           - br_temp*bz_temp*gradbz(3) ) + vll*bphi_temp*binv*rinv
             
            qdot(3) = A*( br_temp*rinv*gradb(2)-bphi_temp*gradb(1) )  - (br_temp*rinv*Efield(2)-bphi_temp*Efield(1))*binv*binv +  &
                      B*( br_temp*br_temp*gradbphi(1) + br_temp*bphi_temp*rinv*gradbphi(2) + br_temp*bz_temp*gradbphi(3) &
                          + br_temp*br_temp*bphi_temp*rinv - bphi_temp*br_temp*gradbr(1) - bphi_temp*bphi_temp*rinv*gradbr(2) &
                          - bphi_temp*bz_temp*gradbr(3) + bphi_temp*bphi_temp*bphi_temp*rinv ) + vll*bz_temp*binv
             
            qdot(4) = -moment*binv*( br_temp*gradb(1) + bphi_temp*rinv*gradb(2) + bz_temp*gradb(3) )/mymass &
                     +mycharge*binv*(br_temp*Efield(1) + bphi_temp*rinv*Efield(2) + bz_temp*Efield(3) )/mymass
         ELSE
            qdot(1:4) = 0
         END IF
      !ELSE
      !   qdot(1) = myv_neut(1)*cos(q(2)) + myv_neut(2)*sin(q(2))
      !   qdot(2) = rinv*( myv_neut(2)*cos(q(2)) - myv_neut(1)*sin(q(2)) )
      !   qdot(3) = myv_neut(3)
      !   qdot(4) = 0
      !   IF (r_temp <= rmin) qdot(1:4) = 0
      !   IF (r_temp > 5*rmax) qdot(1:4) = 0
      !END IF
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fpart_nag
