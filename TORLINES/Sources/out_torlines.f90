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
      USE torlines_runtime, ONLY: dphi, pi2, npoinc, lverb, phi_end
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
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------
      rho_local = q(1)
      theta_local = q(2)
      phi_local = phi
      IF (rho_local >= 1) THEN
         rho_local = 1
         goodline(myline) = .false.
      END IF
      theta_local = MOD(theta_local,thmx)
      IF (theta_local < 0) theta_local = theta_local + pi2
      phi_local = MOD(phi,phmx)
      IF (phi_local < 0) phi_local = phi_local + pi2
      ! Spline R
      CALL EZspline_interp(R_spl,rho_local,theta_local,phi_local,r_temp,ier)
      ! Spline Z
      CALL EZspline_interp(Z_spl,rho_local,theta_local,phi_local,z_temp,ier)
      ! Spline B
      CALL EZspline_interp(B_spl,rho_local,theta_local,phi_local,b_temp,ier)

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
      
      phi = phi + dphi

      myldex = myldex + 1
      IF (myldex == nsteps) phi = phi_end(myline)
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------
      END SUBROUTINE out_torlines
