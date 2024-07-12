      SUBROUTINE obtain_theta(alpha, iotac, thetacn, nsurf, npt,
     1   zetang, thetang, lambdath)
      USE stel_kinds
      USE normalize_data, ONLY: lasym_v                                     ! 110909 RS: logical for Asymmetric input (if TRUE)
      USE ballooning_data
      USE general_dimensions
      USE fmesh_quantities
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN):: nsurf, npt
      REAL(rprec), INTENT(IN) :: iotac, alpha, thetacn
      REAL(rprec), INTENT(IN), DIMENSION(npt):: zetang
      REAL(rprec), INTENT(OUT), DIMENSION(npt):: thetang,
     1   lambdath
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: guard_band = 0.95_dp
      INTEGER :: j, lj, l, kj
      REAL(rprec):: lambda0, dlambda0, ccosi, ssine, arg,
     1   aux1, theta0, dtheta, fun0, dfun0
!-----------------------------------------------

      theta_loop: DO l = 1, npt
         aux1 = iotac*zetang(l) + alpha
         IF(l == 1)  THEN
           theta0 = thetacn
         ELSE
          theta0 = thetang(l-1)
         ENDIF

!
!        Use Newton-Raphson to find angle theta that is the zero of
!        F(theta) = theta + lambda(theta) - (alpha + iota*zeta)
!
         newton_raphson: DO kj = 1, jnewton_max
           lambda0 = 0
           dlambda0 = 0

           fourier: DO j = 1, mnmax_v

             arg = xm_v(j)*theta0-xn_v(j)*zetang(l)
             ssine = SIN(arg)
             ccosi = COS(arg)
             lj = mnmax_v*(nsurf-1)+j
             lambda0 = lambda0 + lmnsf(lj)*ssine                            ! Fourier invert lambda and d(lambda)/dtheta
             IF (lasym_v) lambda0 = lambda0 + lmncf(lj)*ccosi               ! 110909 RS: Asymmetric input
             dlambda0 = dlambda0 + xm_v(j)*lmnsf(lj)*ccosi
             IF (lasym_v) dlambda0 = dlambda0 - xm_v(j)*lmncf(lj)*ssine     ! 110909 RS: Asymmetric input

           ENDDO fourier

           fun0 = theta0 - aux1 + lambda0
           dfun0 = 1 + dlambda0

           IF (kj .ge. jnewton_max/2) dfun0 = MAX(dfun0, guard_band)        ! Added by SPH (10/00)

           dtheta = fun0/dfun0
           theta0 = theta0 - dtheta

           IF (ABS(dtheta) < newton_tol) THEN
             thetang(l) = theta0
             lambdath(l) = dlambda0
             EXIT
           ELSE IF (kj == jnewton_max) THEN
             !STOP 'COBRA: NEWTON fails!'
             lfail_balloon = .true.
             RETURN
           ENDIF

         ENDDO newton_raphson

      ENDDO theta_loop

      END SUBROUTINE obtain_theta
