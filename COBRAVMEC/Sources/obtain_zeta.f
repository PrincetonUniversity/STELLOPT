
       subroutine obtain_zeta(alpha, iotac, zetacn, nsurf, npt,
     1   zetang, thetang, lambdazh)
       use stel_kinds
       use ballooning_data
       use general_dimensions
       use fmesh_quantities
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in):: nsurf, npt
      REAL(rprec), INTENT(in) :: iotac, alpha, zetacn
      REAL(rprec), INTENT(in), dimension(npt):: thetang
      REAL(rprec), INTENT(out), dimension(npt):: zetang, lambdazh
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0._dp, one = 1.0_dp
      REAL(rprec), PARAMETER :: guard_band = 0.95_dp, two = 2.0_dp
      INTEGER :: l, kj
      REAL(rprec):: dlamb0,  aux1, fun0, dfun0,
     1   z0, zl, zh, dzold, dz, fl, fh, two pi
!-----------------------------------------------
!
!  Have switched this from just Newton-Raphson to NR with bisect - ASW (6/05)
!
!  Start the loop over the numbering of starting points
!
       twopi   = 8*atan(one)
       zeta_loop: DO l = 1, npt
         aux1  = thetang(l) - alpha
         zl    =  -twopi+0.01d0                                               ! At the moment, assuming the solution lies within 0 < zeta < 2*pi
         zh    = twopi

         CALL get_lambda(nsurf,zl,thetang(l),aux1,iotac,fl,
     1     dfun0,dlamb0)
         CALL get_lambda(nsurf,zh,thetang(l),aux1,iotac,fh,
     1     dfun0,dlamb0)

         IF (fl*fh.ge.zero) STOP ' BISECT failed! '                           ! Making sure f(0) < 0 or f(2*pi) < 0 (but not both)
         IF (fl.gt.zero) THEN
            z0 = zl
            zl = zh
            zh = z0
            z0 = fl
            fl = fh
            fh = z0
         END IF
         z0 = (zh + zl)/two
         dzold = ABS(zh - zl)
         dz = dzold
         CALL get_lambda(nsurf,z0,thetang(l),aux1,iotac,fun0,                 ! Evaluate the function at the midpoint
     1     dfun0,dlamb0)

         newton_raphson: DO kj = 1, jnewton_max                               ! Use Newton-Raphson with bisect to find zero of
                                                                              ! F(zeta) = theta + lambda(zeta) - alpha - iota*zeta
            IF ( ((z0-zh)*dfun0-fun0)*((z0-zl)*dfun0-fun0).ge. zero           ! Bisect if Newton out of range
     1          .or. abs(two*fun0) .gt. abs(dzold*dfun0) ) THEN               ! or not decreasing fast enough
               dzold = dz
               dz = (zh - zl)/two
               z0 = zl + dz
            ELSE                   ! Take the Newton step
               dzold = dz
               dz = fun0/dfun0
               z0 = z0 - dz
            END IF
            IF (abs(dz) < newton_tol) THEN                                    ! We have found the root!
               zetang(l) = z0
               lambdazh(l) = dlamb0
               EXIT
            ELSE IF (kj == jnewton_max) THEN                                  ! Too many attempts
               WRITE(*,*) '|dz| = ',abs(dz),', zeta=',z0
               STOP 'COBRA: NEWTON fails!'
            END IF
            CALL get_lambda(nsurf,z0,thetang(l),aux1,iotac,fun0,
     1        dfun0,dlamb0)
!
!  Narrow the bracket on the root
!
            IF (fun0.lt.zero) THEN
               zl = z0
               fl = fun0
            ELSE
               zh = z0
               fh = fun0
            END IF
         END DO newton_raphson
       END DO zeta_loop

!      END SUBROUTINE obtain_zeta
       END
