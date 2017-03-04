      SUBROUTINE surnorm(isurf, u, v, xnorm, ynorm, znorm)
      USE stel_kinds
      USE cmnf1
      USE cmnf2
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER isurf
      REAL(rprec) u, v, xnorm, ynorm, znorm
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: si, coh, sih, co, cm, cn, xv, yv, crmag, yu, x, y,
     1     xu, pi2, zu, rv, zv, ru, alp, r, z
!-----------------------------------------------
!
!     compute the normal to the surface at (u,v)
!
!
      pi2 = 4*ASIN(1._dp)
      IF (isurf == 1) THEN
         alp = pi2/q1
      ELSE IF (isurf == 2) THEN
         alp = pi2/q1
      ELSE
         STOP 'ISURF must be surface 1 or 2'
      ENDIF
!
      r = 0;  z = 0;  ru = 0;  zu = 0;  rv = 0;  zv = 0
!
      IF (isurf == 1) THEN
         DO i = 1, nmn1
            cm = xm1(i)*pi2
            cn = xn1(i)*pi2
            co = COS(cm*u + cn*v)
            si = SIN(cm*u + cn*v)
            r = r + rmn1(i)*co
            z = z + zmn1(i)*si
            ru = ru - cm*rmn1(i)*si
            rv = rv - cn*rmn1(i)*si
            zu = zu + cm*zmn1(i)*co
            zv = zv + cn*zmn1(i)*co
         END DO
!
         coh = COS(alp*v)
         sih = SIN(alp*v)
!
         x = coh*r
         y = sih*r
         xu = coh*ru
         yu = sih*ru
         xv = coh*rv - alp*y
         yv = sih*rv + alp*x
!
      ELSE IF (isurf == 2) THEN
!
         DO i = 1, nmn2
            cm = xm2(i)*pi2
            cn = xn2(i)*pi2
            co = COS(cm*u + cn*v)
            si = SIN(cm*u + cn*v)
            r = r + rmn2(i)*co
            z = z + zmn2(i)*si
            ru = ru - cm*rmn2(i)*si
            rv = rv - cn*rmn2(i)*si
            zu = zu + cm*zmn2(i)*co
            zv = zv + cn*zmn2(i)*co
         END DO
!
         coh = COS(alp*v)
         sih = SIN(alp*v)
!
         x = coh*r
         y = sih*r
         xu = coh*ru
         yu = sih*ru
         xv = coh*rv - alp*y
         yv = sih*rv + alp*x
!
      ELSE
!
         STOP 'ISURF must be surface 1 or 2'
!
      ENDIF
!
!     compute the components of (Xu x Xv)
      xnorm = yu*zv - yv*zu
      ynorm = -(xu*zv - xv*zu)
      znorm = xu*yv - xv*yu
!     ! magnitude of (Xu x Xv)
      crmag = SQRT(xnorm*xnorm + ynorm*ynorm + znorm*znorm)
!     ! normalize the components and multiply by -1
!     ! to get right sign of the normal
      xnorm = -xnorm/crmag
      ynorm = -ynorm/crmag
      znorm = -znorm/crmag

      END SUBROUTINE surnorm
