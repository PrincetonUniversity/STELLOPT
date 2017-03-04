      SUBROUTINE funcv(n, xx, fvec)

      USE stel_kinds
      USE cmnf1
      USE geom
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      REAL(rprec), DIMENSION(n) :: xx, fvec
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: zv, cm, rv, u, v,
     1   ru,  zu, r, z, cn, dy, dz, dx, 
     2   coh, sih, y, xv, yu, x, si, co, yv, xu
!-----------------------------------------------
!
!
!
      u = xx(1)
      v = xx(2)
!
!     initializations
      r = 0; z = 0; ru = 0; zu = 0; rv = 0; zv = 0

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
      x = coh*r
      y = sih*r
      xu = coh*ru
      yu = sih*ru
      xv = coh*rv - alp*y
      yv = sih*rv + alp*x
!
      dx = x - xp
      dy = y - yp
      dz = z - zp
!
      fvec(1) = 2*(dx*xu + dy*yu + dz*zu)
      fvec(2) = 2*(dx*xv + dy*yv + dz*zv)
!
      END SUBROUTINE funcv
