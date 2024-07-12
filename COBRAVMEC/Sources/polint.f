      SUBROUTINE polINT(xa, ya, n, x, y, dy)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,INTENT(IN) :: n
      REAL(rprec),INTENT(IN), DIMENSION(n) :: xa, ya
      REAL(rprec),INTENT(IN) :: x
      REAL(rprec),INTENT(OUT) :: y, dy
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER,PARAMETER :: nmax=10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nmax) :: c, d
      REAL(rprec) ::  dif, dift, h0, hp, w, den
      INTEGER :: ns, i, m
!-----------------------------------------------
!      Given arrays XA and YA of length N, and given X,
!      returns Y, and an error estimate DY.
!      If P(x) is the polynomial of degree N-1 such that
!      P(AXi)=AYi, i=1,...,n, then the returned y=p(x).
!      Interpolation is done using Neville's algorithm.
!      (NUMERICAL RECIPES, p. 82).
!
       ns=1
       dif=ABS(x-xa(1))
       DO i=1,n
       dift=ABS(x-xa(i))
         IF(dift.lt.dif)then
           ns=i
           dif=dift
         ENDIF
         c(i)=ya(i)
         d(i)=ya(i)
       ENDDO
       y=ya(ns)
       ns=ns-1
       DO m=1,n-1
         DO i=1,n-m
           h0=xa(i)-x
           hp=xa(i+m)-x
           w=c(i+1)-d(i)
           den=h0-hp
           IF (den .eq. 0._dp)then
             WRITE(*,*)'Two identical xa in input'
             RETURN
           ENDIF
           den=w/den
           d(i)=hp*den
           c(i)=h0*den
         ENDDO
         IF(2*ns.lt.n-m)then
           dy=c(ns+1)
         ELSE
           dy=d(ns)
           ns=ns-1
         ENDIF
         y=y+dy
       ENDDO

      END SUBROUTINE polint
