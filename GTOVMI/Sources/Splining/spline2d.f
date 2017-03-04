      SUBROUTINE splnyr(x,y,n,yp,wk)
      USE precision
      IMPLICIT NONE
      INTEGER :: n
      REAL(rprec) :: x(n),y(n),yp(n),wk(n,5)
c...........................................................................
      REAL(rprec) :: t,v,btmp
      REAL(rprec) :: a2,b2,c2,an2,bn2,cn2
      INTEGER :: j
c
      t=(x(2)-x(1))/(x(3)-x(2))
      v=(x(n)-x(n-2))/(x(n-1)-x(n-2))
c
      DO j=1,n-1
         yp(j)=1./(x(j+1)-x(j))
      ENDdo
c
c     set up the maxtrix elements for y'
      wk(2,2)=(2.+t+1./t)*yp(2)
      wk(2,3)=(1.+t)*yp(2)
      wk(2,4)=(y(2)-y(1))*yp(1)**2+(y(3)-y(2))*(2.+3.*t)*yp(2)**2
      wk(n-1,1)=(v+1./v-2.)*yp(n-2)
      wk(n-1,2)=(v-1.)*yp(n-2)
      wk(n-1,4)=(1./v**2+2.*v-3.)*yp(n-2)**2*(y(n-1)-y(n-2))+
     $     (y(n)-y(n-1))*yp(n-2)**2/v**2
      IF (n .gt. 4)then
         DO j=3,n-2
            wk(j,1)=yp(j-1)
            wk(j,2)=2.*(yp(j-1)+yp(j))
            wk(j,3)=yp(j)
            wk(j,4)=3.*((y(j)-y(j-1))*yp(j-1)**2+(y(j+1)-y(j))*yp(j)**2)
         ENDdo
      ENDIF
c
      btmp=wk(2,2)
      yp(2)=wk(2,4)/btmp
      DO j=3,n-1
         wk(j,5)=wk(j-1,3)/btmp
         btmp=wk(j,2)-wk(j,1)*wk(j,5)
         yp(j)=(wk(j,4)-wk(j,1)*yp(j-1))/btmp
      ENDdo
      DO j=n-2,2,-1
         yp(j)=yp(j)-yp(j+1)*wk(j+1,5)
      ENDdo  
c
      a2=yp(2)*(x(3)-x(2))
      b2=3.*(y(3)-y(2))-(x(3)-x(2))*(2.*yp(2)+yp(3))
      c2=-2.*(y(3)-y(2))+(x(3)-x(2))*(yp(2)+yp(3))
      yp(1)=(a2+t*(-2.*b2+t*3.*c2))/(x(3)-x(2))
c
      an2=yp(n-2)*(x(n-1)-x(n-2))
      bn2=3.*(y(n-1)-y(n-2))-(x(n-1)-x(n-2))*(2.*yp(n-2)+yp(n-1))
      cn2=-2.*(y(n-1)-y(n-2))+(x(n-1)-x(n-2))*(yp(n-2)+yp(n-1))
      yp(n)=(an2+v*(2.*bn2+v*3.*cn2))/(x(n-1)-x(n-2))
c
      END SUBROUTINE splnyr
      SUBROUTINE bcspline(x,y,f,nx,ny,nxp,c)
      USE precision
      IMPLICIT NONE
      INTEGER :: nx,ny,nxp
      REAL(rprec) :: f(nxp,ny),x(nx),y(ny),c(4,nxp,ny)
c..........................................................................
      INTEGER :: kx
      PARAMETER (kx=515)
      REAL(rprec) :: g(kx),gp(kx),wk(5*kx)
      INTEGER :: i,j
c
      IF (nx .lt. 4 .or. nx .gt. kx)then
         WRITE (6,'(a)')'bcspline:   nx < 4 .or. nx > kx = 515'
         STOP
      ENDIF
      IF (ny .lt. 4 .or. ny .gt. kx)then
         WRITE (6,'(a)')'bcspline:   ny < 4 .or. ny > kx = 515'
         STOP
      ENDIF
      IF (nx .gt. nxp)then
         WRITE (6,'(a)')'bcspline:   nx > nxp (storage error)'
         STOP
      ENDIF
c
      IF (x(nx)-x(1) .lt. 0.)then
         WRITE (6,'(a)')'bcspline:  x not monotonically increasing'
         STOP
      ENDIF
      IF (y(ny)-y(1) .lt. 0.)then
         WRITE (6,'(a)')'bcspline:  y not monotonically increasing'
         STOP
      ENDIF
c
      DO j=1,ny
         DO i=1,nx
            c(1,i,j)=f(i,j)
         ENDdo
      ENDdo
c
      DO j=1,ny
         DO i=1,nx
            g(i)=f(i,j)
         ENDdo
         CALL splnyr(x,g,nx,gp,wk)
         DO i=1,nx
            c(2,i,j)=gp(i)
         ENDdo
      ENDdo
c
      DO i=1,nx
         DO j=1,ny
            g(j)=f(i,j)
         ENDdo
         CALL splnyr(y,g,ny,gp,wk)
         DO j=1,ny
            c(3,i,j)=gp(j)
         ENDdo
      ENDdo
c
      DO i=1,nx
         DO j=1,ny
            g(j)=c(2,i,j)
         ENDdo
         CALL splnyr(y,g,ny,gp,wk)
         DO j=1,ny
            c(4,i,j)=gp(j)
         ENDdo
      ENDdo
c
      END SUBROUTINE bcspline
      SUBROUTINE bcsplINT(x,y,c,nx,ny,nxp,xl,yl,pds,icalc,ier)
      USE precision
      IMPLICIT NONE
      INTEGER :: nx,ny,nxp,icalc,ier
      REAL(rprec) :: x(nx),y(ny),c(4,nxp,ny),xl,yl,pds(6)
c...........................................................................
      INTEGER :: i,j
      REAL(rprec) :: hx, hy, sux(2), suy(2), 
     $     su(2), svx(2), sv(2), sxy(2),
     $     u, v, spln0, spln1, spln2, s0, sh, sp0, sph, h, d
c
      spln0(s0,sh,sp0,sph,h,d) = s0+d*(h*sp0+d*(3.*(sh-s0)-
     $  (sph+2.*sp0)*h+d*(2.*(s0-sh)+(sph+sp0)*h)))
      spln1(s0,sh,sp0,sph,h,d) = sp0+d*(6.*(sh-s0)/h-2.0*
     $  (sph+2.*sp0)+3.*d*(2.*(s0-sh)/h+(sph+sp0)))
      spln2(s0,sh,sp0,sph,h,d) = 6.*(sh-s0)/h**2-2.0*
     $  (sph+2.*sp0)/h+d*(2.*(s0-sh)/h**2+(sph+sp0)/h)*6.0
c
      SAVE i,j
      data i,j/0,0/
c
      ier = 0
c --- correlated table search for xl
      CALL tableintrp (x, nx, xl, i)
      IF ( i .eq. 0 )then
         IF (ABS(xl-x(1)) .le. ABS(x(2)-x(1)))then
            i=1
            ier = 3
            WRITE (6,'(a,i2)')'bcsplint: x < x(1), ier = ',ier
            go to 10
         ELSE
            ier=33
            WRITE (6,'(a,i2)')'bscplint: x < x(1), ier = ',ier
            RETURN
         ENDIF
      ENDIF
      IF ( i .eq. nx)then
         IF (ABS(xl-x(nx)) .le. ABS(x(nx-1)-x(nx)))then
            i=nx-1
            ier=5
            WRITE (6,'(a,i2)')'bcsplint: x > x(nx), ier = ',ier
         ELSE
            ier=35
            WRITE (6,'(a,i2)')'bcsplint: x > x(nx), ier = ',ier
            RETURN
         ENDIF
      ENDIF
 10   CONTINUE
c --- correlated table search for yl
      CALL tableintrp (y, ny, yl, j)
      IF ( j .eq. 0 )then
         IF (ABS(yl-y(1)) .le. ABS(y(2)-y(1)))then
            j=1
            ier=4
            WRITE (6,'(a,i2)')'bcsplint: y < y(1), ier = ',ier
            go to 20
         ELSE
            ier=34
            WRITE (6,'(a,i2)')'bcsplint: y < y(1), ier = ',ier
            RETURN
         ENDIF
      ENDIF
      IF ( j .eq. ny)then
         IF (ABS(yl-y(ny)) .le. ABS(y(ny-1)-y(ny)))then
            j=ny-1
            ier=6
            WRITE (6,'(a,i2)')'bcsplint: y > y(ny), ier = ',ier
         ELSE
            ier=36
            WRITE (6,'(a,i2)')'bcsplint: y > y(ny), ier = ',ier
            RETURN
         ENDIF
      ENDIF
 20   CONTINUE
c
      hx   = x(i+1)-x(i)
      hy   = y(j+1)-y(j)
      u    = (xl-x(i))/hx
      v    = (yl-y(j))/hy
c
      sv(1)= spln0(c(1,i,j),c(1,i,j+1),c(3,i,j),c(3,i,j+1),hy,v)
      svx(1)=spln0(c(2,i,j),c(2,i,j+1),c(4,i,j),c(4,i,j+1),hy,v)
      sv(2)= spln0(c(1,i+1,j),c(1,i+1,j+1),c(3,i+1,j),c(3,i+1,j+1),hy,v)
      svx(2)=spln0(c(2,i+1,j),c(2,i+1,j+1),c(4,i+1,j),c(4,i+1,j+1),hy,v)      
      pds(1)=spln0(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF (icalc .le. 1)return
c
      su(1)= spln0(c(1,i,j),c(1,i+1,j),c(2,i,j),c(2,i+1,j),hx,u)
      suy(1)=spln0(c(3,i,j),c(3,i+1,j),c(4,i,j),c(4,i+1,j),hx,u)
      su(2)= spln0(c(1,i,j+1),c(1,i+1,j+1),c(2,i,j+1),c(2,i+1,j+1),hx,u)
      suy(2)=spln0(c(3,i,j+1),c(3,i+1,j+1),c(4,i,j+1),c(4,i+1,j+1),hx,u)
      pds(2)=spln1(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF (icalc .eq. 2)return
      pds(3)=spln1(su(1),su(2),suy(1),suy(2),hy,v)
      IF (icalc .eq. 3)return
c
      sux(1)=spln1(c(1,i,j),c(1,i+1,j),c(2,i,j),c(2,i+1,j),hx,u)
      sxy(1)=spln1(c(3,i,j),c(3,i+1,j),c(4,i,j),c(4,i+1,j),hx,u)
      sux(2)=spln1(c(1,i,j+1),c(1,i+1,j+1),c(2,i,j+1),c(2,i+1,j+1),hx,u)
      sxy(2)=spln1(c(3,i,j+1),c(3,i+1,j+1),c(4,i,j+1),c(4,i+1,j+1),hx,u)
      pds(4) = spln1(sux(1),sux(2),sxy(1),sxy(2),hy,v)
      IF (icalc .eq. 4)  RETURN
      pds(5) = spln2(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF (icalc .eq. 5)  RETURN
      pds(6) = spln2(su(1),su(2),suy(1),suy(2),hy,v)
c
      END SUBROUTINE bcsplint


      SUBROUTINE tableintrp (xx, n, x, jlo)
c
c ----------------------------------------------------------------------
c --- correlated table search routine. USE jlo from previous CALL to
c --- get jlow for current value of x. IF jlo from previous CALL is
c --- no good(jlo=0 or jlo=n) THEN DO a binary search.
c --- this routine RETURNs jlo so that
c   xx(jlo) .le. x  .le. xx(jlo+1) IF jlo=1
c   xx(jlo) .lt. x  .le. xx(jlo+1) IF jlo =2,3...n-1
c   it is assumed that xx(j),j=1,2..n is monotonically
c   increasing,this is NOT checked for.
c  this is a MODified version of SUBROUTINE HUNT (Numerical Recipes)
c ------------------------------------------------------------------ HSJ
c
      USE precision
      IMPLICIT NONE
c
      INTEGER :: n,jlo,jhi,jmid,inc
      REAL(rprec) ::    x,xx(n)
c
      IF (x .lt. xx(1)) THEN
        jlo=0                       ! indicates out of range below xx(1)
      ELSE IF (x .le. xx(2)) THEN
        jlo=1                       ! xx(1) .le. x .le. xx(2)
      ELSE IF (x .le. xx(n)) THEN   ! xx(2) .lt. x .le. xx(n)
c
c         check IF jlo from previous CALL is usable
c
          IF (jlo .eq. 0 .or. jlo .eq. n) THEN
              jlo=2
              jhi=n
              go to 15  ! no correlation go directly to bisection
          ELSE          ! 1 .le. jlo .lt. n
c
c             bracket x, THEN USE bisection
c             start with jlo from previous CALL
c
              inc=1
              IF (x .gt. xx(jlo)) THEN ! search up
    4             jhi=jlo+inc
                  IF (jhi .ge. n) THEN
                      jhi=n
                  ELSE IF (x .gt. xx(jhi)) THEN
                      inc=inc+inc
                      jlo=jhi
                      go to 4
                  END IF
              ELSE
    5             jhi=jlo
                  jlo=jlo-inc
                  IF (jlo .le. 1) THEN
                      jlo=1
                  ELSE IF (x .le. xx(jlo)) THEN
                      inc=inc+inc
                      go to 5
                  END IF
              END IF
c
          END IF
c
c         bisection
c
   10     IF (jhi-jlo .eq. 1)  RETURN
   15     jmid=(jhi+jlo)/2
          IF (x .gt. xx(jmid)) THEN
              jlo=jmid
          ELSE  ! x .le. xx(jmid)
              jhi=jmid
          END IF
          go to 10
c
      ELSE  ! x .gt. xx(n)
          jlo=n
      END IF
c
      END SUBROUTINE tableintrp


