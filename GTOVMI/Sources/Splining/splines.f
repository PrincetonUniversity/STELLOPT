      SUBROUTINE spline2d(fun,x,y,nx,ny,kx,coef)
c-----------------------------------------------------------------------
c     setup routine for bicubic spline of fun[x,y]
c     output of this routine is nx*ny spline coefficients 
c     stored in coef(kx,ny)
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: fun(*),x(*),y(*),coef(*)
      INTEGER :: kx,nx,ny
      INTEGER :: ind,j
      DO j=1,ny
         ind=(j-1)*kx+1
         CALL spline(x,fun(ind),nx,-1.e30_dbl,-1.e30_dbl,coef(ind))
      END DO
      RETURN
      END
c=======================================================================
      SUBROUTINE spline2dt(fun_new,x_new,y_new,kx_new,nx_new,ny_new,
     $     fun_old,x_old,y_old,kx_old,nx_old,ny_old,coef)
c-----------------------------------------------------------------------
c     evaluates bicubic spline: determines fun_new[x_new,y_new]
c     after spline2d has been CALLed to determine coef
c     in the CALLing PROGRAM fun_new,run_old, and coef are 2d arrays:
c     fun_old(kx_old,ky_old),fun_new(kx_new,ky_new),coef(kx_old,ky_old)
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: fun_new(*),x_new(*),y_new(*)
      REAL(rprec) :: fun_old(*),x_old(*),y_old(*)
      REAL(rprec) :: coef(*)
      INTEGER :: kx_new,nx_new,ny_new
      INTEGER :: kx_old,nx_old,ny_old
      INTEGER :: i,j,ind,nlow
      INTEGER, PARAMETER :: ky_old=2000
      REAL(rprec) :: ftemp(ky_old),ctemp(ky_old)
      IF(ny_old.gt.ky_old) THEN
         WRITE(6,'("dimensioning problem in spline2dt")')
         STOP
      END IF
      DO i=1,nx_new
         DO j=1,ny_old
            ind=(j-1)*kx_old+1
            nlow=0
            CALL splINT(x_old,fun_old(ind),coef(ind),nx_old,x_new(i),
     $           ftemp(j),nlow)
         END DO
         CALL  spline(y_old,ftemp,ny_old,-1.e30_dbl,-1.e30_dbl,ctemp)
         DO j=1,ny_new
            nlow=0
            CALL splINT(y_old,ftemp,ctemp,ny_old,y_new(j),
     $           fun_new((j-1)*kx_new+i),nlow)
         END DO
      END DO
      RETURN
      END
c=======================================================================
      SUBROUTINE spline1d(ynew,xnew,nnew,yold,xold,nold,y2old)
c-----------------------------------------------------------------------
c     USE 1d cubic spline on yold[xold] to produce ynew[xnew]
c     y2old(1:nold) is a work array
c     ynew(1:nnew) is the output
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: ynew(*),yold(*),xnew(*),xold(*)
      REAL(rprec) :: y2old(*)
      REAL(rprec) :: yp1,ypn
      INTEGER :: nnew,nold,i,nlow
      yp1=-1.e30_dbl
      ypn=-1.e30_dbl
      CALL spline(xold,yold,nold,yp1,ypn,y2old)
      DO i=1,nnew
         nlow=0
         CALL splINT(xold,yold,y2old,nold,xnew(i),ynew(i),nlow)
      END DO
      RETURN
      END
c=======================================================================
      SUBROUTINE spline1dp(ynew,xnew,nnew,yold,xold,nold,yp1,ypn,y2old)
c-----------------------------------------------------------------------
c     USE 1d cubic spline on yold[xold] to produce ynew[xnew]
c     y2old(1:nold) is a work array
c     ynew(1:nnew) is the output
c     This routine ALLows different boundary conditions
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: ynew(*),yold(*),xnew(*),xold(*)
      REAL(rprec) :: y2old(*)
      REAL(rprec) :: yp1,ypn
      INTEGER :: nnew,nold,i,nlow
c      yp1=-1.e30_dbl
c      ypn=-1.e30_dbl
      CALL spline(xold,yold,nold,yp1,ypn,y2old)
      DO i=1,nnew
         nlow=0
         CALL splINT(xold,yold,y2old,nold,xnew(i),ynew(i),nlow)
      END DO
      RETURN
      END
c=======================================================================
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
c-----------------------------------------------------------------------
c     spline routine based upon numerical recipes
c     this is the setup routine which needs to be CALLed ONLY once
c     splines y as a FUNCTION of x--both arrays have n elements
c     yp1 and ypn are boundary conditions on the spline
c     yp1=y'[x] at x=x[1]
c     IF yp1>=1.e30 THEN y''[x]=0 at x=x[1] is USEd
c     IF yp1<=-1.e30_dbl THEN y'[x[1]] is calculated from first four points
c     ypn=y'[x] at x=x[n]
c     IF ypn>=1.e30 THEN y''[x]=0 at x=x[n] is USEd
c     IF ypn<=-1.e30_dbl THEN y'[x[n]] is calculated from last four points
c     y2[1:n] is calculated array of the second derivatives of the
c     INTerpolating FUNCTION at the x[i]
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      INTEGER :: n,k,i
      INTEGER, PARAMETER :: nmax=2000
      REAL(rprec) :: x(*),y(*),y2(*)
      REAL(rprec) :: u(nmax)
      REAL(rprec) :: yp1,ypn
      REAL(rprec) :: yp1t,ypnt,p,sig,qn,un
      IF (n .gt. nmax)then
         WRITE(6,'("spline;  DIMENSIONal error")')
         STOP
      ENDIF
      IF (yp1.gt..99e30) THEN
        y2(1)=0.
        u(1)=0.
      ELSE IF(yp1.lt.-.99e30) THEN
         yp1t=(3*x(1)**2+x(2)*x(3)+
     $        x(2)*x(4)+x(3)*x(4)-2*x(1)*(x(2)+x(3)+x(4)))*
     $ y(1)/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4)))+
     $ (-x(1)**2+x(1)*x(3)+x(1)*x(4)-x(3)*x(4))*y(2)/
     $ ((x(1)-x(2))*(x(2)-x(3))*(x(2)-x(4)))+
     $ (x(1)**2-x(1)*x(2)-x(1)*x(4)+x(2)*x(4))*y(3)/
     $ ((x(1)-x(3))*(x(2)-x(3))*(x(3)-x(4)))+
     $ (-x(1)**2+x(1)*x(2)+x(1)*x(3)-x(2)*x(3))*y(4)/
     $ ((x(1)-x(4))*(x(2)-x(4))*(x(3)-x(4)))
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1t)
      ELSE
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF
      DO  i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      END DO
      IF (ypn.gt..99e30) THEN
        qn=0.
        un=0.
      ELSE IF(ypn.lt.-.99e30) THEN
         ypnt=(-(x(-2+n)*x(-1+n))+x(-2+n)*x(n)+x(-1+n)*x(n)-x(n)**2)*
     $  y(-3+n)/
     $  ((-x(-3+n)+x(-2+n))*(-x(-3+n)+x(-1+n))*
     $  (-x(-3+n)+x(n)))+
     $  (x(-3+n)*x(-1+n)-x(-3+n)*x(n)-x(-1+n)*x(n)+x(n)**2)*
     $  y(-2+n)/
     $  ((-x(-3+n)+x(-2+n))*(-x(-2+n)+x(-1+n))*
     $  (-x(-2+n)+x(n)))+
     $  (-(x(-3+n)*x(-2+n))+x(-3+n)*x(n)+x(-2+n)*x(n)-x(n)**2)*
     $  y(-1+n)/
     $  ((-x(-3+n)+x(-1+n))*(-x(-2+n)+x(-1+n))*
     $  (-x(-1+n)+x(n)))+
     $        (x(-3+n)*x(-2+n)+x(-3+n)*x(-1+n)+x(-2+n)*x(-1+n)-
     $        2*(x(-3+n)+x(-2+n)+x(-1+n))*x(n)+3*x(n)**2)*y(n)/
     $        ((-x(-3+n)+x(n))*(-x(-2+n)+x(n))*(-x(-1+n)+x(n)))
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypnt-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ELSE
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      DO  k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      END DO
      RETURN
      END
c=======================================================================
      SUBROUTINE splINT(xa,ya,y2a,n,x,y,nlow)
c-----------------------------------------------------------------------
c     cubic spline evaluator--spline must be CALLed first to evaluate
c     y2a
c     ya is a FUNCTION of xa--both are arrays of LENgth n
c     ya2[1:n] contains spline coefficients calculated in spline
c     x is the argument of y[x] WHERE y is to be evaluated
c     y=y[x] is the RETURNed value
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      INTEGER :: n,nlow,nhigh,nave
      REAL(rprec) ::  xa(*),ya(*),y2a(*)
      REAL(rprec) :: x,y
      REAL(rprec) :: a,b,h

      IF(nlow.eq.0) THEN
         nlow=1
         nhigh=n
 100     IF (nhigh-nlow.gt.1) THEN
            nave=(nhigh+nlow)/2
            IF(xa(nave).gt.x)then
               nhigh=nave
            ELSE
               nlow=nave
            ENDIF
            go to 100
         ENDIF
      ELSE
         nhigh=nlow+1
      ENDIF
      h=xa(nhigh)-xa(nlow)
      IF (h.eq.0.) THEN
         WRITE(6,'( "two xa = in  splinet")')
         STOP
      ENDIF
      a=(xa(nhigh)-x)/h
      b=1.-a
      y=a*ya(nlow)+b*ya(nhigh)+
     $     ((a-2.)*y2a(nlow)-(a+1.)*y2a(nhigh))*a*b*h*h/6.
      RETURN
      END
c=======================================================================
      SUBROUTINE zspline(xa,ya,y2a,n,zg,ng,za)
c-----------------------------------------------------------------------
c     zspline INTegrates the cubic spline of ya[xa]
c     it assumes that spline has already been CALLed to evaluate y2a
c     xa[1:n], ya[1:n], y2a[1:n]
c     In Mathematica notation z[i]=zg+Integrate[y[x],{x,x[ng],x[i]}]
c     WHERE both zg and ng are input quantities,
c     za is calculated here. 
c     za can THEN be USEd in zsplINT to determine z at a 
c     specific x
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
c
      REAL(rprec) :: xa(n),ya(n),y2a(n),zg,za(n)
      INTEGER :: n,ng
c
      INTEGER :: j
      REAL(rprec) :: const
c
      IF (ng .lt. 0 .or. ng .gt. n)stop 'zspline: wrong ng'
c
      za(1)=0.
      DO j=2,n
         za(j)=za(j-1)+0.5*(xa(j)-xa(j-1))*(ya(j)+ya(j-1))-
     >        (xa(j)-xa(j-1))**3*(y2a(j)+y2a(j-1))/24.0
      ENDdo
c
      const=zg-za(ng)
      DO j=1,n
         za(j)=za(j)+const
      ENDdo
c
      RETURN
      END
c=======================================================================
      SUBROUTINE zsplINT(xa,ya,y2a,za,n,x,y,yp,z)
c-----------------------------------------------------------------------
c     evaluate cubic spline to determine FUNCTION (y),
c     derivative (yp) and INTegral (z) at location x.
c     first spline must be CALLed to obtain y2a and
c     zspline must be CALLed to obtain za
c-----------------------------------------------------------------------
      DIMENSION xa(n),ya(n),y2a(n),za(n)
      klo=1
      khi=n
1     IF (khi-klo.gt.1) THEN
        k=(khi+klo)/2
        IF(xa(k).gt.x)then
          khi=k
        ELSE
          klo=k
        ENDIF
      GOTO 1
      ENDIF
      h=xa(khi)-xa(klo)
      IF (h.eq.0.) THEN
         WRITE(6,'(a)')"bad xa input"
         STOP
      ENDIF
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      yp=(ya(khi)-ya(klo))/h-(3.*a**2-1.)/6.*h*y2a(klo)
     >     +(3.*b**2-1.)/6.*h*y2a(khi)
      z=za(klo)+0.5*h*ya(klo)*(1.-a**2)+0.5*ya(khi)*h*b**2+
     >     y2a(klo)*h**3*(2.*a**2-a**4-1.)/24.-
     >     y2a(khi)*h**3*(2.*b**2-b**4)/24.
      RETURN
      END
