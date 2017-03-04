      SUBROUTINE rayatt(x,nw,y,nh,cspln,nxd,dx,dy,xaxd,yaxd,psiaxd,
     $     xmin,xmax,ymin,ymax,psivl,thet,xc,yc,pds,ifail)
      USE precision
      IMPLICIT NONE
      INTEGER :: nw,nh,nxd
      REAL(rprec) :: x(nw),y(nh),cspln(2,nxd,1),dx,dy,xaxd,yaxd,psiaxd
      REAL(rprec) :: xmin,xmax,ymin,ymax
      REAL(rprec) :: psivl,thet,xc,yc,pds(6)
      INTEGER IFail
c...........................................................................
      INTEGER, PARAMETER :: kmax=129, nmax=32
      REAL(rprec), PARAMETER :: rndoff0=1.e-9
      REAL(rprec), PARAMETER :: serrt=1.e-9 ,serrs=10.*serrt
      INTEGER :: izone,iflag,kountr,newti,ier
      REAL(rprec) :: pi, COSt, SINt, sgn
      REAL(rprec) :: a, bincp, x1, y1, x2, y2, psi1, psi2
      REAL(rprec) :: xn,yn,delx,dely,serr,xnorm,pnorm,dpsi,dpsids,rndoff
      REAL(rprec) :: rangex,rangey,dum
      data dum/0./
c
      pi = 4*ATAN(1._dbl)
      COSt=COS(thet)
      SINt=SIN(thet)
      izone=INT(2*thet/pi+SIGN(0.5_dbl,thet)) 
      IFlag=MOD(izone,2)
c
      IF (iflag .eq. 0)then
         a=sint/cost
         bincp=yaxd-a*xaxd
         sgn=SIGN(1._dbl,cost)
      ELSE
         a=cost/sint
         bincp=xaxd-a*yaxd
         sgn=SIGN(1._dbl,sint)
      ENDIF
c
c     DO the sliding window search
      x1=xaxd
      y1=yaxd
      psi1=psiaxd
      DO kountr=1,kmax
         IF (iflag .eq. 0)then
            x2=x1+sgn*dx
            y2=a*x2+bincp
         ELSE
            y2=y1+sgn*dy
            x2=a*y2+bincp
         ENDIF
         rangex=(x2-xmin)*(x2-xmax)
         rangey=(y2-ymin)*(y2-ymax)
         IF (rangex .gt. 0. .or. rangey .gt. 0.)then
            IFail=1
            CALL errray1(1,psivl,thet,x2,y2,dum,dum,dum) ! terminates execution
         ENDIF
c         CALL dbcevl2(x,nw,y,nh,cspln,nxd,x2,y2,pds,ier,1)
c         CALL bcsplINT(x,nw,y,nh,nxd,cspln,x2,y2,pds,1,ier)
         CALL bcsplINT(x,y,cspln,nw,nh,nxd,x2,y2,pds,1,ier)
         psi2=pds(1)
         IF ((psivl-psi1)*(psivl-psi2) .le. 0.)go to 10
         x1=x2
         y1=y2
         psi1=psi2
      ENDdo
      IFail=2
      CALL errray1(2,psivl,thet,x2,y2,dum,dum,dum) ! terminates execution
 10   CONTINUE                  ! DOne with sliding search
c     
      xnorm=SQRT(xmax*xmax+xmin*xmin+ymax*ymax+ymin*ymin)
      pnorm=ABS(psi1)+ABS(psi2)
      rndoff=rndoff0*(pnorm/xnorm)
      IF (iflag .eq. 0)then
         xn=x1+0.5*sgn*dx
         yn=a*xn+bincp
      ELSE
         yn=y1+0.5*sgn*dy
         xn=a*yn+bincp
      ENDIF
      DO newti=1,nmax
c         CALL dbcevl2(x,nw,y,nh,cspln,nxd,xn,yn,pds,ier,3)
         CALL bcsplINT(x,y,cspln,nw,nh,nxd,xn,yn,pds,3,ier)
         dpsids=pds(2)*cost+pds(3)*sint
         dpsi=pds(1)-psivl
         IF (ABS(dpsids) .gt. rndoff)then
            serr=-dpsi/dpsids
         ELSE
            serr=-dpsi/rndoff
         ENDIF
         IF (ABS(serr) .lt. serrt)go to 20
         delx=serr*cost
         dely=serr*sint
         xn=xn+delx
         yn=yn+dely
      ENDdo
      IFail=3
      CALL errray1(3,psivl,thet,xn,yn,serr,serrt,serrs)
 20   CONTINUE
      xc=xn
      yc=yn
      IFail=0
c      WRITE (3,'(2e15.6,2i5)')psivl,thet,kountr,newti
c
      RETURN
      END
