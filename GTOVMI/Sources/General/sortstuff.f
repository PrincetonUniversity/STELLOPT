      SUBROUTINE sortr(xp,zp,gsq,npc,xmx,zmx,xset,zset,j,jsep)
      USE precision
      IMPLICIT NONE
      INTEGER :: j, jsep
      REAL(rprec) :: xp(*),zp(*),gsq(*), xmx, zmx, xset, zset
      REAL(rprec) :: dr1, dr2, dz1, dz2, DOt, cross
      INTEGER :: i, npc
c ---------------------------------------
c order points counterclockwise IF needed
c ---------------------------------------
      i=0
   10 i=i+1
      dr1=xp(i)-xmx
      dr2=xp(i+1)-xmx
      dz1=zp(i)-zmx
      dz2=zp(i+1)-zmx
      DOt=dr1*dr2+dz1*dz2
      cross=dr1*dz2-dr2*dz1
      IF(cross.eq.0.) go to 10
      IF(dot.eq.0.) go to 10
      IF(cross/dot.gt.0.) go to 20
c ------------------------------------------------
c reorder the points to make them counterclockwise
c ------------------------------------------------
      CALL sorter(xp,zp,gsq,npc)
   20 CONTINUE
c ---------------------
c delete the last point
c ---------------------
      npc=npc-1
c --------------------------------------------------------------------
c sort so first poINT is one with smallest positive angle about center
c xmx, zmx for j .lt. jsep
c pick poINT CLOSEst to xset, zset for j .ge. jsep
c --------------------------------------------------------------------
      CALL sortog(xp,zp,gsq,npc,xmx,zmx,xset,zset,j,jsep)
c -----------------------------------------
c reset new last poINT to equal first point
c -----------------------------------------
      npc=npc+1
      xp(npc)=xp(1)
      zp(npc)=zp(1)
      gsq(npc)=gsq(1)
      END SUBROUTINE sortr
      SUBROUTINE sorter(x,y,d,n)
      USE precision
      IMPLICIT NONE
cray  lcm (x),(y),(d)
        REAL(rprec) :: x(*),y(*),d(*)
        REAL(rprec) :: t
        INTEGER :: n, nsort, nl, i
        nsort=n/2
        nl=n+1
        DO 100 i=1,nsort
        nl=nl-1
        t=x(i)
        x(i)=x(nl)
        x(nl)=t
        t=y(i)
        y(i)=y(nl)
        y(nl)=t
        t=d(i)
        d(i)=d(nl)
        d(nl)=t
  100   CONTINUE
      END SUBROUTINE sorter
      SUBROUTINE sortog(x,y,d,n,xs,ys,xset,zset,j,jsep)
      USE precision
      IMPLICIT NONE
      INTEGER :: j, jsep,n 
      REAL(rprec) ::  xs, ys, xset, zset
      REAL(rprec) :: pi, ang, angp, tx, ty, td, dist, dt
      INTEGER :: i, npc, isp, im1, nm1, jj
      REAL(rprec) :: x(*),y(*),d(*)
      pi=4*ATAN(1._dbl)
      IF(j.ge.jsep) go to 110
      ang=1.e10
        DO 50 i=1,n
      angp=atan2(y(i)-ys,x(i)-xs)
      angp=ABS(angp)
      IF(angp.gt.ang) go to 50
   30   isp=i
      ang=angp
   50   CONTINUE
   60 IF(isp.eq.1) RETURN
        im1=isp-1
        nm1=n-1
        DO 100 i=1,im1
        tx=x(1)
        ty=y(1)
        td=d(1)
        DO 90 jj=1,nm1
        x(jj)=x(jj+1)
        y(jj)=y(jj+1)
        d(jj)=d(jj+1)
   90   CONTINUE
        x(n)=tx
        y(n)=ty
        d(n)=td
  100   CONTINUE
        RETURN
c ----------------------------------------------
c sorting PROCEDURE for points inside separatrix
c find points CLOSEst to xset, zset
c ----------------------------------------------
  110 dist=1.e10
      DO 130 i=1,n
      dt=(y(i)-zset)**2+(x(i)-xset)**2
      IF(dt.gt.dist) go to 130
      isp=i
      dist=dt
  130 CONTINUE
      go to 60
      END SUBROUTINE sortog

